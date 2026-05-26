#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# Written by Ryo Harada


from pathlib import Path
import os
import subprocess
import sys
sys.dont_write_bytecode = True
import datetime
import shutil
import json

# my packages
from gtrspmix_utils.logger import setup_logger
from gtrspmix_utils.constants import __version__,__ascii__
from gtrspmix_utils.cli import get_options,check_args
from gtrspmix_utils.model_handle import parse_model,write_model,init_model,extract_LL_rate_modelw, extract_GTR,extract_Profile,make_fmix,get_Profile_F
from gtrspmix_utils.clustering import kmeans_profile,sppc_profile
from gtrspmix_utils.iqtree_run import iqtree_est,iqtree_opt_GTR,iqtree_opt_Profile
from gtrspmix_utils.soft_partition import soft_part
from gtrspmix_utils.weighted_F import get_site_F,average_F


def resolve_model_path(models_dir, target):
    if target == None: 
        return target

    if target.exists():
        return target

    model_path = models_dir / f'{str(target)}.nex'
    if model_path.exists():
        return model_path

    return target


def estimate_theta(args, num_it, pre_run, nex, l_fmix, l_mix, pre_target):
    def make_fmix_w(l_fmix, d_fmix_comp):
        l_fmix_w = []
        for i in l_fmix:
            new_comp = []
            comp = i[1].replace('FMIX{', '').replace('}', '')
            for j in comp.split(','):
                comp_name = j.split(':')[0]
                new_comp.append(comp_name + d_fmix_comp[comp_name])
            l_fmix_w.append([i[0], 'FMIX{' + ','.join(new_comp) + '}'])
        return l_fmix_w

    # run on train data
    if pre_target == None:
        est_dir = args.outdir / f'train_est_{num_it}'
    else:
        est_dir = args.outdir / f'train_est_{num_it}_{pre_target}'
    est_dir.mkdir(parents=True, exist_ok=True)

    # model = l_mix[0][0]+'+'+args.mrate
    model = l_mix[-1][0]+'+'+args.mrate
    prefix = est_dir / 'train_est'

    if pre_run == None:
        pre_tree = args.outdir / 'input.tree'
    else:
        pre_tree = args.outdir / f'train_est_{pre_run}' / 'train_est.treefile'

    iqtree_est(args, (args.outdir / 'input_train.fasta'), pre_tree, nex, model, prefix)
    train_LL, rate_param, l_mix_w, d_fmix_comp = extract_LL_rate_modelw(prefix.with_suffix('.iqtree'))
    l_fmix = make_fmix_w(l_fmix, d_fmix_comp)
    l_mix = [l_mix[0], [l_mix[0][0][:-3], l_mix_w]]

    return train_LL, rate_param, l_fmix, l_mix



def main():
    start_time = datetime.datetime.now()
    args = get_options(start_time)
    args.iqtree = args.iqtree.split(' ')
    args.gotree = args.gotree.split(' ')
    if args.overwrite:
        shutil.rmtree(args.outdir)
    (args.outdir).mkdir(parents=True, exist_ok=True)
    logger = setup_logger(args.outdir / 'GTRspmix_maker.log')

    print(f'{__ascii__}v{__version__} \n# ', file=sys.stderr)
    logger.info(f'GTRspmix v{__version__} started')

    logger.info('Command: '+' '.join(sys.argv))
    mode = check_args(args)
    logger.info(f'{mode} mode')
    logger.info('Random Seed: '+args.seed)

    shutil.copy(args.fasta, args.outdir / 'input_train.fasta')
    shutil.copy(args.tree, args.outdir / 'input.tree')

    # empirical models
    # ['C10', 'C20', 'C30', 'C40', 'C50', 'C60', 'S10C10', 'S20C20', 'S30C30', 'S10C60', 'S28C59', 'S28C60']:
    models_dir = Path(__file__).parent.resolve() / 'models'
    args.nexus = resolve_model_path(models_dir, args.nexus)
    args.nexus_few = resolve_model_path(models_dir, args.nexus_few)
    args.nexus_many = resolve_model_path(models_dir, args.nexus_many)

    # Initializing
    num_it = 0
    nex = args.outdir / f'model_{num_it}.nex'
    ## optimize profile only
    if mode == 'ProfileOnly':
        l_freq, _, l_exch, _ = parse_model(args.nexus)
        d_cluster = {str(1): list(range(1, len(l_freq) + 1))}
        l_fmix, l_exch, l_mix = init_model(args, d_cluster, l_freq, l_exch)
        if l_exch[0][1] == '': # dummy l_exch
            write_model(nex, l_freq, l_fmix, [], l_mix)
        else: # user defined l_exch
            write_model(nex, l_freq, l_fmix, l_exch, l_mix)

    ## start with empirical model
    elif mode == 'PreDefined':
        nex_src = models_dir / f'{args.model}.nex'
        json_src = models_dir / f'{args.model}.json'

        if not (nex_src.exists() and json_src.exists()):
            logger.error(f'Pre-defined model files not found.')
            logger.error(f'Check {nex_src.resolve()} and {json_src.resolve()}.')
            sys.exit(1)

        else:
            l_freq, _, l_exch, l_mix = parse_model(nex_src)
            d_cluster = json.loads(json_src.read_text())
            l_fmix = make_fmix(d_cluster, l_freq)
            write_model(nex, l_freq, l_fmix, l_exch, l_mix)

    ## re-start or user defined d_cluster
    elif mode == 'ReStart':
        l_freq, l_fmix, l_exch, l_mix = parse_model(args.nexus)
        d_cluster = json.loads(args.json.read_text())
        if l_fmix == [] or l_mix == []: # user defined d_cluster
            l_fmix, l_exch, l_mix = init_model(args, d_cluster, l_freq, l_exch)
            if l_exch[0][1] == '':
                write_model(nex, l_freq, l_fmix, [], l_mix)
            else:
                write_model(nex, l_freq, l_fmix, l_exch, l_mix)
        else: # real restart
            write_model(nex, l_freq, l_fmix, l_exch, l_mix)


        # d_cluster = json.loads(args.json.read_text())
        # l_fmix, l_exch, l_mix = init_model(args, d_cluster, l_freq, l_exch)
        # if l_exch[0][1] == '':
        #     write_model(nex, l_freq, l_fmix, [], l_mix)
        # else:
        #     write_model(nex, l_freq, l_fmix, l_exch, l_mix)

        # l_fmix, l_exch_new, l_mix_new = init_model(args, d_cluster, l_freq)
        # if l_mix == []: # user defined d_cluster
        #     l_mix = l_mix_new
        #     l_exch = l_exch_new
        #     write_model(nex, l_freq, l_fmix, [], l_mix)
        # else: # normal restart
        #     write_model(nex, l_freq, l_fmix, l_exch, l_mix)

    ## from scratch
    elif mode.split(' ')[0] == 'FromScratch':
        if mode.split(' ')[1] == 'SPPC':
            l_freq, _, l_exch, _ = parse_model(args.nexus_many)
            d_cluster = sppc_profile(args)
        elif mode.split(' ')[1] == 'Kmeans':
            l_freq, _, l_exch, _ = parse_model(args.nexus)
            d_cluster = kmeans_profile(args, l_freq)
        l_fmix, l_exch, l_mix = init_model(args, d_cluster, l_freq, l_exch)
        if l_exch[0][1] == '':
            write_model(nex, l_freq, l_fmix, [], l_mix)
        else:
            write_model(nex, l_freq, l_fmix, l_exch, l_mix)


    # output cluster info
    (args.outdir / 'd_cluster.json').write_text(json.dumps(d_cluster))
    num_exch = len(d_cluster)
    num_prof = 0
    for l in d_cluster.values():
        num_prof += len(l)

    if args.opt_gtr and (args.opt_profile_FO or args.opt_profile_F):
        logger.info(f'Optimizing {num_exch} GTRs AND {num_prof} profiles')
    elif args.opt_gtr and not (args.opt_profile_FO or args.opt_profile_F):
        logger.info(f'Optimizing {num_exch} GTRs with fixed {num_prof} profiles')
    elif not args.opt_gtr and (args.opt_profile_FO or args.opt_profile_F):
        logger.info(f'Optimizing {num_prof} profiles with fixed {num_exch} GTRs')

    pre_run = None
    train_LL, rate_param, l_fmix, l_mix = estimate_theta(args, num_it, pre_run, nex, l_fmix, l_mix, None)
    pre_run = 0
    _, _, tmp_l_exch, _ = parse_model(nex)
    write_model(nex, l_freq, l_fmix, tmp_l_exch, l_mix)

    logger.info(f'  it: {num_it}            \tlnL: {train_LL:.4f}')

    best_train_LL = train_LL
    best_it = num_it

    sitefreq = None

    # Optimizing
    continue_or_not = True
    while continue_or_not:
    # while True:
        # new iteration starts
        num_it += 1
        nex = args.outdir / f'model_{num_it-1}.nex' # previous nexus

        if args.opt_gtr: # opt GTRs
            soft_part(args, d_cluster, num_it, pre_run, 'GTR')
            iqtree_opt_GTR(args, l_fmix, l_exch, num_it, rate_param)
            pre_target = 'GTR'
            l_exch, l_mix = extract_GTR(args, d_cluster, l_freq, num_it)

            nex = args.outdir / f'model_{num_it}.nex'
            write_model(nex, l_freq, l_fmix, l_exch, l_mix)

            # opt theta params after GTR opt
            train_LL, rate_param, l_fmix, l_mix = estimate_theta(args, num_it, pre_run, nex, l_fmix, l_mix, pre_target)
            pre_run = f'{num_it}_GTR'
            write_model(nex, l_freq, l_fmix, l_exch, l_mix)
            logger.info(f'  it: {num_it}   (GTR)    \tlnL: {train_LL:.4f}')

            # Delete soft-partition fasta files
            if not args.keep_tmp:
                for c in d_cluster:
                    fasta_path = args.outdir / f'GTR_opt_{num_it}' / f'cluster{c}' / 'input_train.fasta'
                    if fasta_path.exists():
                        fasta_path.unlink()

        if args.opt_profile_FO: # opt Profiles FO
            soft_part(args, d_cluster, num_it, pre_run, 'Profile')
            iqtree_opt_Profile(args, nex, d_cluster, l_freq, l_exch, num_it, rate_param)
            pre_target = 'ProfileFO'
            l_freq, l_fmix, l_mix = extract_Profile(args, num_it, d_cluster, l_freq, l_exch)

            nex = args.outdir / f'model_{num_it}.nex'
            if l_exch[0][1] == '':
                write_model(nex, l_freq, l_fmix, [], l_mix)
            else:
                write_model(nex, l_freq, l_fmix, l_exch, l_mix)

            # opt theta params after profile opt FO
            train_LL, rate_param, l_fmix, l_mix = estimate_theta(args, num_it, pre_run, nex, l_fmix, l_mix, pre_target)
            pre_run = f'{num_it}_ProfileFO'
            if l_exch[0][1] == '':
                write_model(nex, l_freq, l_fmix, [], l_mix)
            else:
                write_model(nex, l_freq, l_fmix, l_exch, l_mix)
            logger.info(f'  it: {num_it} (ProfileFO)\tlnL: {train_LL:.4f}')

            # Delete soft-partition fasta files
            if not args.keep_tmp:
                for c in range(len(l_freq)):
                    fasta_path = args.outdir / f'Profile_opt_{num_it}' / f'p{c+1}' / 'input_train.fasta'
                    if fasta_path.exists():
                        fasta_path.unlink()

        if args.opt_profile_F: # opt profile F
            if sitefreq is None:
                sitefreq = get_site_F(args.outdir / 'input_train.fasta')
            siteprob = args.outdir / f'train_est_{pre_run}' / 'train_est.siteprob'
            profiles = average_F(sitefreq, siteprob)
            pre_target = 'ProfileF'

            l_freq, l_fmix, l_mix = get_Profile_F(profiles, d_cluster, l_freq, l_exch)

            nex = args.outdir / f'model_{num_it}.nex'
            if l_exch[0][1] == '':
                write_model(nex, l_freq, l_fmix, [], l_mix)
            else:
                write_model(nex, l_freq, l_fmix, l_exch, l_mix)

            # opt theta params after profile opt F
            train_LL, rate_param, l_fmix, l_mix = estimate_theta(args, num_it, pre_run, nex, l_fmix, l_mix, pre_target)
            pre_run = f'{num_it}_ProfileF'
            if l_exch[0][1] == '':
                write_model(nex, l_freq, l_fmix, [], l_mix)
            else:
                write_model(nex, l_freq, l_fmix, l_exch, l_mix)
            logger.info(f'  it: {num_it} (ProfileF) \tlnL: {train_LL:.4f}')


        # continue or stop
        train_dLL = train_LL - best_train_LL
        continue_or_not = (train_dLL > args.epsilon)
        if best_train_LL < train_LL:
            best_it = num_it
            best_train_LL = train_LL
        if args.forever:
            continue_or_not = True

    # finalizing
    shutil.copy(args.outdir / f'model_{best_it}.nex', args.outdir / f'model_best.nex')

    end_time = datetime.datetime.now()
    total_time = (end_time - start_time).total_seconds()
    total_time = f'{int(total_time // 3600)}h:{int((total_time % 3600) // 60)}m'
    logger.info(f'Total Run Time: {total_time}')


if __name__ == "__main__":
    main()
