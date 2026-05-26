
import logging
import subprocess

def iqtree_est(args, fasta, tree, nex, model, prefix):
    l_cmd = args.iqtree + [\
    '--seqtype', 'AA', \
    '--seed', args.seed, \
    '-s', str(fasta), \
    '-te', str(tree), \
    '-mdef', str(nex), \
    '-m', model, \
    '-mwopt', \
    '-wspm', \
    '-nt', str(args.nt), \
    '-mem', args.mem, \
    '--epsilon', str(args.epsilon_theta), \
    '--prefix', str(prefix), \
    '-prec', '10', \
    '--quiet' ]
    subprocess.run(l_cmd, capture_output=True, text=True)
    # try:
    #     subprocess.run(l_cmd, capture_output=True, text=True, check=True)
    # except subprocess.CalledProcessError as e:
    #     logger = logging.getLogger('gtrspmix')
    #     logger.error('IQ-TREE Failed')
    #     logger.error(f'Command: {e.cmd}')
    #     logger.error(f'Stderr: {e.stderr}')


def iqtree_opt_GTR(args, l_fmix, l_exch, num_it, rate_param):
    for i in range(len(l_fmix)):
        l_cmd = args.iqtree + [\
        '--seed', args.seed, \
        '--seqtype', 'AA', \
        '-s', str(args.outdir / f'GTR_opt_{num_it}' / f'cluster{i+1}' / 'input_train.fasta'), \
        '-te', str(args.outdir / f'GTR_opt_{num_it}' / f'cluster{i+1}' / 'guide.treefile'), \
        '-m', f'GTR20+{l_fmix[i][0]}+{args.mrate}{{{rate_param}}}', \
        '--link-exchange-rates', \
        '--gtr20-model', f'{l_exch[i][0]}', \
        '-mdef', str(args.outdir / f'model_{num_it-1}.nex'), \
        '--epsilon', str(args.epsilon_gtr), \
        '-blfix', \
        '-safe', \
        '-nt', str(args.nt), \
        '-mem', args.mem, \
        '-prec', '10', \
        '--quiet'
        ]
        subprocess.run(l_cmd, capture_output=True, text=True)
        # try:
        #     subprocess.run(l_cmd, capture_output=True, text=True, check=True)
        # except subprocess.CalledProcessError as e:
        #     logger = logging.getLogger('gtrspmix')
        #     logger.error('IQ-TREE Failed')
        #     logger.error(f'Command: {e.cmd}')
        #     logger.error(f'Stderr: {e.stderr}')



def iqtree_opt_Profile(args, nex, d_cluster, l_freq, l_exch, num_it, rate_param):
    rev = {x: k for k, lst in d_cluster.items() for x in lst}
    for i in range(len(l_freq)):
        l_cmd = args.iqtree + [\
        '--seed', args.seed, \
        '--seqtype', 'AA', \
        '-s', str(args.outdir / f'Profile_opt_{num_it}' / f'p{i+1}' / 'input_train.fasta'), \
        '-te', str(args.outdir / f'Profile_opt_{num_it}' / f'p{i+1}' / 'guide.treefile'), \
        '-m', f'{l_exch[int(rev[i+1])-1][0]}+FO+{args.mrate}{{{rate_param}}}', \
        '-mdef', str(nex), \
        '--epsilon', str(args.epsilon_profile), \
        '-blfix', \
        '-safe', \
        '-nt', str(args.nt), \
        '-mem', args.mem, \
        '-prec', '10', \
        '--quiet'
        ]
        subprocess.run(l_cmd, capture_output=True, text=True)
        # try:
        #     subprocess.run(l_cmd, capture_output=True, text=True, check=True)
        # except subprocess.CalledProcessError as e:
        #     logger = logging.getLogger('gtrspmix')
        #     logger.error('IQ-TREE Failed')
        #     logger.error(f'Command: {e.cmd}')
        #     logger.error(f'Stderr: {e.stderr}')

