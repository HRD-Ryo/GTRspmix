
import sys
import subprocess
import logging
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from decimal import Decimal, ROUND_HALF_UP

def soft_part(args, d_cluster, num_it, pre_run, target):
    def make_clustered_weights(d_cluster, df):
        df_clu = pd.DataFrame()
        for i in d_cluster:
            columns = ['p' + str(num) for num in d_cluster[i]]
            df_clu[f'cluster{i}'] = df[columns].sum(axis=1)
        return df_clu

    def extract_fasta(fasta_file, df):
        new_records = []
        gap_records = []
        for record in SeqIO.parse(fasta_file, 'fasta'):
            new_seq = ''
            for site, weight in df.items():
                num = Decimal(str(weight))
                num = int(num.quantize(Decimal('1'), rounding=ROUND_HALF_UP))
                new_seq = new_seq + str(record.seq)[site] * num
            if set(new_seq) <= {'-', 'X', 'x'}:
                gap_records.append(record.id)
            else:
                record.seq = Seq(new_seq)
                new_records.append(record)
        return new_records, gap_records

    # input files
    fasta = (args.outdir / 'input_train.fasta')
    pre_run = args.outdir / f'train_est_{pre_run}' / 'train_est'
    siteprob = pre_run.with_suffix('.siteprob')
    tree = pre_run.with_suffix('.treefile')

    df = pd.read_csv(siteprob, sep='\t', comment='#')
    df = df.iloc[:, 1:]
    if target == 'GTR':
        df_clu = make_clustered_weights(d_cluster, df)
        df_clu = df_clu * args.scale_gtr
        work_dir = args.outdir / f'GTR_opt_{num_it}'
    elif target == 'Profile':
        df_clu = df
        df_clu = df_clu * args.scale_profile
        work_dir = args.outdir / f'Profile_opt_{num_it}'

    # split alignment
    for c in df_clu.columns:
        cluster_dir = work_dir / c
        cluster_dir.mkdir(parents=True, exist_ok=True)
        df_f = df_clu[df_clu[c] >= 0.5][c]
        if len(df_f) > 5:
            new_records, gap_records = extract_fasta(fasta, df_f)
            SeqIO.write(new_records, (cluster_dir / 'input_train.fasta'), 'fasta')

            # remove gap seq from tree
            with (cluster_dir / 'guide.treefile').open(mode='w') as fo:
                l_cmd = args.gotree + ['prune', '-i', str(tree)] + gap_records
                subprocess.run(l_cmd, stdout=fo)

        else:
            logger = logging.getLogger('gtrspmix')
            logger.error(f'Soft-partitioning Failed: {c} has too few sites ({len(df_f)}).')
            logger.error(f'  - use larger scale value or remove this cluster/profile.')
            sys.exit(1)
