
import shutil
import numpy as np
import pandas as pd

from gtrspmix_utils.model_handle import parse_model,make_fmix
from gtrspmix_utils.iqtree_run import iqtree_est

def kmeans_profile(args, l_freq):
    from sklearn.cluster import KMeans
    def parse_freq(l_freq):
        l_profiles = []
        for i in l_freq:
            l = [float(j) for j in i[1].split(' ')]
            l_profiles.append(l)
        return np.array(l_profiles)

    np_freq = parse_freq(l_freq)
    kmeans_result = KMeans(n_clusters=args.kmeans, random_state=0).fit(np_freq)
    labels = list(kmeans_result.labels_)

    d_cluster = {}
    for i in range(args.kmeans):
        d_cluster[i+1] = []
    c = 1
    for i in labels:
        d_cluster[i+1].append(c)
        c += 1
    return d_cluster


def sppc_profile(args):
    # run iqtree with few profiles
    l_freq, _, _, _ = parse_model(args.nexus_few)
    l_fmix = make_fmix(None, l_freq)
    num_freq_few = len(l_freq)

    model = f'{args.m_gtr20}+{l_fmix[0][1]}+{args.mrate}'
    dir_few = args.outdir / f'{args.m_gtr20}+Finput{num_freq_few}+{args.mrate}'
    dir_few.mkdir(parents=True, exist_ok=True)
    prefix_few = dir_few / 'train_est'
    iqtree_est(args, (args.outdir / 'input_train.fasta'), (args.outdir / 'input.tree'), args.nexus_few, model, prefix_few)
    siteprob_few = prefix_few.with_suffix('.siteprob')

    # run iqtree with many profiles
    l_freq, _, _, _ = parse_model(args.nexus_many)
    l_fmix = make_fmix(None, l_freq)
    num_freq_many = len(l_freq)

    model = f'{args.m_gtr20}+{l_fmix[0][1]}+{args.mrate}'
    dir_many = args.outdir / f'{args.m_gtr20}+Finput{num_freq_many}+{args.mrate}'
    dir_many.mkdir(parents=True, exist_ok=True)
    prefix_many = dir_many / 'train_est'
    iqtree_est(args, (args.outdir / 'input_train.fasta'), (args.outdir / 'input.tree'), args.nexus_many, model, prefix_many)
    siteprob_many = prefix_many.with_suffix('.siteprob')

    dir_0 = args.outdir / 'train_est_0'
    dir_0.mkdir(parents=True, exist_ok=True)
    l_extensions = ['.ckp.gz', '.iqtree', '.log', '.siteprob', '.treefile']
    for ext in l_extensions:
        src = prefix_many.with_suffix(ext)
        if src.exists():
            shutil.copy(src, dir_0 / src.name)

    df_many = pd.read_csv(siteprob_many, sep='\t', comment='#')
    df_few = pd.read_csv(siteprob_few, sep='\t', comment='#')
    df_many = df_many.iloc[:, 1:]
    df_few = df_few.iloc[:, 1:]

    max_df_many = df_many.idxmax(axis=1).rename('max_df_many')
    max_df_few = df_few.idxmax(axis=1).rename('max_df_few')
    siteprob_all = pd.concat([max_df_many, max_df_few], axis=1)

    df_co = siteprob_all.groupby('max_df_many')['max_df_few'].agg(lambda x: x.mode().iloc[0])
    df_co = df_co.reset_index()
    df_co.columns = ['max_df_many', 'max_df_few']
    df_co['max_df_many'] = df_co['max_df_many'].str.replace(r"p", "", regex=True).astype(int)
    df_co['max_df_few'] = df_co['max_df_few'].str.replace(r"p", "", regex=True).astype(int)

    l_cluster = df_co.groupby('max_df_few')['max_df_many'].apply(list).tolist()

    d_cluster = {}
    for i in range(len(l_cluster)):
        d_cluster[i+1] = sorted(l_cluster[i])

    return d_cluster
