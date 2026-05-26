
import numpy as np
import pandas as pd
from Bio import AlignIO


def get_site_F(fasta_file):
    ali = AlignIO.read(fasta_file, 'fasta')
    np_ali = np.array([list(rec) for rec in ali])
    aa = np.array(list('ARNDCQEGHILKMFPSTWYV'))
    
    # Count amino acids using broadcasting
    # Shape: (20, N_sequences, L_sites) -> Sum over sequences (axis=1)
    counts = (np_ali[None, :, :] == aa[:, None, None]).sum(axis=1)
    
    # # Total counts
    # column_totals = counts.sum(axis=0)
    
    # # avoiding devide by 0
    # sitefreq = np.divide(
    #     counts, 
    #     column_totals, 
    #     out=np.zeros_like(counts, dtype=float), 
    #     where=column_totals != 0
    # )
    # return sitefreq
    return counts


def average_F(sitefreq, siteprob):
    df = pd.read_csv(siteprob, sep="\t")
    # Shape: (L_sites, K_components)
    siteprob = df.iloc[:, 1:].values
    
    # Compute weighted sum using Dot Product
    # (20, L) @ (L, K) -> (20, K)
    # This sums (freq * prob) over all sites (L)
    weighted_sum = np.dot(sitefreq, siteprob)

    # Normalize
    profiles = weighted_sum / weighted_sum.sum(axis=0)

    # Set lower bound & re-normalize
    lower = 1e-6
    profiles = np.maximum(profiles, lower)
    profiles = profiles / profiles.sum(axis=0)
    
    return profiles

