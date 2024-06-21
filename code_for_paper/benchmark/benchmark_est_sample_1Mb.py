import time
import numpy as np
import pandas as pd
import math
from random import sample
from hicpep import peptools

def flip_tracks(track1_np: np.ndarray, track2_np: np.ndarray):
    if np.corrcoef(track1_np[~np.isnan(track1_np)], track2_np[~np.isnan(track2_np)])[0][1] < 0:
        track2_np = -track2_np
    return track1_np, track2_np

def read_file(pearson_path):
    pearson_np = peptools.read_pearson(pearson=pearson_path)

    if len(pearson_np) != len(pearson_np[0]):
        print("Pearson matrix has a different number of rows and columns")

    pearson_np = pearson_np.astype('float64')
    diag = np.diag(pearson_np)
    diag_valid = ~np.isnan(diag)
    ixgrid = np.ix_(diag_valid, diag_valid) # Record the position of the valid sub-matrix. 
    pearson_np = pearson_np[ixgrid]
    
    has_nan = np.isnan(pearson_np).any()
    if has_nan:
        print("NaN entries still exist in the Pearson matrix.")

    return pearson_np, diag_valid

def hicpep_est_sample(pearson_path, proportion):
    pearson_np, diag_valid = read_file(pearson_path) 
    pearson_np -= pearson_np.mean(axis=1, keepdims=True)
    start = time.time()
    n = len(pearson_np[0])
    sample_indexes = sample(list(range(n)), math.floor(n * proportion))
    partial_cov_np = pearson_np @ pearson_np[sample_indexes].T / n
    partial_cov_abs_sum = [(index, np.sum(np.abs(row))) for index, row in enumerate(partial_cov_np.T)] 
    sorted_partial_cov_abs_sum = sorted(partial_cov_abs_sum, key=lambda x: x[1], reverse=True) # Sorted from the maximum to the minimum 
    est_np = np.full(len(diag_valid), np.nan)
    est_np[diag_valid] = partial_cov_np.T[sorted_partial_cov_abs_sum[0][0]]
    end = time.time()
    print(f"Time spent for creating the Estimated PC1-pattern by finding the CxMax in the partial-covariance matrix through sampling (seconds): {end - start}")

    return est_np

if __name__ == "__main__":
    pearson_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/rao_2014/juicer_outputs/gm12878/1000000/pearsons/pearson_chr2.txt"
    juicer_pc1_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/rao_2014/juicer_outputs/gm12878/1000000/eigenvector/pc1_chr2.txt"

    est_np_2 = hicpep_est_sample(pearson_path, proportion=0.1)

    juicer_pc1_df = pd.read_table(juicer_pc1_path, header=None)
    juicer_pc1_np = juicer_pc1_df.values.flatten()
    juicer_pc1_np, est_np_2  = flip_tracks(track1_np=juicer_pc1_np, track2_np=est_np_2)
    similarity_info = peptools.calc_similarity(track1_np=juicer_pc1_np, track2_np=est_np_2)

    print(similarity_info)