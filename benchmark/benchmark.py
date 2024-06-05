import time
import numpy as np
import pandas as pd
import math
from random import sample
from sklearn.decomposition import PCA
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

@profile
def scikit_pc1(pearson_path):
    pearson_np, diag_valid = read_file(pearson_path) 
    start = time.time()
    pca = PCA(n_components=1)
    pca.fit(pearson_np)
    pc1_np = np.full(len(diag_valid), np.nan)
    pc1_np[diag_valid] = pca.components_[0]
    end = time.time()
    print(f"Time spent for performing PCA through Scikit-learn to get the PC1 (seconds): {end - start}")
    return pc1_np

@profile
def hicpep_est_all(pearson_path):
    pearson_np, diag_valid = read_file(pearson_path) 
    start = time.time()
    pearson_np -= pearson_np.mean(axis=1, keepdims=True)
    n = len(pearson_np[0])
    cov_np = pearson_np @ pearson_np.T / n
    cov_abs_sum = [(index, np.sum(np.abs(row))) for index, row in enumerate(cov_np)] 
    sorted_cov_abs_sum = sorted(cov_abs_sum, key=lambda x: x[1], reverse=True)
    est_np = np.full(len(diag_valid), np.nan)
    est_np[diag_valid] = cov_np[sorted_cov_abs_sum[0][0]]
    end = time.time()
    print(f"Time spent for creating the Estimated PC1-pattern by finding the CxMax in the full-covariance matrix (seconds): {end - start}")
    return est_np

@profile
def hicpep_est_sample(pearson_path, proportion):
    pearson_np, diag_valid = read_file(pearson_path) 
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
    pearson_path = "/home/jordan990301/Projects/HiCPEP/code_for_paper/notebooks/data/gm12878_pearson_25000_chr2.txt"
    juicer_pc1_path = "/home/jordan990301/Projects/HiCPEP/code_for_paper/notebooks/data/gm12878_pc1_25000_chr2.txt" # Ground Truth.

    # pearson_path = "/home/jordan990301/Projects/HiCPEP/test/gm12878_1000000_pearson_chr1.txt"
    # juicer_pc1_path = "/home/jordan990301/Projects/HiCPEP/test/gm12878_1000000_pc1_chr1.txt" # Ground Truth.

    juicer_pc1_df = pd.read_table(juicer_pc1_path, header=None)
    juicer_pc1_np = juicer_pc1_df.values.flatten()

    pc1_np = scikit_pc1(pearson_path)
    juicer_pc1_np, pc1_np  = flip_tracks(track1_np=juicer_pc1_np, track2_np=pc1_np)
    similarity_info_1 = peptools.calc_similarity(track1_np=juicer_pc1_np, track2_np=pc1_np)

    est_np_1 = hicpep_est_all(pearson_path)
    juicer_pc1_np, est_np_1  = flip_tracks(track1_np=juicer_pc1_np, track2_np=est_np_1)
    similarity_info_2 = peptools.calc_similarity(track1_np=juicer_pc1_np, track2_np=est_np_1)

    est_np_2 = hicpep_est_sample(pearson_path, proportion=0.1)
    juicer_pc1_np, est_np_2  = flip_tracks(track1_np=juicer_pc1_np, track2_np=est_np_2)
    similarity_info_3 = peptools.calc_similarity(track1_np=juicer_pc1_np, track2_np=est_np_2)

    print(similarity_info_1)
    print(similarity_info_2)
    print(similarity_info_3)