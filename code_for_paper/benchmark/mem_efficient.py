import time
import sys
import numpy as np
import pandas as pd
import math
from random import sample
from hicpep import peptools
from scipy import sparse
np.set_printoptions(threshold=sys.maxsize)

def flip_tracks(track1_np: np.ndarray, track2_np: np.ndarray):
    if np.corrcoef(track1_np[~np.isnan(track1_np)], track2_np[~np.isnan(track2_np)])[0][1] < 0:
        track2_np = -track2_np
    return track1_np, track2_np

### `store_oe_sparse` is not included in benchmark.
def store_oe_sparse(oe_path):
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values
    del oe_df
    oe_np = oe_np.astype('float64')
    valid = oe_np.any(axis=1)
    ixgrid = np.ix_(valid, valid) # Record the position of the valid sub-matrix.
    oe_np = oe_np[ixgrid]
    x = sparse.csr_matrix(oe_np)
    sparse.save_npz('/tmp/oe_sparse.npz', x)
    return

def load_oe_sparse():
    x = sparse.load_npz('/tmp/oe_sparse.npz')
    print(x)
    return

def mem_efficient_sampling(proportion=0.1):
    x = sparse.load_npz('/tmp/oe_sparse.npz')
    start = time.time()
    index_s = 0
    std = []
    c = []
    n = x.get_shape()[0]

    for i in range(n):
        tmp = np.std(x[i].toarray())
        std.append(tmp)

    for i in range(n):
        tmp = x[i].mean()
        c.append(tmp)

    std = np.array(std)
    c = np.array(c)
    I = np.full(n, 1)

    max = 0
    est_np = np.array([]) 
    sample_indexes = sample(list(range(n)), math.floor(n * proportion))

    for k in sample_indexes:
        index_s = k # The target column (or row) in the Pearson's covariance matrix.
        x_s = x[index_s].toarray()[0]
        x_s -= x_s.mean()
        alpha = I @ x_s
        corr_s = (x.dot(x_s) - c * alpha) / ((std[index_s] * std) * n)
        corr_s -= corr_s.mean()
        cov_s = []

        for i in range(n):
            x_i = x[i].toarray()[0]
            x_i -= x_i.mean()
            alpha = I @ x_i
            corr_i = (x.dot(x_i) - c * alpha) / ((std[i] * std) * n)
            corr_i -= corr_i.mean()
            entry = corr_i @ corr_s / n # Calculate each entry of the target column (or row) in the Pearson's covariance matrix.
            cov_s.append(entry)

        cov_s = np.array(cov_s)
        track_sum_abs = np.sum(np.abs(cov_s))

        if track_sum_abs > max:
            max = track_sum_abs
            est_np = cov_s.copy()

    end = time.time()
    print(f"Time spent for creating the Estimated PC1-pattern by sampling and new algo (seconds): {end - start}")
    print(f"mem_efficient cov: {est_np[:5]}")

    return est_np

if __name__ == '__main__':
    oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr2_chr2_100000_obsexp.txt"
    # store_oe_sparse(oe_path) # Not include in benchmarking 
    est_np = mem_efficient_sampling(proportion=0.1)

    pc1_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/lieberman_2009/eigenvectors/GM-combined.ctg2.ctg2.1000000bp.hm.eigenvector.tab"
    pc1_df = pd.read_table(pc1_path, header=None, sep="\s+")
    pc1_df = pc1_df.iloc[:, [2]]
    pc1_np = pc1_df.values # Turn into numpy format
    pc1_np = pc1_np.flatten() # Turn into 1D vector
    pc1_np = pc1_np[pc1_np != 0]

    print(len(est_np))
    print(len(pc1_np))

    pc1_np, est_np  = flip_tracks(track1_np=pc1_np, track2_np=est_np)
    similarity_info = peptools.calc_similarity(track1_np=pc1_np, track2_np=est_np)
    print(similarity_info)