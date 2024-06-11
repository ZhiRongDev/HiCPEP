import sys
import numpy as np
import pandas as pd
import math
from random import sample
from hicpep import peptools
from sklearn.decomposition import PCA
from scipy import sparse
np.set_printoptions(threshold=sys.maxsize)

def store_oe_sparse(oe_path):
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values
    del oe_df
    oe_np = oe_np.astype('float64')
    diag = np.diag(oe_np)
    diag_valid = diag != 0 
    ixgrid = np.ix_(diag_valid, diag_valid) # Record the position of the valid sub-matrix.
    oe_np = oe_np[ixgrid]
    x = sparse.csr_matrix(oe_np)
    sparse.save_npz('/tmp/oe_sparse.npz', x)
    return

def load_oe_sparse():
    x = sparse.load_npz('/tmp/oe_sparse.npz')
    print(x)
    return

def sum_zero_percent(oe_path):
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values 
    size = len(oe_np) * len(oe_np)
    zero = size - np.count_nonzero(oe_np)
    print(zero)
    print(size)
    print(zero / size)
    return

@profile
def test(oe_path):
    index_s = 10
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values
    del oe_df
    oe_np = oe_np.astype('float64')
    diag = np.diag(oe_np)
    diag_valid = diag != 0 
    ixgrid = np.ix_(diag_valid, diag_valid) # Record the position of the valid sub-matrix.
    oe_np = oe_np[ixgrid]
    corr = np.corrcoef(oe_np)
    cov = np.cov(corr, bias=True)
    corr[index_s] -= corr[index_s].mean()

    print(f"test corr: {corr[index_s][:5]}")
    print(f"test cov: {cov[index_s][:5]}")
    return

@profile
def mem_efficient_random1():
    index_s = 10
    x = sparse.load_npz('/tmp/oe_sparse.npz')
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
        entry = corr_i @ corr_s / n 
        cov_s.append(entry)

    cov_s = np.array(cov_s)

    print(f"mem_efficient corr: {corr_s[:5]}")
    print(f"mem_efficient cov: {cov_s[:5]}")
    return

@profile
def normal(oe_path):
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values
    del oe_df
    oe_np = oe_np.astype('float64')
    diag = np.diag(oe_np)
    diag_valid = diag != 0 
    ixgrid = np.ix_(diag_valid, diag_valid) # Record the position of the valid sub-matrix.
    oe_np = oe_np[ixgrid]
    pearson_np = np.corrcoef(oe_np)
    est_np = peptools.create_est(pearson_np=pearson_np)

    print(f"normal cov: {est_np[:5]}")
    return

@profile
def mem_efficient_find_best():
    x = sparse.load_npz('/tmp/oe_sparse.npz')
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
    proportion = 0.1
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

    print(f"mem_efficient cov: {est_np[:5]}")
    return

if __name__ == '__main__':
    oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr1_chr1_1000000_obsexp.txt"
    # load_oe_sparse()
    # sum_zero_percent(oe_path)

    store_oe_sparse(oe_path)
    ####
    # test(oe_path)
    # mem_efficient_random1()
    ####
    normal(oe_path)
    mem_efficient_find_best()