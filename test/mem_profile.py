import sys
import numpy as np
import pandas as pd
import hicstraw
from hicpep import peptools
from sklearn.decomposition import PCA
from scipy import sparse
np.set_printoptions(threshold=sys.maxsize)

oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr1_chr1_100000_obsexp.txt"

# @profile
def test():
    A = np.array([
        [1, 0, 2],
        [0, 0, 3],
        [4, 5, 6]
    ])
    # row = np.array([0, 0, 1, 2, 2, 2])
    # col = np.array([0, 2, 2, 0, 1, 2])
    # data = np.array([1, 2, 3, 4, 5, 6])
    # sA = sparse.csr_matrix((data, (row, col)), shape=(3, 3))

    nonzeros = A.nonzero()
    sA = sparse.csr_matrix((A[nonzeros], nonzeros), shape=A.shape)
    print(sA, '\n')
    print(sA.toarray())

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

# @profile
def normal(oe_path):
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values 
    oe_np = oe_np.astype('float64')
    del oe_df
    pearson_np = np.corrcoef(oe_np)
    peptools.create_est(pearson_np=pearson_np)
    return

# @profile
def mem_efficient(oe_path):
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values
    del oe_df
    oe_np = oe_np.astype('float64')
    diag = np.diag(oe_np)
    diag_valid = diag != 0 
    ixgrid = np.ix_(diag_valid, diag_valid) # Record the position of the valid sub-matrix.
    oe_np = oe_np[ixgrid] 
    n = len(oe_np[0])
    X = sparse.csr_matrix(oe_np)
    cov = np.cov(oe_np, bias=True)
    del oe_np
    index_s = 0
    I = np.full(n, 1)
    std = []
    C = []

    for i in range(n):
        tmp = np.std(X[i].toarray())
        std.append(tmp)

    for i in range(n):
        tmp = X[i].mean()
        C.append(tmp)

    v = (X[index_s].toarray() - X[index_s].mean())[0] # X selected.
    std = np.array(std)
    C = np.array(C)
    S = sparse.csr_matrix((1/std, np.arange(len(1/std)), np.arange(len(1/std) + 1)), shape=(len(1/std), len(1/std)))
    X = X.dot(S)
    r = C / std
    corr_s = X.dot(v).flatten()
    alpha_s = r @ v 
    corr_s = (corr_s - I * alpha_s) / n
    print(corr_s[:5])

    return

if __name__ == '__main__':
    # test()
    # sum_zero_percent(oe_path)
    # normal(oe_path)
    mem_efficient(oe_path)