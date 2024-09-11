# This is the implementation of the PC1-pattern Estimation using sparse O/E matrix. The algorithm is written in `mem_efficient_sampling` function.
import time
import sys
import numpy as np
from numpy import dot
from numpy.linalg import norm
import pandas as pd
import math
from random import sample
from hicpep import peptools
from scipy import sparse
np.set_printoptions(threshold=sys.maxsize)

'''
`store_oe_sparse` is not included in the benchmark calculation, otherwise the RAM usage of the dense O/E matrix will be recorded. 
We only used this function for storing the O/E matrix as the .npz file. 
'''
def store_oe_sparse(oe_path, name):
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values
    del oe_df
    oe_np = oe_np.astype('float64')
    valid = oe_np.any(axis=1)
    ixgrid = np.ix_(valid, valid) # Record the position of the valid sub-matrix.
    oe_np = oe_np[ixgrid]
    x = sparse.csr_matrix(oe_np)
    sparse.save_npz(f'/tmp/{name}', x)
    return

if __name__ == '__main__':
    oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpep/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr2_chr2_1000000_obsexp.txt"
    store_oe_sparse(oe_path, name="oe_sparse_1Mb.npz") # Not included in benchmark, we comment this line after storing the O/E matrix to .npz file.

    oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpep/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr2_chr2_100000_obsexp.txt"
    store_oe_sparse(oe_path, name="oe_sparse_100Kb.npz") # Not included in benchmark, we comment this line after storing the O/E matrix to .npz file.