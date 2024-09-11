# This Python script is used for calculating the RAM usage fo GM12878, chromosome 2 at 1Mb, 100Kb, 25Kb resolution using the NumPy matrix format.
import numpy as np
import pandas as pd
import gc
import sys

def read_file():
    pearson_path = "/home/jordan990301/Projects/HiCPEP/code_for_paper/notebooks/data/gm12878_pearson_25000_chr2.txt"
    pearson_df = pd.read_table(pearson_path, header=None, sep="\s+").fillna(0)
    pearson_np = pearson_df.values # Turn into numpy.ndarray.
    pearson_np = pearson_np.astype('float64')
    valid = pearson_np.any(axis=1)
    ixgrid = np.ix_(valid, valid) # Record the position of the valid sub-matrix.
    pearson_np = pearson_np[ixgrid]
    print(f"Length: {len(pearson_np)}, size: {sys.getsizeof(pearson_np) / (1024 * 1024)} Mb")
    del pearson_df, valid, ixgrid, pearson_path, pearson_np
    gc.collect()

    pearson_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpep/data_store/data/rao_2014/juicer_outputs/gm12878/100000/pearsons/pearson_chr2.txt"
    pearson_df = pd.read_table(pearson_path, header=None, sep="\s+").fillna(0)
    pearson_np = pearson_df.values # Turn into numpy.ndarray.
    pearson_np = pearson_np.astype('float64')
    valid = pearson_np.any(axis=1)
    ixgrid = np.ix_(valid, valid) # Record the position of the valid sub-matrix.
    pearson_np = pearson_np[ixgrid]
    print(f"Length: {len(pearson_np)}, size: {sys.getsizeof(pearson_np) / (1024 * 1024)} Mb")
    del pearson_df, valid, ixgrid, pearson_path, pearson_np
    gc.collect()

    pearson_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpep/data_store/data/rao_2014/juicer_outputs/gm12878/1000000/pearsons/pearson_chr2.txt"
    pearson_df = pd.read_table(pearson_path, header=None, sep="\s+").fillna(0)
    pearson_np = pearson_df.values # Turn into numpy.ndarray.
    pearson_np = pearson_np.astype('float64')
    valid = pearson_np.any(axis=1)
    ixgrid = np.ix_(valid, valid) # Record the position of the valid sub-matrix.
    pearson_np = pearson_np[ixgrid]
    print(f"Length: {len(pearson_np)}, size: {sys.getsizeof(pearson_np) / (1024 * 1024)} Mb")
    del pearson_df, valid, ixgrid, pearson_path, pearson_np
    gc.collect()

    oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpep/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr2_chr2_100000_obsexp.txt"
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values
    oe_np = oe_np.astype('float64')
    valid = oe_np.any(axis=1)
    ixgrid = np.ix_(valid, valid) # Record the position of the valid sub-matrix.
    oe_np = oe_np[ixgrid].copy()
    print(f"Length: {len(oe_np)}, size: {sys.getsizeof(oe_np) / (1024 * 1024)} Mb")
    del oe_df, valid, ixgrid, oe_path, oe_np
    gc.collect()

    oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpep/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr2_chr2_1000000_obsexp.txt"
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values
    oe_np = oe_np.astype('float64')
    valid = oe_np.any(axis=1)
    ixgrid = np.ix_(valid, valid) # Record the position of the valid sub-matrix.
    oe_np = oe_np[ixgrid]
    print(f"Length: {len(oe_np)}, size: {sys.getsizeof(oe_np) / (1024 * 1024)} Mb")
    del oe_df, valid, ixgrid, oe_path, oe_np
    gc.collect()

if __name__ == "__main__":
    read_file()