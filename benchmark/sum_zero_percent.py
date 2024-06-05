import numpy as np
import pandas as pd

if __name__ == '__main__':
    oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr2_chr2_100000_obsexp.txt"
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values 
    size = len(oe_np) * len(oe_np)
    zero = size - np.count_nonzero(oe_np)
    print(zero)
    print(size)
    print(zero / size)