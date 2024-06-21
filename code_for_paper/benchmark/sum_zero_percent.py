# This Python script is used for calculating the sparsity of a matrix (i.e., The percentage of zero accounts for the matrix).
import numpy as np
import pandas as pd
import hicstraw

if __name__ == '__main__':
    print("GM06690")
    oe_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/lieberman_2009/heatmaps/HIC_gm06690_chr2_chr2_100000_obsexp.txt"
    oe_df = pd.read_table(oe_path, index_col=0, header=1, sep="\s+")
    oe_np = oe_df.values 
    size = len(oe_np) * len(oe_np)
    zero = size - np.count_nonzero(oe_np)
    print(zero)
    print(size)
    print(zero / size)

    print("\nGM12878")
    hic_path = "/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpap/data_store/data/rao_2014/hic/GSE63525_GM12878_insitu_primary_replicate_combined_30.hic"
    hic = hicstraw.HiCFile(hic_path)
    chrom = "2"
    normalization = "KR"
    resolution = 100000 

    for chromosome in hic.getChromosomes():
        if chromosome.name == chrom:
            chrom_size = int(chromosome.length)

    oe_mat = hic.getMatrixZoomData(chrom, chrom, "oe", normalization, "BP", resolution)
    oe_mat_np = oe_mat.getRecordsAsMatrix(0, chrom_size, 0, chrom_size)
    size = len(oe_mat_np) * len(oe_mat_np)
    zero = size - np.count_nonzero(oe_mat_np)
    print(zero)
    print(size)
    print(zero / size)