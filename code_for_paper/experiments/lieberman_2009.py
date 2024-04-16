"""
The dataset used for this package can be downloaded from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199
"""
import os
import datetime
import numpy as np
import pandas as pd
from hicpap.paptools import create_approx, calc_similarity, plot_comparison
from experiments.utils import read_pearson, flip_track_gc, pca_on_pearson

import logging
logging.basicConfig(format='%(message)s', level=logging.INFO)

def data_prepare(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 data_prepare start")

    data_path = f"{data_store}/data/lieberman_2009"
    output_path=f"{data_store}/outputs/approx_pc1_pattern/lieberman_2009"
    resolutions = [1000000, 100000]
    cell_lines = ["gm06690", "k562"]
    methods = ["cxmax", "cxmin"]

    for resolution in resolutions:
        for cell_line in cell_lines:
            for method in methods:
                logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {resolution} {cell_line} {method} start")
                
                # There are some missing chromosomes in Lieberman's dataset.
                if cell_line == "k562" and resolution == 100000:
                    chroms = [str(i) for i in range(1, 23)]
                else:
                    chroms = [str(i) for i in range(1, 23)]
                    chroms.extend(["X"])
                
                for chrom in chroms:
                    pearson=f"{data_path}/heatmaps/HIC_{cell_line}_chr{chrom}_chr{chrom}_{resolution}_pearson.txt"
                    output=f"{output_path}/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chr{chrom}.txt"
                    pearson_np = read_pearson(pearson=pearson, format="aiden_2009")
                    create_approx(pearson_np=pearson_np, output=output, method=method)

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 data_prepare end")
    return

def summary_similarity(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 summary_similarity start")
    cxmax_df = pd.DataFrame(columns = {
        'cell_line': [], 
        'resolution': [], 
        'chrom':[], 
        "method":[], 
        "total_entry_num":[], 
        "valid_entry_num":[], 
        "similar_num":[], 
        "similar_rate": []
    })
    cxmin_df = pd.DataFrame(columns = {
        'cell_line': [], 
        'resolution': [], 
        'chrom':[], 
        "method":[], 
        "total_entry_num":[], 
        "valid_entry_num":[], 
        "similar_num":[], 
        "similar_rate": []
    })
    pc1_path = f"{data_store}/data/lieberman_2009/eigenvectors"
    approx_path = f"{data_store}/outputs/approx_pc1_pattern/lieberman_2009"

    cell_lines = ["gm06690", "k562"]
    methods = ["cxmax", "cxmin"]

    '''
        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Freadme%5Fv4.txt
        "chromosome" nomenclature: 
        0 Mitochondrial
        23 X
        24 Y

        Notes: 
        In the dataset provided by Lieberman, 
        1. The chromosome 23 are called the chromosome X in the`heatmaps` diractory.
        2. There's no chromosome Y (24) in the datasets. 
        3. There's no Pearson matrix nor PC1 for the K562 in chrX, resolution 100000 .
    '''
    
    for cell_line in cell_lines:
        if cell_line == "k562":
            resolutions = [1000000]
            chroms = [str(i) for i in range(1, 24)]
        else:
            resolutions = [1000000, 100000]
            chroms = [str(i) for i in range(1, 24)]

        for resolution in resolutions:
            for method in methods:
                for chrom in chroms:
                    if cell_line == "gm06690":
                        if chrom == "23":
                            pc1 = f"{pc1_path}/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{approx_path}/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chrX.txt"
                        else:
                            pc1 = f"{pc1_path}/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{approx_path}/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chr{chrom}.txt"
                    elif cell_line == "k562":
                        if chrom == "23":
                            pc1 = f"{pc1_path}/K562-HindIII.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{approx_path}/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chrX.txt"
                        else:
                            pc1 = f"{pc1_path}/K562-HindIII.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{approx_path}/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chr{chrom}.txt"
                    
                    pc1_df = pd.read_table(pc1, header=None, sep="\s+")
                    pc1_df = pc1_df.iloc[:, [2]]
                    pc1_np = pc1_df.values # Turn into numpy format
                    pc1_np = pc1_np.flatten() # Turn into 1D vector
                    tmp = np.full(len(pc1_np), np.nan)
                    tmp[pc1_np != 0] = pc1_np[pc1_np != 0]
                    pc1_np = tmp

                    approx_df = pd.read_table(approx, header=None)
                    approx_np = approx_df.values # Turn into numpy format
                    approx_np = approx_np.flatten() # Turn into 1D vector
                    
                    del pc1_df, approx_df

                    # Flip according to the GC content.
                    if chrom == "23":
                        chrom_name = "X"
                    else:
                        chrom_name = chrom

                    gc_df = pd.read_table(f"./reference_gc/hg18/hg18_gc{resolution}_chr{chrom_name}.txt", skiprows=[0], names=["bin", "GC"])
                    gc_np = gc_df["GC"].values.flatten()
                    pc1_np = flip_track_gc(track_np=pc1_np, gc_np=gc_np)
                    approx_np = flip_track_gc(track_np=approx_np, gc_np=gc_np)

                    similarity_info = calc_similarity(track1_np=pc1_np, track2_np=approx_np)

                    if method == "cxmax":
                        cxmax_df.loc[len(cxmax_df)] = [cell_line, resolution, f"chr{chrom}", method, similarity_info["total_entry_num"], similarity_info["valid_entry_num"], similarity_info["similar_num"], similarity_info["similar_rate"]] 
                    elif method == "cxmin":
                        cxmin_df.loc[len(cxmin_df)] = [cell_line, resolution, f"chr{chrom}", method, similarity_info["total_entry_num"], similarity_info["valid_entry_num"], similarity_info["similar_num"], similarity_info["similar_rate"]] 

    output_df = pd.concat([cxmax_df, cxmin_df], ignore_index=True)

    filename = f"{data_store}/outputs/summary/summary_similarity_2009.xlsx"
    
    if os.path.dirname(filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)

    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_similarity_2009")

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 summary_similarity end")
    return

def plot_all_comparisons(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 plot_comparison start")
    data_path = f"{data_store}/data"
    cell_lines = ["gm06690", "k562"]
    methods = ["cxmax", "cxmin"]

    for cell_line in cell_lines:
        if cell_line == "k562":
            resolutions = [1000000]
            chroms = [str(i) for i in range(1, 23)]
        else:
            resolutions = [1000000, 100000]
            chroms = [str(i) for i in range(1, 24)]

        for resolution in resolutions:
            for method in methods:
                logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 plot_comparison {resolution} {cell_line} {method}")

                for chrom in chroms:
                    if resolution == 1000000:
                        figsize = 20
                    elif resolution == 100000:
                        figsize = 40

                    if cell_line == "gm06690":
                        if chrom == "23":
                            pc1 = f"{data_path}/lieberman_2009/eigenvectors/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{data_store}/outputs/approx_pc1_pattern/lieberman_2009/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chrX.txt"
                        else:
                            pc1 = f"{data_path}/lieberman_2009/eigenvectors/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{data_store}/outputs/approx_pc1_pattern/lieberman_2009/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chr{chrom}.txt"
                    elif cell_line == "k562":
                        if chrom == "23":
                            pc1 = f"{data_path}/lieberman_2009/eigenvectors/K562-HindIII.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{data_store}/outputs/approx_pc1_pattern/lieberman_2009/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chrX.txt"
                        else:
                            pc1 = f"{data_path}/lieberman_2009/eigenvectors/K562-HindIII.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{data_store}/outputs/approx_pc1_pattern/lieberman_2009/{cell_line}/{resolution}/{method}/approx_PC1_pattern_chr{chrom}.txt"

                    output_path = f"{data_store}/outputs/plots/lieberman_2009/{cell_line}/{resolution}/{method}"
                    os.makedirs(f"{output_path}/scatter", exist_ok=True)
                    os.makedirs(f"{output_path}/relative_magnitude", exist_ok=True)

                    pc1_df = pd.read_table(pc1, header=None, sep="\s+")
                    pc1_df = pc1_df.iloc[:, [2]]
                    pc1_np = pc1_df.values # Turn into numpy format
                    pc1_np = pc1_np.flatten() # Turn into 1D vector
                    tmp = np.full(len(pc1_np), np.nan)
                    tmp[pc1_np != 0] = pc1_np[pc1_np != 0]
                    pc1_np = tmp

                    approx_df = pd.read_table(approx, header=None)
                    approx_np = approx_df.values # Turn into numpy format
                    approx_np = approx_np.flatten() # Turn into 1D vector
                    del pc1_df, approx_df
                    relative_magnitude = f"{output_path}/relative_magnitude/relative_magnitude_chr{chrom}.png"
                    scatter=f"{output_path}/scatter/scatter_chr{chrom}.png"

                    # Flip according to the GC content.
                    if chrom == "23":
                        chrom_name = "X"
                    else:
                        chrom_name = chrom

                    gc_df = pd.read_table(f"./reference_gc/hg18/hg18_gc{resolution}_chr{chrom_name}.txt", skiprows=[0], names=["bin", "GC"])
                    gc_np = gc_df["GC"].values.flatten()

                    pc1_np = flip_track_gc(track_np=pc1_np, gc_np=gc_np)
                    approx_np = flip_track_gc(track_np=approx_np, gc_np=gc_np)

                    plot_comparison(pc1_np=pc1_np, approx_np=approx_np, figsize=figsize, scatter=scatter, relative_magnitude=relative_magnitude)

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 plot_comparison end")
    return

def summary_self_pca(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 summary_self_pca start")
    output_df = pd.DataFrame(columns = {
        "cell_line": [],
        "resolution": [],
        "chrom": [],
        "total_entry_num": [],
        "valid_entry_num": [],
        "exp_var_ratio_pc1": [],
        "exp_var_ratio_pc2": [],
        "exp_var_ratio_pc3": [],
        "similar_rate_cxmax": [],
        "similar_rate_cxmin": [],
    })

    cell_lines = ["gm06690", "k562"]

    for cell_line in cell_lines:
        if cell_line == "k562":
            resolutions = [1000000]
        else:
            resolutions = [1000000, 100000]

        for resolution in resolutions:
            chroms = [str(i) for i in range(1, 23)]
            chroms.extend(["X"])

            logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 summary_self_pca_2009 {resolution} {cell_line}")
            for chrom in chroms:
                pearson = f"{data_store}/data/lieberman_2009/heatmaps/HIC_{cell_line}_chr{chrom}_chr{chrom}_{resolution}_pearson.txt"
                pearson_np = read_pearson(pearson=pearson, format="aiden_2009")
                Vh, explained_variances_ratio, total_entry_num, valid_entry_num = pca_on_pearson(pearson_np=pearson_np)

                ## Compute similar_rate for cxmax
                pc1_np = Vh[0].copy()
                approx_np = create_approx(pearson_np=pearson_np, method="cxmax")

                # Flip according to the GC content.
                gc_df = pd.read_table(f"./reference_gc/hg18/hg18_gc{resolution}_chr{chrom}.txt", skiprows=[0], names=["bin", "GC"])
                gc_np = gc_df["GC"].values.flatten()

                pc1_np = flip_track_gc(track_np=pc1_np, gc_np=gc_np)

                similarity_info_cxmax = calc_similarity(track1_np=pc1_np, track2_np=approx_np)
                similar_rate_cxmax = similarity_info_cxmax["similar_rate"]

                ## Compute similar_rate for cxmin
                pc1_np = Vh[0].copy()
                approx_np = create_approx(pearson_np=pearson_np, method="cxmin")

                # Flip according to the GC content.
                gc_df = pd.read_table(f"./reference_gc/hg18/hg18_gc{resolution}_chr{chrom}.txt", skiprows=[0], names=["bin", "GC"])
                gc_np = gc_df["GC"].values.flatten()

                pc1_np = flip_track_gc(track_np=pc1_np, gc_np=gc_np)

                similarity_info_cxmin = calc_similarity(track1_np=pc1_np, track2_np=approx_np)
                similar_rate_cxmin = similarity_info_cxmin["similar_rate"]
                
                output_df.loc[len(output_df)] = [
                    cell_line, 
                    resolution, 
                    chrom, 
                    total_entry_num, 
                    valid_entry_num,
                    explained_variances_ratio[0],
                    explained_variances_ratio[1],
                    explained_variances_ratio[2],
                    similar_rate_cxmax,
                    similar_rate_cxmin
                ] 

    filename = f"{data_store}/outputs/summary/summary_self_pca_2009.xlsx"

    if os.path.dirname(filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)

    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_self_pca_2009")

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 summary_self_pca end")
    return

def run_all(data_store):
    data_prepare(data_store) # Create the Approximated PC1-pattern .txt files.
    summary_similarity(data_store) # Compare the similarity difference with the PC1 and the Approximated PC1-pattern.
    plot_all_comparisons(data_store) # Plot the scatter and relative-magnitude chart.
    summary_self_pca(data_store) # Performing the PCA by self and get the information of the explained variances.
    return