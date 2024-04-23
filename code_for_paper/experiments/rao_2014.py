"""
The dataset used for this package can be downloaded from:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525
"""
import os
import datetime
import pandas as pd
import numpy as np
from hicpap.paptools import create_approx, calc_similarity, plot_comparison
from experiments.utils import read_pearson, flip_tracks, pca_on_pearson

import logging
logging.basicConfig(format='%(message)s', level=logging.INFO)

def data_prepare(data_store):
    """
    Note that GSE63525 doesn't provide the `.hic` files for HeLa, we skip this cell line.
    """
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 data_prepare start")
    data_path = f"{data_store}/data/rao_2014/juicer_outputs"
    output_path=f"{data_store}/outputs/approx_pc1_pattern/rao_2014"
    resolutions = [1000000, 100000]
    cell_lines = ["gm12878", "imr90", "hmec", "nhek", "k562", "kbm7", "huvec", "ch12-lx"]
    methods = ["cxmax", "cxmin"]

    for resolution in resolutions:
        for cell_line in cell_lines:
            for method in methods:
                logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {resolution} {cell_line} {method} start")
                
                # There are only 19 chromosome in CH12-LX, 
                if cell_line == "ch12-lx":
                    chroms = [str(i) for i in range(1, 20)]
                    chroms.extend(["X", "Y"])
                else:
                    chroms = [str(i) for i in range(1, 23)]
                    chroms.extend(["X", "Y"])
                
                for chrom in chroms:
                    pearson=f"{data_path}/{cell_line}/{resolution}/pearsons/pearson_chr{chrom}.txt"
                    output=f"{output_path}/{cell_line}/{resolution}/{method}/approx_pc1_pattern_chr{chrom}.txt"
                    pearson_np = read_pearson(pearson=pearson, format="juicer")
                    create_approx(pearson_np=pearson_np, output=output, method=method)

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 data_prepare end")
    return

def summary_similarity(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_similarity start")
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
    juicer_outputs_path = f"{data_store}/data/rao_2014/juicer_outputs"
    approx_path = f"{data_store}/outputs/approx_pc1_pattern/rao_2014"

    for species in ["human", "mouse"]:
        if species == "human":
            cell_lines = ["gm12878", "imr90", "hmec", "nhek", "k562", "kbm7", "huvec"]
            chroms = [str(i) for i in range(1, 23)]
        elif species == "mouse":
            cell_lines = ["ch12-lx"]
            chroms = [str(i) for i in range(1, 20)]
    
        chroms.extend(["X", "Y"])
        resolutions = [1000000, 100000]
        methods = ["cxmax", "cxmin"]

        for resolution in resolutions:
            for cell_line in cell_lines:
                for chrom in chroms:
                    for method in methods:
                        pc1 = f"{juicer_outputs_path}/{cell_line}/{resolution}/eigenvector/pc1_chr{chrom}.txt"
                        approx = f"{approx_path}/{cell_line}/{resolution}/{method}/approx_pc1_pattern_chr{chrom}.txt"

                        pc1_df = pd.read_table(pc1, header=None)
                        pc1_np = pc1_df.values # Turn into numpy format
                        pc1_np = pc1_np.flatten() # Turn into 1D vector
                        approx_df = pd.read_table(approx, header=None)
                        approx_np = approx_df.values # Turn into numpy format
                        approx_np = approx_np.flatten() # Turn into 1D vector

                        del pc1_df, approx_df

                        ### Flip tracks according to the similarity between track1 and track2
                        pc1_np, approx_np  = flip_tracks(track1_np=pc1_np, track2_np=approx_np)

                        similarity_info = calc_similarity(track1_np=pc1_np, track2_np=approx_np)

                        if method == "cxmax":
                            cxmax_df.loc[len(cxmax_df)] = [cell_line, resolution, f"chr{chrom}", method, similarity_info["total_entry_num"], similarity_info["valid_entry_num"], similarity_info["similar_num"], similarity_info["similar_rate"]] 
                        elif method == "cxmin":
                            cxmin_df.loc[len(cxmin_df)] = [cell_line, resolution, f"chr{chrom}", method, similarity_info["total_entry_num"], similarity_info["valid_entry_num"], similarity_info["similar_num"], similarity_info["similar_rate"]] 

    output_df = pd.concat([cxmax_df, cxmin_df], ignore_index=True)

    filename = f"{data_store}/outputs/summary/summary_similarity_2014.xlsx"

    if os.path.dirname(filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)

    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_similarity_2014")

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_similarity end")
    return

def summary_similar_rate_percentage(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_similar_rate_percentage start")
    output_df = pd.DataFrame(columns = {
        'cell_line': [], 
        'resolution': [], 
        'chrom':[], 
        "method":[], 
        "total_entry_num":[], 
        "valid_entry_num":[], 
        "similar_rate_over90": [],
        "similar_rate_over95": [],
        "similar_rate_over99": [],
    })
    juicer_outputs_path = f"{data_store}/data/rao_2014/juicer_outputs"

    for species in ["human", "mouse"]:
        if species == "human":
            cell_lines = ["gm12878", "imr90", "hmec", "nhek", "k562", "kbm7", "huvec"]
            chroms = [str(i) for i in range(1, 23)]
        elif species == "mouse":
            cell_lines = ["ch12-lx"]
            chroms = [str(i) for i in range(1, 20)]
    
        chroms.extend(["X", "Y"])
        resolutions = [1000000, 100000]

        for resolution in resolutions:
            for cell_line in cell_lines:
                for chrom in chroms:

                    pc1 = f"{juicer_outputs_path}/{cell_line}/{resolution}/eigenvector/pc1_chr{chrom}.txt"
                    pc1_df = pd.read_table(pc1, header=None)
                    pc1_np = pc1_df.values # Turn into numpy format
                    pc1_np = pc1_np.flatten() # Turn into 1D vector

                    del pc1_df

                    pearson = f"{juicer_outputs_path}/{cell_line}/{resolution}/pearsons/pearson_chr{chrom}.txt"
                    pearson_np = read_pearson(pearson=pearson, format="juicer")
                    # Extract the valid sub-matrix and record the position of NaN entries (The calculation will only be performed on the not-NaN entries).
                    diag = np.diag(pearson_np)
                    diag_valid = ~np.isnan(diag)
                    ixgrid = np.ix_(diag_valid, diag_valid) # Record the position of the valid sub-matrix.

                    pearson_np = pearson_np[ixgrid] 
                    pearson_np -= pearson_np.mean(axis=1, keepdims=True)
                    n = len(pearson_np[0])
                    cov_np = pearson_np @ pearson_np.T / n

                    del pearson_np
                    # Start to accumulate the similar_rate percentages. 
                    similar_rates_over90_count = 0
                    similar_rates_over95_count = 0
                    similar_rates_over99_count = 0

                    for i in range(len(cov_np)):
                        # Place back the valid entries to it's origin position in the chromosome.  
                        approx_np = np.full(len(diag_valid), np.nan)
                        approx_np[diag_valid] = cov_np[i]

                        ### Flip tracks according to the similarity between track1 and track2
                        pc1_np, approx_np  = flip_tracks(track1_np=pc1_np, track2_np=approx_np)

                        similar_info = calc_similarity(track1_np=pc1_np, track2_np=approx_np)

                        if similar_info["similar_rate"] >= 0.9:
                            similar_rates_over90_count += 1
                        
                        if similar_info["similar_rate"] >= 0.95:
                            similar_rates_over95_count += 1

                        if similar_info["similar_rate"] >= 0.99:
                            similar_rates_over99_count += 1
                    
                    similar_rates_over90 = float(similar_rates_over90_count / len(cov_np[0]))
                    similar_rates_over95 = float(similar_rates_over95_count / len(cov_np[0]))
                    similar_rates_over99 = float(similar_rates_over99_count / len(cov_np[0]))
                    output_df.loc[len(output_df)] = [cell_line, resolution, f"chr{chrom}", "cxmax", similar_info["total_entry_num"], similar_info["valid_entry_num"], similar_rates_over90, similar_rates_over95, similar_rates_over99] 

    filename = f"{data_store}/outputs/summary/summary_similar_rate_percentage_2014.xlsx"

    if os.path.dirname(filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)

    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="similar_rate_percentage")

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_similar_rate_percentage end")
    return

def plot_all_comparisons(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 plot_comparison start")
    data_path = f"{data_store}/data"

    resolutions = [1000000, 100000]
    cell_lines = ["gm12878", "imr90", "hmec", "nhek", "k562", "kbm7", "huvec", "ch12-lx"]
    methods = ["cxmax", "cxmin"]

    for resolution in resolutions:
        for cell_line in cell_lines:
            for method in methods:
                if resolution == 1000000:
                    figsize = 20
                elif resolution == 100000:
                    figsize = 40

                if cell_line == "ch12-lx":
                    chroms = [str(i) for i in range(1, 20)]
                else:
                    chroms = [str(i) for i in range(1, 23)]
                chroms.extend(["X", "Y"])

                logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 plot_comparison {resolution} {cell_line} {method}")
                for chrom in chroms:
                    pc1 = f"{data_path}/rao_2014/juicer_outputs/{cell_line}/{resolution}/eigenvector/pc1_chr{chrom}.txt"
                    approx=f"{data_store}/outputs/approx_pc1_pattern/rao_2014/{cell_line}/{resolution}/{method}/approx_pc1_pattern_chr{chrom}.txt"
                    output_path = f"{data_store}/outputs/plots/rao_2014/{cell_line}/{resolution}/{method}"
                    os.makedirs(f"{output_path}/scatter", exist_ok=True)
                    os.makedirs(f"{output_path}/relative_magnitude", exist_ok=True)

                    pc1_df = pd.read_table(pc1, header=None)
                    pc1_np = pc1_df.values # Turn into numpy format
                    pc1_np = pc1_np.flatten() # Turn into 1D vector
                    approx_df = pd.read_table(approx, header=None)
                    approx_np = approx_df.values # Turn into numpy format
                    approx_np = approx_np.flatten() # Turn into 1D vector
                    relative_magnitude = f"{output_path}/relative_magnitude/relative_magnitude_chr{chrom}.png"
                    scatter = f"{output_path}/scatter/scatter_chr{chrom}.png"

                    del pc1_df, approx_df

                    ### Flip tracks according to the similarity between track1 and track2
                    pc1_np, approx_np  = flip_tracks(track1_np=pc1_np, track2_np=approx_np)

                    plot_comparison(pc1_np=pc1_np, approx_np=approx_np, figsize=figsize, scatter=scatter, relative_magnitude=relative_magnitude)

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 plot_comparison end")
    return

def summary_self_pca(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_self_pca start")
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

    resolutions = [1000000, 100000]
    cell_lines = ["gm12878", "imr90", "hmec", "nhek", "k562", "kbm7", "huvec", "ch12-lx"]

    for resolution in resolutions:
        for cell_line in cell_lines:
            if cell_line == "ch12-lx":
                chroms = [str(i) for i in range(1, 20)]
            else:
                chroms = [str(i) for i in range(1, 23)]
            chroms.extend(["X", "Y"])

            logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_self_pca_2014 {resolution} {cell_line}")
            for chrom in chroms:
                pearson = f"{data_store}/data/rao_2014/juicer_outputs/{cell_line}/{resolution}/pearsons/pearson_chr{chrom}.txt"
                pearson_np = read_pearson(pearson=pearson, format="juicer")
                Vh, explained_variances_ratio, total_entry_num, valid_entry_num = pca_on_pearson(pearson_np=pearson_np)

                ## Compute similar_rate for cxmax
                pc1_np = Vh[0].copy()
                approx=f"{data_store}/outputs/approx_pc1_pattern/rao_2014/{cell_line}/{resolution}/cxmax/approx_pc1_pattern_chr{chrom}.txt"
                approx_df = pd.read_table(approx, header=None)
                approx_np = approx_df.values # Turn into numpy format
                approx_np = approx_np.flatten() # Turn into 1D vector

                ### Flip tracks according to the similarity between track1 and track2
                pc1_np, approx_np  = flip_tracks(track1_np=pc1_np, track2_np=approx_np)

                similarity_info = calc_similarity(track1_np=pc1_np, track2_np=approx_np)
                similar_rate_cxmax = similarity_info["similar_rate"]

                ## Compute similar_rate for cxmin
                pc1_np = Vh[0].copy()
                approx=f"{data_store}/outputs/approx_pc1_pattern/rao_2014/{cell_line}/{resolution}/cxmin/approx_pc1_pattern_chr{chrom}.txt"
                approx_df = pd.read_table(approx, header=None)
                approx_np = approx_df.values # Turn into numpy format
                approx_np = approx_np.flatten() # Turn into 1D vector

                ### Flip tracks according to the similarity between track1 and track2
                pc1_np, approx_np  = flip_tracks(track1_np=pc1_np, track2_np=approx_np)

                similarity_info = calc_similarity(track1_np=pc1_np, track2_np=approx_np)
                similar_rate_cxmin = similarity_info["similar_rate"]
                
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

    filename = f"{data_store}/outputs/summary/summary_self_pca_2014.xlsx"

    if os.path.dirname(filename):
        os.makedirs(os.path.dirname(filename), exist_ok=True)

    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_self_pca_2014")

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_self_pca end")
    return

def run_all(data_store):
    # data_prepare(data_store) # Create the Approximated PC1-pattern .txt files.
    summary_similarity(data_store) # Compare the similarity difference with the PC1 and the Approximated PC1-pattern.
    plot_all_comparisons(data_store) # Plot the scatter and relative-magnitude chart.
    summary_self_pca(data_store) # Performing the PCA by self and get the information of the explained variance ratios.
    summary_similar_rate_percentage(data_store) # Summarize the percentage of columns in the covariance matrix that has a similar_rate over 90%, 95% or 99%. 
    return
