import os
import datetime
import pandas as pd
import numpy as np
from hicpap.paptools import read_pearson, create_approx, calc_correctness, plot_comparison, pca_on_pearson, flip_tracks

import logging
logging.basicConfig(format='%(message)s', level=logging.INFO)

def data_prepare(data_store):
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
                    pearson_np = read_pearson(pearson=pearson, zero_mean=True, format="rao_2014")
                    create_approx(pearson_np=pearson_np, output=output, method=method)

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 data_prepare end")
    return

def summary_correctness(data_store):
    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_correctness start")
    cxmax_df = pd.DataFrame(columns = {
        'cell_line': [], 
        'resolution': [], 
        'chrom':[], 
        "method":[], 
        "total_entry_num":[], 
        "valid_entry_num":[], 
        "correct_num":[], 
        "correct_rate": []
    })
    cxmin_df = pd.DataFrame(columns = {
        'cell_line': [], 
        'resolution': [], 
        'chrom':[], 
        "method":[], 
        "total_entry_num":[], 
        "valid_entry_num":[], 
        "correct_num":[], 
        "correct_rate": []
    })
    pc1_path = f"{data_store}/data/rao_2014/juicer_outputs"
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
                        pc1 = f"{pc1_path}/{cell_line}/{resolution}/eigenvector/pc1_chr{chrom}.txt"
                        approx = f"{approx_path}/{cell_line}/{resolution}/{method}/approx_pc1_pattern_chr{chrom}.txt"

                        pc1_df = pd.read_table(pc1, header=None)
                        pc1_df = pc1_df.fillna(0)
                        pc1_np = pc1_df.values # Turn into numpy format
                        pc1_np = pc1_np.flatten() # Turn into 1D vector
                        approx_df = pd.read_table(approx, header=None)
                        approx_np = approx_df.values # Turn into numpy format
                        approx_np = approx_np.flatten() # Turn into 1D vector

                        del pc1_df, approx_df
                        pc1_np, approx_np = flip_tracks(track1_np=pc1_np, track2_np=approx_np)
                        correctness_info = calc_correctness(pc1_np=pc1_np, approx_np=approx_np)

                        if method == "cxmax":
                            cxmax_df.loc[len(cxmax_df)] = [cell_line, resolution, f"chr{chrom}", method, correctness_info["total_entry_num"], correctness_info["valid_entry_num"], correctness_info["correct_num"], correctness_info["correct_rate"]] 
                        elif method == "cxmin":
                            cxmin_df.loc[len(cxmin_df)] = [cell_line, resolution, f"chr{chrom}", method, correctness_info["total_entry_num"], correctness_info["valid_entry_num"], correctness_info["correct_num"], correctness_info["correct_rate"]] 

    output_df = pd.concat([cxmax_df, cxmin_df], ignore_index=True)

    filename = f"{data_store}/outputs/summary/summary_correctness_2014.xlsx"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_correctness_2014")

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_correctness end")
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
                    pc1_df = pc1_df.fillna(0)
                    pc1_np = pc1_df.values # Turn into numpy format
                    pc1_np = pc1_np.flatten() # Turn into 1D vector
                    approx_df = pd.read_table(approx, header=None)
                    approx_np = approx_df.values # Turn into numpy format
                    approx_np = approx_np.flatten() # Turn into 1D vector
                    relative_magnitude = f"{output_path}/relative_magnitude/relative_magnitude_chr{chrom}.png"
                    scatter = f"{output_path}/scatter/scatter_chr{chrom}.png"

                    del pc1_df, approx_df
                    pc1_np, approx_np = flip_tracks(track1_np=pc1_np, track2_np=approx_np)
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
        "exp_var_pc1": [],
        "exp_var_pc2": [],
        "exp_var_pc3": [],
        "correct_rate_cxmax": [],
        "correct_rate_cxmin": [],
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
                pearson_np = read_pearson(pearson=pearson, zero_mean=True, format="rao_2014")
                Vh, explained_variances, total_entry_num, valid_entry_num = pca_on_pearson(pearson_np=pearson_np)

                ## Compute correct_rate for cxmax
                pc1_np = Vh[0].copy()
                approx_np = create_approx(pearson_np=pearson_np, method="cxmax")
                pc1_np, approx_np = flip_tracks(track1_np=pc1_np, track2_np=approx_np)
                correctness_info = calc_correctness(pc1_np=pc1_np, approx_np=approx_np)
                correct_rate_cxmax = correctness_info["correct_rate"]

                ## Compute correct_rate for cxmax
                pc1_np = Vh[0].copy()
                approx_np = create_approx(pearson_np=pearson_np, method="cxmin")
                pc1_np, approx_np = flip_tracks(track1_np=pc1_np, track2_np=approx_np)
                correctness_info = calc_correctness(pc1_np=pc1_np, approx_np=approx_np)
                correct_rate_cxmin = correctness_info["correct_rate"]
                
                output_df.loc[len(output_df)] = [
                    cell_line, 
                    resolution, 
                    chrom, 
                    total_entry_num, 
                    valid_entry_num,
                    explained_variances[0],
                    explained_variances[1],
                    explained_variances[2],
                    correct_rate_cxmax,
                    correct_rate_cxmin
                ] 

    filename = f"{data_store}/outputs/summary/summary_self_pca_2014.xlsx"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_self_pca_2014")

    logging.info(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_self_pca end")
    return

def run_all(data_store):
    data_prepare(data_store)
    summary_correctness(data_store)
    plot_all_comparisons(data_store)
    summary_self_pca(data_store)
    return
