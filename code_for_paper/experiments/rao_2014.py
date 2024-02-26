import os
import datetime
import pandas as pd
import numpy as np
from experiments.process import create_approx, calc_correctness, plot_comparison, calc_explained_variance

def data_prepare(docker_volume_path):
    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 data_prepare start")
    data_path = f"{docker_volume_path}/data/rao_2014/juicer_outputs"
    output_path=f"{docker_volume_path}/outputs/approx_pc1_pattern/rao_2014"
    resolutions = [1000000, 100000]
    cell_lines = ["gm12878", "imr90", "hmec", "nhek", "k562", "kbm7", "huvec", "hela", "ch12-lx"]
    methods = ["cxmax", "cxmin"]

    for resolution in resolutions:
        for cell_line in cell_lines:
            for method in methods:
                print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {resolution} {cell_line} {method} start")
                
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
                    create_approx(pearson, output, method, source="2014")

    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 data_prepare end")
    return

def summary_correctness(docker_volume_path):
    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_correctness start")
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
    pc1_path = f"{docker_volume_path}/data/rao_2014/juicer_outputs"
    approx_path = f"{docker_volume_path}/outputs/approx_pc1_pattern/rao_2014"

    for species in ["human", "mouse"]:
        if species == "human":
            cell_lines = ["gm12878", "hela", "hmec", "huvec", "imr90", "k562", "kbm7", "nhek"]
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
                        correctness_info = calc_correctness(pc1, approx, source="2014")

                        if method == "cxmax":
                            cxmax_df.loc[len(cxmax_df)] = [cell_line, resolution, f"chr{chrom}", method, correctness_info["total_entry_num"], correctness_info["valid_entry_num"], correctness_info["correct_num"], correctness_info["correct_rate"]] 
                        elif method == "cxmin":
                            cxmin_df.loc[len(cxmin_df)] = [cell_line, resolution, f"chr{chrom}", method, correctness_info["total_entry_num"], correctness_info["valid_entry_num"], correctness_info["correct_num"], correctness_info["correct_rate"]] 

    output_df = pd.concat([cxmax_df, cxmin_df], ignore_index=True)

    filename = f"{docker_volume_path}/outputs/summary/summary_correctness_2014.xlsx"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_correctness_2014")

    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_correctness end")
    return

def plot_all_comparisons(docker_volume_path):
    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 plot_comparison start")
    data_path = f"{docker_volume_path}/data"

    resolutions = [1000000, 100000]
    cell_lines = ["gm12878", "hela", "hmec", "huvec", "imr90", "k562", "kbm7", "nhek", "ch12-lx"]
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

                print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 plot_comparison {resolution} {cell_line} {method}")
                for chrom in chroms:
                    pc1 = f"{data_path}/rao_2014/juicer_outputs/{cell_line}/{resolution}/eigenvector/pc1_chr{chrom}.txt"
                    approx=f"{docker_volume_path}/outputs/approx_pc1_pattern/rao_2014/{cell_line}/{resolution}/{method}/approx_pc1_pattern_chr{chrom}.txt"
                    output_path = f"{docker_volume_path}/outputs/plots/rao_2014/{cell_line}/{resolution}/{method}"
                    os.makedirs(f"{output_path}/scatter", exist_ok=True)
                    os.makedirs(f"{output_path}/relative_magnitude", exist_ok=True)
                    plot_comparison(pc1, approx, relative_magnitude=f"{output_path}/relative_magnitude/relative_magnitude_chr{chrom}.png", scatter=f"{output_path}/scatter/scatter_chr{chrom}.png", figsize=figsize, source="2014")

    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 plot_comparison end")
    return

def summary_explained_variance(docker_volume_path):
    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_explained_variance start")
    output_df = pd.DataFrame(columns = {
        "cell_line": [],
        "resolution": [],
        "chrom": [],
        "total_entry_num": [],
        "valid_entry_num": [],
        "exp_var_pc1": [],
        "exp_var_pc2": [],
        "exp_var_pc3": [],
        "cos_sim_pc1": [],
        "corr_pc1": []
    })

    resolutions = [1000000, 100000]
    cell_lines = ["gm12878", "hela", "hmec", "huvec", "imr90", "k562", "kbm7", "nhek", "ch12-lx"]

    for resolution in resolutions:
        for cell_line in cell_lines:
            if cell_line == "ch12-lx":
                chroms = [str(i) for i in range(1, 20)]
            else:
                chroms = [str(i) for i in range(1, 23)]
            chroms.extend(["X", "Y"])

            print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_explained_variance_2014 {resolution} {cell_line}")
            for chrom in chroms:
                pearson = f"{docker_volume_path}/data/rao_2014/juicer_outputs/{cell_line}/{resolution}/pearsons/pearson_chr{chrom}.txt"
                pc1 = f"{docker_volume_path}/data/rao_2014/juicer_outputs/{cell_line}/{resolution}/eigenvector/pc1_chr{chrom}.txt"
                pc1_df = pd.read_table(pc1, header=None)
                pc1_df = pc1_df.fillna(0)
                pc1_np = pc1_df.values # Turn into numpy format
                pc1_np = pc1_np.flatten() # Turn into 1D vector
                pc1_np = pc1_np[pc1_np != 0] # Remove 0

                Vh, explained_variances, total_entry_num, valid_entry_num = calc_explained_variance(pearson, source="2014")

                if len(pc1_np) != len(Vh[0]):
                    print("Juicer PC1 and self calculation PC1 has a different valid_entry_num")
                    return

                # Compare the pc1 calculated by numpy with the Juicer's pc1. 
                cos_sim = np.dot(Vh[0], pc1_np) / (np.linalg.norm(Vh[0]) * np.linalg.norm(pc1_np))
                corr = np.corrcoef(Vh[0], pc1_np)[0][1]
                
                output_df.loc[len(output_df)] = [
                    cell_line, 
                    resolution, 
                    chrom, 
                    total_entry_num, 
                    valid_entry_num,
                    explained_variances[0],
                    explained_variances[1],
                    explained_variances[2],
                    cos_sim,
                    corr,
                ] 

    filename = f"{docker_volume_path}/outputs/summary/summary_explained_variance_2014.xlsx"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_explained_variance_2014")

    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} rao_2014 summary_explained_variance end")
    return

def run_all(docker_volume_path):
    data_prepare(docker_volume_path)
    summary_correctness(docker_volume_path)
    plot_all_comparisons(docker_volume_path)
    summary_explained_variance(docker_volume_path)
    return
