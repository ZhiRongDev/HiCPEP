import os
import datetime
import pandas as pd
from experiments.process import create_approx, calc_correctness

def data_prepare(docker_volume_path):
    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 create_approx start")

    data_path = f"{docker_volume_path}/data/lieberman_2009"
    output_path=f"{docker_volume_path}/outputs/approx_pc1_pattern/lieberman_2009"
    resolutions = [1000000, 100000]
    cell_lines = ["gm06690", "k562"]
    types = ["cxmax", "cxmin"]

    for resolution in resolutions:
        for cell_line in cell_lines:
            for type in types:
                print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {resolution} {cell_line} {type} start")
                
                # There are some missing chromosomes in Lieberman's dataset.
                if cell_line == "k562" and resolution == 100000:
                    chrom_list = [str(i) for i in range(1, 23)]
                else:
                    chrom_list = [str(i) for i in range(1, 23)]
                    chrom_list.extend(["X"])
                
                for chrom in chrom_list:
                    pearson=f"{data_path}/heatmaps/HIC_{cell_line}_chr{chrom}_chr{chrom}_{resolution}_pearson.txt"
                    output=f"{output_path}/{cell_line}/{resolution}/{type}/approx_PC1_pattern_chr{chrom}.txt"
                    create_approx(pearson, output, type, source="2009")
                
                print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} {resolution} {cell_line} {type} end")

    print(f"{datetime.datetime.now():%Y-%m-%d %H:%M:%S} lieberman_2009 create_approx end")
    return

def summary_correctness(docker_volume_path):
    output = f"{docker_volume_path}/outputs/summary/summary_2009.xlsx"
    cxmax_df = pd.DataFrame()
    cxmin_df = pd.DataFrame()
    pc1_path = f"{docker_volume_path}/data/lieberman_2009/eigenvectors"
    approx_path = f"{docker_volume_path}/outputs/approx_pc1_pattern/lieberman_2009"

    cell_lines = ["gm06690", "k562"]
    types = ["cxmax", "cxmin"]

    '''
        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Freadme%5Fv4.txt
        "chromosome" nomenclature: 
        0 Mitochondrial
        23 X
        24 Y

        Notes: 
        In Lieberman's dataset, 
        1. The symbol used for the 23 chromosome are not the same in the `eigenvectors` and `heatmaps` diractory. (eigenvectors: 23, heatmaps: X)
        2. There's no chromosome Y(24) in the datasets, 
        3. There's 100000 resolution pearson matrix for k562 (chr1 - chr22, "no" chrX), 
        4. There's "no" 100000 resolution PC1 for k562.
    '''
    
    for cell_line in cell_lines:
        if cell_line == "k562":
            resolutions = ["1000000"]
            chroms = [str(i) for i in range(1, 23)]
        else:
            resolutions = ["1000000", "100000"]
            chroms = [str(i) for i in range(1, 24)]

        for resolution in resolutions:
            for chrom in chroms:
                for type in types:
                    if cell_line == "gm06690":
                        if chrom == "23":
                            pc1 = f"{pc1_path}/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{approx_path}/{cell_line}/{resolution}/{type}/approx_PC1_pattern_chrX.txt"
                        else:
                            pc1 = f"{pc1_path}/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                            approx = f"{approx_path}/{cell_line}/{resolution}/{type}/approx_PC1_pattern_chr{chrom}.txt"
                    elif cell_line == "k562":
                        pc1 = f"{pc1_path}/K562-HindIII.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab"
                        approx = f"{approx_path}/{cell_line}/{resolution}/{type}/approx_PC1_pattern_chr{chrom}.txt"

                    pc1_df = pd.read_table(pc1, header=None)
                    pc1_df = pc1_df.iloc[:, [2]]
                    approx_df = pd.read_table(approx, header=None)

                    correctness_info = calc_correctness(pc1_df, approx_df)

                    if type == "cxmax":
                        if cxmax_df.empty:
                            cxmax_df = pd.DataFrame(
                                [[cell_line, resolution, f"chr{chrom}", type, correctness_info["valid_entry_num"], correctness_info["correct_num"], correctness_info["correct_rate"]]],
                                columns=['cell', 'resolution', 'chromosome', "type", "valid_entry_num", "correct_num", "correct_rate"]
                            )
                        else:
                            new_row_df = pd.DataFrame(
                                [[cell_line, resolution, f"chr{chrom}", type, correctness_info["valid_entry_num"], correctness_info["correct_num"], correctness_info["correct_rate"]]],
                                columns=['cell', 'resolution', 'chromosome', "type", "valid_entry_num", "correct_num", "correct_rate"]
                            )
                            cxmax_df = pd.concat([cxmax_df, new_row_df], ignore_index=True)
                    elif type == "cxmin":
                        if cxmin_df.empty:
                            cxmin_df = pd.DataFrame(
                                [[cell_line, resolution, f"chr{chrom}", type, correctness_info["valid_entry_num"], correctness_info["correct_num"], correctness_info["correct_rate"]]],
                                columns=['cell', 'resolution', 'chromosome', "type", "valid_entry_num", "correct_num", "correct_rate"]
                            )
                        else:
                            new_row_df = pd.DataFrame(
                                [[cell_line, resolution, f"chr{chrom}", type, correctness_info["valid_entry_num"], correctness_info["correct_num"], correctness_info["correct_rate"]]],
                                columns=['cell', 'resolution', 'chromosome', "type", "valid_entry_num", "correct_num", "correct_rate"]
                            )
                            cxmin_df = pd.concat([cxmin_df, new_row_df], ignore_index=True)

    output_df = pd.concat([cxmax_df, cxmin_df], ignore_index=True)

    filename = output
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_2009")
    return

def run_all(docker_volume_path):
    # data_prepare(docker_volume_path)
    summary_correctness(docker_volume_path)
