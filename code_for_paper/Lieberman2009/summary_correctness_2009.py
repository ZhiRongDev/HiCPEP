import numpy as np
import pandas as pd

def summary_correctness(Lieberman_PC1, approx):
    # Read in the Eigenvector 1
    PC1_df = pd.read_table(Lieberman_PC1, header=None)
    PC1_df = PC1_df.iloc[:, [2]]
    PC1_np = PC1_df.values # Turn into numpy format
    PC1_np = PC1_np.flatten() # Turn into 1D vector
    PC1_np = PC1_np[PC1_np != 0] # Remove 0

    approx_df = pd.read_table(approx, header=None)
    approx_np = approx_df.values # Turn into numpy format
    approx_np = approx_np.flatten() # Turn into 1D vector

    del PC1_df, approx_df

    if len(PC1_np) != len(approx_np): 
        print("juicer_PC1 and approx has a different number of elements")
        return
    
    entryNum = len(PC1_np)

    if np.corrcoef(PC1_np, approx_np)[0][1] < 0:
        approx_np = -approx_np

    PC1_pos_np = PC1_np > 0
    approx_pos_np = approx_np > 0
    PC1_pos_VS_approx_pos_np = PC1_pos_np == approx_pos_np 
    
    correctNum = list(PC1_pos_VS_approx_pos_np).count(True)
    correctRate = correctNum / entryNum

    return {
        "entryNum": entryNum,
        "correctNum": correctNum,
        "correctRate": correctRate,
    }

def write_to_excel(output_file):
    CxMax_df = pd.DataFrame(columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"])
    CxMin_df = pd.DataFrame(columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"])
    Lieberman_PC1_path = "/home/jordan990301/Projects/HiC-PC1_approx/data/Lieberman_Aiden-datasets/eigenvectors"
    approx_path = "/home/jordan990301/Projects/HiC-PC1_approx/outputs/approx_PC1_pattern"
    cell = "GM06690"
    resolution_names = ["1Mb", "100Kb"]
    types = ["CxMax", "CxMin"]

    '''
        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Freadme%5Fv4.txt
        "chromosome" nomenclature: 
        0 Mitochondrial
        23 X
        24 Y

        Notes: There's no chromosome Y(24) in the datasets
    '''
    chroms = [str(i) for i in range(1, 24)]
    
    for resolution_name in resolution_names:
        for chrom in chroms:
            for type in types:
                if resolution_name == "1Mb":
                    resolution = "1000000"
                elif resolution_name == "100Kb":
                    resolution = "100000"

                if chrom == "23":
                    kwargs = {
                        "Lieberman_PC1": f"{Lieberman_PC1_path}/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab",
                        "approx": f"{approx_path}/{cell}/{resolution_name}/{type}/approx_PC1_pattern_chrX.txt",
                    }
                else:
                    kwargs = {
                        "Lieberman_PC1": f"{Lieberman_PC1_path}/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab",
                        "approx": f"{approx_path}/{cell}/{resolution_name}/{type}/approx_PC1_pattern_chr{chrom}.txt",
                    }

                correctness_info = summary_correctness(**kwargs)

                new_row_df = pd.DataFrame(
                    [[cell, resolution_name, f"chr{chrom}", type, correctness_info["entryNum"], correctness_info["correctNum"], correctness_info["correctRate"]]],
                    columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"]
                )

                if type == "CxMax":
                    CxMax_df = pd.concat([CxMax_df, new_row_df], ignore_index=True)
                elif type == "CxMin":
                    CxMin_df = pd.concat([CxMin_df, new_row_df], ignore_index=True)

    output_df = pd.concat([CxMax_df, CxMin_df], ignore_index=True)
    with pd.ExcelWriter(output_file, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="Lieberman2009_summary")

if __name__ == "__main__":
    output_file = "/home/jordan990301/Projects/HiC-PC1_approx/outputs/summary_2009.xlsx"
    write_to_excel(output_file)