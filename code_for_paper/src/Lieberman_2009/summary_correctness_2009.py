import numpy as np
import pandas as pd
import argparse
import os

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
        print("Lieberman_PC1 and approx has a different number of elements")
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

def write_to_excel(docker_volume_path, output):
    CxMax_df = pd.DataFrame()
    CxMin_df = pd.DataFrame()
    Lieberman_PC1_path = f"{docker_volume_path}/data/Lieberman_2009/eigenvectors"
    approx_path = f"{docker_volume_path}/outputs/approx_PC1_pattern/Lieberman_2009"
    
    cells = ["gm06690", "k562"]
    types = ["CxMax", "CxMin"]

    '''
        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Freadme%5Fv4.txt
        "chromosome" nomenclature: 
        0 Mitochondrial
        23 X
        24 Y

        Notes: There's no chromosome Y(24) in the datasets, and there's some missing data in k562
    '''
    
    for cell in cells:
        if cell == "k562":
            resolutions = ["1000000"]
            chroms = [str(i) for i in range(1, 23)]
        else:
            resolutions = ["1000000", "100000"]
            chroms = [str(i) for i in range(1, 24)]

        for resolution in resolutions:
            for chrom in chroms:
                for type in types:
                    if cell == "gm06690":
                        if chrom == "23":
                            kwargs = {
                                "Lieberman_PC1": f"{Lieberman_PC1_path}/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab",
                                "approx": f"{approx_path}/{cell}/{resolution}/{type}/approx_PC1_pattern_chrX.txt",
                            }
                        else:
                            kwargs = {
                                "Lieberman_PC1": f"{Lieberman_PC1_path}/GM-combined.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab",
                                "approx": f"{approx_path}/{cell}/{resolution}/{type}/approx_PC1_pattern_chr{chrom}.txt",
                            }
                    elif cell == "k562":
                        kwargs = {
                            "Lieberman_PC1": f"{Lieberman_PC1_path}/K562-HindIII.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab",
                            "approx": f"{approx_path}/{cell}/{resolution}/{type}/approx_PC1_pattern_chr{chrom}.txt",
                        }
                        # if chrom == "23":
                        #     kwargs = {
                        #         "Lieberman_PC1": f"{Lieberman_PC1_path}/K562-HindIII.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab",
                        #         "approx": f"{approx_path}/{cell}/{resolution}/{type}/approx_PC1_pattern_chrX.txt",
                        #     }
                        # else:
                        #     kwargs = {
                        #         "Lieberman_PC1": f"{Lieberman_PC1_path}/K562-HindIII.ctg{chrom}.ctg{chrom}.{resolution}bp.hm.eigenvector.tab",
                        #         "approx": f"{approx_path}/{cell}/{resolution}/{type}/approx_PC1_pattern_chr{chrom}.txt",
                        #     }

                    correctness_info = summary_correctness(**kwargs)

                    if type == "CxMax":
                        if CxMax_df.empty:
                            CxMax_df = pd.DataFrame(
                                [[cell, resolution, f"chr{chrom}", type, correctness_info["entryNum"], correctness_info["correctNum"], correctness_info["correctRate"]]],
                                columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"]
                            )
                        else:
                            new_row_df = pd.DataFrame(
                                [[cell, resolution, f"chr{chrom}", type, correctness_info["entryNum"], correctness_info["correctNum"], correctness_info["correctRate"]]],
                                columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"]
                            )
                            CxMax_df = pd.concat([CxMax_df, new_row_df], ignore_index=True)
                    elif type == "CxMin":
                        if CxMin_df.empty:
                            CxMin_df = pd.DataFrame(
                                [[cell, resolution, f"chr{chrom}", type, correctness_info["entryNum"], correctness_info["correctNum"], correctness_info["correctRate"]]],
                                columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"]
                            )
                        else:
                            new_row_df = pd.DataFrame(
                                [[cell, resolution, f"chr{chrom}", type, correctness_info["entryNum"], correctness_info["correctNum"], correctness_info["correctRate"]]],
                                columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"]
                            )
                            CxMin_df = pd.concat([CxMin_df, new_row_df], ignore_index=True)

    output_df = pd.concat([CxMax_df, CxMin_df], ignore_index=True)

    filename = output
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    with pd.ExcelWriter(filename, mode="w") as writer:
        output_df.to_excel(writer, sheet_name="summary_2009")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='summary_correctness_2009.py',
        allow_abbrev=False,
        description='What the program does', epilog='Text at the bottom of help'
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        help="Input blablabla"
    )
    parser.add_argument(
        "--docker_volume_path",
        type=str,
        required=True,
        help="Input blablabla"
    )

    args = parser.parse_args()
    output = args.output
    docker_volume_path = args.docker_volume_path

    write_to_excel(docker_volume_path=docker_volume_path, output=output)