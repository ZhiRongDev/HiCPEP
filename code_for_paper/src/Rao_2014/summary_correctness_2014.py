import numpy as np
import pandas as pd
import argparse
import os

def summary_correctness(Juicer_PC1, approx):
    # Read in the Eigenvector 1
    PC1_df = pd.read_table(Juicer_PC1, header=None)
    PC1_df = PC1_df.dropna(axis=0, how="all").reset_index(drop=True)
    PC1_np = PC1_df.values # Turn into numpy format
    PC1_np = PC1_np.flatten() # Turn into 1D vector

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

def write_to_excel(docker_volume_path, output):
    CxMax_df = pd.DataFrame()
    CxMin_df = pd.DataFrame()
    Juicer_PC1_path = f"{docker_volume_path}/data/Rao_2014/juicer_outputs"
    approx_path = f"{docker_volume_path}/outputs/approx_PC1_pattern/Rao_2014"

    for species in ["Human", "Mouse"]:
        if species == "Human":
            cells = ["GM12878", "HeLa", "HMEC", "HUVEC", "IMR90", "K562", "KBM7", "NHEK"]
            chroms = [str(i) for i in range(1, 23)]
        elif species == "Mouse":
            cells = ["CH12-LX"]
            chroms = [str(i) for i in range(1, 20)]
    
        chroms.extend(["X", "Y"])
        resolutions = ["1000000", "100000"]
        types = ["CxMax", "CxMin"]

        for resolution in resolutions:
            for cell in cells:
                for chrom in chroms:
                    for type in types:
                        kwargs = {
                            "Juicer_PC1": f"{Juicer_PC1_path}/{cell}/{resolution}/eigenvector/pc1_chr{chrom}.txt",
                            "approx": f"{approx_path}/{cell}/{resolution}/{type}/approx_PC1_pattern_chr{chrom}.txt",
                        }

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
        output_df.to_excel(writer, sheet_name="summary_2014")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='summary_correctness_2014.py',
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