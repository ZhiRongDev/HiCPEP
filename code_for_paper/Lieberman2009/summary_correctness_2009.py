# import numpy as np
# import pandas as pd

# def summary_correctness(Juicer_PC1, approx):
#     # Read in the Eigenvector 1
#     EV1_df = pd.read_table(Juicer_PC1, header=None)
#     EV1_df = EV1_df.dropna(axis=0, how="all").reset_index(drop=True)
#     EV1_np = EV1_df.values # Turn into numpy format
#     EV1_np = EV1_np.flatten() # Turn into 1D vector

#     approx_df = pd.read_table(approx, header=None)
#     approx_np = approx_df.values # Turn into numpy format
#     approx_np = approx_np.flatten() # Turn into 1D vector

#     del EV1_df, approx_df

#     if len(EV1_np) != len(approx_np): 
#         print("juicer_PC1 and approx has a different number of elements")
#         return
    
#     entryNum = len(EV1_np)

#     if np.corrcoef(EV1_np, approx_np)[0][1] < 0:
#         approx_np = -approx_np

#     EV1_pos_np = EV1_np > 0
#     approx_pos_np = approx_np > 0
#     EV1_pos_VS_approx_pos_np = EV1_pos_np == approx_pos_np 
    
#     correctNum = list(EV1_pos_VS_approx_pos_np).count(True)
#     correctRate = correctNum / entryNum

#     return {
#         "entryNum": entryNum,
#         "correctNum": correctNum,
#         "correctRate": correctRate,
#     }

# def write_to_excel(output_file):
#     CxMax_df = pd.DataFrame(columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"])
#     CxMin_df = pd.DataFrame(columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"])
#     Juicer_PC1_path = "/home/jordan990301/Projects/HiC-PC1_approx/data/Rao2014/juicer_outputs"
#     approx_path = "/home/jordan990301/Projects/HiC-PC1_approx/outputs/approx_PC1_pattern"

#     for species in ["Human", "Mouse"]:
#         if species == "Human":
#             cells = ["GM12878", "HeLa", "HMEC", "HUVEC", "IMR90", "K562", "KBM7", "NHEK"]
#             chroms = [str(i) for i in range(1, 23)]
#         elif species == "Mouse":
#             cells = ["CH12-LX"]
#             chroms = [str(i) for i in range(1, 20)]
    
#         chroms.extend(["X", "Y"])
#         resolutions = ["1Mb", "100Kb"]
#         types = ["CxMax", "CxMin"]

#         for resolution in resolutions:
#             for cell in cells:
#                 for chrom in chroms:
#                     for type in types:
#                         kwargs = {
#                             "Juicer_PC1": f"{Juicer_PC1_path}/{cell}/{resolution}/data/eigenvector/pc1_chr{chrom}.txt",
#                             "approx": f"{approx_path}/{cell}/{resolution}/{type}/approx_PC1_pattern_chr{chrom}.txt",
#                         }

#                         correctness_info = summary_correctness(**kwargs)

#                         new_row_df = pd.DataFrame(
#                             [[cell, resolution, f"chr{chrom}", type, correctness_info["entryNum"], correctness_info["correctNum"], correctness_info["correctRate"]]],
#                             columns=['cell', 'resolution', 'chromosome', "type", "entryNum", "correctNum", "correctRate"]
#                         )

#                         if type == "CxMax":
#                             CxMax_df = pd.concat([CxMax_df, new_row_df], ignore_index=True)
#                         elif type == "CxMin":
#                             CxMin_df = pd.concat([CxMin_df, new_row_df], ignore_index=True)

#     output_df = pd.concat([CxMax_df, CxMin_df], ignore_index=True)
#     with pd.ExcelWriter(output_file, mode="w") as writer:
#         output_df.to_excel(writer, sheet_name="summary")


# if __name__ == "__main__":
#     output_file = "/home/jordan990301/Projects/HiC-PC1_approx/outputs/summary.xlsx"
#     write_to_excel(output_file)