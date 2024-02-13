import argparse
import numpy as np
import pandas as pd
from copy import deepcopy
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
np.set_printoptions(suppress=True)

parser = argparse.ArgumentParser(
    prog='Rao2014_approximation.py',
    allow_abbrev=False,
    description='What the program does',
    epilog='Text at the bottom of help'
)

parser.add_argument(
    "--pearsons",
    type=str,
    required=True,
    help="Input blablabla"
)
parser.add_argument(
    "--eigenvector",
    type=str,
    required=True,
    help="Input int"
)
parser.add_argument(
    "--norm",
    type=str,
    default="NONE",
    help="Input int"
)
parser.add_argument(
    "--binsize",
    type=str,
    default="NONE",
    help="Input int"
)
parser.add_argument(
    "--approx",
    type=str,
    default="None",
    help="Input int"
)
parser.add_argument(
    "--scatterPlot",
    type=str,
    default="None",
    help="Input int"
)
parser.add_argument(
    "--linePlot",
    type=str,
    default="None",
    help="Input int"
)
parser.add_argument(
    "--log",
    type=str,
    default="None",
    help="Input int"
)

args = parser.parse_args()


# Read in the Pearson correlatin matrix
pearson_df = pd.read_table(f"{input_path}/pearsons/pearson_chr{chrom}.txt", header=None, sep=" ")
pearson_df.pop(pearson_df.columns[-1])
pearson_df = pearson_df.dropna(axis=0, how="all").reset_index(drop=True)
pearson_df = pearson_df.dropna(axis=1, how="all")
pearson_np = deepcopy(pearson_df.values) # Turn into numpy format
pearson_np = pearson_np - pearson_np.mean(axis=1, keepdims=True) # Zero mean of Pearson correlaton matrix

# Read in the Eigenvector 1
EV1_df = pd.read_table(f"{input_path}/eigenvector/pc1_chr{chrom}.txt", header=None, sep=" ")
EV1_df = EV1_df.dropna(axis=0, how="all").reset_index(drop=True)
EV1_np = deepcopy(EV1_df.values) # Turn into numpy format
EV1_np = EV1_np.flatten() # Turn into 1D vector

del pearson_df, EV1_df

if len(pearson_np) == len(EV1_np) and len(pearson_np) == len(pearson_np[0]):
    # According the steps in SVD, set x as pearson_df, and set y as x' / np.sqrt(n)
    n = len(pearson_np)
    pearson_T_np = deepcopy(np.transpose(pearson_np))
    y_np = deepcopy(pearson_T_np / np.sqrt(n))
    cov_np = deepcopy(np.matmul(np.transpose(y_np), y_np)) # Covariance matrix of pearson_df

    # Main idea
    cov_absSum = [np.sum(np.abs(row)) for row in cov_np] 
    cov_absSum = list(enumerate(cov_absSum)) # Turn list into tuple with index, ex: (index, absSum)
    sorted_cov_absSum = sorted(cov_absSum, key=lambda x: x[1], reverse=True) 

    for sorted_index in [0, -1]:
        # The sign of the pearson with the largest absSum in cov_pearson_np_absSum should correspond with the patterns of EV1.
        cov_selected_np = cov_np[sorted_cov_absSum[sorted_index][0]]

        # Flip the sign if the corrcoef of cov_pearson_np_Selected and EV1_np is negative.
        if np.corrcoef(cov_selected_np, EV1_np)[0][1] < 0:
            cov_selected_np = -cov_selected_np

        EV1_pos_np = EV1_np > 0
        cov_selected_pos_np = cov_selected_np > 0
        EV1_pos_VS_cov_selected_pos_np = EV1_pos_np == cov_selected_pos_np

        output_path="/home/jordan990301/Projects/HiC-PC1_Approximation/outputs"

        if sorted_index == 0:
            cov_selected_type = "Max"
            log_path = f"{output_path}/logs/EV1-covD_absSumMax/{cell_line}/{resolution}" 
            linePlot_path = f"{output_path}/plots/EV1-covD_absSumMax/{cell_line}/{resolution}/line" 
            scatterPlot_path = f"{output_path}/plots/EV1-covD_absSumMax/{cell_line}/{resolution}/scatter" 
        elif sorted_index == -1:
            cov_selected_type = "Min"
            log_path = f"{output_path}/logs/EV1-covD_absSumMin/{cell_line}/{resolution}" 
            linePlot_path = f"{output_path}/plots/EV1-covD_absSumMin/{cell_line}/{resolution}/line" 
            scatterPlot_path = f"{output_path}/plots/EV1-covD_absSumMin/{cell_line}/{resolution}/scatter" 
        
        # with open(f"{log_path}/chr{chrom}/patterns.txt", "w+") as f:
        #     for i in EV1_pos_VS_cov_selected_pos_np:
        #         f.write(str(i))
        #         f.write('\n')
    
        correctNum = list(EV1_pos_VS_cov_selected_pos_np).count(True)
        correctRate = correctNum / len(pearson_np)

        new_row_df = pd.DataFrame(
            [[cell_line, resolution, f"chr{chrom}", cov_selected_type, len(pearson_np), correctNum, correctRate]],
            columns=['cellLine', 'resolution', 'chromosome', "cov_selected_type(absSum)", "binsNum", "correctNum", "correctRate"]
        )

        if sorted_index == 0:
            outputMax_df = pd.concat([outputMax_df, new_row_df], ignore_index=True)
        elif sorted_index == -1:
            outputMin_df = pd.concat([outputMin_df, new_row_df], ignore_index=True)
    
        # Visualization
        plot_x_axis = [i + 1 for i in range(len(pearson_np))]
        cov_selected_Dots = [1 if i else -1 for i in cov_selected_pos_np]
        EV1_colors_values = [1 if i else 0 for i in EV1_pos_np]
        EV1_colors = ListedColormap(['r', 'b'])
        scatter_labels = ["Juicer's PC1 < 0", "Juicer's PC1 > 0"]

        plt.xticks(np.arange(0, len(pearson_np), 50)) 
        plt.rcParams["figure.figsize"] = [figsize, 5]
        plt.rcParams["figure.autolayout"] = True
        scatter =  plt.scatter(plot_x_axis, cov_selected_Dots, c=EV1_colors_values, cmap=EV1_colors)
        plt.legend(handles=scatter.legend_elements()[0], labels=scatter_labels, fontsize="20", loc="center left")
        # print(scatter.legend_elements()[0])
        plt.title(f"chromosome: {cell_line}_chromosome{chrom}, resolution: {resolution}, entryNum: {len(pearson_np)}, correctNum = {correctNum}, correctRate={np.round(correctRate, 2)}", fontsize=20, loc="left")
        plt.savefig(f'{scatterPlot_path}/{cell_line}_chr{chrom}.png')
        plt.clf() 

        cov_pearson_np_Selected_Norm = deepcopy((cov_selected_np - np.mean(cov_selected_np)) / np.std(cov_selected_np))
        EV1_np_Norm = deepcopy((EV1_np - np.mean(EV1_np)) / np.std(EV1_np))
        
        plt.xticks(np.arange(0, len(pearson_np), 50)) 
        plt.rcParams["figure.figsize"] = [figsize, 5]
        plt.rcParams["figure.autolayout"] = True
        plt.plot(EV1_np_Norm, c='r')
        plt.plot(cov_pearson_np_Selected_Norm, c='b')
        plt.legend(["Juicer's PC1", "approximated PC1-pattern"], fontsize="20", loc ="upper left")
        plt.title(f"chromosome: {cell_line}_chromosome{chrom}, resolution: {resolution}, entryNum: {len(pearson_np)}", fontsize=20, loc="left")
        plt.savefig(f'{linePlot_path}/{cell_line}_chr{chrom}.png')
        plt.clf()