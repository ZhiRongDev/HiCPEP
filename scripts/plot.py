import argparse
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

parser = argparse.ArgumentParser(
    prog='plot.py',
    allow_abbrev=False,
    description='What the program does',
    epilog='Text at the bottom of help'
)
parser.add_argument(
    "--Juicer_PC1",
    type=str,
    required=True,
    help="Input blablabla"
)
parser.add_argument(
    "--approx",
    type=str,
    required=True,
    help="Input int"
)
parser.add_argument(
    "--relative_magnitude",
    type=str,
    default="None",
    help="Input int"
)
parser.add_argument(
    "--scatter",
    type=str,
    default="None",
    help="Input int"
)
parser.add_argument(
    "--resolution",
    type=str,
    default="1Mb",
    help="Input int"
)

args = parser.parse_args()

def plot(**kwargs):
    # Read in the Eigenvector 1
    EV1_df = pd.read_table(kwargs["Juicer_PC1"], header=None)
    EV1_df = EV1_df.dropna(axis=0, how="all").reset_index(drop=True)
    EV1_np = EV1_df.values # Turn into numpy format
    EV1_np = EV1_np.flatten() # Turn into 1D vector

    approx_df = pd.read_table(kwargs["approx"], header=None)
    approx_np = approx_df.values # Turn into numpy format
    approx_np = approx_np.flatten() # Turn into 1D vector

    del EV1_df, approx_df

    if len(EV1_np) != len(approx_np): 
        print("juicer_PC1 and approx has a different number of elements")
        return
    
    entryNum = len(EV1_np)

    if np.corrcoef(EV1_np, approx_np)[0][1] < 0:
        approx_np = -approx_np

    EV1_pos_np = EV1_np > 0
    approx_pos_np = approx_np > 0
    EV1_pos_VS_approx_pos_np = EV1_pos_np == approx_pos_np 
    
    correctNum = list(EV1_pos_VS_approx_pos_np).count(True)
    correctRate = correctNum / entryNum

    # Visualization
    if kwargs["resolution"] == "25Kb":
        figsize = 50
    else:
        figsize = 30
    
    if kwargs["scatter"] != "None":
        plot_x_axis = [i + 1 for i in range(entryNum)]
        approx_Dots = [1 if i else -1 for i in approx_pos_np]
        EV1_colors_values = [1 if i else 0 for i in EV1_pos_np]
        EV1_colors = ListedColormap(['r', 'b'])
        scatter_labels = ["Juicer's PC1 < 0", "Juicer's PC1 > 0"]

        plt.xticks(np.arange(0, entryNum, 50)) 
        # plt.rcParams["figure.figsize"] = [figsize, 5]
        # plt.rcParams["figure.autolayout"] = True
        plt.figure(figsize=(figsize, 6))
        scatter =  plt.scatter(plot_x_axis, approx_Dots, c=EV1_colors_values, cmap=EV1_colors)
        plt.legend(handles=scatter.legend_elements()[0], labels=scatter_labels, fontsize="20", loc="center left")
        plt.title(f"entryNum: {entryNum}, correctNum = {correctNum}, correctRate={np.round(correctRate, 2)}", fontsize=20, loc="left")
        plt.savefig(f'{kwargs["scatter"]}/1_scatter.png')
        plt.clf() 
    
    if kwargs["relative_magnitude"] != "None":
        approx_np_Norm = (approx_np - np.mean(approx_np)) / np.std(approx_np)
        EV1_np_Norm = (EV1_np - np.mean(EV1_np)) / np.std(EV1_np)
        
        plt.xticks(np.arange(0, entryNum, 50)) 
        # plt.rcParams["figure.figsize"] = [figsize, 5]
        # plt.rcParams["figure.autolayout"] = True
        plt.figure(figsize=(figsize, 6))
        plt.plot(EV1_np_Norm, c='r')
        plt.plot(approx_np_Norm, c='b')
        plt.legend(["Juicer's PC1", "approximated PC1-pattern"], fontsize="20", loc ="upper left")
        plt.title(f"entryNum: {entryNum}", fontsize=20, loc="left")
        plt.savefig(f'{kwargs["relative_magnitude"]}/1_line.png')
        plt.clf()

plot(Juicer_PC1=args.Juicer_PC1, approx=args.approx, relative_magnitude=args.relative_magnitude, scatter=args.scatter, resolution=args.resolution)