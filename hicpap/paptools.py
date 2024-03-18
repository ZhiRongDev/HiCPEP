import os
import hicstraw
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import logging
logging.basicConfig(format='%(message)s', level=logging.INFO)

def read_pearson(pearson: str, zero_mean: bool=True) -> np.ndarray:
    """_summary_
    Note that the return pearson matrix will be symmetric.
    """
    pearson_df = pd.read_table(pearson, header=None, sep=" ")
    pearson_df.pop(pearson_df.columns[-1])
    pearson_df = pearson_df.fillna(0)
    pearson_np = pearson_df.values # Turn into numpy.ndarray

    # Zero mean the Pearson correlaton matrix but still keep the all 0 rows/columns.
    if zero_mean:
        diag = np.diag(pearson_np)
        diag_valid = diag != 0 
        ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.
        pearson_np[ixgrid] -= pearson_np[ixgrid].mean(axis=1, keepdims=True)

    return pearson_np

def straw_to_pearson(hic_path: str, chrom_x: str, chrom_y: str, resolution: int, norm: str="KR", method: str="oe", zero_mean: bool=True) -> np.ndarray:
    """_summary_

            # # # #
    chrom_y # # # #
            # # # #
            chrom_x

    Note that the return pearson matrix will be symmetric.
    """
    hic = hicstraw.HiCFile(hic_path)

    for chrom in hic.getChromosomes():
        if chrom.name == chrom_x:
            chrom_x_size = int(chrom.length)
        if chrom.name == chrom_y:
            chrom_y_size = int(chrom.length)

    matrix = hic.getMatrixZoomData(chrom_y, chrom_x, method, norm, "BP", resolution)
    matrix_np = matrix.getRecordsAsMatrix(0, chrom_y_size, 0, chrom_x_size)

    pearson_np = np.corrcoef(matrix_np)
    pearson_df = pd.DataFrame(pearson_np)
    pearson_df = pearson_df.fillna(0)
    pearson_np = pearson_df.values # Turn into numpy.ndarray
    del pearson_df

    # Zero mean the Pearson correlaton matrix but still keep the all 0 rows/columns.
    if zero_mean:
        diag = np.diag(pearson_np)
        diag_valid = diag != 0
        ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.
        pearson_np[ixgrid] -= pearson_np[ixgrid].mean(axis=1, keepdims=True)

    return pearson_np

# I should add an parameter which accepts the fasta file as an argument(to calc the GC content)
def create_approx(pearson_np: np.ndarray, output: str="None", method: str="cxmax") -> np.ndarray:
    if len(pearson_np) != len(pearson_np[0]): 
        logging.info("Pearson matrix given has a different number of rows and columns")
        return

    """_summary_
    Calaulate the covariance matrix of the zero-means pearson matrix, according to the steps in PCA.
    Note that we set the degree of freedom as n.
    """
    diag = np.diag(pearson_np)
    diag_valid = diag != 0
    ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.

    n = len(pearson_np[ixgrid])
    cov_np = np.zeros((len(pearson_np), len(pearson_np))) 
    cov_np[ixgrid] = np.matmul(pearson_np[ixgrid], pearson_np[ixgrid].T) / n # covariance matrix

    ### Core idea, note that the covariance matrix is symmetric
    cov_abs_sum = [np.sum(np.abs(row)) for row in cov_np] 
    cov_abs_sum = list(enumerate(cov_abs_sum)) # Turn list into tuple with index, ex: (index, absSum)
    sorted_cov_abs_sum = sorted(cov_abs_sum, key=lambda x: x[1], reverse=True) # Sorted from the maximum to the minimum 

    if method == "cxmax":
        sorted_index = 0
        cov_selected_np = cov_np[sorted_cov_abs_sum[sorted_index][0]]
    elif method == "cxmin":
        sorted_index = -1
        cov_selected_np = cov_np[sorted_cov_abs_sum[sorted_index][0]]
        # To avoid selecting the column with all 0 values. 
        while cov_selected_np.sum() == 0:
            sorted_index -= 1
            cov_selected_np = cov_np[sorted_cov_abs_sum[sorted_index][0]]

    if output == "None":
        return cov_selected_np
    
    os.makedirs(os.path.dirname(output), exist_ok=True)
    with open(output, 'w') as f:
        tmp_str = ''

        for i in cov_selected_np:
            tmp_str += f"{str(i)}\n"

        tmp_str = tmp_str[:-1]
        f.write(tmp_str)

    return cov_selected_np

def pca_on_pearson(pearson_np: np.ndarray):
    if len(pearson_np) != len(pearson_np[0]): 
        logging.info("Pearson matrix has a different number of rows and columns")
        return
    
    total_entry_num = len(pearson_np)
    diag = np.diag(pearson_np)
    diag_valid = diag != 0
    ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.

    n = valid_entry_num = len(pearson_np[ixgrid])
    y = pearson_np[ixgrid].T / np.sqrt(n)

    """_summary_
    These two lines of code will both calculate the covariance matrix of the pearson matrix, and is confirmed to have the same results.
    logging.info(np.matmul(y.T, y), '\n')
    logging.info(np.cov(pearson_np[ixgrid], bias=True)) # `bias=True` will set the degree of freedom as n.
    """
    U, S, Vh = np.linalg.svd(y, full_matrices=True)
    eigenvalues = S * S
    sum_eigenvalues = np.sum(eigenvalues)
    explained_variances = eigenvalues / sum_eigenvalues

    tmp = np.zeros((len(pearson_np), len(pearson_np))) 
    tmp[ixgrid] = Vh
    Vh = tmp

    # Return Principal components(Vector) and the explained variances of PC1(Vector).
    return Vh[diag_valid], explained_variances, total_entry_num, valid_entry_num

def calc_correctness(pc1_np: np.ndarray, approx_np: np.ndarray):
    total_entry_num = len(pc1_np)
    if total_entry_num != len(approx_np): 
        logging.info("pc1_np and approx_np has a different total_entry_num")
        return
    
    pc1_np = pc1_np[pc1_np != 0] # Remove 0
    approx_np = approx_np[approx_np != 0] # Remove 0
    valid_entry_num = len(pc1_np)
    if valid_entry_num != len(approx_np): 
        logging.info("pc1_np and approx_np has a different valid_entry_num")
        return

    if np.corrcoef(pc1_np, approx_np)[0][1] < 0:
        approx_np = -approx_np

    pc1_pos_np = pc1_np > 0
    approx_pos_np = approx_np > 0
    pc1_pos_vs_approx_pos_np = pc1_pos_np == approx_pos_np 
    correct_num = list(pc1_pos_vs_approx_pos_np).count(True)
    correct_rate = correct_num / valid_entry_num

    return {
        "total_entry_num": total_entry_num,
        "valid_entry_num": valid_entry_num,
        "correct_num": correct_num,
        "correct_rate": correct_rate,
    }

def plot_comparison(pc1_np: np.ndarray, approx_np: np.ndarray, figsize: int=20, scatter: str="None", relative_magnitude: str="None"):
    correctness_info = calc_correctness(pc1_np, approx_np)
    total_entry_num = correctness_info["total_entry_num"]
    valid_entry_num = correctness_info["valid_entry_num"]
    correct_num = correctness_info["correct_num"]
    correct_rate = correctness_info["correct_rate"]

    if np.corrcoef(pc1_np, approx_np)[0][1] < 0:
        approx_np = -approx_np

    if scatter != "None":
        plot_x_axis = [i + 1 for i in range(total_entry_num)]
        approx_dots = [1 if i > 0 else -1 if i < 0 else 0 for i in approx_np]
        pc1_colors_values = [2 if i > 0 else 0 if i < 0 else 1 for i in pc1_np]
        pc1_colors = ListedColormap(['r', 'g', 'b'])
        scatter_labels = ["PC1 < 0", "PC1 == 0", "PC1 > 0"]

        plt.figure(figsize=(figsize, 6))
        plt.xticks(np.arange(0, valid_entry_num, 50)) 
        scatter_config =  plt.scatter(plot_x_axis, approx_dots, c=pc1_colors_values, cmap=pc1_colors)
        plt.legend(handles=scatter_config.legend_elements()[0], labels=scatter_labels, fontsize="20", loc="center left")
        plt.title(f"total_entry_num: {total_entry_num}, valid_entry_num: {valid_entry_num}, correct_num = {correct_num}, correct_rate={np.round(correct_rate, 2)}", fontsize=20, loc="left")
        plt.savefig(scatter)
        plt.clf() 

    if relative_magnitude != "None":
        approx_np_norm = (approx_np - np.mean(approx_np)) / np.std(approx_np)
        pc1_np_norm = (pc1_np - np.mean(pc1_np)) / np.std(pc1_np)
        
        plt.figure(figsize=(figsize, 6))
        plt.xticks(np.arange(0, valid_entry_num, 50)) 
        plt.plot(pc1_np_norm, c='r')
        plt.plot(approx_np_norm, c='b')
        plt.legend(["PC1", "approximated PC1-pattern"], fontsize="20", loc ="upper left")
        plt.title(f"total_entry_num: {total_entry_num}, valid_entry_num: {valid_entry_num}", fontsize=20, loc="left")
        plt.savefig(relative_magnitude)        
        plt.clf()
    
    plt.close('all')
    return

def flip_track(track_np: np.ndarray, gc_df: pd.DataFrame, chrom: str) -> np.ndarray:
    tmp_df = gc_df.loc[gc_df['chrom'] == chrom]
    gc_np = tmp_df["GC"].fillna(0).values

    if np.mean(gc_np[track_np > 0]) < np.mean(gc_np[track_np < 0]):
        track_np = -track_np

    return track_np