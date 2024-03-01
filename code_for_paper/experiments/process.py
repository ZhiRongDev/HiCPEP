import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

# np.set_printoptions(suppress=True)
# np.set_printoptions(precision=2)

def create_approx(pearson, output, method, source="2014"):
    # Read in the Pearson correlatin matrix
    if source == "2014":
        pearson_df = pd.read_table(pearson, header=None, sep=" ")
    elif source == "2009":
        pearson_df = pd.read_table(pearson, index_col=0, header=1, sep="\t")

    pearson_df.pop(pearson_df.columns[-1])
    pearson_df = pearson_df.fillna(0)
    pearson_np = pearson_df.values # Turn into numpy.ndarray
    
    diag = np.diag(pearson_np)
    diag_valid = diag != 0 
    ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.
    pearson_np[ixgrid] -= pearson_np[ixgrid].mean(axis=1, keepdims=True) # Zero mean the Pearson correlaton matrix but still keep the all 0 rows/columns.

    del pearson_df

    if len(pearson_np) != len(pearson_np[0]): 
        print("Pearson matrix has a different number of rows and columns")
        return

    """_summary_
    Calaulate the covariance matrix of the zero-means pearson matrix, according to the steps in PCA.
    Note that we set the degree of freedom as n.
    """
    n = len(pearson_np[ixgrid])
    cov_np = np.zeros((len(pearson_np), len(pearson_np))) 
    cov_np[ixgrid] = np.matmul(pearson_np[ixgrid], pearson_np[ixgrid].T) / n # covariance matrix

    # Core idea, note that the covariance matrix is symmetric
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
        for val in cov_selected_np:
            print(val)
        return
    
    os.makedirs(os.path.dirname(output), exist_ok=True)

    with open(output, 'w') as f:
        tmp_str = ''

        for i in cov_selected_np:
            tmp_str += f"{str(i)}\n"

        tmp_str = tmp_str[:-1]
        f.write(tmp_str)
    return

def calc_correctness(pc1, approx, source="2014"):
    if source == "2014":
        pc1_df = pd.read_table(pc1, header=None)
        pc1_df = pc1_df.fillna(0)
    elif source == "2009":
        pc1_df = pd.read_table(pc1, header=None, sep="\t")
        pc1_df = pc1_df.iloc[:, [2]]

    pc1_np = pc1_df.values # Turn into numpy format
    pc1_np = pc1_np.flatten() # Turn into 1D vector
    approx_df = pd.read_table(approx, header=None)
    approx_np = approx_df.values # Turn into numpy format
    approx_np = approx_np.flatten() # Turn into 1D vector

    total_entry_num = len(pc1_np)
    if total_entry_num != len(approx_np): 
        print("PC1 and approx has a different total_entry_num")
        return
    
    pc1_np = pc1_np[pc1_np != 0] # Remove 0
    approx_np = approx_np[approx_np != 0] # Remove 0
    valid_entry_num = len(pc1_np)
    if valid_entry_num != len(approx_np): 
        print("PC1 and approx has a different valid_entry_num")
        return

    del pc1_df, approx_df

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

def plot_comparison(pc1, approx, figsize, scatter, relative_magnitude, source="2014"):
    if source == "2014":
        pc1_df = pd.read_table(pc1, header=None)
        pc1_df = pc1_df.fillna(0)
    elif source == "2009":
        pc1_df = pd.read_table(pc1, header=None, sep="\t")
        pc1_df = pc1_df.iloc[:, [2]]
    
    pc1_np = pc1_df.values # Turn into numpy format
    pc1_np = pc1_np.flatten() # Turn into 1D vector

    approx_df = pd.read_table(approx, header=None)
    approx_np = approx_df.values # Turn into numpy format
    approx_np = approx_np.flatten() # Turn into 1D vector

    del pc1_df, approx_df

    correctness_info = calc_correctness(pc1, approx, source)
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

def calc_explained_variance(pearson, source="2014"):
    # Read in the Pearson correlatin matrix
    if source == "2014":
        pearson_df = pd.read_table(pearson, header=None, sep=" ")
    elif source == "2009":
        pearson_df = pd.read_table(pearson, index_col=0, header=1, sep="\t")

    pearson_df.pop(pearson_df.columns[-1])
    pearson_df = pearson_df.fillna(0)
    pearson_np = pearson_df.values # Turn into numpy.ndarray
    total_entry_num = len(pearson_np)
    
    diag = np.diag(pearson_np)
    diag_valid = diag != 0 
    ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.
    pearson_np[ixgrid] -= pearson_np[ixgrid].mean(axis=1, keepdims=True) # Zero mean the Pearson correlaton matrix but still keep the all 0 rows/columns.

    n = valid_entry_num = len(pearson_np[ixgrid])
    y = pearson_np[ixgrid].T / np.sqrt(n)

    """_summary_
    These two lines of code will both calculate the covariance matrix of the pearson matrix, and is confirmed to have the same results.
    print(np.matmul(y.T, y), '\n')
    print(np.cov(pearson_np[ixgrid], bias=True)) # `bias=True` will set the degree of freedom as n.
    """

    U, S, Vh = np.linalg.svd(y, full_matrices=True)
    eigenvalues = S * S
    sum_eigenvalues = np.sum(eigenvalues)
    explained_variances = eigenvalues / sum_eigenvalues

    tmp = np.zeros((len(pearson_np), len(pearson_np))) 
    tmp[ixgrid] = Vh
    Vh = tmp

    # Return Principal components(Vector) and the explained variances of PC1(Vector).
    return Vh, explained_variances, total_entry_num, valid_entry_num