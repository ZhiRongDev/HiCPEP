import os
import numpy as np
import pandas as pd

def create_approx(pearson, output, type, source="2014"):
    # Read in the Pearson correlatin matrix
    if source == "2014":
        pearson_df = pd.read_table(pearson, header=None, sep=" ")
    elif source == "2009":
        pearson_df = pd.read_table(pearson, index_col=0, header=1, sep="\t")

    pearson_df.pop(pearson_df.columns[-1])
    pearson_df = pearson_df.fillna(0)
    pearson_np = pearson_df.values # Turn into numpy.ndarray
    pearson_np = pearson_np - pearson_np.mean(axis=1, keepdims=True) # Zero mean of Pearson correlaton matrix

    del pearson_df

    if len(pearson_np) != len(pearson_np[0]): 
        print("Pearson matrix has a different number of rows and columns")
        return

    # According the steps in SVD, here we set the degree of freedom as n   
    n = len(pearson_np[0])
    cov_np = np.matmul(pearson_np, pearson_np.T) / n

    # Main idea, note that the covariance matrix is symmetric
    cov_abs_sum = [np.sum(np.abs(row)) for row in cov_np] 
    cov_abs_sum = list(enumerate(cov_abs_sum)) # Turn list into tuple with index, ex: (index, absSum)
    sorted_cov_abs_sum = sorted(cov_abs_sum, key=lambda x: x[1], reverse=True) # Sorted from the maximum to the minimum 

    if type == "cxmax":
        sorted_index = 0
        cov_selected_np = cov_np[sorted_cov_abs_sum[sorted_index][0]]
    elif type == "cxmin":
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
    
    filename = output
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    with open(filename, 'w') as f:
        tmp_str = ''

        for i in cov_selected_np:
            tmp_str += f"{str(i)}\n"

        tmp_str = tmp_str[:-1]
        f.write(tmp_str)

    return

# Accept pandas Dataframe as the parameters.
def calc_correctness(pc1_df, approx_df):
    pc1_np = pc1_df.values # Turn into numpy format
    pc1_np = pc1_np.flatten() # Turn into 1D vector
    pc1_np = pc1_np[pc1_np != 0] # Remove 0

    approx_np = approx_df.values # Turn into numpy format
    approx_np = approx_np.flatten() # Turn into 1D vector
    approx_np = approx_np[approx_np != 0] # Remove 0

    if len(pc1_np) != len(approx_np): 
        print("PC1 and approx has a different number of elements")
        return
    
    valid_entry_num = len(pc1_np)

    if np.corrcoef(pc1_np, approx_np)[0][1] < 0:
        approx_np = -approx_np

    pc1_pos_np = pc1_np > 0
    approx_pos_np = approx_np > 0
    pc1_pos_vs_approx_pos_np = pc1_pos_np == approx_pos_np 
    
    correct_num = list(pc1_pos_vs_approx_pos_np).count(True)
    correct_rate = correct_num / valid_entry_num

    return {
        "valid_entry_num": valid_entry_num,
        "correct_num": correct_num,
        "correct_rate": correct_rate,
    }

def plot_comparison():
    return
