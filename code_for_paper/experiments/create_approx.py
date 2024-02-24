# This code is for creating the approximated PC1-pattern. 
import os
import numpy as np
import pandas as pd

def create_approx(pearson, output, type, source):
    # Read in the Pearson correlatin matrix
    if source == "Rao_2014":
        pearson_df = pd.read_table(pearson, header=None, sep=" ")
    elif source == "Lieberman_2009":
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
    cov_absSum = [np.sum(np.abs(row)) for row in cov_np] 
    cov_absSum = list(enumerate(cov_absSum)) # Turn list into tuple with index, ex: (index, absSum)
    sorted_cov_absSum = sorted(cov_absSum, key=lambda x: x[1], reverse=True) # Sorted from the maximum to the minimum 

    if type == "CxMax":
        sorted_index = 0
        cov_selected_np = cov_np[sorted_cov_absSum[sorted_index][0]]
    elif type == "CxMin":
        sorted_index = -1
        cov_selected_np = cov_np[sorted_cov_absSum[sorted_index][0]]
        # To avoid selecting the column with all 0 values. 
        while cov_selected_np.sum() == 0:
            sorted_index -= 1
            cov_selected_np = cov_np[sorted_cov_absSum[sorted_index][0]]

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