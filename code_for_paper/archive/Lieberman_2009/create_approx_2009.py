# This code is for creating the approximated PC1-pattern. 
import os
import argparse
import numpy as np
import pandas as pd

# np.set_printoptions(threshold = np.inf)

def create_approx(pearson, output, type):
    # Read in the Pearson correlatin matrix
    pearson_df = pd.read_table(pearson, index_col=0, header=1, sep="\t")
    pearson_df.pop(pearson_df.columns[-1])
    pearson_np = pearson_df.values # Turn into numpy.ndarray

    # In Lieberman's datasets, the NaN value are stored as zero.
    # https://stackoverflow.com/questions/11188364/remove-zero-lines-2-d-numpy-array
    pearson_np = pearson_np[~np.all(pearson_np == 0, axis=1)] # remove rows with all 0 
    pearson_np = pearson_np[:, ~np.all(pearson_np == 0, axis=0)] # remove columns with all 0

    # Zero means
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

    # Decide to choose CxMax or CxMin
    if type == "CxMax":
        sorted_index = 0
    elif type == "CxMin":
        sorted_index = -1

    cov_selected_np = cov_np[sorted_cov_absSum[sorted_index][0]]

    if output == "None":
        return
    
    filename = output
    os.makedirs(os.path.dirname(filename), exist_ok=True)

    with open(filename, 'w') as f:
        tmp_str = ''
        for i in cov_selected_np:
            tmp_str += f"{str(i)}\n"

        tmp_str = tmp_str[:-1]
        
        f.write(tmp_str)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='xxx.py',
        allow_abbrev=False,
        description='What the program does', epilog='Text at the bottom of help'
    )
    parser.add_argument(
        "--pearson",
        type=str,
        required=True,
        help="Input blablabla"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="None",
        help="Input int"
    )
    parser.add_argument(
        "--type",
        type=str,
        default="CxMax",
        help="Input int"
    )

    args = parser.parse_args()
    kwargs = vars(args) # Turn into dict

    create_approx(**kwargs)