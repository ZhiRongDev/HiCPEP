import os
import hicstraw
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import logging
logging.basicConfig(format='%(message)s', level=logging.INFO)

# Please make sure the pearson_np is already fillna with 0.
def zero_means_pearson(pearson_np: np.ndarray) -> np.ndarray:
    """
    Perform zero means on the given Hi-C Pearson matrix. 
    The calculation is only performed on the valid sub-matrix (rows/columns with all ``NaN`` will be excluded, however these all ``NaN`` rows/columns will not be removed in the Pearson matrix returned).

    :param pearson_np: The origin Hi-C Pearson matrix, which typically includes some all ``NaN`` rows/columns.
    :type pearson_np: numpy.ndarray

    :return: Zero-means Hi-C Pearson matrix. 
    :rtype: numpy.ndarray
    """
    pearson_np = pearson_np.astype('float64')
    diag = np.diag(pearson_np)
    diag_valid = ~np.isnan(diag)
    ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.
    pearson_np[ixgrid] -= pearson_np[ixgrid].mean(axis=1, keepdims=True)

    return pearson_np

def read_pearson(pearson: str) -> np.ndarray:
    """
    Read a intra-chromosomal Hi-C Pearson matrix's ``.txt`` file created by `juicer_tools <https://github.com/aidenlab/juicer/wiki/Pearsons>`_ 
    and return a Pearson matrix in NumPy format.

    :param pearson: Path of the juicer_tools created intra-chromosomal Hi-C Pearson matrix ``.txt`` file.
    :type pearson: str

    :return: Intra-chromosomal Hi-C Pearson matrix in NumPy format.
    :rtype: numpy.ndarray
    """
    pearson_df = pd.read_table(pearson, header=None, sep="\s+")
    pearson_np = pearson_df.values # Turn into numpy.ndarray

    return pearson_np

def straw_to_pearson(hic_path: str, chrom: str, resolution: int, normalization: str="KR") -> np.ndarray:
    """
    Read a ``.hic`` file created by `Juicer <https://github.com/aidenlab/juicer>`_ and return the Pearson matrix in NumPy format. 
    We only support creating the intra-chromosomal Hi-C Pearson matrix from O/E matrix, please check the `Straw API reference <https://pypi.org/project/hic-straw/>`_ if this function doesn't meet your needs.

    :param hic_path: Path of the Juicer created ``.hic`` file.
    :type hic_path: str

    :param chrom: chromosome name in string, from chromosome 1 to chromosome 22 and chromosome X, chromosome Y. (e.g. ``'1'``, ``'22'``, ``'X'``, ``'Y'``, etc.)
    :type chrom: str

    :param resolution: Typically ``2500000``, ``1000000``, ``500000``, ``100000``, ``50000``, ``25000``, ``10000``, ``5000``, etc. 
    :type resolution: int

    :param normalization: ``'NONE'``, ``'VC'``, ``'VC_SQRT'``, ``'KR'``, ``'SCALE'``, etc.
    :type normalization: str

    :return: Intra-chromosomal Pearson matrix in NumPy format.
    :rtype: numpy.ndarray
    """
    hic = hicstraw.HiCFile(hic_path)

    for chromosome in hic.getChromosomes():
        if chromosome.name == chrom:
            chrom_size = int(chromosome.length)

    matrix = hic.getMatrixZoomData(chrom, chrom, "oe", normalization, "BP", resolution)
    matrix_np = matrix.getRecordsAsMatrix(0, chrom_size, 0, chrom_size)
    pearson_np = np.corrcoef(matrix_np)

    return pearson_np

def create_approx(pearson_np: np.ndarray, output: str | None = None, method: str="cxmax") -> np.ndarray:
    """
    Create the Approximated PC1-pattern of the given Hi-C Pearson matrix. 
    The calculation is only performed on the valid sub-matrix. (rows/columns with all ``NaN`` will be excluded, however these all ``NaN`` rows/columns will not be removed in the Approximated PC1-pattern returned)

    :param pearson_np: Hi-C Pearson matrix in NumPy format.
    :type pearson_np: numpy.ndarray
    
    :param output: (Optional) If the file path is specified, the Approximated PC1-pattern will be stored (e.g. ``output="./test/approx_pc1.txt"``).
    :type output: str
    
    :param method: ``cxmax`` or ``cxmin``.
    :type method: str

    :return: Approximated PC1-pattern in NumPy format.
    :rtype: numpy.ndarray.
    """
    if len(pearson_np) != len(pearson_np[0]): 
        logging.info("Pearson matrix given has a different number of rows and columns")
        return

    # Zero means pearson_np
    pearson_np = zero_means_pearson(pearson_np=pearson_np)

    diag = np.diag(pearson_np)
    diag_valid = ~np.isnan(diag)
    ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.

    # Core idea, note that the covariance matrix is symmetric
    n = len(pearson_np[ixgrid])
    cov_np = np.zeros((len(pearson_np), len(pearson_np))) 
    cov_np[ixgrid] = np.matmul(pearson_np[ixgrid], pearson_np[ixgrid].T) / n # covariance matrix

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

    # Add back the origin `NaN` values to cov_selected_np 
    tmp = np.full(len(cov_selected_np), np.nan)
    tmp[diag_valid] = cov_selected_np[diag_valid]
    cov_selected_np = tmp

    if output == None:
        return cov_selected_np
    
    if os.path.dirname(output):
        os.makedirs(os.path.dirname(output), exist_ok=True)

    with open(output, 'w') as f:
        tmp_str = ''

        for i in cov_selected_np:
            tmp_str += f"{str(i)}\n"

        tmp_str = tmp_str[:-1]
        f.write(tmp_str)

    return cov_selected_np

def pca_on_pearson(pearson_np: np.ndarray):
    """
    Perform PCA on the given Hi-C Pearson matrix.
    The calculation is only performed on the valid sub-matrix. (rows/columns with all ``NaN`` will be excluded, however these all ``NaN`` rows/columns will not be removed in the Principal component vectors returned)

    Note that we set the degree of freedom as n (length of Pearson matrix's x-axis).

    :param pearson_np: Hi-C Pearson matrix in NumPy format.
    :type pearson_np: numpy.ndarray

    :return:

    .. code::
        
        Vh (numpy.ndarray) : Principal Component vectors of the valid sub-matrix, start from the PC1 in index 0, PC2 in index 1, PC3 in index 2, etc.
        explained_variances (numpy.ndarray): A vector contains the explained variance of PC1, PC2, PC3 etc.
        total_entry_num (int): Entry numbers including NaN.
        valid_entry_num (int): Entry numbers excluding NaN.

    :rtype: numpy.ndarray, numpy.ndarray, int, int
    """
    if len(pearson_np) != len(pearson_np[0]): 
        logging.info("Pearson matrix has a different number of rows and columns")
        return

    # Zero means pearson_np
    pearson_np = zero_means_pearson(pearson_np=pearson_np)
    
    total_entry_num = len(pearson_np)
    diag = np.diag(pearson_np)
    diag_valid = ~np.isnan(diag)
    ixgrid = np.ix_(diag_valid, diag_valid) # Extract the submatrix.

    n = valid_entry_num = len(pearson_np[ixgrid])
    y = pearson_np[ixgrid].T / np.sqrt(n)

    U, S, Vh = np.linalg.svd(y, full_matrices=True)
    eigenvalues = S * S
    sum_eigenvalues = np.sum(eigenvalues)
    explained_variances = eigenvalues / sum_eigenvalues

    tmp = np.full((len(pearson_np), len(pearson_np)), np.nan) 
    tmp[ixgrid] = Vh
    Vh = tmp

    # Return Principal components(Vector) and the explained variances of PC1(Vector).
    return Vh[diag_valid], explained_variances, total_entry_num, valid_entry_num

def calc_similarity(pc1_np: np.ndarray, approx_np: np.ndarray):
    """
    Compare the similarity information between the given PC1 and Approximated PC1-pattern.
    (The ``NaN`` value entries will be excluded).

    :param pc1_np: PC1 in NumPy format. 
    :type pc1_np: numpy.ndarray.

    :param approx_np: Approximated PC1-pattern in NumPy format. 
    :type approx_np: numpy.ndarray.

    :return:

    .. code::

        {
            total_entry_num (int): Entry numbers including NaN.
            valid_entry_num (int): Entry numbers excluding NaN.
            similar_num (int): Number of entries that the PC1 and Approximated PC1-pattern have the some sign (positive/negative).
            similar_rate (float): similar_num / valid_entry_num.
        }

    :rtype: dict
    """
    total_entry_num = len(pc1_np)
    
    if total_entry_num != len(approx_np): 
        logging.info("pc1_np and approx_np has a different total_entry_num")
        return
    
    pc1_np = pc1_np[~np.isnan(pc1_np)]    
    approx_np = approx_np[~np.isnan(approx_np)]    
    valid_entry_num = len(pc1_np)
    
    if valid_entry_num != len(approx_np): 
        logging.info("pc1_np and approx_np has a different valid_entry_num")
        return

    pc1_pos_np = pc1_np > 0
    approx_pos_np = approx_np > 0
    pc1_pos_vs_approx_pos_np = pc1_pos_np == approx_pos_np 
    similar_num = list(pc1_pos_vs_approx_pos_np).count(True)
    similar_rate = similar_num / valid_entry_num

    return {
        "total_entry_num": total_entry_num,
        "valid_entry_num": valid_entry_num,
        "similar_num": similar_num,
        "similar_rate": similar_rate,
    }

def plot_comparison(pc1_np: np.ndarray, approx_np: np.ndarray, figsize: int=20, scatter: str | None = None, relative_magnitude: str | None = None):
    """
    Plot the comparison figure between PC1 and Approximated PC1-pattern we used in our paper.
    Note that for the plot of relative_magnitude, ``NaN`` value entries will be replaced with ``0`` in advance, 
    and both the PC1 and Approximated PC1-pattern will be Z-score normalized.   

    :param pc1_np: PC1 in NumPy format. 
    :type pc1_np: numpy.ndarray.

    :param approx_np: Approximated PC1-pattern in NumPy format. 
    :type approx_np: numpy.ndarray.

    :param scatter: (Optional) If the file path is specified, the scatter plot will be stored (e.g. ``scatter="./test/scatter.png"``).
    :type scatter: str

    :param relative_magnitude: (Optional) If the file path is specified, the relative_magnitude plot will be stored (e.g. ``relative_magnitude="./test/scatter.png"``).
    :type relative_magnitude: str
    """

    if os.path.dirname(scatter):
        os.makedirs(os.path.dirname(scatter), exist_ok=True)

    if os.path.dirname(relative_magnitude):
        os.makedirs(os.path.dirname(relative_magnitude), exist_ok=True)

    similarity_info = calc_similarity(pc1_np, approx_np)
    total_entry_num = similarity_info["total_entry_num"]
    valid_entry_num = similarity_info["valid_entry_num"]
    similar_num = similarity_info["similar_num"]
    similar_rate = similarity_info["similar_rate"]

    if scatter != None:
        plot_x_axis = [i + 1 for i in range(total_entry_num)]
        approx_dots = [1 if i > 0 else -1 if i < 0 else 0 for i in approx_np]
        pc1_colors_values = [2 if i > 0 else 0 if i < 0 else 1 for i in pc1_np]
        pc1_colors = ListedColormap(['r', 'g', 'b'])
        scatter_labels = ["PC1 < 0", "PC1 == NaN", "PC1 > 0"]

        plt.figure(figsize=(figsize, 6))
        plt.xticks(np.arange(0, total_entry_num, 50)) 
        scatter_config =  plt.scatter(plot_x_axis, approx_dots, c=pc1_colors_values, cmap=pc1_colors)
        plt.legend(handles=scatter_config.legend_elements()[0], labels=scatter_labels, fontsize="20", loc="center left")
        plt.title(f"total_entry_num: {total_entry_num}, valid_entry_num: {valid_entry_num}, similar_num = {similar_num}, similar_rate={np.round(similar_rate, 2)}", fontsize=20, loc="left")
        plt.savefig(scatter)
        plt.clf() 

    if relative_magnitude != None:
        # Fill NaN with 0 for plotting
        pc1_np[np.isnan(pc1_np)] = float(0) 
        approx_np[np.isnan(approx_np)] = float(0)

        # Z-score Normalization
        approx_np_norm = (approx_np - np.mean(approx_np)) / np.std(approx_np)
        pc1_np_norm = (pc1_np - np.mean(pc1_np)) / np.std(pc1_np)
        
        plt.figure(figsize=(figsize, 6))
        plt.xticks(np.arange(0, total_entry_num, 50)) 
        plt.plot(pc1_np_norm, c='r')
        plt.plot(approx_np_norm, c='b')
        plt.legend(["PC1", "approximated PC1-pattern"], fontsize="20", loc ="upper left")
        plt.title(f"total_entry_num: {total_entry_num}, valid_entry_num: {valid_entry_num}", fontsize=20, loc="left")
        plt.savefig(relative_magnitude)        
        plt.clf()
    
    plt.close('all')
    return