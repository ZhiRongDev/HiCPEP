import os
import random
import math
import numpy as np
import pandas as pd
from copy import deepcopy
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import logging
logging.basicConfig(format='%(message)s', level=logging.INFO)

def read_pearson(pearson: str) -> np.ndarray:
    """
    Read a ``.txt`` file of the intra-chromosomal Hi-C Pearson matrix created by `juicer_tools <https://github.com/aidenlab/juicer/wiki/Pearsons>`_ 
    and return the Pearson matrix in NumPy format.

    :param pearson: Path of the juicer_tools created intra-chromosomal Hi-C Pearson matrix ``.txt`` file.
    :type pearson: ``str``

    :return: Intra-chromosomal Hi-C Pearson matrix in NumPy format.
    :rtype: ``numpy.ndarray``
    """

    pearson_df = pd.read_table(pearson, header=None, sep="\s+")
    pearson_np = pearson_df.values # Turn into numpy.ndarray.
    pearson_np = pearson_np.astype('float64')
    return pearson_np

def create_est(pearson_np: np.ndarray, output: str | None = None, method: str="cxmax", sampling_proportion: float=1.0) -> np.ndarray:
    """
    Create the Estimated PC1-pattern of the given Hi-C Pearson matrix. 
    The calculation is only performed on the valid sub-matrix. 
    (We exclude the rows and columns which the corresponding diagonal value is ``NaN``, implies that these rows and columns are all ``NaN``.    
    However these all ``NaN`` rows or columns will not be removed in the Estimated PC1-pattern returned)

    :param pearson_np: Hi-C Pearson matrix in NumPy format.
    :type pearson_np: ``numpy.ndarray``

    :param output: (Optional) If the file path is specified, the Estimated PC1-pattern will be stored (e.g. ``output="./test/est_pc1.txt"``).
    :type output: ``str``

    :param method: ``cxmax`` or ``cxmin``.
    :type method: ``str``

    :param sampling_proportion: If this parameter is specified (e.g. 0.1), than the function will randomly sample the given percentage of rows in the Pearson matrix to create a partial covariance matrix, 
                                and select the ``cxmax`` in this partial covariance matrix as the Estimated PC1-pattern.
    :type sampling_proportion: ``float``

    :return: Estimated PC1-pattern in NumPy format.
    :rtype: ``numpy.ndarray``
    """

    if len(pearson_np) != len(pearson_np[0]): 
        logging.info("Pearson matrix given has a different number of rows and columns")
        return
    
    # Extract the valid sub-matrix and record the position of NaN entries (The calculation will only be performed on the not-NaN entries).
    pearson_np = pearson_np.astype('float64')
    diag = np.diag(pearson_np)
    diag_valid = ~np.isnan(diag)
    ixgrid = np.ix_(diag_valid, diag_valid) # Record the position of the valid sub-matrix.
    pearson_np = pearson_np[ixgrid] 
    has_nan = np.isnan(pearson_np).any()

    if has_nan:
        print("NaN entries still exist in the Pearson matrix.")
        return

    # Center the pearson_np
    pearson_np -= pearson_np.mean(axis=1, keepdims=True)

    # Core idea, note that the covariance matrix is symmetric.
    if isinstance(sampling_proportion, float) and sampling_proportion > 0 and sampling_proportion <= 1.0:
        n = len(pearson_np)
        sample_indexes = random.sample(list(range(n)), math.floor(n * float(sampling_proportion)))

        if len(sample_indexes) < 1:
            sample_indexes = [i for i in range(n)]

        partial_cov_np = pearson_np @ pearson_np[sample_indexes].T / n
        partial_cov_abs_sum = [(index, np.sum(np.abs(row))) for index, row in enumerate(partial_cov_np.T)] 
        sorted_partial_cov_abs_sum = sorted(partial_cov_abs_sum, key=lambda x: x[1], reverse=True) # Sorted from the maximum to the minimum 

        if method == "cxmax":
            sorted_index = 0
            est_np = partial_cov_np.T[sorted_partial_cov_abs_sum[sorted_index][0]]
        elif method == "cxmin":
            sorted_index = -1
            est_np = partial_cov_np.T[sorted_partial_cov_abs_sum[sorted_index][0]]

    # Place back the valid entries to it's origin position in the chromosome.
    tmp = np.full(len(diag_valid), np.nan)
    tmp[diag_valid] = est_np
    est_np = tmp

    if output == None:
        return est_np

    if os.path.dirname(output):
        os.makedirs(os.path.dirname(output), exist_ok=True)

    with open(output, 'w') as f:
        tmp_str = ''

        for i in est_np:
            tmp_str += f"{str(i)}\n"

        tmp_str = tmp_str[:-1]
        f.write(tmp_str)

    return est_np

def calc_similarity(track1_np: np.ndarray, track2_np: np.ndarray):
    """
    Compare the similarity information between the given track1 and track2, 
    The ``similar_rate`` is defined as the proportion of the entries in track1 (e.g. ``pc1_np``) 
    that have a same positive/negative sign as the track2 (e.g. Estimated PC1-pattern) entries compared.

    Note that the ``NaN`` value entries will be excluded in advance.

    :param track1_np: PC1 or Estimated PC1-pattern in ``numpy.ndarray`` format. 
    :type track1_np: ``numpy.ndarray``

    :param track2_np: PC1 or Estimated PC1-pattern in ``numpy.ndarray`` format. 
    :type track2_np: ``numpy.ndarray``

    :return:

    .. code::

        {
            total_entry_num (int): Entry numbers including ``NaN``.
            valid_entry_num (int): Entry numbers excluding ``NaN``.
            similar_num (int): Number of entries that the track1_np and track2_np have the same positive or negative sign.
            similar_rate (float): similar_num divide by valid_entry_num.
        }

    :rtype: ``dict``
    """

    track1_np = deepcopy(track1_np)
    track2_np = deepcopy(track2_np)
    total_entry_num = len(track1_np)
    
    if total_entry_num != len(track2_np): 
        logging.info("track1_np and track2_np has a different total_entry_num")
        return
    
    track1_np = track1_np[~np.isnan(track1_np)]
    track2_np = track2_np[~np.isnan(track2_np)]
    valid_entry_num = len(track1_np)
    
    if valid_entry_num != len(track2_np): 
        logging.info("track1_np and track2_np has a different valid_entry_num")
        return

    track1_pos_np = track1_np > 0
    track2_pos_np = track2_np > 0
    track1_pos_vs_track2_pos_np = track1_pos_np == track2_pos_np 
    similar_num = list(track1_pos_vs_track2_pos_np).count(True)
    similar_rate = float(similar_num / valid_entry_num)

    return {
        "total_entry_num": total_entry_num,
        "valid_entry_num": valid_entry_num,
        "similar_num": similar_num,
        "similar_rate": similar_rate,
    }

### If there's no `NaN` value, then the legend will have some mistakes. 
def plot_comparison(pc1_np: np.ndarray, est_np: np.ndarray, figsize: int=20, scatter: str | None = None, relative_magnitude: str | None = None, xticks: int=50):
    """
    Plot the scatter or relative-magnitude comparison figure between the PC1 and Estimated PC1-pattern. 
    Please specified at least one of the figure storing path among the scatter plot or the relative_magnitude plot.

    Note that for the plot of relative_magnitude, all the ``NaN`` value entries will be replaced with ``0`` in advance, 
    and both the PC1 and Estimated PC1-pattern will be Z-score normalized.

    :param pc1_np: PC1 in ``numpy.ndarray`` format. 
    :type pc1_np: ``numpy.ndarray``

    :param est_np: Estimated PC1-pattern in ``numpy.ndarray`` format. 
    :type est_np: ``numpy.ndarray``

    :param figsize: Scaling the figure size. 
    :type figsize: ``int``

    :param scatter: (Optional) If the file path is specified, the scatter plot will be stored (e.g. ``scatter="./test/scatter.png"``).
    :type scatter: ``str``

    :param relative_magnitude: (Optional) If the file path is specified, the relative_magnitude plot will be stored (e.g. ``relative_magnitude="./test/scatter.png"``).
    :type relative_magnitude: ``str``
    """

    pc1_np = deepcopy(pc1_np)
    est_np = deepcopy(est_np)

    if os.path.dirname(scatter):
        os.makedirs(os.path.dirname(scatter), exist_ok=True)

    if os.path.dirname(relative_magnitude):
        os.makedirs(os.path.dirname(relative_magnitude), exist_ok=True)

    similarity_info = calc_similarity(track1_np=pc1_np, track2_np=est_np)
    total_entry_num = similarity_info["total_entry_num"]
    valid_entry_num = similarity_info["valid_entry_num"]
    similar_num = similarity_info["similar_num"]
    similar_rate = similarity_info["similar_rate"]

    if scatter != None:
        plot_x_axis = [i + 1 for i in range(total_entry_num)]
        est_dots = [1 if i > 0 else -1 if i < 0 else 0 for i in est_np]
        pc1_colors_values = ['b' if i > 0 else 'r' if i < 0 else 'g' for i in pc1_np]

        # https://matplotlib.org/stable/users/explain/axes/legend_guide.html
        red_patch = mpatches.Patch(color='r', label='PC1 < 0')
        green_patch = mpatches.Patch(color='g', label='PC1 == NaN')
        blue_patch = mpatches.Patch(color='b', label='PC1 > 0')

        plt.figure(figsize=(figsize, 6))
        plt.xticks(np.arange(0, total_entry_num, xticks)) 
        plt.scatter(plot_x_axis, est_dots, c=pc1_colors_values)
        plt.legend(handles=[red_patch, green_patch, blue_patch], fontsize="20", loc="center left")
        plt.title(f"total_entry_num: {total_entry_num}, valid_entry_num: {valid_entry_num}, similar_num = {similar_num}, similar_rate={np.round(similar_rate, 2)}", fontsize=20, loc="left")
        plt.savefig(scatter)
        plt.clf() 

    if relative_magnitude != None:
        # Fill NaN with 0 for plotting
        pc1_np[np.isnan(pc1_np)] = float(0) 
        est_np[np.isnan(est_np)] = float(0)

        # Z-score Normalization
        est_np_norm = (est_np - np.mean(est_np)) / np.std(est_np)
        pc1_np_norm = (pc1_np - np.mean(pc1_np)) / np.std(pc1_np)
        
        plt.figure(figsize=(figsize, 6))
        plt.xticks(np.arange(0, total_entry_num, xticks)) 
        plt.plot(pc1_np_norm, c='r')
        plt.plot(est_np_norm, c='b')
        plt.legend(["PC1", "Estimated PC1-pattern"], fontsize="20", loc ="upper left")
        plt.title(f"total_entry_num: {total_entry_num}, valid_entry_num: {valid_entry_num}", fontsize=20, loc="left")
        plt.savefig(relative_magnitude)        
        plt.clf()

    if scatter is None and relative_magnitude is None:
        print("Please specified at least one of the figure storing path among the scatter plot or the relative_magnitude plot.")
    
    plt.close('all')
    return