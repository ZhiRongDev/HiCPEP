import numpy as np
import pandas as pd
import logging
logging.basicConfig(format='%(message)s', level=logging.INFO)

def read_pearson(pearson: str, format="juicer") -> np.ndarray:
    """
    :param pearson: The text file path of `juicer_tools <https://github.com/aidenlab/juicer/wiki/Pearsons>`_  created Pearson matrix.
    :type pearson: str
    :return All `0` rows/columns removed pearson_np.
    :rtype: numpy.ndarray
    """
    if format == "juicer":
        pearson_df = pd.read_table(pearson, header=None, sep="\s+")
        pearson_np = pearson_df.values # Turn into numpy.ndarray
        pearson_np = pearson_np.astype('float64')
    elif format == "aiden_2009":
        pearson_df = pd.read_table(pearson, index_col=0, header=1, sep="\s+")
        pearson_np = pearson_df.values # Turn into numpy.ndarray
        pearson_np = pearson_np.astype('float64')
        diag = np.diag(pearson_np)
        diag_valid = diag != 0
        ixgrid = np.ix_(diag_valid, diag_valid) # Extract the valid sub-matrix.

        # Fill the all-zero rows/columns with `NaN`.
        length = len(pearson_np)
        tmp = np.full((length, length), np.nan) 
        tmp[ixgrid] = pearson_np[ixgrid] 
        pearson_np = tmp 

    return pearson_np

def flip_tracks(track1_np: np.ndarray, track2_np: np.ndarray):
    if len(track1_np) != len(track2_np):
        logging.info("The length of track1_np is different with track2_np")
        logging.info(f"Length of track1_np: {len(track1_np)}")
        logging.info(f"Length of track2_np: {len(track2_np)}")

    if np.corrcoef(track1_np[~np.isnan(track1_np)], track2_np[~np.isnan(track2_np)])[0][1] < 0:
        track2_np = -track2_np

    return track1_np, track2_np

def flip_track_gc(track_np: np.ndarray, gc_np: np.ndarray) -> np.ndarray:
    """
    Note that the GC content information files created by the UCSC Genome Browser tool missed the last bin of each chromosome.
    """
    if len(track_np[:-1]) != len(gc_np):
        logging.info("The length of track_np[:-1] is different with gc_np")
        logging.info(f"Length of track_np[:-1]: {len(track_np[:-1])}")
        logging.info(f"Length of gc_np: {len(gc_np)}")

    if np.nanmean(gc_np[track_np[:-1] > 0]) < np.nanmean(gc_np[track_np[:-1] < 0]):
        track_np = -track_np

    return track_np

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
    
    pearson_np = pearson_np.astype('float64')
    diag = np.diag(pearson_np)
    diag_valid = ~np.isnan(diag)
    ixgrid = np.ix_(diag_valid, diag_valid) # Record the position of the valid sub-matrix.

    total_entry_num = len(pearson_np)
    pearson_np = pearson_np[ixgrid] # Remove NaN 
    pearson_np -= pearson_np.mean(axis=1, keepdims=True) # Center the pearson_np

    # PCA
    n = valid_entry_num = len(pearson_np)
    y = pearson_np.T / np.sqrt(n)
    U, S, Vh = np.linalg.svd(y, full_matrices=True)
    eigenvalues = S * S
    sum_eigenvalues = np.sum(eigenvalues)
    explained_variances_ratio = eigenvalues / sum_eigenvalues

    # Place back the valid entries to it's origin position in the chromosome.  
    tmp = np.full((len(diag_valid), len(diag_valid)), np.nan) 
    tmp[ixgrid] = Vh
    Vh = tmp

    # Return Principal components(Vector) and the explained variances of PC1(Vector).
    return Vh[diag_valid], explained_variances_ratio, total_entry_num, valid_entry_num