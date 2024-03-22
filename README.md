# Hi-C Pearson matrix's Approximated-PC1-pattern (HiCPAP)

HiCPAP is a python package used for creating the Approximated-PC1-pattern, which can be used to identify the A/B compartment.

## Features

## Prerequisites

To execute this example successfully, Please `cd` to the root directory:

1. `python -m pip install -e .` (This will install hicpap and all the dependency packages through the setup.py).
2. `pip install -r requirements.txt`.

## Installation

## Quick start

For instance:
>hicmaptools -in_map examples/fly_30k.binmap -in_bin examples/fly_30k.bins -bait examples/bait.bed -output temp.tsv

## References

How the Juicer implemented there eigenvector-calculation:
https://github.com/aidenlab/Juicebox/blob/12bc67454c460159c758606ff05eaecaf6601d1b/src/juicebox/data/MatrixZoomData.java#L708

* Fill the matrix with 0 first, and remove rows/columns with all 0 entries.
* Record the 0 (NaN) index position with numpy.diag, which should be 0.
* Take advantage of numpy.ix_() to get the submatrix.
