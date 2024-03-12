# Hi-C Pearson matrix's Approximated-PC1-pattern (HiCPAP)

HiCPAP is a python package used for creating the Approximated-PC1-pattern, which can be used to identify the A/B compartment.

## Features

## Prerequisites

## Installation

## Quick start

For instance:
>hicmaptools -in_map examples/fly_30k.binmap -in_bin examples/fly_30k.bins -bait examples/bait.bed -output temp.tsv

## Todo

* The flipped of approx_pc1 should be decided by GC content, but not the correlation coefficient between PC1.
https://github.com/vaquerizaslab/fanc/blob/55ec7376df1ec0cb6df51cfe908b15883a4045da/fanc/architecture/compartments.py#L155

* Biopython gc_fraction
https://biopython.org/docs/dev/api/Bio.SeqUtils.html#Bio.SeqUtils.gc_fraction

* Flip the approx_pc1 with the support of cooler, provide the 1Mb and 100Kb resolution GC_content tsv for both 
human and mouse cells.(The total_entry_num are all the same for gm12878, k562....., since they all use the same reference genome.)

* The reference genome use for rao_2014 are the b37(human) and mm9(mouse).

## References

How the Juicer implemented there eigenvector-calculation:
https://github.com/aidenlab/Juicebox/blob/12bc67454c460159c758606ff05eaecaf6601d1b/src/juicebox/data/MatrixZoomData.java#L708

* Fill the matrix with 0 first, and remove rows/columns with all 0 entries.
* Record the 0 (NaN) index position with numpy.diag, which should be 0.
* Take advantage of numpy.ix_() to get the submatrix.
