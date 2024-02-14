#!/bin/bash

# PEARSON_PATH="/home/jordan990301/Projects/HiC-PC1_approx/data/Rao2014/juicer_outputs/GM12878/1Mb/data/pearsons/pearson_chr1.txt"
# OUTPUT_FILE="/home/jordan990301/Projects/HiC-PC1_approx/outputs/test_scripts/approxPC1.txt"
# python create_approxPattern.py --pearson $PEARSON_PATH --output $OUTPUT_FILE 

EIGEN_PATH="/home/jordan990301/Projects/HiC-PC1_approx/data/Rao2014/juicer_outputs/GM12878/1Mb/data/eigenvector/pc1_chr1.txt"
APPROX_PATH="/home/jordan990301/Projects/HiC-PC1_approx/outputs/test_scripts/approxPC1.txt"
LINE_PLOT="/home/jordan990301/Projects/HiC-PC1_approx/outputs/test_scripts/"
SCATTER_PLOT="/home/jordan990301/Projects/HiC-PC1_approx/outputs/test_scripts/"

python plot.py --Juicer_PC1 $EIGEN_PATH --approx $APPROX_PATH --relative_magnitude $LINE_PLOT --scatter $SCATTER_PLOT 