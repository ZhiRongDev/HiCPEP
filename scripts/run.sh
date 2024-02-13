#!/bin/bash

EIGEN_PATH="/home/jordan990301/Projects/HiC-PC1_approx/data/Rao2014/juicer_outputs/GM12878/1Mb/data/pearsons/pearson_chr1.txt"
OUTPUT_FILE="/home/jordan990301/Projects/HiC-PC1_approx/outputs/test_scripts/approxPC1.txt"

python create_approxPattern.py --pearson $EIGEN_PATH --output $OUTPUT_FILE 