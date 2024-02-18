#!/bin/bash 

# Usage:
#    ./build_PearsonsAndEigenvector1_humanCells.sh
#
# This script is for building the Pearson correlation matrix and Eigenvector1 with the support of juicer_tools_1.22.01.jar (Or the other version of juicer_tools).
# Note that you should set up your own RESOLUTION_NAME, RESOLUTION, JUICER_TOOLS_PATH, OUTPUT_PATH.
# 
# Here we list the example outputs' directory structure:
# .
# └── $OUTPUT_PATH
#     ├── CH12-LX
#     │   ├── $RESOLUTION_NAME
#     │   │   ├── $RESOLUTION_NAME.log
#     │   │   └── data
#     │   │       ├── eigenvector
#     │   │       │   └── pc1_chr*.txt
#     │   │       └── pearsons
#     │   │           └── pearson_chr*.txt

RESOLUTION_NAME="100Kb" # Set the matrix resolution (This is used for the directory name in the $OUTPUT_PATH)
RESOLUTION=100000 # Set the matrix resolution
JUICER_TOOLS_PATH="/content/drive/MyDrive/HiC_experiments/juicer_tools/juicer_tools.jar" # Customize your path for the juicer_tools

cell_line_names=(\
    "CH12-LX"\
)

declare -a cell_line_links
cell_line_links=(\
    ["CH12-LX"]="https://hicfiles.s3.amazonaws.com/hiseq/ch12-lx-b-lymphoblasts/in-situ/combined.hic" \
)

for CELL_LINE in "${cell_line_names[@]}"
do
    HIC_PATH="${cell_line_links[$CELL_LINE]}"
    OUTPUT_PATH="/content/drive/MyDrive/HiC_experiments/juicer_outputs/$CELL_LINE/$RESOLUTION_NAME/" # Customize your path for the outputs

    mkdir -p "$OUTPUT_PATH/data/pearsons"
    mkdir -p "$OUTPUT_PATH/data/eigenvector"
    touch $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log

    for CHROM in chrY chrX chr19 chr18 chr17 chr16 chr15 chr14 chr13 chr12 chr11 chr10 chr9 chr8 chr7 chr6 chr5 chr4 chr3 chr2 chr1
    do

        echo "====== $CELL_LINE $RESOLUTION_NAME $CHROM Pearson start ======" 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        java -Xmx13004m -Xms13004m -jar $JUICER_TOOLS_PATH pearsons KR $HIC_PATH $CHROM BP $RESOLUTION $OUTPUT_PATH/data/pearsons/pearson_$CHROM.txt -p 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        echo "====== Pearson End ======" 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        echo " " 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log

        echo "====== $CELL_LINE $RESOLUTION_NAME $CHROM PC1 start ======" 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        java -Xmx13004m -Xms13004m -jar $JUICER_TOOLS_PATH eigenvector KR $HIC_PATH $CHROM BP $RESOLUTION $OUTPUT_PATH/data/eigenvector/pc1_$CHROM.txt -p 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        echo "====== PC1 End ======" 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        echo " " 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log

    done

done