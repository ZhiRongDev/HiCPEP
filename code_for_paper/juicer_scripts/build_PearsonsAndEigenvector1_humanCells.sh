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
#     ├── GM12878
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
    "GM12878"\
    "IMR90" \
    "HMEC" \
    "NHEK" \
    "K562" \
    "KBM7" \
    "HUVEC" \
    "HeLa" \
)

declare -a cell_line_links
cell_line_links=(\
    ["GM12878"]="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic" \
    ["IMR90"]="https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic" \
    ["HMEC"]="https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic" \
    ["NHEK"]="https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic" \
    ["K562"]="https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic" \
    ["KBM7"]="https://hicfiles.s3.amazonaws.com/hiseq/kbm7/in-situ/combined.hic" \
    ["HUVEC"]="https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic" \
    ["HeLa"]="https://hicfiles.s3.amazonaws.com/hiseq/hela/in-situ/combined.hic" \
)

for CELL_LINE in "${cell_line_names[@]}"
do
    HIC_PATH="${cell_line_links[$CELL_LINE]}"
    OUTPUT_PATH="/content/drive/MyDrive/HiC_experiments/juicer_outputs/$CELL_LINE/$RESOLUTION_NAME/" # Customize your path for the outputs

    mkdir -p "$OUTPUT_PATH/data/pearsons"
    mkdir -p "$OUTPUT_PATH/data/eigenvector"
    touch $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log

    for CHROM in Y X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
    do

        echo "====== $CELL_LINE $RESOLUTION_NAME Chr$CHROM Pearson start ======" 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        java -Xmx13004m -Xms13004m -jar $JUICER_TOOLS_PATH pearsons KR $HIC_PATH $CHROM BP $RESOLUTION $OUTPUT_PATH/data/pearsons/pearson_chr$CHROM.txt -p 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        echo "====== Pearson End ======" 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        echo " " 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log

        echo "====== $CELL_LINE $RESOLUTION_NAME Chr$CHROM PC1 start ======" 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        java -Xmx13004m -Xms13004m -jar $JUICER_TOOLS_PATH eigenvector KR $HIC_PATH $CHROM BP $RESOLUTION $OUTPUT_PATH/data/eigenvector/pc1_chr$CHROM.txt -p 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        echo "====== PC1 End ======" 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log
        echo " " 2>&1 | tee -a $OUTPUT_PATH/$CELL_LINE_$RESOLUTION_NAME.log

    done

done