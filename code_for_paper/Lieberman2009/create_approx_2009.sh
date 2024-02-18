#!/bin/bash

# mkdir -p "/home/jordan990301/Projects/HiC-PC1_approx/outputs/approx_PC1_pattern/"

LIEBERMAN_DATA_PATH="/home/jordan990301/Projects/HiC-PC1_approx/data/Lieberman_Aiden-datasets/"
OUTPUT_PATH="/home/jordan990301/Projects/HiC-PC1_approx/outputs/approx_PC1_pattern"
PY_FILE="/home/jordan990301/Projects/HiC-PC1_approx/code_for_paper/Lieberman2009/create_approx_2009.py"

for RESOLUTION in "1000000" "100000"
do
    if [ $RESOLUTION = "1000000" ]
    then
        RESOLUTION_NAME="1Mb"
    elif [ $RESOLUTION = "100000" ]
    then
        RESOLUTION_NAME="100Kb"
    fi

    echo "==== ${RESOLUTION_NAME} ====="

    for TYPE in "CxMax" "CxMin"
    do
        echo "==== ${TYPE} ====="
        echo "GM06690 start"
        # Lieberman's dataset doesn't provide chromY
        for CHROM in X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
        do
            PEARSON_FILE="${LIEBERMAN_DATA_PATH}/heatmaps/HIC_gm06690_chr${CHROM}_chr${CHROM}_${RESOLUTION}_pearson.txt"
            OUTPUT_FILE="${OUTPUT_PATH}/GM06690/${RESOLUTION_NAME}/${TYPE}/approx_PC1_pattern_chr${CHROM}.txt"
            python $PY_FILE --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE
        done
        echo "GM06690 end"
        echo " "
    done
done