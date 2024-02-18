#!/bin/bash

mkdir -p "/home/jordan990301/Projects/HiC-PC1_approx/outputs/approx_PC1_pattern/"

JUICER_OUTPUTS_PATH="/home/jordan990301/Projects/HiC-PC1_approx/data/Rao2014/juicer_outputs"
OUTPUT_PATH="/home/jordan990301/Projects/HiC-PC1_approx/outputs/approx_PC1_pattern"
PY_FILE="/home/jordan990301/Projects/HiC-PC1_approx/scripts/create_approx.py"

for RESOLUTION_NAME in "1Mb" "100Kb" 
do
    echo "==== ${RESOLUTION_NAME} ====="
    ## Human cells (has chrom 1-22, x, y)
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

    for TYPE in "CxMax" "CxMin"
    do
        echo "==== ${TYPE} ====="
        for CELL_LINE in "${cell_line_names[@]}"
        do
            echo "${CELL_LINE} start"
            for CHROM in Y X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
            do
                PEARSON_FILE="${JUICER_OUTPUTS_PATH}/${CELL_LINE}/${RESOLUTION_NAME}/data/pearsons/pearson_chr${CHROM}.txt"
                OUTPUT_FILE="${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION_NAME}/${TYPE}/approx_PC1_pattern_chr${CHROM}.txt"
                python $PY_FILE --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE
            done
            echo "${CELL_LINE} end"
            echo " "
        done

        ## Mouse cells (has chrom 1-19, x, y)
        echo "CH12-LX start"
        for CHROM in Y X 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
        do
            PEARSON_FILE="${JUICER_OUTPUTS_PATH}/CH12-LX/${RESOLUTION_NAME}/data/pearsons/pearson_chr${CHROM}.txt"
            OUTPUT_FILE="${OUTPUT_PATH}/CH12-LX/${RESOLUTION_NAME}/${TYPE}/approx_PC1_pattern_chr${CHROM}.txt"
            python $PY_FILE --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE
        done
        echo "CH12-LX end"
        echo " "
    done
done