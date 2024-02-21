
#!/bin/bash

JUICER_OUTPUTS_PATH="/home/jordan990301/Projects/HiC-PC1_approx/data/Rao2014/juicer_outputs"
OUTPUT_PATH="/home/jordan990301/Projects/HiC-PC1_approx/outputs/approx_PC1_pattern"
PY_FILE="/home/jordan990301/Projects/HiC-PC1_approx/scripts/create_approx.py"

mkdir -p "${OUTPUT_PATH}"

for RESOLUTION in "1000000" "100000" 
do
    echo "==== ${RESOLUTION} ====="
    cell_line_names=(\
        "GM12878"\
        "IMR90" \
        "HMEC" \
        "NHEK" \
        "K562" \
        "KBM7" \
        "HUVEC" \
        "HeLa" \
        "CH12-LX" \
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
        ["CH12-LX"]="https://hicfiles.s3.amazonaws.com/hiseq/ch12-lx-b-lymphoblasts/in-situ/combined.hic" \
    )
###

    for TYPE in "CxMax" "CxMin"
    do
        echo "==== ${TYPE} ====="
        for CELL_LINE in "${cell_line_names[@]}"
        do
            echo "${CELL_LINE} start"
            for CHROM in Y X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
            do
                PEARSON_FILE="${JUICER_OUTPUTS_PATH}/${CELL_LINE}/${RESOLUTION}/data/pearsons/pearson_chr${CHROM}.txt"
                OUTPUT_FILE="${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/${TYPE}/approx_PC1_pattern_chr${CHROM}.txt"
                python $PY_FILE --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE
            done
            echo "${CELL_LINE} end"
            echo " "
        done

        ## Mouse cells (has chrom 1-19, x, y)
        echo "CH12-LX start"
        for CHROM in Y X 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1
        do
            PEARSON_FILE="${JUICER_OUTPUTS_PATH}/CH12-LX/${RESOLUTION}/data/pearsons/pearson_chr${CHROM}.txt"
            OUTPUT_FILE="${OUTPUT_PATH}/CH12-LX/${RESOLUTION}/${TYPE}/approx_PC1_pattern_chr${CHROM}.txt"
            python $PY_FILE --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE
        done
        echo "CH12-LX end"
        echo " "
    done
done