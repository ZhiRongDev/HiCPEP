#!/bin/bash
PY_FILE="src/common/create_approx.py"

### Rao 2014
DATA_PATH="${DOCKER_VOLUME_PATH}/data/Rao_2014/juicer_outputs"
OUTPUT_PATH="${DOCKER_VOLUME_PATH}/outputs/approx_PC1_pattern/Rao_2014"
mkdir -p "${OUTPUT_PATH}"

echo "$(date '+%Y-%m-%d %H:%M:%S') Rao_2014 create_approx start"

for RESOLUTION in "1000000" "100000" 
do
    CELL_LINE_NAMES=(\
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

    declare -A CELL_LINE_LINKS
    CELL_LINE_LINKS=(\
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

    for TYPE in "CxMax" "CxMin"
    do
        for CELL_LINE in "${CELL_LINE_NAMES[@]}"
        do
            echo "$(date '+%Y-%m-%d %H:%M:%S') [${RESOLUTION}] [${TYPE}] [${CELL_LINE}] start"

            HIC_PATH="${CELL_LINE_LINKS[$CELL_LINE]}"
            CHROM_LIST=$(python scripts/straw.py --hic_path "${HIC_PATH}")

            for CHROM in $CHROM_LIST
            do
                CHROM_ID="${CHROM/chr/""}" 
                CHROM_ID="${CHROM_ID/MT/"M"}" 

                PEARSON_FILE="${DATA_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons/pearson_chr${CHROM_ID}.txt"
                OUTPUT_FILE="${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/${TYPE}/approx_PC1_pattern_chr${CHROM_ID}.txt"
                python "${PY_FILE}" --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE --source "Rao_2014"
            done

            echo "$(date '+%Y-%m-%d %H:%M:%S') [${RESOLUTION}] [${TYPE}] [${CELL_LINE}] end"
        done
    done
done

echo "$(date '+%Y-%m-%d %H:%M:%S') Rao_2014 create_approx end"

### Lieberman 2009 
DATA_PATH="${DOCKER_VOLUME_PATH}/data/Lieberman_2009"
OUTPUT_PATH="${DOCKER_VOLUME_PATH}/outputs/approx_PC1_pattern/Lieberman_2009"
mkdir -p "${OUTPUT_PATH}"

echo "$(date '+%Y-%m-%d %H:%M:%S') Lieberman_2009 create_approx start"
for RESOLUTION in "1000000" "100000" 
do
    CELL_LINE_NAMES=(\
        "gm06690"\
        "k562" \
    )

    for TYPE in "CxMax" "CxMin"
    do
        for CELL_LINE in "${CELL_LINE_NAMES[@]}"
        do
            echo "$(date '+%Y-%m-%d %H:%M:%S') [${RESOLUTION}] [${TYPE}] [${CELL_LINE}] start"
            
            # There are some missing chromosomes in Lieberman's dataset.
            if [[ "${CELL_LINE}" == "gm06690" ]]
            then
                CHROM_LIST="X 22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1"
            else
                CHROM_LIST="22 21 20 19 18 17 16 15 14 13 12 11 10 9 8 7 6 5 4 3 2 1"
            fi

            for CHROM_ID in $CHROM_LIST
            do
                PEARSON_FILE="${DATA_PATH}/heatmaps/HIC_${CELL_LINE}_chr${CHROM_ID}_chr${CHROM_ID}_${RESOLUTION}_pearson.txt"
                OUTPUT_FILE="${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/${TYPE}/approx_PC1_pattern_chr${CHROM_ID}.txt"
                python "${PY_FILE}" --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE --source "Lieberman_2009"
            done

            echo "$(date '+%Y-%m-%d %H:%M:%S') [${RESOLUTION}] [${TYPE}] [${CELL_LINE}] end"
        done
    done
done

echo "$(date '+%Y-%m-%d %H:%M:%S') Lieberman_2009 create_approx end"