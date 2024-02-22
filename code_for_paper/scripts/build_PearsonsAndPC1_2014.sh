#!/bin/bash 

JUICER_TOOLS="${DOCKER_VOLUME_PATH}/juicer_tools/juicer_tools.jar" # Customize your path for the juicer_tools
OUTPUT_PATH="${DOCKER_VOLUME_PATH}/data/Rao_2014/juicer_outputs"

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

    ### Note that `-A` should be set (Associative Array: key-value), not `-a` (Index Array)
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

    for CELL_LINE in "${CELL_LINE_NAMES[@]}"
    do
        LOG="${OUTPUT_PATH}/juicer_outputs.log"
        HIC_PATH="${CELL_LINE_LINKS[$CELL_LINE]}"
        CHROM_LIST=$(python scripts/straw.py --hic_path "${HIC_PATH}")

        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons"
        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector"
        touch "${LOG}"

        for CHROM in $CHROM_LIST
        do
            CHROM_ID="${CHROM/chr/""}" 
            CHROM_ID="${CHROM_ID/MT/"M"}" 

            echo "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} Pearson start" 2>&1 | tee -a "${LOG}"
            java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" pearsons KR "${HIC_PATH}" "${CHROM}" BP "${RESOLUTION}" "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons/pearson_chr${CHROM_ID}.txt" -p 2>&1 | tee -a "${LOG}"
            echo "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} Pearson End" 2>&1 | tee -a "${LOG}"

            echo "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} PC1 start" 2>&1 | tee -a "${LOG}"
            java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" eigenvector KR "${HIC_PATH}" "${CHROM}" BP "${RESOLUTION}" "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector/pc1_chr${CHROM_ID}.txt" -p 2>&1 | tee -a "${LOG}"
            echo "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} PC1 End" 2>&1 | tee -a "${LOG}"
        done
    done
done
