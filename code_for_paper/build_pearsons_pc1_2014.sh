#!/bin/bash 
JUICER_TOOLS="${DOCKER_VOLUME_PATH}/juicer_tools/juicer_tools.jar" 
OUTPUT_PATH="${DOCKER_VOLUME_PATH}/data/rao_2014/juicer_outputs"

for RESOLUTION in "1000000" "100000"
do
    CELL_LINES=(\
        "gm12878"\
        "imr90" \
        "hmec" \
        "nhek" \
        "k562" \
        "kbm7" \
        "huvec" \
        "hela" \
        "ch12-lx" \
    )

    ### Note that `-A` should be set (Associative Array: key-value), not `-a` (Index Array)
    declare -A CELL_LINE_LINKS
    CELL_LINE_LINKS=(\
        ["gm12878"]="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic" \
        ["imr90"]="https://hicfiles.s3.amazonaws.com/hiseq/imr90/in-situ/combined.hic" \
        ["hmec"]="https://hicfiles.s3.amazonaws.com/hiseq/hmec/in-situ/combined.hic" \
        ["nhek"]="https://hicfiles.s3.amazonaws.com/hiseq/nhek/in-situ/combined.hic" \
        ["k562"]="https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic" \
        ["kbm7"]="https://hicfiles.s3.amazonaws.com/hiseq/kbm7/in-situ/combined.hic" \
        ["huvec"]="https://hicfiles.s3.amazonaws.com/hiseq/huvec/in-situ/combined.hic" \
        ["hela"]="https://hicfiles.s3.amazonaws.com/hiseq/hela/in-situ/combined.hic" \
        ["ch12-lx"]="https://hicfiles.s3.amazonaws.com/hiseq/ch12-lx-b-lymphoblasts/in-situ/combined.hic" \
    )

    for CELL_LINE in "${CELL_LINES[@]}"
    do
        LOG="${OUTPUT_PATH}/juicer_outputs.log"
        HIC_PATH="${CELL_LINE_LINKS[$CELL_LINE]}"

        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons"
        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector"
        touch "${LOG}"

        # There are only 19 chromosome in ch12-lx, 
        if [[ "${CELL_LINE}" == "ch12-lx" ]]
        then
            CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY)
        else
            CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
        fi

        for CHROM in "${CHROMS[@]}"
        do
            CHROM_ID="${CHROM/chr/""}" 

            echo "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} Pearson start" 2>&1 | tee -a "${LOG}"
            java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" pearsons KR "${HIC_PATH}" "${CHROM}" BP "${RESOLUTION}" "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons/pearson_chr${CHROM_ID}.txt" -p 2>&1 | tee -a "${LOG}"
            echo "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} Pearson End" 2>&1 | tee -a "${LOG}"

            echo "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} PC1 start" 2>&1 | tee -a "${LOG}"
            java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" eigenvector KR "${HIC_PATH}" "${CHROM}" BP "${RESOLUTION}" "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector/pc1_chr${CHROM_ID}.txt" -p 2>&1 | tee -a "${LOG}"
            echo "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} PC1 End" 2>&1 | tee -a "${LOG}"
        done
    done
done
