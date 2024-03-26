#!/bin/bash 

# This script is for creating the Pearson and PC1 ground truth files for the experiments of Rao, 2014.
# Note that the logs will be stored in `juicer_outputs.log`.

JUICER_TOOLS="${DATA_STORE}/juicer_tools/juicer_tools_1.22.01.jar" 
OUTPUT_PATH="${DATA_STORE}/data/rao_2014/juicer_outputs"
RAO_HIC_PATH="${DATA_STORE}/data/rao_2014/hic"

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
        "ch12-lx" \
    )

    ### Note that `-A` should be set (Associative Array: key-value), not `-a` (Index Array)
    declare -A CELL_LINE_LINKS
    CELL_LINE_LINKS=(\
        ["gm12878"]="${RAO_HIC_PATH}/GSE63525_GM12878_insitu_primary_replicate_combined_30.hic" \
        ["imr90"]="${RAO_HIC_PATH}/GSE63525_IMR90_combined_30.hic" \
        ["hmec"]="${RAO_HIC_PATH}/GSE63525_HMEC_combined_30.hic" \
        ["nhek"]="${RAO_HIC_PATH}/GSE63525_NHEK_combined_30.hic" \
        ["k562"]="${RAO_HIC_PATH}/GSE63525_K562_combined_30.hic" \
        ["kbm7"]="${RAO_HIC_PATH}/GSE63525_KBM7_combined_30.hic" \
        ["huvec"]="${RAO_HIC_PATH}/GSE63525_HUVEC_combined_30.hic" \
        ["ch12-lx"]="${RAO_HIC_PATH}/GSE63525_CH12_LX_combined_30.hic" \
    )

    for CELL_LINE in "${CELL_LINES[@]}"
    do
        LOG="${OUTPUT_PATH}/juicer_outputs.log"
        HIC_PATH="${CELL_LINE_LINKS[$CELL_LINE]}"

        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons"
        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector"
        touch "${LOG}"

        # There are only 19 chromosome in ch12-lx (Mouse), and the naming format is different with the other Human cell lines.
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
