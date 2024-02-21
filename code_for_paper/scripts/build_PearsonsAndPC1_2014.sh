#!/bin/bash 

helpFunction()
{
    echo ""
    echo "Usage: $0 -p DOCKER_VOLUME_PATH"
    echo -e "\t-p Description of what is parameterP"
    exit 1 # Exit script after printing help
}

while getopts "p:" opt
do
    case "${opt}" in
        p ) DOCKER_VOLUME_PATH="${OPTARG}" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "${DOCKER_VOLUME_PATH}" ]
then
    echo "The -p parameter is empty.";
    helpFunction
fi

JUICER_TOOLS="${DOCKER_VOLUME_PATH}/juicer_tools/juicer_tools.jar" # Customize your path for the juicer_tools
OUTPUT_PATH="${DOCKER_VOLUME_PATH}/data/Rao_2014/juicer_outputs"

for RESOLUTION in "1000000" "100000"
do
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

    ### Note that `-A` should be set, not `-a` !
    declare -A cell_line_links
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

    for CELL_LINE in "${cell_line_names[@]}"
    do
        LOG="${OUTPUT_PATH}/juicer_outputs.log"
        HIC_PATH="${cell_line_links[$CELL_LINE]}"
        CHROM_LIST=$(python scripts/straw.py --hic_path "${HIC_PATH}")

        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons"
        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector"
        touch "${LOG}"

        for CHROM in $CHROM_LIST
        do
            CHROM_ID="${CHROM/chr/""}" 
            CHROM_ID="${CHROM_ID/MT/"M"}" 

            echo -e "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} Pearson start" 2>&1 | tee -a "${LOG}"
            java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" pearsons KR "${HIC_PATH}" "${CHROM}" BP "${RESOLUTION}" "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons/pearson_chr${CHROM_ID}.txt" -p 2>&1 | tee -a "${LOG}"
            echo -e "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} Pearson End \n" 2>&1 | tee -a "${LOG}"

            echo -e "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} PC1 start" 2>&1 | tee -a "${LOG}"
            java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" eigenvector KR "${HIC_PATH}" "${CHROM}" BP "${RESOLUTION}" "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector/pc1_chr${CHROM_ID}.txt" -p 2>&1 | tee -a "${LOG}"
            echo -e "$(date '+%Y-%m-%d %H:%M:%S') ${CELL_LINE} ${RESOLUTION} chr${CHROM_ID} PC1 End \n" 2>&1 | tee -a "${LOG}"
        done
    done
done
