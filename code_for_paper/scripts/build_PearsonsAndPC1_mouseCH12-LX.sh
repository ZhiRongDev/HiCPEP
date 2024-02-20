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
        "CH12-LX"\
    )

    declare -a cell_line_links
    cell_line_links=(\
        ["CH12-LX"]="https://hicfiles.s3.amazonaws.com/hiseq/ch12-lx-b-lymphoblasts/in-situ/combined.hic" \
    )

    for CELL_LINE in "${cell_line_names[@]}"
    do
        HIC_PATH="${cell_line_links[$CELL_LINE]}"
        LOG="${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/${RESOLUTION}.log"

        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons"
        mkdir -p "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector"
        touch "${LOG}"

        for CHROM in chrY chrX chr19 chr18 chr17 chr16 chr15 chr14 chr13 chr12 chr11 chr10 chr9 chr8 chr7 chr6 chr5 chr4 chr3 chr2 chr1
        do
            echo -e ">>> ${CELL_LINE} ${RESOLUTION} chr${CHROM} Pearson start" 2>&1 | tee -a "${LOG}"
            java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" pearsons KR "${HIC_PATH}" "${CHROM}" BP "${RESOLUTION}" "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons/pearson_chr${CHROM}.txt" -p 2>&1 | tee -a "${LOG}"
            echo -e ">>> Pearson End \n" 2>&1 | tee -a "${LOG}"

            echo -e ">>> ${CELL_LINE} ${RESOLUTION} ${CHROM} PC1 start" 2>&1 | tee -a "${LOG}"
            java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" eigenvector KR "${HIC_PATH}" "${CHROM}" BP "${RESOLUTION}" "${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/eigenvector/pc1_chr${CHROM}.txt" -p 2>&1 | tee -a "${LOG}"
            echo -e ">>> PC1 End \n" 2>&1 | tee -a "${LOG}"
        done

    done
done
