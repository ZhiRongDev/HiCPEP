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

DATA_PATH="${DOCKER_VOLUME_PATH}/data"
OUTPUT_PATH="${DOCKER_VOLUME_PATH}/outputs"
JUICER_TOOLS_PATH="${DOCKER_VOLUME_PATH}/juicer_tools"

echo "${DATA_PATH}"
echo "${OUTPUT_PATH}"
echo "${JUICER_TOOLS_PATH}"

mkdir -p "${DATA_PATH}" "${OUTPUT_PATH}" "${JUICER_TOOLS_PATH}"

echo -e ">>> Download Juicer tools Version 1.22.01 \n"
wget "https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar" -P "${JUICER_TOOLS_PATH}" 
ln -s "${JUICER_TOOLS_PATH}/juicer_tools_1.22.01.jar" "${JUICER_TOOLS_PATH}/juicer_tools.jar" 

mkdir "${DATA_PATH}/Lieberman_2009/zip"

echo -e ">>> Download Lieberman 2009's GSE18199_binned_heatmaps.zip.gz \n"
mkdir -p "${DATA_PATH}/Lieberman_2009/heatmaps"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Fbinned%5Fheatmaps.zip.gz" -P "${DATA_PATH}/Lieberman_2009/zip"
gzip -d "${DATA_PATH}/Lieberman_2009/zip/GSE18199_binned_heatmaps.zip.gz"
unzip "${DATA_PATH}/Lieberman_2009/zip/GSE18199_binned_heatmaps.zip" -d "${DATA_PATH}/Lieberman_2009/heatmaps"

echo -e ">>> Download Lieberman 2009's GSE18199_eigenvector_files.zip.gz \n"
mkdir -p "${DATA_PATH}/Lieberman_2009/eigenvectors"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Feigenvector%5Ffiles.zip.gz" -P "${DATA_PATH}/Lieberman_2009/zip"
gzip -d "${DATA_PATH}/Lieberman_2009/zip/GSE18199_eigenvector_files.zip.gz"
unzip "${DATA_PATH}/Lieberman_2009/zip/GSE18199_eigenvector_files.zip" -d "${DATA_PATH}/Lieberman_2009/eigenvectors"

rm -rf "${DATA_PATH}/Lieberman_2009/zip" 