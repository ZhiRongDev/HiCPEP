#!/bin/bash
DATA_PATH="${DOCKER_VOLUME_PATH}/data"
OUTPUT_PATH="${DOCKER_VOLUME_PATH}/outputs"
JUICER_TOOLS_PATH="${DOCKER_VOLUME_PATH}/juicer_tools"

mkdir -p "${DATA_PATH}" "${OUTPUT_PATH}" "${JUICER_TOOLS_PATH}"

echo "(date '+%Y-%m-%d %H:%M:%S') Download Juicer tools Version 1.22.01"
wget "https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar" -P "${JUICER_TOOLS_PATH}" 
ln -s "${JUICER_TOOLS_PATH}/juicer_tools_1.22.01.jar" "${JUICER_TOOLS_PATH}/juicer_tools.jar" 

mkdir -p "${DATA_PATH}/lieberman_2009/zip"

echo "(date '+%Y-%m-%d %H:%M:%S') Download Lieberman 2009's GSE18199_binned_heatmaps.zip.gz"
mkdir -p "${DATA_PATH}/lieberman_2009/heatmaps"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Fbinned%5Fheatmaps.zip.gz" -P "${DATA_PATH}/lieberman_2009/zip"
gzip -d "${DATA_PATH}/lieberman_2009/zip/GSE18199_binned_heatmaps.zip.gz"
unzip "${DATA_PATH}/lieberman_2009/zip/GSE18199_binned_heatmaps.zip" -d "${DATA_PATH}/lieberman_2009/heatmaps"

echo "(date '+%Y-%m-%d %H:%M:%S') Download Lieberman 2009's GSE18199_eigenvector_files.zip.gz"
mkdir -p "${DATA_PATH}/lieberman_2009/eigenvectors"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Feigenvector%5Ffiles.zip.gz" -P "${DATA_PATH}/lieberman_2009/zip"
gzip -d "${DATA_PATH}/lieberman_2009/zip/GSE18199_eigenvector_files.zip.gz"
unzip "${DATA_PATH}/lieberman_2009/zip/GSE18199_eigenvector_files.zip" -d "${DATA_PATH}/lieberman_2009/eigenvectors"

rm -rf "${DATA_PATH}/lieberman_2009/zip" 