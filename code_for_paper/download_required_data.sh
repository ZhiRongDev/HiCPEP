#!/bin/bash

# This script is used for downloading the required data for the experiments of Lieberman, 2009 and Rao, 2014.
# The data for the experiments of Lieberman, 2009 are downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199.
# The data for the experiments of Rao, 2014 are downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525.
# The juicer_tools_1.22.01.jar is downloaded from https://github.com/aidenlab/juicer/wiki/Download.

DATA_PATH="${DATA_STORE}/data"
OUTPUT_PATH="${DATA_STORE}/outputs"
JUICER_TOOLS_PATH="${DATA_STORE}/juicer_tools"

mkdir -p "${DATA_PATH}" "${OUTPUT_PATH}" "${JUICER_TOOLS_PATH}"
mkdir -p "${DATA_PATH}/rao_2014/hic"
mkdir -p "${DATA_PATH}/lieberman_2009/zip"
mkdir -p "${DATA_PATH}/lieberman_2009/heatmaps"
mkdir -p "${DATA_PATH}/lieberman_2009/eigenvectors"

echo "(date '+%Y-%m-%d %H:%M:%S') Download required .hic files from GEO GSE63525"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%5F30.hic" -P "${DATA_PATH}/rao_2014/hic"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FIMR90%5Fcombined%5F30.hic" -P "${DATA_PATH}/rao_2014/hic"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FHMEC%5Fcombined%5F30.hic" -P "${DATA_PATH}/rao_2014/hic"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FNHEK%5Fcombined%5F30.hic" -P "${DATA_PATH}/rao_2014/hic"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FK562%5Fcombined%5F30.hic" -P "${DATA_PATH}/rao_2014/hic"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FKBM7%5Fcombined%5F30.hic" -P "${DATA_PATH}/rao_2014/hic"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FHUVEC%5Fcombined%5F30.hic" -P "${DATA_PATH}/rao_2014/hic"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FCH12%2DLX%5Fcombined%5F30.hic" -P "${DATA_PATH}/rao_2014/hic"

# Rename the .hic file to avoid error in the juicer_tools. 
mv "${DATA_PATH}/rao_2014/hic/GSE63525_GM12878_insitu_primary+replicate_combined_30.hic" "${DATA_PATH}/rao_2014/hic/GSE63525_GM12878_insitu_primary_replicate_combined_30.hic"
mv "${DATA_PATH}/rao_2014/hic/GSE63525_CH12-LX_combined_30.hic" "${DATA_PATH}/rao_2014/hic/GSE63525_CH12_LX_combined_30.hic"

echo "(date '+%Y-%m-%d %H:%M:%S') Download Juicer tools Version 1.22.01"
wget "https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar" -P "${JUICER_TOOLS_PATH}" 

echo "(date '+%Y-%m-%d %H:%M:%S') Download Lieberman 2009's GSE18199_binned_heatmaps.zip.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Fbinned%5Fheatmaps.zip.gz" -P "${DATA_PATH}/lieberman_2009/zip"
gzip -d "${DATA_PATH}/lieberman_2009/zip/GSE18199_binned_heatmaps.zip.gz"
unzip "${DATA_PATH}/lieberman_2009/zip/GSE18199_binned_heatmaps.zip" -d "${DATA_PATH}/lieberman_2009/heatmaps"

echo "(date '+%Y-%m-%d %H:%M:%S') Download Lieberman 2009's GSE18199_eigenvector_files.zip.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18199/suppl/GSE18199%5Feigenvector%5Ffiles.zip.gz" -P "${DATA_PATH}/lieberman_2009/zip"
gzip -d "${DATA_PATH}/lieberman_2009/zip/GSE18199_eigenvector_files.zip.gz"
unzip "${DATA_PATH}/lieberman_2009/zip/GSE18199_eigenvector_files.zip" -d "${DATA_PATH}/lieberman_2009/eigenvectors"

### Optional
rm -rf "${DATA_PATH}/lieberman_2009/zip" 