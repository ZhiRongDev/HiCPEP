#!/bin/bash
JUICER_TOOLS="/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpep/data_store/juicer_tools/juicer_tools_1.22.01.jar" 
HIC_PATH="/media/jordan990301/Samsung_T5/HiC_Datasets/data_for_hicpep/data_store/data/rao_2014/hic/GSE63525_GM12878_insitu_primary_replicate_combined_30.hic"
OUTPUT_PATH="/home/jordan990301/Projects/HiCPEP/benchmark/juicer/output"
LOG="${OUTPUT_PATH}/juicer.log"

java -Xmx13004m -Xms13004m -jar "${JUICER_TOOLS}" eigenvector KR "${HIC_PATH}" "2" BP 100000 "${OUTPUT_PATH}/pc1_chr2_100000.txt" -p 2>&1 | tee -a "${LOG}"
