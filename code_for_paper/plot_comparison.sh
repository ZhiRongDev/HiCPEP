#!/bin/bash
DATA_PATH="${DOCKER_VOLUME_PATH}/data"
OUTPUT_PATH="${DOCKER_VOLUME_PATH}/outputs/plots/Rao_2014"
mkdir -p "${OUTPUT_PATH}"

PC1=$DATA_PATH/Rao_2014/juicer_outputs/GM12878/1000000/eigenvector/pc1_chr1.txt
APPROX=${DOCKER_VOLUME_PATH}/outputs/approx_PC1_pattern/Rao_2014/GM12878/1000000/CxMax/approx_PC1_pattern_chr1.txt

python src/common/plot_comparison.py \
    --PC1 $PC1 \
    --approx $APPROX \
    --relative_magnitude "$OUTPUT_PATH/line.png" \
    --scatter "$OUTPUT_PATH/scatter.png" \
    --figsize 50

# #!/bin/bash
# PY_FILE="src/common/plot_comparison.py"

# ### Rao 2014
# DATA_PATH="${DOCKER_VOLUME_PATH}/data/Rao_2014/juicer_outputs"
# OUTPUT_PATH="${DOCKER_VOLUME_PATH}/outputs/approx_PC1_pattern/Rao_2014"
# mkdir -p "${OUTPUT_PATH}"

# echo "$(date '+%Y-%m-%d %H:%M:%S') Rao_2014 plot_comparison start"

# for RESOLUTION in "1000000" "100000"
# do
#     CELL_LINES=(\
#         "GM12878"\
#         "IMR90" \
#         "HMEC" \
#         "NHEK" \
#         "K562" \
#         "KBM7" \
#         "HUVEC" \
#         "HeLa" \
#         "CH12-LX" \
#     )

#     for TYPE in "CxMax" "CxMin"
#     do
#         for CELL_LINE in "${CELL_LINES[@]}"
#         do
#             echo "$(date '+%Y-%m-%d %H:%M:%S') [${RESOLUTION}] [${TYPE}] [${CELL_LINE}] start"

#             # There are only 19 chromosome in CH12-LX, 
#             if [[ "${CELL_LINE}" == "CH12-LX" ]]
#             then
#                 CHROM_LIST=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y)
#             else
#                 CHROM_LIST=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
#             fi

#             for CHROM in "${CHROM_LIST[@]}"
#             do
#                 PEARSON_FILE="${DATA_PATH}/${CELL_LINE}/${RESOLUTION}/pearsons/pearson_chr${CHROM}.txt"
#                 OUTPUT_FILE="${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/${TYPE}/approx_PC1_pattern_chr${CHROM}.txt"
#                 python "${PY_FILE}" --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE --source "Rao_2014"
#             done

#             echo "$(date '+%Y-%m-%d %H:%M:%S') [${RESOLUTION}] [${TYPE}] [${CELL_LINE}] end"
#         done
#     done
# done

# echo "$(date '+%Y-%m-%d %H:%M:%S') Rao_2014 plot_comparison end"

# ### Lieberman 2009 
# DATA_PATH="${DOCKER_VOLUME_PATH}/data/Lieberman_2009"
# OUTPUT_PATH="${DOCKER_VOLUME_PATH}/outputs/approx_PC1_pattern/Lieberman_2009"
# mkdir -p "${OUTPUT_PATH}"

# echo "$(date '+%Y-%m-%d %H:%M:%S') Lieberman_2009 plot_comparison start"
# for RESOLUTION in "1000000" "100000" 
# do
#     CELL_LINES=(\
#         "gm06690"\
#         "k562" \
#     )

#     for TYPE in "CxMax" "CxMin"
#     do
#         for CELL_LINE in "${CELL_LINES[@]}"
#         do
#             echo "$(date '+%Y-%m-%d %H:%M:%S') [${RESOLUTION}] [${TYPE}] [${CELL_LINE}] start"
            
#             # There are some missing chromosomes in Lieberman's dataset.
#             if [[ "${CELL_LINE}" == "k562" && "${RESOLUTION}" == "100000" ]]
#             then
#                 CHROM_LIST=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
#             else
#                 CHROM_LIST=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#             fi

#             for CHROM_ID in "${CHROM_LIST[@]}"
#             do
#                 PEARSON_FILE="${DATA_PATH}/heatmaps/HIC_${CELL_LINE}_chr${CHROM_ID}_chr${CHROM_ID}_${RESOLUTION}_pearson.txt"
#                 OUTPUT_FILE="${OUTPUT_PATH}/${CELL_LINE}/${RESOLUTION}/${TYPE}/approx_PC1_pattern_chr${CHROM_ID}.txt"
#                 python "${PY_FILE}" --pearson "${PEARSON_FILE}" --output "${OUTPUT_FILE}" --type $TYPE --source "Lieberman_2009"
#             done

#             echo "$(date '+%Y-%m-%d %H:%M:%S') [${RESOLUTION}] [${TYPE}] [${CELL_LINE}] end"
#         done
#     done
# done

# echo "$(date '+%Y-%m-%d %H:%M:%S') Lieberman_2009 plot_comparison end"