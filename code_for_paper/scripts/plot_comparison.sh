#!/bin/bash
# DATA_PATH="${DOCKER_VOLUME_PATH}/data"
# OUTPUT_PATH="${DOCKER_VOLUME_PATH}/outputs/plots/Rao_2014"
# mkdir -p "${OUTPUT_PATH}"

# PC1=$DATA_PATH/Rao_2014/juicer_outputs/GM12878/1000000/eigenvector/pc1_chr1.txt
# APPROX=${DOCKER_VOLUME_PATH}/outputs/approx_PC1_pattern/Rao_2014/GM12878/1000000/CxMax/approx_PC1_pattern_chr1.txt

# python src/common/plot_comparison.py \
#     --PC1 $PC1 \
#     --approx $APPROX \
#     --relative_magnitude "$OUTPUT_PATH/line.png" \
#     --scatter "$OUTPUT_PATH/scatter.png" \
#     --figsize 50
