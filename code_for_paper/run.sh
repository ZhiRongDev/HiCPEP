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

export DOCKER_VOLUME_PATH # Set as the environment variable

### Begin 
echo "$(date '+%Y-%m-%d %H:%M:%S') The docker volume path is set as ${DOCKER_VOLUME_PATH}"
echo "$(date '+%Y-%m-%d %H:%M:%S') Program start."
# sleep 3

### Prepare for the data required.
# echo "$(date '+%Y-%m-%d %H:%M:%S') Start to download the required data."
# bash download_required_data.sh

# echo "$(date '+%Y-%m-%d %H:%M:%S') Create the required juicer's pearsons and PC1 for the Rao 2014 experiments"
# mkdir -p "${DOCKER_VOLUME_PATH}/data/rao_2014/juicer_outputs"
# bash build_pearsons_pc1_2014.sh

### Start the program
python main.py --docker_volume_path "${DOCKER_VOLUME_PATH}"
echo "$(date '+%Y-%m-%d %H:%M:%S') Program end."