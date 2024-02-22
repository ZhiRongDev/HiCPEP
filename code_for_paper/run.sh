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

export DOCKER_VOLUME_PATH # Set as the global variable

# Begin script in case all parameters are correct
echo "$(date '+%Y-%m-%d %H:%M:%S') The docker volume path is set as ${DOCKER_VOLUME_PATH}"
echo "$(date '+%Y-%m-%d %H:%M:%S') Start the process......"
# sleep 3

### Prepare for the data required ###
# echo "$(date '+%Y-%m-%d %H:%M:%S') Start to download the required data."
# bash scripts/download_required_data.sh

# echo "$(date '+%Y-%m-%d %H:%M:%S') Create the required juicer's pearsons and PC1 for the Rao 2014 experiments"
# mkdir -p "${DOCKER_VOLUME_PATH}/data/Rao_2014/juicer_outputs"
# bash scripts/build_PearsonsAndPC1_2014.sh

### Create the approximation ###
echo "$(date '+%Y-%m-%d %H:%M:%S') Create the approximated PC1-pattern for the experiments of Rao 2014 and Lieberman 2009"
bash scripts/create_approx.sh