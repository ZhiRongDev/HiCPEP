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

# Begin script in case all parameters are correct
echo -e ">>> The docker volume path is set as ${DOCKER_VOLUME_PATH}"
echo -e ">>> Start the process...... \n"
sleep 3

### Prepare for data required ###
echo -e ">>> Start to download the required data. \n"
bash scripts/download_required_data.sh -p "${DOCKER_VOLUME_PATH}"

echo -e ">>> Create the required juicer's pearsons and PC1 for Rao 2014 experiments \n"
mkdir -p "${DOCKER_VOLUME_PATH}/data/Rao_2014/juicer_outputs"
bash scripts/build_PearsonsAndPC1_humanCells.sh -p "${DOCKER_VOLUME_PATH}"
bash scripts/build_PearsonsAndPC1_mouseCH12-LX.sh -p "${DOCKER_VOLUME_PATH}"