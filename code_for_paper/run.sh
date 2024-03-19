#!/bin/bash

helpFunction()
{
    echo ""
    echo "Usage: $0 -p DATA_STORE"
    echo -e "\t-p DATA_STORE is the directory used for saving the required data and the experiment results."
    exit 1 # Exit script after printing help
}

while getopts "p:" opt
do
    case "${opt}" in
        p ) DATA_STORE="${OPTARG}" ;;
        ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
    esac
done

# Print helpFunction in case parameters are empty
if [ -z "${DATA_STORE}" ]
then
    echo "The -p parameter is empty.";
    helpFunction
fi

export DATA_STORE # Set as the environment variable

### Begin 
echo "$(date '+%Y-%m-%d %H:%M:%S') The DATA_STORE path is set as ${DATA_STORE}"
echo "$(date '+%Y-%m-%d %H:%M:%S') Program start."
sleep 3

### Prepare for the data required.
# echo "$(date '+%Y-%m-%d %H:%M:%S') Start to download the required data."
# bash download_required_data.sh

# echo "$(date '+%Y-%m-%d %H:%M:%S') Create the required juicer's pearsons and PC1 for the Rao 2014 experiments"
# mkdir -p "${DATA_STORE}/data/rao_2014/juicer_outputs"
# bash build_pearsons_pc1_2014.sh

### Start the program
python main.py --data_store "${DATA_STORE}"
echo "$(date '+%Y-%m-%d %H:%M:%S') Program end."