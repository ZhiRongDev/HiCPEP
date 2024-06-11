#!/bin/bash

LOG="20240609.log"

start_time=$(date +%s.%N)
echo "$(date '+%Y-%m-%d %H:%M:%S') 1Mb Chr2 start" 2>&1 | tee -a "${LOG}"
memusage python benchmark_scikit.py
end_time=$(date +%s.%N)
elapsed_time=$(awk -v start="$start_time" -v end="$end_time" 'BEGIN{print (end - start)}')
echo "Elapsed time in seconds: $elapsed_time" 2>&1 | tee -a "${LOG}"