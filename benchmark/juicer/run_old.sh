#!/bin/bash

start_time=$(date +%s.%N)
echo "$(date '+%Y-%m-%d %H:%M:%S') 1Mb Chr2 start" 2>&1 | tee -a "${LOG}"
memusage bash run_1Mb_pearson.sh
end_time=$(date +%s.%N)
elapsed_time=$(awk -v start="$start_time" -v end="$end_time" 'BEGIN{print (end - start)}')
echo "Elapsed time in seconds: $elapsed_time" 2>&1 | tee -a "${LOG}"
echo "1Mb pearson end ======================="

# start_time=$(date +%s.%N)
# echo "$(date '+%Y-%m-%d %H:%M:%S') 1Mb Chr2 start" 2>&1 | tee -a "${LOG}"
# memusage bash run_1Mb_pc1.sh
# end_time=$(date +%s.%N)
# elapsed_time=$(awk -v start="$start_time" -v end="$end_time" 'BEGIN{print (end - start)}')
# echo "Elapsed time in seconds: $elapsed_time" 2>&1 | tee -a "${LOG}"
# echo "1Mb pc1 end ======================="

# start_time=$(date +%s.%N)
# echo "$(date '+%Y-%m-%d %H:%M:%S') 1Mb Chr2 start" 2>&1 | tee -a "${LOG}"
# memusage bash run_100Kb_pearson.sh
# end_time=$(date +%s.%N)
# elapsed_time=$(awk -v start="$start_time" -v end="$end_time" 'BEGIN{print (end - start)}')
# echo "Elapsed time in seconds: $elapsed_time" 2>&1 | tee -a "${LOG}"
# echo "100Kb pearson end ======================="

# start_time=$(date +%s.%N)
# echo "$(date '+%Y-%m-%d %H:%M:%S') 1Mb Chr2 start" 2>&1 | tee -a "${LOG}"
# memusage bash run_100Kb_pc1.sh
# end_time=$(date +%s.%N)
# elapsed_time=$(awk -v start="$start_time" -v end="$end_time" 'BEGIN{print (end - start)}')
# echo "Elapsed time in seconds: $elapsed_time" 2>&1 | tee -a "${LOG}"
# echo "100Kb pc1 end ======================="

# start_time=$(date +%s.%N)
# echo "$(date '+%Y-%m-%d %H:%M:%S') 1Mb Chr2 start" 2>&1 | tee -a "${LOG}"
# memusage bash run_25Kb_pearson.sh
# end_time=$(date +%s.%N)
# elapsed_time=$(awk -v start="$start_time" -v end="$end_time" 'BEGIN{print (end - start)}')
# echo "Elapsed time in seconds: $elapsed_time" 2>&1 | tee -a "${LOG}"
# echo "25Kb pearson end ======================="

# start_time=$(date +%s.%N)
# echo "$(date '+%Y-%m-%d %H:%M:%S') 1Mb Chr2 start" 2>&1 | tee -a "${LOG}"
# memusage bash run_25Kb_pc1.sh
# end_time=$(date +%s.%N)
# elapsed_time=$(awk -v start="$start_time" -v end="$end_time" 'BEGIN{print (end - start)}')
# echo "Elapsed time in seconds: $elapsed_time" 2>&1 | tee -a "${LOG}"
# echo "25Kb pc1 end ======================="
