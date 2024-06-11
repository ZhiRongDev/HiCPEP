#!/bin/bash
rm -rf /home/jordan990301/Projects/HiCPEP/benchmark/juicer/output/run_all.log
touch /home/jordan990301/Projects/HiCPEP/benchmark/juicer/output/run_all.log
LOG="/home/jordan990301/Projects/HiCPEP/benchmark/juicer/output/run_all.log"

/usr/bin/time -v bash run_1Mb_pearson.sh -p 2>&1 | tee -a "${LOG}"
echo "1Mb pearson end =======================" -p 2>&1 | tee -a "${LOG}" 

/usr/bin/time -v bash run_1Mb_pc1.sh -p 2>&1 | tee -a "${LOG}"
echo "1Mb pc1 end =======================" -p 2>&1 | tee -a "${LOG}"

/usr/bin/time -v bash run_100Kb_pearson.sh -p 2>&1 | tee -a "${LOG}"
echo "100Kb pearson end =======================" -p 2>&1 | tee -a "${LOG}"

/usr/bin/time -v bash run_100Kb_pc1.sh -p 2>&1 | tee -a "${LOG}"
echo "100Kb pc1 end =======================" -p 2>&1 | tee -a "${LOG}"

/usr/bin/time -v bash run_25Kb_pearson.sh -p 2>&1 | tee -a "${LOG}"
echo "25Kb pearson end =======================" -p 2>&1 | tee -a "${LOG}"

/usr/bin/time -v bash run_25Kb_pc1.sh -p 2>&1 | tee -a "${LOG}"
echo "25Kb pc1 end =======================" -p 2>&1 | tee -a "${LOG}"