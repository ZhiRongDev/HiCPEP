#!/bin/bash

echo "benchmark_scikit"
python -m memory_profiler benchmark_scikit.py
echo "benchmark_est_all"
python -m memory_profiler benchmark_est_all.py
echo "benchmark_est_sample"
python -m memory_profiler benchmark_est_sample.py

###
# echo "mem_efficient"
# python -m memory_profiler ./new_algo/mem_efficient.py