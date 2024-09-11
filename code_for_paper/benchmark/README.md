# Benchmark the time and RAM requirement

* `1_calc_matrix_size.py` is for calculating the matrix size of GM12878 chromosome 2 and GM06690 chromosome 2 mentioned in thesis page 31, which is not calcuated from GNU time tool
* `2_sum_zero_percent.py` is for calculating the sparsity of matrices mentioned in thesis page 40.

---

> We manually tested the Python script and record the Elapsed (wall clock) time, Maximum resident set size and similar_rate.

For benchmarking the time and RAM requirement of the hicpep PC1-pattern estimation start from dense Pearson matrix without random sampling techniques at 1Mb 100Kb or 25Kb resolution, use the following command (The only different in these scripts is the parameter setting):

```shell
/usr/bin/time -v python benchmark_est_all_1Mb.py 
/usr/bin/time -v python benchmark_est_all_100Kb.py 
/usr/bin/time -v python benchmark_est_all_25Kb.py 
```

For benchmarking the time and RAM requirement of the hicpep PC1-pattern estimation start from dense Pearson matrix with random sampling techniques at 1Mb 100Kb or 25Kb resolution, use the following command (The only different in these scripts is the parameter setting):

```shell
/usr/bin/time -v python benchmark_est_sample_1Mb.py 
/usr/bin/time -v python benchmark_est_sample_100Kb.py 
/usr/bin/time -v python benchmark_est_sample_25Kb.py 
```

For benchmarking the time and RAM requirement of the Power iteration PC1 approximation using Scikit-learn start from dense Pearson matrix  at 1Mb 100Kb or 25Kb resolution, use the following command (The only different in these scripts is the parameter setting):

```shell
/usr/bin/time -v python benchmark_scikit_1Mb.py 
/usr/bin/time -v python benchmark_scikit_100Kb.py 
/usr/bin/time -v python benchmark_scikit_25Kb.py 
```

For benchmarking the time and RAM requirement of the hicpep PC1-pattern estimation start from sparse O/E matrix with random sampling techniques at 1Mb or 100Kb resolution, use the following command:

```shell
python 3_store_oe_sparse.py # Store the sparse O/E matrix as the .npz file, which is not included in benchmark.
/usr/bin/time -v python benchmark_est_mem_efficient_1Mb.py
/usr/bin/time -v python benchmark_est_mem_efficient_100Kb.py
```

---
For benchmarking the time and RAM requirement of Juicer to create the Pearson matrix and PC1 txt files at 1Mb 100Kb or 25Kb resolution, use the following command:

```shell
cd juicer
/usr/bin/time -v run_all.sh
```

For benchmarking the time and RAM requirement of POSSUMM to create the PC1 txt files at 1Mb 100Kb or 25Kb resolution. The code in the `code_for_paper/benchmark/POSSUMM/R` folder are downloaded from https://github.com/aidenlab/EigenVector.

**References**  
* Fine-mapping of nuclear compartments using ultra-deep Hi-C shows that active promoter and enhancer elements localize in the active A compartment even when adjacent sequences do not (https://www.biorxiv.org/content/10.1101/2021.10.03.462599v2)
* Harris, H.L., Gu, H., Olshansky, M. et al. Chromatin alternates between A and B compartments at kilobase scale for subgenic organization. Nat Commun 14, 3303 (2023). https://doi.org/10.1038/s41467-023-38429-1

After all the requirements installed, use the following command:

```shell
 /usr/bin/time -v Rscript eigFromHicRscript.R -v TRUE -n KR GSE63525_GM12878_insitu_primary_replicate_combined_30.hic 2 possumm_pc1.txt 1000000
 /usr/bin/time -v Rscript eigFromHicRscript.R -v TRUE -n KR GSE63525_GM12878_insitu_primary_replicate_combined_30.hic 2 possumm_pc1.txt 100000
 /usr/bin/time -v Rscript eigFromHicRscript.R -v TRUE -n KR GSE63525_GM12878_insitu_primary_replicate_combined_30.hic 2 possumm_pc1.txt 25000
```
It will shows the time and RAM usage of POSSUMM. The `similar_rate` compared with Juicer's PC1 will be calculated using `code_for_paper/benchmark/POSSUMM/R/benchmark_similarity_POSSUMM.ipynb`.

---
**Importent**
* Note that the 'time' we mentioned in our thesis is the wall clock time calculated by `/usr/bin/time`, but not the Python time package imported in the Python script. 
For the example showing below, "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:27.54" will be recorded (And the "Maximum resident set size (kbytes): 2528488").

```code
(env) (base) ➜  benchmark git:(main) ✗ /usr/bin/time -v python benchmark_est_all_25Kb.py
Time spent for creating the Estimated PC1-pattern by finding the CxMax in the full-covariance matrix (seconds): 6.320481300354004
{'total_entry_num': 9728, 'valid_entry_num': 9519, 'similar_num': 9234, 'similar_rate': 0.9700598802395209}
        Command being timed: "python benchmark_est_all_25Kb.py"
        User time (seconds): 86.22
        System time (seconds): 5.04
        Percent of CPU this job got: 331%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 0:27.54
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        y: 2528488
        Average resident set size (kbytes): 0
        Major (requiring I/O) page faults: 0
        Minor (reclaiming a frame) page faults: 339242
        Voluntary context switches: 26
        Involuntary context switches: 7897
        Swaps: 0
        File system inputs: 0
        File system outputs: 0
        Socket messages sent: 0
        Socket messages received: 0
        Signals delivered: 0
```

* Besides, we mention that the time and RAM usage for calculating `similar_rate` can be ignored since they are both too small.