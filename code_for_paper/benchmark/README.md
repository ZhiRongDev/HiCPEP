# Benchmark the time and RAM requirement

For benchmarking the time and RAM requirement of the hicpep PC1-pattern estimation without random sampling techniques at 1Mb 100Kb or 25Kb resolution, use the following command:

```shell
/usr/bin/time -v python benchmark_est_all_1Mb.py 
```

For benchmarking the time and RAM requirement of the hicpep PC1-pattern estimation with random sampling techniques at 1Mb 100Kb or 25Kb resolution, use the following command:

```shell
/usr/bin/time -v python benchmark_est_sample_1Mb.py 
```

For benchmarking the time and RAM requirement of the Power iteration PC1 approximation using Scikit-learn at 1Mb 100Kb or 25Kb resolution, use the following command:

```shell
/usr/bin/time -v python benchmark_scikit_1Mb.py 
```

For benchmarking the time and RAM requirement of Juicer to create the Pearson matrix and PC1 txt files at 1Mb 100Kb or 25Kb resolution, use the following command:

```shell
cd juicer
/usr/bin/time -v run_all.sh
```

For benchmarking the time and RAM requirement of POSSUMM to create the PC1 txt files at 1Mb 100Kb or 25Kb resolution, you should download the Rscript from https://github.com/aidenlab/EigenVector, after all the requirement installed, use the following command:

```shell
 /usr/bin/time -v Rscript eigFromHicRscript.R -v TRUE -n KR GSE63525_GM12878_insitu_primary_replicate_combined_30.hic 2 possumm_pc1.txt 1000000
```