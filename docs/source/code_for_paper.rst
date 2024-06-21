Code for paper
==============

This page explains how to reproduce the experiment results of our research. 
We store all the programs in the `code_for_paper <https://github.com/ZhiRongDev/HiCPEP/blob/main/code_for_paper>`_ directory.

Prerequisites
-------------
* Make sure you have enough disk size for at least 120 GB, and the wifi connection is stable. 
* Make sure you have at least 16 GB memory.
* Make sure you have installed HiCPEP in your system.

Quick start
-----------

First you have to specify the path for storing all the data required and experiment results, we highly recommand to created an empty directory for this purpose. 
Here we assume the data will all be stored in the ``/tmp/data_store`` directory, all you need to do is to paste the commands below and wait for it finish (It will takes few hours).

.. code:: bash

    mkdir /tmp/data_store
    git clone git@github.com:ZhiRongDev/HiCPEP.git
    cd HiCPEP/code_for_paper
    bash run.sh -p /tmp/data_store

We explain the directory structure in the `data_store.rst <https://github.com/ZhiRongDev/HiCPEP/blob/main/docs/source/data_store.rst>`_.

Guidance
--------

The following is the content of the ``code_for_paper`` directory:

.. code:: bash

    code_for_paper
    ├── benchmark
    │   ├── juicer
    │   │   ├── run_100Kb_pc1.sh
    │   │   ├── run_100Kb_pearson.sh
    │   │   ├── run_1Mb_pc1.sh
    │   │   ├── run_1Mb_pearson.sh
    │   │   ├── run_25Kb_pc1.sh
    │   │   ├── run_25Kb_pearson.sh
    │   │   └── run_all.sh
    │   ├── benchmark_est_all_100kb.py
    │   ├── benchmark_est_all_1Mb.py
    │   ├── benchmark_est_all_25Kb.py
    │   ├── benchmark_est_sample_100kb.py
    │   ├── benchmark_est_sample_1Mb.py
    │   ├── benchmark_est_sample_25kb.py
    │   ├── benchmark_matrix_size.py
    │   ├── benchmark_scikit_100Kb.py
    │   ├── benchmark_scikit_1Mb.py
    │   ├── benchmark_scikit_25Kb.py
    │   ├── mem_efficient.py
    │   └── sum_zero_percent.py
    ├── experiments
    │   ├── __init__.py
    │   ├── lieberman_2009.py
    │   ├── rao_2014.py
    │   └── utils.py
    ├── notebooks
    │   ├── compare_explained_variance.ipynb
    │   ├── compare_similarity.ipynb
    │   ├── compare_time.ipynb
    │   ├── heatmaps.ipynb
    │   └── plots_for_explaining_formula.ipynb
    ├── reference_gc
    │   ├── hg18
    │   ├── hg19
    │   └── mm9
    ├── build_pearsons_pc1_2014.sh
    ├── create_ref_gc.sh
    ├── download_required_data.sh
    ├── main.py
    └── run.sh

* The ``run.sh`` is the entry point to carry out all the experiments except for ``benchmark``, please start the code tracing from this script if you're interested in how the entire programs work. 
* The ``benchmark`` directory contains Python script we used for benchmarking the time and RAM requirement for our estimation algorithm, Juicer and POSSUMM. Please check `README.md <https://github.com/ZhiRongDev/HiCPEP/blob/main/code_for_paper/benchmark/README.md>`_ for more information.
* The ``notebooks`` directory contains some of the examples we explained in our paper, all the details are written in the markdown of these notebooks.
* The ``experiments`` directory is a Python package used for the experiments, including the experiment process for rao_2014 (GSE63525) and lieberman_2009 (GSE18199).
* The ``reference_gc`` directory contains the chromosome GC content references created by UCSC tools, please read `create_ref_gc.sh <https://github.com/ZhiRongDev/HiCPEP/blob/main/code_for_paper/create_ref_gc.sh>`_ for more information.