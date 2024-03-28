data_store directory structure
==============================

This page explains the output directory structure (i.e. ``data_store``) of the ``code_for_paper``, 
please read the `code_for_paper.rst <https://github.com/ZhiRongDev/HiCPAP/blob/main/docs/code_for_paper.rst>`_ before proceeding.

In the first layer of ``data_store`` there are three directories:

1. data
2. juicer_tools
3. outputs

The ``data`` directory is used for storing the data required for the experiments of Lieberman, 2009 and Rao, 2014, 
including the Pearson matrices and PC1s for each cell line at the resolution of 1Mb and 100Kb.
For the experiments of Lieberman, 2009, the data are directly downloaded from `GSE18199 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199>`_; 
For the experiments of Rao, 2014, the ``.hic`` data are downloaded from `GSE63525 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525>`_, 
and processed with `juicer_tools 1.22.01 <https://github.com/aidenlab/juicer/wiki/Download>`_ for creating the Pearsons and PC1s.

The ``outputs`` directory is used for storing the experiment results, including the Approximated PC1-pattern ``.txt`` files, scatter & relative_magnitude plots and the summary informations of all the experiments.

Here we list the full directory structure, which might be informative before/after you execute the `run.sh <https://github.com/ZhiRongDev/HiCPAP/blob/main/code_for_paper/run.sh>`_ 
(We only list the directories and skip the files):

.. code:: bash

    data_store
    ├── data
    │   ├── lieberman_2009
    │   │   ├── eigenvectors
    │   │   ├── heatmaps
    │   │   └── zip
    │   ├── rao_2014
    │   │   ├── hic
    │   │   └── juicer_outputs
    │   │       ├── ch12-lx
    │   │       │   ├── 100000
    │   │       │   │   ├── eigenvector
    │   │       │   │   └── pearsons
    │   │       │   └── 1000000
    │   │       │       ├── eigenvector
    │   │       │       └── pearsons
    │   │       ├── gm12878
    │   │       │   ├── 100000
    │   │       │   │   ├── eigenvector
    │   │       │   │   └── pearsons
    │   │       │   └── 1000000
    │   │       │       ├── eigenvector
    │   │       │       └── pearsons
    │   │       ├── hmec
    │   │       │   ├── 100000
    │   │       │   │   ├── eigenvector
    │   │       │   │   └── pearsons
    │   │       │   └── 1000000
    │   │       │       ├── eigenvector
    │   │       │       └── pearsons
    │   │       ├── huvec
    │   │       │   ├── 100000
    │   │       │   │   ├── eigenvector
    │   │       │   │   └── pearsons
    │   │       │   └── 1000000
    │   │       │       ├── eigenvector
    │   │       │       └── pearsons
    │   │       ├── imr90
    │   │       │   ├── 100000
    │   │       │   │   ├── eigenvector
    │   │       │   │   └── pearsons
    │   │       │   └── 1000000
    │   │       │       ├── eigenvector
    │   │       │       └── pearsons
    │   │       ├── k562
    │   │       │   ├── 100000
    │   │       │   │   ├── eigenvector
    │   │       │   │   └── pearsons
    │   │       │   └── 1000000
    │   │       │       ├── eigenvector
    │   │       │       └── pearsons
    │   │       ├── kbm7
    │   │       │   ├── 100000
    │   │       │   │   ├── eigenvector
    │   │       │   │   └── pearsons
    │   │       │   └── 1000000
    │   │       │       ├── eigenvector
    │   │       │       └── pearsons
    │   │       └── nhek
    │   │           ├── 100000
    │   │           │   ├── eigenvector
    │   │           │   └── pearsons
    │   │           └── 1000000
    │   │               ├── eigenvector
    │   │               └── pearsons
    │   └── ucsc
    ├── juicer_tools
    └── outputs
        ├── approx_pc1_pattern
        │   ├── lieberman_2009
        │   │   ├── gm06690
        │   │   │   ├── 100000
        │   │   │   │   ├── cxmax
        │   │   │   │   └── cxmin
        │   │   │   └── 1000000
        │   │   │       ├── cxmax
        │   │   │       └── cxmin
        │   │   └── k562
        │   │       ├── 100000
        │   │       │   ├── cxmax
        │   │       │   └── cxmin
        │   │       └── 1000000
        │   │           ├── cxmax
        │   │           └── cxmin
        │   └── rao_2014
        │       ├── ch12-lx
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   └── cxmin
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       └── cxmin
        │       ├── gm12878
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   └── cxmin
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       └── cxmin
        │       ├── hmec
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   └── cxmin
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       └── cxmin
        │       ├── huvec
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   └── cxmin
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       └── cxmin
        │       ├── imr90
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   └── cxmin
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       └── cxmin
        │       ├── k562
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   └── cxmin
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       └── cxmin
        │       ├── kbm7
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   └── cxmin
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       └── cxmin
        │       └── nhek
        │           ├── 100000
        │           │   ├── cxmax
        │           │   └── cxmin
        │           └── 1000000
        │               ├── cxmax
        │               └── cxmin
        ├── plots
        │   ├── lieberman_2009
        │   │   ├── gm06690
        │   │   │   ├── 100000
        │   │   │   │   ├── cxmax
        │   │   │   │   │   ├── relative_magnitude
        │   │   │   │   │   └── scatter
        │   │   │   │   └── cxmin
        │   │   │   │       ├── relative_magnitude
        │   │   │   │       └── scatter
        │   │   │   └── 1000000
        │   │   │       ├── cxmax
        │   │   │       │   ├── relative_magnitude
        │   │   │       │   └── scatter
        │   │   │       └── cxmin
        │   │   │           ├── relative_magnitude
        │   │   │           └── scatter
        │   │   └── k562
        │   │       └── 1000000
        │   │           ├── cxmax
        │   │           │   ├── relative_magnitude
        │   │           │   └── scatter
        │   │           └── cxmin
        │   │               ├── relative_magnitude
        │   │               └── scatter
        │   └── rao_2014
        │       ├── ch12-lx
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   │   ├── relative_magnitude
        │       │   │   │   └── scatter
        │       │   │   └── cxmin
        │       │   │       ├── relative_magnitude
        │       │   │       └── scatter
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       │   ├── relative_magnitude
        │       │       │   └── scatter
        │       │       └── cxmin
        │       │           ├── relative_magnitude
        │       │           └── scatter
        │       ├── gm12878
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   │   ├── relative_magnitude
        │       │   │   │   └── scatter
        │       │   │   └── cxmin
        │       │   │       ├── relative_magnitude
        │       │   │       └── scatter
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       │   ├── relative_magnitude
        │       │       │   └── scatter
        │       │       └── cxmin
        │       │           ├── relative_magnitude
        │       │           └── scatter
        │       ├── hmec
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   │   ├── relative_magnitude
        │       │   │   │   └── scatter
        │       │   │   └── cxmin
        │       │   │       ├── relative_magnitude
        │       │   │       └── scatter
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       │   ├── relative_magnitude
        │       │       │   └── scatter
        │       │       └── cxmin
        │       │           ├── relative_magnitude
        │       │           └── scatter
        │       ├── huvec
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   │   ├── relative_magnitude
        │       │   │   │   └── scatter
        │       │   │   └── cxmin
        │       │   │       ├── relative_magnitude
        │       │   │       └── scatter
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       │   ├── relative_magnitude
        │       │       │   └── scatter
        │       │       └── cxmin
        │       │           ├── relative_magnitude
        │       │           └── scatter
        │       ├── imr90
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   │   ├── relative_magnitude
        │       │   │   │   └── scatter
        │       │   │   └── cxmin
        │       │   │       ├── relative_magnitude
        │       │   │       └── scatter
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       │   ├── relative_magnitude
        │       │       │   └── scatter
        │       │       └── cxmin
        │       │           ├── relative_magnitude
        │       │           └── scatter
        │       ├── k562
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   │   ├── relative_magnitude
        │       │   │   │   └── scatter
        │       │   │   └── cxmin
        │       │   │       ├── relative_magnitude
        │       │   │       └── scatter
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       │   ├── relative_magnitude
        │       │       │   └── scatter
        │       │       └── cxmin
        │       │           ├── relative_magnitude
        │       │           └── scatter
        │       ├── kbm7
        │       │   ├── 100000
        │       │   │   ├── cxmax
        │       │   │   │   ├── relative_magnitude
        │       │   │   │   └── scatter
        │       │   │   └── cxmin
        │       │   │       ├── relative_magnitude
        │       │   │       └── scatter
        │       │   └── 1000000
        │       │       ├── cxmax
        │       │       │   ├── relative_magnitude
        │       │       │   └── scatter
        │       │       └── cxmin
        │       │           ├── relative_magnitude
        │       │           └── scatter
        │       └── nhek
        │           ├── 100000
        │           │   ├── cxmax
        │           │   │   ├── relative_magnitude
        │           │   │   └── scatter
        │           │   └── cxmin
        │           │       ├── relative_magnitude
        │           │       └── scatter
        │           └── 1000000
        │               ├── cxmax
        │               │   ├── relative_magnitude
        │               │   └── scatter
        │               └── cxmin
        │                   ├── relative_magnitude
        │                   └── scatter
        └── summary
