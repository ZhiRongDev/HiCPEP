Code for paper
==============

This page explains how to reproduce the experiment results of our research. 
We store all the programs in the `code_for_paper <https://github.com/ZhiRongDev/HiCPAP/blob/main/code_for_paper>`_ directory.

Prerequisites
-------------
* Make sure you have enough disk size for at least 120 GB, and the wifi connection is stable. 
* Make sure you have at least 16 GB memory.
* Make sure you have installed HiCPAP in your system.

Quick start
-----------

First you have to specify the path for storing all the data required and experiment results, we highly recommand to created an empty directory for this purpose. 
Here we assume the data will all be stored in the ``/tmp/data_store`` directory, all you need to do is to paste the commands below and wait for it finish (It will takes few hours).

.. code:: bash

    mkdir /tmp/data_store
    git clone git@github.com:ZhiRongDev/HiCPAP.git
    cd HiCPAP/code_for_paper
    bash run.sh -p /tmp/data_store

We explain the directory structure in the `data_store.rst <https://github.com/ZhiRongDev/HiCPAP/blob/main/docs/source/data_store.rst>`_.

Guidance
--------

The following is the content of the ``code_for_paper`` directory:

.. code:: bash

    code_for_paper
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


* The ``run.sh`` is the entry point to carry out all the experiments, please start the code tracing from this script if you're interested in how the entire programs work. 
* The ``notebooks`` directory contains some of the examples we explained in our paper, all the details are written in the markdown of these notebooks.
* The ``experiments`` directory is a Python package used for the experiments, including the experiment process for rao_2014 and lieberman_2009.
* The ``reference_gc`` directory contains the chromosome GC content references created by UCSC tools, please read `create_ref_gc.sh <https://github.com/ZhiRongDev/HiCPAP/blob/main/code_for_paper/create_ref_gc.sh>`_ for more information.