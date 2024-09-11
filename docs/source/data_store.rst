data_store directory structure
==============================

This page explains the output directory structure (i.e. ``data_store``) of the ``code_for_paper``, 
please read the `code_for_paper.rst <https://github.com/ZhiRongDev/HiCPEP/blob/main/docs/code_for_paper.rst>`_ before proceeding.

In the first layer of ``data_store`` there are three directories:

1. data
2. juicer_tools
3. outputs

The ``data`` directory is used for storing the data required for the experiments of Lieberman, 2009 (GSE18199) and Rao, 2014 (GSE63525), 
including the Pearson matrices and PC1s for each cell line at the resolution of 1Mb and 100Kb.
For the experiments of Lieberman, 2009, the data are directly downloaded from `GSE18199 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18199>`_; 
For the experiments of Rao, 2014, the ``.hic`` data are downloaded from `GSE63525 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525>`_, 
and processed with `juicer_tools 1.22.01 <https://github.com/aidenlab/juicer/wiki/Download>`_ for creating the Pearsons and PC1s.

The ``outputs`` directory is used for storing the experiment results, including the Estimated PC1-pattern ``.txt`` files, scatter & relative_magnitude plots and the summary informations of all the experiments.
Note that we created the Estimated PC1-pattern by selecting the ``cxmax`` or ``cxmin`` of the Pearson's covariance matrix, at the resolution of 1Mb and 100Kb; 
Besides, since GSE63525 doesn't provide the `.hic` files for HeLa, we skip this cell line. 

Here we further explain the details of the ``outputs`` directory structure, in the first layer of ``outputs`` there are three directories:

1. est_pc1_pattern
2. plots
3. summary

The ``est_pc1_pattern`` directory contains the text files of the Estimated PC1-pattern.
The ``plots`` directories including the scatter and the relative_magnitude plots. 
In the ``summary`` directory there are (2009 means using the data from GSE18199; 2014 means using the data from GSE63525):

1. ``summary_similarity_2009.xlsx`` and ``summary_similarity_2014.xlsx``, `summary_similarity_2014_sample10.xlsx`, which is for comparing the ``similar_rate`` between the juicer_tools calculated PC1 and the Estimated PC1-pattern, with and  without using sampling method.
2. ``summary_self_pca_2009.xlsx`` and ``summary_self_pca_2014.xlsx``, which is for recording the explained variance ratio of the first 3 Principal components of the Pearson matrix, and for recording the ``similar_rate`` between the self calaulated PC1 (NOT the juicer_tools calculated PC1) and the Estimated PC1-pattern.
3. ``summary_similar_rate_percentage_2014.xlsx``, which is used for summarizing the percentage of columns in the covariance matrix that has a similar_rate over 90%, 95% or 99%.  