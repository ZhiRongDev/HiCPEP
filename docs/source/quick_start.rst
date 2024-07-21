Quick start
===========

Here we summarize the main usage of HiCPEP:

1. Get the Pearson matrix through the ``peptools.read_pearson()`` or other tools such as `Straw <https://github.com/aidenlab/straw>`_ or `cooltools <https://github.com/open2c/cooltools>`_.
2. Create the Estimated PC1-pattern with ``peptools.create_est()``.

.. code::

    from hicpep import peptools
    import hicstraw

    hic_path="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic", # Path to the Juicer's `.hic` file.
    chrom = "1"
    resolution = 1000000
    normalization = "KR"

    hic = hicstraw.HiCFile(hic_path)

    for chromosome in hic.getChromosomes():
        if chromosome.name == chrom:
            chrom_size = int(chromosome.length)

    matrix = hic.getMatrixZoomData(chrom, chrom, "oe", normalization, "BP", resolution)
    matrix_np = matrix.getRecordsAsMatrix(0, chrom_size, 0, chrom_size)
    pearson_np = np.corrcoef(matrix_np)

    est_np = peptools.create_est(pearson_np=pearson_np)
    print(f"est_np: {est_np}")

For more details, please check the tutorial in the `examples directory <https://github.com/ZhiRongDev/HiCPEP/blob/main/examples/>`_. 
If you are looking for the programs we used in the paper, please check the `code_for_paper <https://github.com/ZhiRongDev/HiCPEP/blob/main/code_for_paper>`_.