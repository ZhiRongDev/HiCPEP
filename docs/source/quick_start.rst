Quick start
===========

Here we summarize the main usage of HiCPEP:

1. Get the Pearson matrix through the ``peptools.read_pearson()`` or ``peptools.straw_to_pearson()``.
2. Create the Estimated PC1-pattern with ``peptools.create_est()``.

.. code::

    from hicpep import peptools

    pearson_np = peptools.straw_to_pearson(
        hic_path="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic", # Path to the Juicer's `.hic` file.
        chrom="1", 
        resolution=1000000,
        normalization="KR",
    )

    est_np = peptools.create_est(pearson_np=pearson_np)
    print(f"est_np: {est_np}")

For more details, please check the tutorial in the `examples directory <https://github.com/ZhiRongDev/HiCPEP/blob/main/examples/>`_. 
If you are looking for the programs we used in the paper, please check the `code_for_paper <https://github.com/ZhiRongDev/HiCPEP/blob/main/code_for_paper>`_.