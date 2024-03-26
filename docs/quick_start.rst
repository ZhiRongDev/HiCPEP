Quick start
===========

Here we summarize the main usage of HiCPAP:

1. Get the Pearson matrix through the ``paptools.read_pearson()`` or ``paptools.straw_to_pearson()``.
2. Create the Approximated PC1-pattern with ``paptools.create_approx()``.

::

    from hicpap import paptools

    pearson_np = paptools.straw_to_pearson(
        hic_path="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic", # Path to the Juicer's `.hic` file.
        chrom_x="1", 
        chrom_y="1",
        resolution=1000000,
        normalization="KR",
        data_type="oe", # Note that the Pearson matrix should be derived from the O/E matrix.
    )

    approx_np = paptools.create_approx(pearson_np=pearson_np)
    print(f"approx_np: {approx_np}")

For more details, please check the `examples <https://github.com/ZhiRongDev/HiCPAP/blob/main/examples/>`_. 
If you would like to check the program we used for the paper, please check the `code_for_paper <https://github.com/ZhiRongDev/HiCPAP/blob/main/code_for_paper>`_.