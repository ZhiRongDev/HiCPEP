# Hi-C Pearson matrix's Approximated PC1-pattern (HiCPAP)

HiCPAP is a Python package for creating the Approximated PC1-pattern of the Hi-C Pearson matrix, which can be used for identifying the A/B compartments.

## Requirements and Installation

All the programs were tested in Ubuntu 22.04.4 LTS, HiCPAP requires `python3`, `pip` and `libcurl4-openssl-dev` installed on your system. 

For example (Paste these commands in Bash or Zsh):

```shell
sudo apt-get update
sudo apt-get install -y libcurl4-openssl-dev
sudo apt-get install -y python3
sudo apt-get install -y pip
sudo apt-get install -y git 
git clone git@github.com:ZhiRongDev/HiCPAP.git
cd HiCPAP
python3 -m pip install -e .
```

If you have already installed the requirements, just paste these commands:

```shell
git clone git@github.com:ZhiRongDev/HiCPAP.git
cd HiCPAP
python3 -m pip install -e .
```

## Quick start

```py
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
```

For more details, please check the [examples](https://github.com/ZhiRongDev/HiCPAP/blob/main/examples/). If you are interested in the programs we used for the paper, please check the [code_for_paper](https://github.com/ZhiRongDev/HiCPAP/blob/main/code_for_paper/).

## References

* Zhi-Rong Cheng, Jia-Ming Chang. The exploration and optimization of the chromatin compartment analysis.
