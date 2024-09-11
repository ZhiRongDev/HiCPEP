# Hi-C Pearson matrix's Estimated PC1-pattern (HiCPEP)

HiCPEP is a Python package for creating the Estimated PC1-pattern of the Hi-C Pearson matrix, which can be used for identifying the A/B compartments.

## Requirements and Installation

All the programs were tested in Ubuntu 22.04.4 LTS, HiCPEP requires `python3`, `pip` and `libcurl4-openssl-dev` installed on your system.

For example (Paste these commands in Bash or Zsh):

```shell
sudo apt-get update
sudo apt-get install -y libcurl4-openssl-dev # For installing hic-straw
sudo apt-get install -y python3
sudo apt-get install -y pip
sudo apt-get install -y git 
git clone git@github.com:ZhiRongDev/HiCPEP.git
cd HiCPEP
python3 -m pip install -e .
```

If you have already installed the requirements, just paste these commands:

```shell
git clone git@github.com:ZhiRongDev/HiCPEP.git
cd HiCPEP
python3 -m pip install -e .
```

## Quick start

Case 1: Using tools such as [Straw](https://github.com/aidenlab/straw) to create the Pearson matrix as input.

```py
# Using hic-straw
from hicpep import peptools
import hicstraw
import numpy as np

hic_path="https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic" # Path to the Juicer's `.hic` file.
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
```

Case 2:Using the [Juicer](https://github.com/aidenlab/juicer/wiki/Pearsons) created Pearson text file as input.

```py
from hicpep import peptools
pearson_np = peptools.read_pearson(
    pearson="gm12878_1000000_pearson_chr1.txt"
)

est_np = peptools.create_est(pearson_np=pearson_np)
print(f"est_np: {est_np}")
```

For more details, please check the [examples](https://github.com/ZhiRongDev/HiCPEP/blob/main/examples/). If you are interested in the programs we used for the paper, please check the [code_for_paper](https://github.com/ZhiRongDev/HiCPEP/blob/main/code_for_paper/). The HiCPEP Python library depends on [NumPy](https://numpy.org/doc/stable/index.html), [pandas](https://pandas.pydata.org/), [SciPy](https://scipy.org/) and [Matplotlib](https://matplotlib.org/).

## References

Zhi-Rong Cheng, Jia-Ming Chang. Decoding the Power of PC1: A Fast and Accurate Covariance-Based Method for A/B Compartment Identification in Hi-C Data.