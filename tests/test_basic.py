import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from hicpapp import papptools

kwargs = {
    "pearson":  "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/data/rao_2014/juicer_outputs/gm12878/1000000/pearsons/pearson_chr1.txt",
    "zero_mean": True 
}
pearson_np = papptools.read_pearson(**kwargs)
approx_np = papptools.create_approx(pearson_np)
Vh, explained_variances, total_entry_num, valid_entry_num = papptools.calc_explained_variance(pearson_np)

print(papptools.calc_correctness(pc1_np=Vh[0], approx_np=approx_np))

kwargs_2 = {
    "pc1_np": Vh[0],
    "approx_np": approx_np,
    "figsize": 20,
    "scatter": "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/outputs/scatter.png",
    "relative_magnitude": "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/outputs/line.png",
}
papptools.plot_comparison(**kwargs_2)

####
kwargs = {
    "hic_path": "https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic",
    "chrom_x": "1",
    "chrom_y": "1",
    "resolution": 1000000,
    "norm": "KR",
    "method": "oe",
    "zero_mean": True 
}
pearson_np = papptools.straw_to_pearson(**kwargs)
approx_np = papptools.create_approx(pearson_np)
Vh, explained_variances, total_entry_num, valid_entry_num = papptools.calc_explained_variance(pearson_np)


### Test
print(papptools.calc_correctness(pc1_np=Vh[0], approx_np=approx_np))

Vh[0] = -Vh[0]
approx_np = -approx_np

kwargs_2 = {
    "pc1_np": Vh[0],
    "approx_np": approx_np,
    "figsize": 20,
    "scatter": "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/outputs/scatter2.png",
    "relative_magnitude": "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/outputs/line2.png",
}
papptools.plot_comparison(**kwargs_2)