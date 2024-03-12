import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import numpy as np
from hicpap import paptools

kwargs = {
    "pearson":  "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/data/rao_2014/juicer_outputs/gm12878/1000000/pearsons/pearson_chr1.txt",
    "zero_mean": True 
}
pearson_np = paptools.read_pearson(**kwargs)
approx_np = paptools.create_approx(pearson_np)
Vh, explained_variances, total_entry_num, valid_entry_num = paptools.pca_on_pearson(pearson_np)

print(paptools.calc_correctness(pc1_np=Vh[0], approx_np=approx_np))

kwargs_2 = {
    "pc1_np": Vh[0],
    "approx_np": approx_np,
    "figsize": 20,
    "scatter": "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/outputs/scatter.png",
    "relative_magnitude": "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/outputs/line.png",
}
paptools.plot_comparison(**kwargs_2)

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
pearson_np = paptools.straw_to_pearson(**kwargs)
approx_np = paptools.create_approx(pearson_np)
Vh, explained_variances, total_entry_num, valid_entry_num = paptools.pca_on_pearson(pearson_np)


### Test
print(paptools.calc_correctness(pc1_np=Vh[0], approx_np=approx_np))

Vh[0] = -Vh[0]
approx_np = -approx_np

kwargs_2 = {
    "pc1_np": Vh[0],
    "approx_np": approx_np,
    "figsize": 20,
    "scatter": "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/outputs/scatter2.png",
    "relative_magnitude": "/home/jordan990301/Projects/HiC-PAPP/docker_volume_test/outputs/line2.png",
}
paptools.plot_comparison(**kwargs_2)