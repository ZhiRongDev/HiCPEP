# Hi-C Pearson matrix Approximated PC1 pattern (HiCPAP)

## Requirements
Create the requirements.txt with the support of `pip-tools`.
`pip-compile requirements.in`

## Notes
* Only provide scripts for calculating the approximation track.
* 用 Docker Volume 存程式跑完的資料，等全部跑完以後再用 docker cp 指令存到本機。
* The NaN value are all fill with 0.
* Teacher asked create_approx function should also acccept the hic contact map, not only the pearson correlation matrix.

## hicapp

* Read in pearson matrix and output approx pc1.
* Take the approx pc1 and the compared pc1, plot scatter or line.

How the Juicer implemented there eigenvector-calculation:
https://github.com/aidenlab/Juicebox/blob/12bc67454c460159c758606ff05eaecaf6601d1b/src/juicebox/data/MatrixZoomData.java#L708

* Fill the matrix with 0 first, and remove rows/columns with all 0 entries.
* Record the 0 (NaN) index position with numpy.diag, which should be 0.
* Take advantage of numpy.ix_() to get the submatrix.
