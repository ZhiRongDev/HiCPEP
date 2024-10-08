# Downloaded from https://github.com/aidenlab/EigenVector.
# **References**  
# * Fine-mapping of nuclear compartments using ultra-deep Hi-C shows that active promoter and enhancer elements localize in the active A compartment even when adjacent sequences do not (https://www.biorxiv.org/content/10.1101/2021.10.03.462599v2)
# * Harris, H.L., Gu, H., Olshansky, M. et al. Chromatin alternates between A and B compartments at kilobase scale for subgenic organization. Nat Commun 14, 3303 (2023). https://doi.org/10.1038/s41467-023-38429-1

library(strawr)
library(optparse,quietly=TRUE)
library(Matrix,quietly=TRUE)
source("bestEigen3.R")

option_list <- list(
  make_option(c("-t","--tolerance"), help = "precision", type="double", default=1.0e-6),
  make_option(c("-m","--maxiter"), help = "maximum iterations", type="integer",default=100),
  make_option(c("-s","--size"),help = "chromosome length (in basepairs)", type="integer",default=NA),
  make_option(c("-v","--verbose"),help = "verbose", type="logical",default=FALSE),
  make_option(c("-n","--norm"),help = "normalization", type="character",default="NONE"),
  make_option(c("-o","--matrix"),help = "use o for observed, othewise o/e", type="character",default="oe")
)
parser <- OptionParser(
   usage = paste("Rscript %prog [OPTIONS] hicFfile chr outFile resolution (basepairs)",
                 "computes the principal eigenvector of the correlation matrix of HiC contacts matrix",
                 "",
                 "hicFile = hic file",
                 "chr = chromosome for which it is requested",
                 "outFile = file to output the eigenvector; outFile.report will also be produced to report run summary",
                 "resolution = resolution in basepairs", sep="\n"),
   epilogue = "hicFile, chr, outFile and resolution are required.\n",
   option_list=option_list)

arguments=NA
tryCatch(
{ arguments = parse_args(parser, positional_arguments=TRUE);},
   error = function(e) { })
if (all(is.na(arguments))) {
     stop (paste("Failed to parse command-line parameters",
                 "Use --help for help"))
}
opts = arguments$options
if (length(arguments$args) < 4) stop("Need inFile, chr, outFile and resolution\nUse --help for help")
hicFile = arguments$args[1]
chr = arguments$args[2]
outFile = arguments$args[3]
binsize = as.numeric(arguments$args[4])

verbose <- opts$verbose
tol <- as.numeric(opts$tolerance)
maxiter <- as.numeric(opts$maxiter)
norm <- opts$norm
matrix <- opts$matrix

print(paste(norm,matrix))
t1 <- system.time(y <- straw(norm,hicFile,chr,chr,"BP",binsize,matrix))
y$counts[is.na(y$counts)] <- 0
y$counts[y$counts == Inf] <- 0
t1 <- t1["elapsed"]
k <- nrow(y)
if (verbose) print(paste("took",t1,"seconds to read",k,"records"),digits=6, quote=FALSE)
i <- y$x/binsize+1
j <- y$y/binsize+1
n <- max(j)
if (!is.na(opts$size)) n <- ceiling(as.numeric(opts$size)/binsize)
t2 <- system.time(x <- sparseMatrix(i=i,j=j,x=y$counts,dims=c(n,n),symmetric=TRUE))
t2 <- t2["elapsed"]
if (verbose) print(paste("took",t2,"seconds to build sparse matrix"),digits=6, quote=FALSE)
rm(y)
rm(i)
rm(j)
for (i in 1:2) gc(verbose=FALSE)
t3 <- system.time(res <- bestEigen(x,maxiter=maxiter,tol=tol))
t3 <- t3["elapsed"]
if (verbose) {
	print(paste("took",t1,"seconds to read",k,"records"),digits=6, quote=FALSE)
	print(paste("took",t2,"seconds to build sparse matrix"),digits=6, quote=FALSE)
	print(paste("eigenvector computation took",t3,"seconds"),digits=6, quote=FALSE)
	print(paste("lam1 =",res$lam1,",er =",res$er,",lam2 =",res$lam2," niter=",res$niter),digits=6, quote=FALSE)
}	

write.table(res$v,outFile,quote=FALSE,row.names=FALSE,col.names=FALSE)
writeLines(
	   c(
	     paste("took",t1,"seconds to read",k,"records"),
	     paste("took",t2,"seconds to build sparse matrix"),
	     paste("eigenvector computation took",t3,"seconds"),
	     paste("lam1 =",res$lam1,",er =",res$er,",lam2 =",res$lam2," niter=",res$niter)
	     ),
	     paste0(outFile,".report")
	 )
