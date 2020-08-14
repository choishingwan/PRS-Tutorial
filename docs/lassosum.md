# Background

`lassosum` is one of the dedicated PRS programs which is an `R` package that uses penalised regression (LASSO) in its approach to PRS calculation.

## Installing lassosum

!!! note
    The script used here is based on lassosum version 0.4.4

!!! note
    For more details, please refer to [lassosum's homepage](https://github.com/tshmak/lassosum/)

You can install `lassosum` and its dependencies in `R` with the following command:

```R
install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
library(devtools)
install_github("tshmak/lassosum")
```

## Required Data

Again, we assume that we have the following files (or you can download it from [here](https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view?usp=sharing)): 

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.cov**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the PCs of the samples |

## Running PRS analysis

We can run lassosum as follows: 

``` R
library(lassosum)
# Prefer to work with data.table as it speeds up file reading
library(data.table)
library(methods)
library(magrittr)
# For multi-threading, you can use the parallel package and 
# invoke cl which is then passed to lassosum.pipeline
library(parallel)
# This will invoke 2 threads. 
cl <- makeCluster(2)

sum.stat <- "Height.QC.gz"
bfile <- "EUR.QC"
# Read in and process the covariates
covariate <- fread("EUR.cov")
pcs <- fread("EUR.eigenvec") %>%
    setnames(., colnames(.), c("FID","IID", paste0("PC",1:6)))
# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well
cov <- merge(covariate, pcs)

# We will need the EUR.hg19 file provided by lassosum 
# which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
ld.file <- "EUR.hg19"
# output prefix
prefix <- "EUR"
# Read in the target phenotype file
target.pheno <- fread("EUR.height")[,c("FID", "IID", "Height")]
# Read in the summary statistics
ss <- fread(sum.stat)
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P == 0]
# Transform the P-values into correlation
cor <- p2cor(p = ss$P,
        n = ss$N,
        sign = log(ss$OR)
        )
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]


# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing
out <- lassosum.pipeline(
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$A1,
    A2 = ss$A2,
    ref.bfile = bfile,
    test.bfile = bfile,
    LDblocks = ld.file, 
    cluster=cl
)
# Store the R2 results
target.res <- validate(out, pheno = target.pheno, covar=cov)
# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
```


??? note "How much phenotypic variation does the "best-fit" PRS explain?"
    0.2395471
