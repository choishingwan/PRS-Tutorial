`lassosum` is an `R` package for PRS calculation. 
It uses LASSO/Elastic Net estimates rather than p-value thresholding to generate PRS and are 
expected to provide a higher $R^2$ when compared to p-value thresholding.

You can install `lassosum` and its dependencies in `R` with the following

```R
install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
install_github("tshmak/lassosum")
```

Assuming we have the following files

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.valid.sample**| This file contains the samples that passed all the QC |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.covariate**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the PCs of the samples |

# Running PRS analysis

We can run lassosum as follow

``` R
library(lassosum)
# Prefer to work with data.table as it speeds up file reading
library(data.table)
library(methods)
# We like to use dplyr for it makes codes much more readable
library(dplyr)
sum.stat <- "Height.QC.gz"
ref.bfile <- "EUR.QC"
bfile <- "EUR.QC"
# Read in and process the covariates
covariate <- fread("EUR.covariate")
pcs <- fread("EUR.eigenvec")
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
cov <- merge(covariate, pcs, by=c("FID", "IID"))

# We will need the EUR.hg19 file provided by lassosum 
# which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
ld.file <-  "EUR.hg19" 
# output prefix
prefix <- "EUR"
# Read in the target phenotype file
target.pheno <- fread("EUR.height")[,c("FID", "IID", "Pheno")]
# Read in samples to include in the analysis
target.keep <- fread("EUR.valid.sample")[,c("FID", "IID")]
# Read in the summary statistics
ss <- fread(sum.stat)
# Number of sample in base
size <- 253288
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P == 0]
# Read in the LD blocks
ld <- fread(ld.file)
# Transform the P-values into correlation
cor <- p2cor(p = ss$P,
        n = size,
        sign = log(ss$OR)
        )
# Run the lassosum pipeline
out <- lassosum.pipeline(
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$A1,
    A2 = ss$A2,
    ref.bfile = ref.bfile,
    keep.ref = target.keep,
    test.bfile = bfile,
    keep.test = target.keep,
    LDblocks = ld,
    trace = 2
)
# Store the R2 results
target.res <- validate(out, pheno = target.pheno, covar=cov)
# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
```