Here we use another PRS program, `lassosum`, which is an `R` package that uses penalised regression (LASSO) in its approach to PRS calculation.

!!! note
    The script used here is based on lassosum version 0.4.4


You can install `lassosum` and its dependencies in `R` with the following command:

```R
install.packages(c("devtools","RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
library(devtools)
install_github("tshmak/lassosum")
```

Again, we assume that we have the following files: 

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.covariate**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the PCs of the samples |

# Running PRS analysis

We can run lassosum as follows: 

``` R
library(lassosum)
# Prefer to work with data.table as it speeds up file reading
library(data.table)
library(methods)
library(magrittr)
sum.stat <- "Height.QC.gz"
bfile <- "EUR.QC"
# Read in and process the covariates
covariate <- fread("EUR.covariate")
pcs <- fread("EUR.eigenvec") %>%
    setnames(., colnames(.), c("FID","IID", paste0("PC",1:6)))
# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well
cov <- merge(covariate, pcs)

# We will need the EUR.hg19 file provided by lassosum 
# which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
ld.file <- system.file("data", "Berisa.EUR.hg19.bed",package="lassosum")
# output prefix
prefix <- "EUR"
# Read in the target phenotype file
target.pheno <- fread("EUR.height")[,c("FID", "IID", "Height")]
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
# Because FID of our samples are all 0, we might encounter problem with lassosum
# we need to provide a T/F vector instead of the target.keep file
target.keep[, ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]

keep <- fam$ID %in% target.keep$ID
# Run the lassosum pipeline
out <- lassosum.pipeline(
    cor = cor,
    chr = ss$CHR,
    pos = ss$BP,
    A1 = ss$A1,
    A2 = ss$A2,
    ref.bfile = bfile,
    keep.ref = keep,
    test.bfile = bfile,
    keep.test = keep,
    LDblocks = ld
)
# Store the R2 results
target.res <- validate(out, pheno = target.pheno, covar=cov)
# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2
```


??? note "How much phenotypic variation does the "best-fit" PRS explain?"
    0.03818004