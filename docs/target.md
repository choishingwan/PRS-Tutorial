# Obtaining the target data
Target data consist of individual-level genotype-phenotype data, usually generated within your lab/department/collaboration. For this tutorial, we have simulated some genotype-phenotype data using the 1000 Genomes Project European samples. 
You can download the data [here](https://github.com/choishingwan/PRS-Tutorial/raw/master/resources/EUR.zip) or you can download the data using the following command:

```bash
curl https://github.com/choishingwan/PRS-Tutorial/raw/master/resources/EUR.zip -L -O
```

Unzip the data as follow:

```bash
unzip EUR.zip
```

!!! note
    Install the program PLINK and include its location in your PATH directory, which allows us to use `plink` instead of `./plink` in the commands below. If PLINK is not in your PATH directory and is instead in your working directory, then replace all instances of `plink` in the tutorial with `./plink`.

# QC checklist: Target data
Below are the QC steps that comprise the QC checklist for the target data.

# \# Sample size
We recommend that users only perform PRS analyses on target data of at least 100 individuals. The sample size of our target data here is 503 individuals. 

# \# File transfer
Usually we do not need to download and transfer the target data file because it is typically generated locally. However, the file should contain an md5sum code in case we send the data file to collaborators who may want to confirm that the file has not changed during the transfer.

??? note "What is the md5sum code for each of the target files?"

    |File|md5sum|
    |:-:|:-:|
    |**EUR.bed**           |940f5a760b41270662eba6264b262a2d|
    |**EUR.bim**           |a528020cc2448aa04a7499f13bf9f16a|
    |**EUR.covariate**     |afff13f8f9e15815f2237a62b8bec00b|
    |**EUR.fam**           |17e8184fb03c690db6980bb7499d4982|
    |**EUR.height**        |052beb4cae32ac7673f1d6b9e854c85b|

# \# Genome build
As stated in the base data section, the genome build for our base and target data is the same, as it should be.

# \# Standard GWAS QC
The target data must be quality controlled to at least the standards 
implemented in GWAS studies, e.g. removing SNPs with low genotyping rate, 
low minor allele frequency, out of Hardy-Weinberg Equilibrium, removing
individuals with low genotyping rate 
(see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).

The following `plink` command applies some of these QC metrics to the target data:

```bash
plink \
    --bfile EUR \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --write-snplist \
    --make-just-fam \
    --out EUR.QC
```
Each of the parameters corresponds to the following

| Paramter | Value | Description|
|:-:|:-:|:-|
| bfile | EUR | Informs `plink` that the input genotype files should have a prefix of `EUR` |
| maf | 0.05 | Removes all SNPs with minor allele frequency less than 0.05. Genotyping errors typically have a larger influence on SNPs with low MAF. Studies with large sample sizes could apply a lower MAF threshold|
| hwe | 1e-6 | Removes SNPs with low P-value from the Hardy-Weinberg Equilibrium Fisher's exact or chi-squared test. SNPs with significant P-values from the HWE test are more likely affected by genotyping error or the effects of natural selection. Filtering should be performed on the control samples to avoid filtering SNPs that are causal (under selection in cases)|
| geno | 0.01 | Excludes SNPs that are missing in a high fraction of subjects. A two-stage filtering process is usually performed (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).|
| mind | 0.01 | Excludes individuals who have a high rate of genotype missingness, since this may indicate problems in the DNA sample or processing. (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/) for more details).|
| make-just-fam | - | Informs `plink` to only generate the QC'ed sample name to avoid generating the .bed file.  |
| write-snplist | - | Informs `plink` to only generate the QC'ed SNP list to avoid generating the .bed file. |
| out | EUR | Informs `plink` that all output should have a prefix of `EUR` |

??? note "How many SNPs and samples were filtered?"
    - `5` samples were removed due to a high rate of genotype missingness
    - `1` SNP were removed due to missing genotype data
    -  `872` SNPs were removed due to being out of Hardy-Weinberg Equilibrium
    - `242,459` SNPs were removed due to low minor allele frequency

 
!!! note
    Normally, we can generate a new genotype file using the new sample list.
    However,  this will use up a lot of storage space. Using `plink`'s
    `--extract`, `--exclude`, `--keep`, `--remove`, `--make-just-fam` and `--write-snplist` functions, we can work 
    solely on the list of samples and SNPs without duplicating the 
    genotype file, reducing the storage space usage.  
    
Very high or low heterozygosity rates in individuals could be due to DNA contamination or to high levels of inbreeding. Therefore, samples with extreme heterozygosity are typically removed prior to downstream analyses. 

First, we perform prunning to remove highly correlated SNPs:
```bash
plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC
```

This will generate two files 1) **EUR.QC.prune.in** and 2) **EUR.QC.prune.out**.All SNPs within **EUR.QC.prune.in** have a pairwise $r^2 < 0.25$. 


Heterozygosity rates can then be computed using `plink`:
```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --het \
    --out EUR.QC
```

This will generate the **EUR.QC.het** file, which contains F coefficient estimates for assessing heterozygosity.
We will remove individuals with F coefficients that are more than 3 standard deviation (SD) units from the mean, which can be performed using the following `R` command (assuming that you have R downloaded, then you can open an `R` session by typing `R` in your terminal):

```R tab="Without library"
dat <- read.table("EUR.QC.het", header=T) # Read in the EUR.het file, specify it has header
m <- mean(dat$F) # Calculate the mean  
s <- sd(dat$F) # Calculate the SD
valid <- subset(dat, F <= m+3*s & F >= m-3*s) # Get any samples with F coefficient within 3 SD of the population mean
write.table(valid[,c(1,2)], "EUR.valid.sample", quote=F, row.names=F) # print FID and IID for valid samples
q() # exit R
```

```R tab="With data.table"
library(data.table)
# Read in file
dat <- fread("EUR.QC.het")
# Get samples with F coefficient within 3 SD of the population mean
valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 
# print FID and IID for valid samples
fwrite(valid[,c("FID","IID")], "EUR.valid.sample", sep="\t") 
q() # exit R
```

??? note "How many samples were excluded due to high heterozygosity rate?"
    - `7` samples were excluded

# \# Ambiguous SNPs
They were removed during the base data QC


# \# Sex chromosomes 

Sometimes sample mislabelling can occur, which may lead to invalid results. 
A good indication of a mislabelled sample is a mismatch between biological sex and reported sex. 
If the biological sex does not match up with the reported sex, then the sample may have been mislabelled.

Before performing a sex check, pruning should be performed (see [here](target.md#35-standard-gwas-qc)).
A sex check can then easily be conducted using `plink`
```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC
```

This will generate a file called **EUR.QC.sexcheck** containing the F-statistics for each individual. Individuals are typically called as being biologically male if the F-statistic is > 0.8 and biologically female if F < 0.2.

```R tab="Without library"
# Read in file
valid <- read.table("EUR.valid.sample", header=T)
dat <- read.table("EUR.QC.sexcheck", header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], "EUR.QC.valid", row.names=F, col.names=F, sep="\t", quote=F) 
q() # exit R
```

```R tab="With data.table"
library(data.table)
# Read in file
valid <- fread("EUR.valid.sample")
dat <- fread("EUR.QC.sexcheck")[FID%in%valid$FID]
fwrite(dat[STATUS=="OK",c("FID","IID")], "EUR.QC.valid", sep="\t") 
q() # exit R
```

??? note "How many samples were excluded due mismatched Sex information?"
    - `2` samples were excluded


# \# Mismatching genotypes
In addition, when there are non-ambiguous mismatches in allele 
coding between the data sets, such as A/C in the base
and G/T in the target data, then this can be resolved by 
‘flipping’ the alleles in the target data to their complementary alleles. 
This can be achieved with the following steps: 

1\. Load the bim file, the GIANT summary statistic and the QC SNP list into R

```R tab="Without data.table"
# Read in bim file
bim <- read.table("EUR.bim")
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
# Read in QCed SNPs
qc <- read.table("EUR.QC.snplist", header = F, stringsAsFactors = F)
# Read in GIANT data
height <-
    read.table(gzfile("Height.QC.gz"),
               header = T,
               stringsAsFactors = F, 
               sep="\t")
# Change all alleles to upper case for easy comparison
height$A1 <- toupper(height$A1)
height$A2 <- toupper(height$A2)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)
```

```R tab="With data.table"
library(data.table)
# Read in bim file 
bim <- fread("EUR.bim")
setnames(bim, colnames(bim), c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2"))
# Read in GIANT data (require data.table v1.12.0+)
height <- fread("Height.QC.gz")
# Change all alleles to upper case for easy comparison
height[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
bim[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]
# Read in QCed SNPs
qc <- fread("EUR.QC.snplist", header=F)
```


2\. Identify SNPs that require strand flipping 

```R tab="Without data.table"
# Merge GIANT with target
info <- merge(bim, height, by = c("SNP", "CHR", "BP"))
# Filter QCed SNPs
info <- info[info$SNP %in% qc$V1,]
# Function for finding the complementary allele
complement <- function(x) {
    switch (
        x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}
# Get SNPs that have the same alleles across base and target
info.match <- subset(info, A1 == B.A1 & A2 == B.A2)
# Identify SNPs that are complementary between base and target
info$C.A1 <- sapply(info$B.A1, complement)
info$C.A2 <- sapply(info$B.A2, complement)
info.complement <- subset(info, A1 == C.A1 & A2 == C.A2)
# Update these allele coding in the bim file
bim[bim$SNP %in% info.complement$SNP,]$B.A1 <-
    sapply(bim[bim$SNP %in% info.complement$SNP,]$B.A1, complement)
bim[bim$SNP %in% info.complement$SNP,]$B.A2 <-
    sapply(bim[bim$SNP %in% info.complement$SNP,]$B.A2, complement)
```

```R tab="With data.table"
# Merge GIANT with target
info <- merge(bim, height, by=c("SNP", "CHR", "BP"))
info <- info[SNP %in% qc$V1]
# Function for calculating the complementary allele
complement <- function(x){
    switch (x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}
# Identify SNPs that are complementary between base and target
com.snps <- info[sapply(B.A1, complement) == A1 &
                     sapply(B.A2, complement) == A2, SNP]
# Now update the bim file
bim[SNP %in% com.snps, c("B.A1", "B.A2") :=
        list(sapply(B.A1, complement),
             sapply(B.A2, complement))]
# And update the info structure
info[SNP %in% com.snps, c("B.A1", "B.A2") :=
        list(sapply(B.A1, complement),
             sapply(B.A2, complement))]
```


3\. Identify SNPs that require recoding in the target (to ensure the coding allele in the target data is the effective allele in the base summary statistic)

```R tab="Without data.table"
# identify SNPs that need recoding
info.recode <- subset(info, A1 == B.A2 & A2 == B.A1)
# identify SNPs that need recoding & complement
info.crecode <- subset(info, A1 == C.A2 & A2 == C.A1)
# Update these allele coding in the bim file
com.snps <- bim$SNP %in% info.crecode$SNP
tmp <- bim[com.snps,]$B.A1
bim[com.snps,]$B.A1 <- as.character(sapply(bim[com.snps,]$B.A2, complement))
bim[com.snps,]$B.A2 <- as.character(sapply(tmp, complement))
# Output updated bim file
write.table(
    bim,
    "EUR.QC.adj.bim",
    quote = F,
    row.names = F,
    col.names = F
)
```

```R tab="With data.table"
# identify SNPs that need recoding & complement
com.recode <- info[sapply(B.A1, complement) == A2 &
                     sapply(B.A2, complement) == A1, SNP]
# Now update the bim file
bim[SNP %in% com.recode, c("B.A1", "B.A2") :=
        list(sapply(B.A2, complement),
             sapply(B.A1, complement))]
# And update the info structure
info[SNP %in% com.recode, c("B.A1", "B.A2") :=
        list(sapply(B.A2, complement),
             sapply(B.A1, complement))]
# Write the updated bim file
fwrite(bim, "EUR.QC.adj.bim", col.names=F, sep="\t")
```


4\. Identify SNPs that have different allele in base and target (usually due to difference in genome build or Indel)

```R tab="Without data.table"
mismatch <-
    bim$SNP[!(bim$SNP %in% info.match$SNP |
                  bim$SNP %in% info.complement$SNP | 
                  bim$SNP %in% info.recode$SNP |
                  bim$SNP %in% info.crecode$SNP)]
write.table(
    mismatch,
    "EUR.mismatch",
    quote = F,
    row.names = F,
    col.names = F
)
q() # exit R
```

``` R tab="With data.table"
matched <- info[(A1 == B.A1 & A2 == B.A2) |
                    (A1 == B.A2 & A2 == B.A1)]
mismatch <- bim[!SNP%in%matched$SNP, SNP]
write.table(mismatch, "EUR.mismatch", quote=F, row.names=F, col.names=F)
q() # exit R
```


5\. Replace **EUR.bim** with **EUR.QC.adj.bim**:

```bash
# Make a back up
mv EUR.bim EUR.bim.bk
ln -s EUR.QC.adj.bim EUR.bim
```

!!! note
    Most PRS software will perform flipping automatically, thus this step is usually not required.


# \# Duplicate SNPs
Make sure to remove any duplicate SNPs in your target data (these target data were simulated and so include no duplicated SNPs)


# \# Sample overlap
Since the target data were simulated there are no overlapping samples between the base and target data here (see the relevant section of [the paper](https://doi.org/10.1101/416545) for discussion of the importance of avoiding sample overlap). 

# \# Relatedness
Closely related individuals in the target data may lead to overfitted results, limiting the generalisability of the results. 

Before calculating the relatedness, pruning should be performed (see [here](target.md#35-standard-gwas-qc)).
Individuals that have a first or second degree relative in the sample ($\hat{\pi} > 0.125$) can be removed with the following command:

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.valid \
    --rel-cutoff 0.125 \
    --out EUR.QC
```

??? note "How many related samples were excluded?"
    - `2` samples were excluded

!!! note
    A greedy algorithm is used to remove closely related individuals in a way that optimises the size of the sample retained.                However, the algorithm is dependent on the random seed used, which can generate different results. Therefore, to reproduce
    the same result, you will need to specify the same random seed. 
    
    PLINK's algorithm for removing related individuals does not account for the phenotype under study. 
    To minimize the removal of cases of a disease, the following algorithm can be used instead: 
    [GreedyRelated](https://github.com/choishingwan/GreedyRelated).

# Generate final QC'ed target data file
After performing the full analysis, you can generate a QC'ed data set with the following command:
```bash
plink \
    --bfile EUR \
    --make-bed \
    --keep EUR.QC.rel.id \
    --out EUR.QC \
    --extract EUR.QC.snplist
```
