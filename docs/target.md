In this section we will perform some quality control (QC) steps on the target data. 

For this tutorial, we have simulated some genotype-phenotype data using the 1000 Genomes Project European samples. 
You can download the data [here](https://github.com/choishingwan/PRS-Tutorial/raw/master/resources/EUR.zip) or you can download the data using the following command:

```bash
curl https://github.com/choishingwan/PRS-Tutorial/raw/master/resources/EUR.zip -L -O
```

Unzip the data as follow:

```bash
unzip EUR.zip
```

??? note "What is the md5sum code for each of the target files?"

    |File|md5sum|
    |:-:|:-:|
    |**EUR.bed**           |940f5a760b41270662eba6264b262a2d|
    |**EUR.bim**           |a528020cc2448aa04a7499f13bf9f16a|
    |**EUR.covariate**     |afff13f8f9e15815f2237a62b8bec00b|
    |**EUR.fam**           |17e8184fb03c690db6980bb7499d4982|
    |**EUR.height**        |052beb4cae32ac7673f1d6b9e854c85b|

!!! note
    Install the program PLINK and include its location in your PATH directory, which allows us to use `plink` instead of `./plink` in the commands below. If PLINK is not in your PATH directory and is instead in your working directory, then replace all instances of `plink` in the tutorial with `./plink`.

# Genotype file format

# Basic filterings
The power and validity of PRS analyses depend on 
the quality of the base and target data. Therefore, 
both data sets must be quality controlled to at least the standards 
implemented in GWAS studies, e.g. removing SNPs with low genotyping rate, 
low minor allele frequency, out of Hardy-Weinberg Equilibrium and removing
individuals with low genotyping rate 
(see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).

The following `plink` command applies some basic filtering:

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

# Filter related samples
Closely related individuals in the target data may lead to overfitted results, limiting the generalisability of the results. 

As a first step to removing related individuals we perform pruning, which removes highly correlated SNPs:
```bash
plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC
```

This will generate two files 1) **EUR.QC.prune.in** and 2) **EUR.QC.prune.out**
All SNPs within **EUR.QC.prune.in** have a pairwise $r^2 < 0.25$

Individuals that have a first or second degree relative in the sample ($\text{pi-hat} > 0.125$) can be removed with the following command:

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
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
    
# Remove samples with extreme heterozygosity rate
Very high or low heterozygosity rates in individuals could be due to DNA contamination or to high levels of inbreeding. Therefore, samples with extreme heterozygosity are typically removed prior to downstream analyses. Heterozygosity rates can be computed using `plink` after performing pruning.
```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.rel.id \
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
```

```R tab="With data.table"
library(data.table)
# Read in file
dat <- fread("EUR.QC.het")
# Get samples with F coefficient within 3 SD of the population mean
valid <- dat[F<=mean(F)+3*sd(F) & F>=mean(F)-3*sd(F)] 
# print FID and IID for valid samples
fwrite(valid[,c("FID","IID")], "EUR.valid.sample", sep="\t") 
```

??? note "How many samples were excluded due to high heterozygosity rate?"
    - `7` samples were excluded

# Check for mismatching sex information

Sometimes sample mislabelling can occur, which may lead to invalid results. 
A good indication of a mislabelled sample is a mismatch between biological sex and reported sex. 
If the biological sex does not match up with the reported sex, then the sample may have been mislabelled.

Before performing a sex check, pruning should be performed (see [here](target.md#filter-related-samples)).
A sex check can then easily be conducted using `plink`
```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC
```

This will generate a file called **EUR.sexcheck** containing the F-statistics for each individual. Individuals are typically called as being biologically male if the F-statistic is > 0.8 and biologically female if F < 0.2.

```R tab="Without library"
# Read in file
valid <- read.table("EUR.valid.sample", header=T)
dat <- read.table("EUR.QC.sexcheck", header=T)
valid <- subset(dat, STATUS=="OK" & FID %in% valid$FID)
write.table(valid[,c("FID", "IID")], "EUR.QC.valid", row.names=F, col.names=F, sep="\t") 
```

```R tab="With data.table"
library(data.table)
# Read in file
valid <- fread("EUR.valid.sample")
dat <- fread("EUR.QC.sexcheck")[FID%in%valid$FID]
fwrite(dat[STATUS=="OK",c("FID","IID")], "EUR.QC.valid", sep="\t") 
```

??? note "How many samples were excluded due mismatched Sex information?"
    - `2` samples were excluded

# Generate final QC'ed target sample
After performing the full analysis, you can generate a QC'ed data set with the following command:
```bash
plink \
    --bfile EUR \
    --make-bed \
    --out EUR.QC \
    --keep EUR.QC.valid \
    --extract EUR.QC.snplist
```

!!! note
    For some software, the **EUR.QC.valid** and **EUR.QC.snplist** can be passed as a parameter to perform the 
    extraction directly. For those software (e.g. PRSice-2, lassosum, etc), this step is not required.
