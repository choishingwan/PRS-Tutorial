Next, we'd like to perform basic quality controls (QC) on the target genotype data. 

In this tutorial, we've simulated some samples using the 1000 genome european genotypes. 
You can download the data [here](https://github.com/choishingwan/PRS-Tutorial/raw/master/resources/EUR.zip). 

Or you can download using the following script:
```bash
curl https://github.com/choishingwan/PRS-Tutorial/raw/master/resources/EUR.zip -L -O
```

Unzip the data as follow:

```bash
unzip EUR.zip
```

??? note "What's the md5sum of the genotype files?"

    |File|md5sum|
    |:-:|:-:|
    |**EUR.bed**           |7d163129c79277ec9858008f10cac0e5|
    |**EUR.bim**           |62b278b3338fc86b89ecc8d622731701|
    |**EUR.covariate**     |afff13f8f9e15815f2237a62b8bec00b|
    |**EUR.fam**           |0189b2f82b5d20f63ab8f667e2feb100|
    |**EUR.height**        |052beb4cae32ac7673f1d6b9e854c85b|

!!! note
    We assume PLINK is installed in your PATH directory, which allow us to use `plink` instead of `./plink`.
    If PLINK is not in your PATH directory, replace all instance of `plink` in the tutorial to `./plink` assuming
    the PLINK executable is located within your working directory

# Genotype file format

# Basic filterings
The power and validity of PRS analyses are highly dependent on 
the quality of the base and target data, therefore 
both data sets must be quality controlled to the high standards 
implemented in GWAS studies, e.g. removing SNPs with low genotyping rate, 
low minor allele frequency, violates the Hardy-Weinberg Equilibrium and
individuals with low genotyping rate 
(see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).

The following `plink` command perform some basic filterings

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
| bfile | EUR | Inform `plink` that the input genotype files should have a prefix of `EUR` |
| maf | 0.05 | Filter out any SNPs with minor allele frequency less than 0.05. Genotype error are more likely to influence SNPs with low MAF. Large sample size can adapt a lower MAF threshold|
| hwe | 1e-6 | Filtering SNPs with low p-value from the Hardy-Weinberg exact test. SNPs with significant p-value from the HWE test are more likely to harbor genotyping error or are under selection. Filtering should be performed on the control samples to avoid filtering SNPs that are causal (under selection in cases)|
| geno | 0.01 | Exclude SNPs that are missing in large proportion of subjects. A two pass filtering is usually performed (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/)).|
| mind | 0.01 | Exclude individual who have a high rate of genotype missingness. This might indicate problems in the DNA sample. (see [Marees et al](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6001694/) for more details).|
| make-just-fam | - | Inform `plink` to only generate the QCed sample name to avoid generating the .bed file.  |
| write-snplist | - | Inform `plink` to only generate the QCed SNP list to avoid generating the .bed file. |
| out | EUR | Inform `plink` that all output should have a prefix of `EUR` |

??? note "How many SNPs were filtered?"
    A total of `5359` SNPs were removed due to Hardy-Weinberg exact test results
    and `241,485` SNPs were removed due to minor allele frequency


!!! note
    Normally, we can generate a new genotype file using the new sample list.
    However,  this will use up a lot of storage space. Using `plink`'s
    `--extract`, `--exclude`, `--keep`, `--remove`, `--make-just-fam` and `--write-snplist` functions, we can work 
    solely on the list of samples and SNPs without duplicating the 
    genotype file, therefore reducing the storage space usage.  

# Filter related samples
Related samples in the target data might lead to overfitted results, 
hampering the generalizability of the results. 

To remove related samples, we first need to perform prunning to remove highly correlated SNPs:
```bash
plink \
    --bfile EUR \
    --keep EUR.QC.fam \
    --extract EUR.QC.snplist \
    --indep-pairwise 200 50 0.25 \
    --out EUR.QC
```

This will generate two files 1) **EUR.QC.prune.in** and 2) **EUR.QC.prune.out**
All SNPs within **EUR.QC.prune.in** has a pairwise $r^2 < 0.25$

Samples with more than third-degree relatedness ($\text{pi-hat} > 0.125$) can then be removed with 

```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.fam \
    --rel-cutoff 0.125 \
    --out EUR.QC
```

!!! note
    A greedy algorithm is used to remove the related samples. Which depending
    on the random seed used, might generate different results. To reproduce
    the same result, you might need to specify the random seed usage. 
    
    PLINK's related sample removal does not take into account of the sample 
    phenotype. If one would like to minimize lost of cases for example, 
    a software called
    [GreedyRelated](https://github.com/choishingwan/GreedyRelated) can be used.
    
# Remove samples with abnormal heterozygosity rate
Individual with high or low heterozygosity rate can be contaminated or are inbreed.
It is therefore a good idea to remove these samples from our dataset before continuing the analyse.
Heterozygosity rate can be calculated using `plink` after performing prunning. 
```bash
plink \
    --bfile EUR.QC \
    --extract EUR.QC.prune.in \
    --keep EUR.QC.rel.id \
    --het \
    --out EUR.QC
```
This will generate the **EUR.QC.het** file which contains the F coefficient estimates.
It will be easier to filter the samples using `R` instead of `awk`:
Open a `R` section by tying `R` in your terminal

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

# Check for mis-matched Sex information
Sometimes, sample mislabeling can occur and will lead to invalid results. 
A good indication of mislabeled sample is a mismatch between the biological sex and the reported sex. 
If the biological sex does not match up with the reported sex, it is likely that the sample has been mislabeled.

Before performing sex check, prunning should be performed (see [here](target.md#filter-related-samples)).
Sex check can then easily be carried out using `plink`
```bash
plink \
    --bfile EUR \
    --extract EUR.QC.prune.in \
    --keep EUR.valid.sample \
    --check-sex \
    --out EUR.QC
```

This will generate a file called **EUR.sexcheck** containing the F-statistics for each individual.
For male, the F-statistic should be > 0.8 and Female should have a value < 0.2.

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

# Generate final QCed sample
After performing the full analysis, you can generate a QCed data set with the following command
```bash
plink \
    --make-bed \
    --out EUR.QC \
    --keep EUR.QC.valid \
    --extract EUR.QC.snplist
```

!!! note
    For some software, the **EUR.QC.valid** and **EUR.QC.snplist** can be passed as a parameter to perform the 
    extraction directly. For those software (e.g. PRSice-2, lassosum, etc), this step is not required