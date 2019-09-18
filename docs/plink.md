# Background
In this section of the tutorial you will use four different software programs to compute PRS from the base and target data that you QC'ed in the previous two sections. On this page, you will compute PRS using the popular genetic analyses tool `plink` - while `plink` is not a dedicated PRS software, each of the steps required to compute PRS using the C+T standard approach can be performed in `plink` and carrying out this multi-step process can be a good way to learn the processes involved in computing PRS (which are typically performed automatically by PRS software). 

# Required Data

In the previous sections, we have generated the following files:

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.valid.sample**| This file contains the samples that passed all the QC |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.covariate**| This file contains the covariates of the samples |


# Remove Ambiguous SNPs
If the base and target data were generated using different genotyping chips and the chromosome strand (+/-) for either is unknown, then it is not possible to match ambiguous SNPs (i.e. those with complementary alleles, either C/G or A/T) across the data sets, because it will be unknown whether the base and target data are referring to the same allele or not. 

Ambiguous SNPs can be obtained by examining the bim file:
```bash
awk '!( ($5=="A" && $6=="T") || \
        ($5=="T" && $6=="A") || \
        ($5=="G" && $6=="C") || \
        ($5=="C" && $6=="G")) {print}' \
        EUR.QC.bim > EUR.unambig.snp 
```

??? note "How many ambiguous SNPs were there?"
    There are `330,818` ambiguous SNPs

# Strand Flipping
In addition, when there are non-ambiguous mismatches in allele 
coding between the data sets, such as A/C in the base
and G/T in the target data, then this can be resolved by 
‘flipping’ the alleles in the target data to their complementary alleles. 
This can be achieved with the following steps: 

1. Get the correct A1 alleles for the bim file

```R tab="Without data.table"
bim <- read.table("EUR.QC.bim", header = F, stringsAsFactors = F)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
height <-
    read.table(gzfile("Height.QC.gz"),
               header = T,
               stringsAsFactors = F)
# Change all alleles to upper case for easy comparison
height$A1 <- toupper(height$A1)
height$A2 <- toupper(height$A2)
bim$B.A1 <- toupper(bim$B.A1)
bim$B.A2 <- toupper(bim$B.A2)
info <- merge(bim, height, by = c("SNP", "CHR", "BP"))

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
# identify SNPs that need flipping
info.flip <- subset(info, A1 == B.A2 & A2 == B.A1)
# identify SNPs that need flipping & complement
info.cflip <- subset(info, A1 == C.A2 & A2 == C.A1)
# Update these allele coding in the bim file
com.snps <- bim$SNP %in% info.cflip$SNP
bim[com.snps,]$B.A1 <- sapply(bim[com.snps,]$B.A1, complement)
bim[com.snps,]$B.A2 <- sapply(bim[com.snps,]$B.A2, complement)
# Get list of SNPs that need to change the A1 encoding
flip <- rbind(info.flip, info.cflip)
flip.snp <- data.frame(SNP = flip$SNP, A1 = flip$A1)
write.table(flip.snp,
            "EUR.update.a1",
            quote = F,
            row.names = F)
write.table(
    bim,
    "EUR.QC.adj.bim",
    quote = F,
    row.names = F,
    col.names = F
)
# And we want to remove any SNPs that do not match with the base data
mismatch <-
    bim$SNP[!(bim$SNP %in% info.match$SNP |
                  bim$SNP %in% info.complement$SNP | 
                  bim$SNP %in% flip$SNP)]
write.table(
    mismatch,
    "EUR.mismatch",
    quote = F,
    row.names = F,
    col.names = F
)
```

```R tab="With data.table"
library(data.table)
bim <- fread("EUR.QC.bim")
bim.col <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
setnames(bim, colnames(bim), bim.col)
height <- fread("Height.QC.gz")
# Change all alleles to upper case for easy comparison
height[,c("A1","A2"):=list(toupper(A1), toupper(A2))]
bim[,c("B.A1","B.A2"):=list(toupper(B.A1), toupper(B.A2))]
info <- merge(bim, height, by=c("SNP", "CHR", "BP"))
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
# identify SNPs that need flipping & complement
com.flip <- info[sapply(B.A1, complement) == A2 &
                     sapply(B.A2, complement) == A1, SNP]
# Now update the bim file
bim[SNP %in% com.flip, c("B.A1", "B.A2") :=
        list(sapply(B.A1, complement),
             sapply(B.A2, complement))]
# Obtain list of SNPs that require flipping
flip <- info[B.A1==A2 & B.A2==A1]
# Now generate file for PLINK 
fwrite(flip[,c("SNP", "A1")], "EUR.update.a1", sep="\t")
# Write the updated bim file
fwrite(bim, "EUR.QC.adj.bim", col.names=F, sep="\t")
# We can then remove all mismatch SNPs
matched <- info[(A1 == B.A1 & A2 == B.A2) |
                    (A1 == B.A2 & A2 == B.A1) |
                    (A1 == sapply(B.A1, complement) &
                         A2 == sapply(B.A2, complement)) |
                    (A1 == sapply(B.A2, complement) &
                         A2 == sapply(B.A1, complement))]
mismatch <- bim[!SNP%in%matched$SNP, SNP]
write.table(mismatch, "EUR.mismatch", quote=F, row.names=F, col.names=F)
```

The above script will generate three files: **EUR.QC.adj.bim**, **EUR.update.a1** and **EUR.mismatch**. 

2. Replace **EUR.QC.bim** with **EUR.QC.adj.bim**:

```bash
# Make a back up
mv EUR.QC.bim EUR.QC.bim.bk
ln -s EUR.QC.adj.bim EUR.QC.bim
```

3. Generate a new genotype file with the alleles flipped so that the alleles in the base and target data match:
```bash
plink \
    --bfile EUR.QC \
    --a1-allele EUR.update.a1 \
    --make-bed \
    --keep EUR.valid.sample \
    --extract EUR.unambig.snp \
    --exclude EUR.mismatch \
    --out EUR.QC.flipped
```

# Update Effect Size
When the effect size relates to disease risk and is thus given as an odd ratio (OR), rather than BETA (for continuous traits), then the PRS is computed as a product of ORs. To simplify this calculation, we usually take the natural logarithm of the OR so that the PRS can be computed using a simple summation instead (which can be back-transformed afterwards). 
We can obtain the transformed summary statistics with `R`:

```R tab="Without data.table"
dat <- read.table(gzfile("Height.QC.gz"), header=T)
dat$OR <- log(dat$OR)
write.table(dat, "Height.QC.Transformed", quote=F, row.names=F)
```

```R tab="With data.table"
library(data.table)
dat <- fread("Height.QC.gz")
fwrite(dat[,OR:=log(OR)], "Height.QC.Transformed", sep="\t")
```


!!! warning
    It may be tempting to perform the log transofrmation using `awk`.
    However, due to rounding of values performed in `awk`, less accurate results
    may be obtained. Therefore, we recommend performing the transformation in `R` or allow the PRS software to perform the transformation directly.

# Clumping
Linkage disequilibrium, which corresponds to the correlation between the genotypes of genetic variants across the genome, makes identifying the contribution from causal independent genetic variants extremely challenging. One way of approximately capturing the right level of causal signal is to perform clumping, which removes SNPs in such a way that only weakly correlated SNPs are retained but preferentially retaining the SNPs most associated with the phenotype under study. Clumping can be performed using the following command in `plink`: 

```bash
plink \
    --bfile EUR.QC.flipped \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Height.QC.transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out EUR
```
Each of the new parameters corresponds to the following

| Paramter | Value | Description|
|:-:|:-:|:-|
| clump-p1 | 1 | P-value threshold for a SNP to be included as an index SNP. 1 is selected such that all SNPs are include for clumping|
| clump-r2 | 0.1 | SNPs having $r^2$ higher than 0.1 with the index SNPs will be removed |
| clump-kb | 250 | SNPs within 250k of the index SNP are considered for clumping|
| clump | Height.QC.transformed | Summary statistic file containing the p-value information|
| clump-snp-field | SNP | Specify that the column `SNP` contains the SNP IDs |
| clump-field | P | Specify that the column `P` contains the P-value information |

A more detailed description of the clumping process can be found [here](https://www.cog-genomics.org/plink/1.9/postproc#clump)

!!! note
    The $r^2$ values computed by `--clump` are based on maximum likelihood haplotype frequency estimates


This will generate **EUR.clumped**, containing the index SNPs after clumping is performed.
We can extract the index SNP ID by performing the following:

```bash
awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp
```

> `$3` because the third column contains the SNP ID


!!! note
    If your target sample is small (e.g. <500), you can use the 1000 Genomes Project samples for the LD calculation.
    Make sure that you use the population that most closely reflects represents the base sample.

# Generate PRS
`plink` provides a convenient function `--score` and `--q-score-range` for calculating polygenic score.

We will need three files:

1. The summary statistic file: **Height.QC.Transformed**
2. A file containing SNP ID and their corresponding p-value (`$1` because SNP ID is located at the first column; `$8` because P-value is located at the eigth column)
```bash
awk '{print $1,$8}' Height.QC.Transformed > SNP.pvalue
```
3. A file containing the P-value thresholds that we are testing. Here we will only test a few thresholds for illustration purposes:
```bash
echo "0.001 0 0.001" > range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list
```
The format of the **range_list** file should be as follows:

|Name of Threshold|Lower bound| Upper Bound|
|:-:|:-:|:-:|

!!! note
    The threshold boundaries are inclusive. For example, for the `0.05` threshold, we include all SNPs with P-value from 
    `0` to `0.05`, **including** any SNPs with P-value equal to `0.05`

We can then calculate the PRS with the following `plink` command:

```bash
plink \
    --bfile EUR.QC.flipped \
    --extract EUR.valid.snp \
    --score Height.QC.Transformed 1 4 11 header \
    --q-score-range range_list SNP.pvalue \
    --out EUR
```
The meaning of the new parameters are as follows:

| Paramter | Value | Description|
|:-:|:-:|:-|
|score|Height.QC.Transformed 1 4 11 header| We read from the **Height.QC.Transformed** file, assuming the `1`st column to be the SNP ID; `4`th column to be the effective allele information; `11`th column to be the effect size estimate; and the file contains a `header`|
|q-score-range| range_test SNP.pvalue| We want to calculate PRS based on the thresholds defined in **range_test**, where the threshold values (p-values) were stored in **SNP.pvalue**|

The above command and range_list will generate 7 files:

1. EUR.0.5.profile
2. EUR.0.4.profile
3. EUR.0.3.profile
4. EUR.0.2.profile
5. EUR.0.1.profile
6. EUR.0.05.profile
7. EUR.0.001.profile

!!! Note
    The default formular for PRS calculation in PLINK is:
    (Assuming the effect size of SNP $i$ is $S_i$;  the number of effective allele observed in sample $j$ is $G_{ij}$; the ploidy of the sample is $P$ (It should be 2 for human); the number of samples included in the PRS be $N$; and the number of non-missing SNPs observed in sample $j$ be $M_j$)
    $$
    PRS_j =\frac{ \sum_i^NS_i*G_{ij}}{P*M_j}
    $$

    If sample has a missing genotype for SNP $i$, the population minor allele frequency times ploidy ($MAF_i*P$) is used inplace of $G_{ij}$

# Accounting for Population Stratification

Population structure is the principal source of confounding in GWAS and are usually accounted for by incorporating the principal components (PCs) as a covariate. 
Similarly, we can incorporate PCs in our PRS analysis to account for population stratification.

Again, we can calculate the PCs using `plink` 
```bash
# First, we need to perform prunning
plink \
    --bfile EUR.QC.flipped \
    --extract EUR.valid.snp \
    --indep-pairwise 200 50 0.25 \
    --out EUR
# Then we calculate the first 6 PCs
plink \
    --bfile EUR.QC.flipped \
    --extract EUR.prune.in \
    --pca 6 \
    --out EUR
```

!!! note
    One way to select the appropriate number of PCs is to perform GWAS on the trait of interest with different number of PCs.
    [LDSC](https://github.com/bulik/ldsc) analysis can then be performed on each of the resulted GWAS summary statistics. 
    By observing the estimated intercept, one can select the number of PCs that provide an intercept estimate closer to 1, which might
    suggest a smaller influence of population stratification.

The eigen-vector (PCs) are stored in **EUR.eigenvec** and can be used as a covariate in the regression model to account for population stratification.

!!! important
    If the base and target samples are collected from different population (e.g. Caucasian vs African ), the results from PRS analysis will be biased (see [Martin et al](https://www.ncbi.nlm.nih.gov/pubmed/28366442)).


# Finding the "Best" P-value threshold
The "best" p-value threshold for PRS construction are usually not known. 
To identify the "best" PRS, we can perform a regression between the calculated PRS and the 
sample phenotype and select the PRS that explains most of the phenotypic variation. 
This can be achieved using `R`.

```R tab="detail"
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
# Read in the phenotype file 
phenotype <- read.table("EUR.height", header=T)
# Read in the PCs
pcs <- read.table("EUR.eigenvec", header=F)
# The default output from plink does not include a header
# To make things simple, we will add the appropriate headers
# (1:6 because there are 6 PCs)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
# Read in the covariates (here, it is sex)
covariate <- read.table("EUR.covariate", header=T)
# Now merge the files
pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
# We can then calculate the null model (model with PRS) using a linear regression 
# (as height is quantitative)
null.model <- lm(Height~., data=pheno[,!colnames(pheno)%in%c("FID","IID")])
# And the R2 of the null model is 
null.r2 <- summary(null.model)$r.squared
prs.result <- NULL
for(i in p.threshold){
    # Go through each p-value threshold
    prs <- read.table(paste0("EUR.",i,".profile"), header=T)
    # Merge the prs with the phenotype matrix
    # We only want the FID, IID and PRS from the PRS file, therefore we only select the 
    # relevant columns
    pheno.prs <- merge(pheno, prs[,c("FID","IID", "SCORE")], by=c("FID", "IID"))
    # Now perform a linear regression on Height with PRS and the covariates
    # ignoring the FID and IID from our model
    model <- lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")])
    # model R2 is obtained as 
    model.r2 <- summary(model)$r.squared
    # R2 of PRS is simply calculated as the model R2 minus the null R2
    prs.r2 <- model.r2-null.r2
    # We can also obtain the coeffcient and p-value of association of PRS as follow
    prs.coef <- summary(model)$coeff["SCORE",]
    prs.beta <- as.numeric(prs.coef[1])
    prs.se <- as.numeric(prs.coef[2])
    prs.p <- as.numeric(prs.coef[4])
    # We can then store the results
    prs.result <- rbind(prs.result, data.frame(Threshold=i, R2=prs.r2, P=prs.p, BETA=prs.beta,SE=prs.se))
}
# Best result is:
prs.result[which.max(prs.result$R2),]
```

```R tab="quick"
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
phenotype <- read.table("EUR.height", header=T)
pcs <- read.table("EUR.eigenvec", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
covariate <- read.table("EUR.covariate", header=T)
pheno <- merge(merge(phenotype, covariate, by=c("FID", "IID")), pcs, by=c("FID","IID"))
null.r2 <- summary(lm(Height~., data=pheno[,!colnames(pheno)%in%c("FID","IID")]))$r.squared
prs.result <- NULL
for(i in p.threshold){
    pheno.prs <- merge(pheno, 
                        read.table(paste0("EUR.",i,".profile"), header=T)[,c("FID","IID", "SCORE")],
                        by=c("FID", "IID"))
    model <- summary(lm(Height~., data=pheno.prs[,!colnames(pheno.prs)%in%c("FID","IID")]))
    model.r2 <- model$r.squared
    prs.r2 <- model.r2-null.r2
    prs.coef <- model$coeff["SCORE",]
    prs.result <- rbind(prs.result, 
        data.frame(Threshold=i, R2=prs.r2, 
                    P=as.numeric(prs.coef[4]), 
                    BETA=as.numeric(prs.coef[1]),
                    SE=as.numeric(prs.coef[2])))
}
print(prs.result[which.max(prs.result$R2),])
```

??? note "Which p-value threshold generate the "best" PRS?"
    0.05

??? note "How much phenotypic variation does the "best" PRS explains?"
    0.03921047
