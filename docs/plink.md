# Background
On this page, you will compute PRS using the popular genetic analyses tool `plink` - while `plink` is not a dedicated PRS software, you can perform every required steps of the C+T approach with `plink`. 
This multi-step process is a good way to learn the processes involved in computing PRS, which are typically performed automatically by PRS software.

## Required Data

In the previous sections, we have generated the following files:

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.cov**| This file contains the covariates of the samples |


## Update Effect Size
When the effect size relates to disease risk and is thus given as an odds ratio (OR), rather than BETA (for continuous traits), then the PRS is computed as a product of ORs. To simplify this calculation, we take the natural logarithm of the OR so that the PRS can be computed using summation instead (which can be back-transformed afterwards). 
We can obtain the transformed summary statistics with `R`:

=== "Without data.table"

    ```R
    dat <- read.table(gzfile("Height.QC.gz"), header=T)
    dat$BETA <- log(dat$OR)
    write.table(dat, "Height.QC.Transformed", quote=F, row.names=F)
    q() # exit R
    ```

=== "With data.table"

    ```R
    library(data.table)
    dat <- fread("Height.QC.gz")
    fwrite(dat[,BETA:=log(OR)], "Height.QC.Transformed", sep="\t")
    q() # exit R
    ```


!!! warning
    Due to rounding of values, using `awk` to log transform OR can lead to less accurate results. Therefore, we recommend performing the transformation in `R` or allow the PRS software to perform the transformation directly.

## Clumping
Linkage disequilibrium, which corresponds to the correlation between the genotypes of genetic variants across the genome, makes identifying the contribution from causal independent genetic variants extremely challenging. 
One way of approximately capturing the right level of causal signal is to perform clumping, which removes SNPs in ways that only weakly correlated SNPs are retained but preferentially retaining the SNPs most associated with the phenotype under study. 
Clumping can be performed using the following command in `plink`: 

```bash
plink \
    --bfile EUR.QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Height.QC.Transformed \
    --clump-snp-field SNP \
    --clump-field P \
    --out EUR
```

Each of the new parameters corresponds to the following

| Parameter | Value | Description|
|:-:|:-:|:-|
| clump-p1 | 1 | P-value threshold for a SNP to be included as an index SNP. 1 is selected such that all SNPs are include for clumping|
| clump-r2 | 0.1 | SNPs having $r^2$ higher than 0.1 with the index SNPs will be removed |
| clump-kb | 250 | SNPs within 250k of the index SNP are considered for clumping|
| clump | Height.QC.Transformed | Base data (summary statistic) file containing the P-value information|
| clump-snp-field | SNP | Specifies that the column `SNP` contains the SNP IDs |
| clump-field | P | Specifies that the column `P` contains the P-value information |

A more detailed description of the clumping process can be found [here](https://www.cog-genomics.org/plink/1.9/postproc#clump)

!!! note
    The $r^2$ values computed by `--clump` are based on maximum likelihood haplotype frequency estimates


This will generate **EUR.clumped**, containing the index SNPs after clumping is performed.
We can extract the index SNP ID by performing the following command:

```awk
awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp
```

> `$3` because the third column contains the SNP ID


!!! note
    If your target data are small (e.g. N < 500) then you can use the 1000 Genomes Project samples for the LD calculation.
    Make sure to use the population that most closely reflects represents the base sample.

## Generate PRS
`plink` provides a convenient function `--score` and `--q-score-range` for calculating polygenic scores.

We will need three files:

1. The base data file: **Height.QC.Transformed**
2. A file containing SNP IDs and their corresponding P-values (`$3` because SNP ID is located in the third column; `$8` because the P-value is located in the eighth column)

```awk
awk '{print $3,$8}' Height.QC.Transformed > SNP.pvalue
```

3. A file containing the different P-value thresholds for inclusion of SNPs in the PRS. Here calculate PRS corresponding to a few thresholds for illustration purposes:

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
    `0` to `0.05`, **including** any SNPs with P-value equal to `0.05`.

We can then calculate the PRS with the following `plink` command:

```bash
plink \
    --bfile EUR.QC \
    --score Height.QC.Transformed 3 4 12 header \
    --q-score-range range_list SNP.pvalue \
    --extract EUR.valid.snp \
    --out EUR
```
The meaning of the new parameters are as follows:

| Paramter | Value | Description|
|:-:|:-:|:-|
|score|Height.QC.Transformed 3 4 12 header| We read from the **Height.QC.Transformed** file, assuming that the `3`st column is the SNP ID; `4`th column is the effective allele information; the `12`th column is the effect size estimate; and that the file contains a `header`|
|q-score-range| range_list SNP.pvalue| We want to calculate PRS based on the thresholds defined in **range_list**, where the threshold values (P-values) were stored in **SNP.pvalue**|

The above command and range_list will generate 7 files:

1. EUR.0.5.profile
2. EUR.0.4.profile
3. EUR.0.3.profile
4. EUR.0.2.profile
5. EUR.0.1.profile
6. EUR.0.05.profile
7. EUR.0.001.profile

!!! Note
    The default formula for PRS calculation in PLINK is:
    
    $$
    PRS_j =\frac{ \sum_i^NS_i*G_{ij}}{P*M_j}
    $$

    where the effect size of SNP $i$ is $S_i$;  the number of effect alleles observed in sample $j$ is $G_{ij}$; the ploidy of the sample is $P$ (is generally 2 for humans); the total number of SNPs included in the PRS is $N$; and the number of non-missing SNPs observed in sample $j$ is $M_j$. If the sample has a missing genotype for SNP $i$, then the population minor allele frequency multiplied by the ploidy ($MAF_i*P$) is used instead of $G_{ij}$.

## Accounting for Population Stratification

Population structure is the principal source of confounding in GWAS and is usually accounted for by incorporating principal components (PCs) as covariates. We can incorporate PCs into our PRS analysis to account for population stratification.

Again, we can calculate the PCs using `plink`: 

```bash
# First, we need to perform prunning
plink \
    --bfile EUR.QC \
    --indep-pairwise 200 50 0.25 \
    --out EUR
# Then we calculate the first 6 PCs
plink \
    --bfile EUR.QC \
    --extract EUR.prune.in \
    --pca 6 \
    --out EUR
```

!!! note
    One way to select the appropriate number of PCs is to perform GWAS on the phenotype under study with different numbers of PCs.
    [LDSC](https://github.com/bulik/ldsc) analysis can then be performed on the set of GWAS summary statistics and the GWAS that used the number of PCs that gave an LDSC intercept closest to 1 should correspond to that for which population structure was most accurately controlled for. 

Here the PCs have been stored in the **EUR.eigenvec** file and can be used as covariates in the regression model to account for population stratification.

!!! important
    If the base and target samples are collected from different worldwide populations then the results from the PRS analysis may be biased (see Section 3.4 of our papper).


## Finding the "best-fit" PRS
The P-value threshold that provides the "best-fit" PRS under the C+T method is usually unknown. 
To approximate the "best-fit" PRS, we can perform a regression between PRS calculated at a range of P-value thresholds and then select the PRS that explains the highest phenotypic variance (please see Section 4.6 of our paper on overfitting issues). 
This can be achieved using `R` as follows:

=== "detail"
    ```R 
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
    covariate <- read.table("EUR.cov", header=T)
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
    q() # exit R
    ```

=== "quick"
    ```R
    p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
    phenotype <- read.table("EUR.height", header=T)
    pcs <- read.table("EUR.eigenvec", header=F)
    colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
    covariate <- read.table("EUR.cov", header=T)
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
    q() # exit R
    ```

=== "with data.table and magrittr"

    ```R
    library(data.table)
    library(magrittr)
    p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
    phenotype <- fread("EUR.height")
    pcs <- fread("EUR.eigenvec", header=F) %>%
        setnames(., colnames(.), c("FID", "IID", paste0("PC",1:6)) )
    covariate <- fread("EUR.cov")
    pheno <- merge(phenotype, covariate) %>%
            merge(., pcs)
    null.r2 <- summary(lm(Height~., data=pheno[,-c("FID", "IID")]))$r.squared
    prs.result <- NULL
    for(i in p.threshold){
        pheno.prs <- paste0("EUR.", i, ".profile") %>%
            fread(.) %>%
            .[,c("FID", "IID", "SCORE")] %>%
            merge(., pheno, by=c("FID", "IID"))

        model <- lm(Height~., data=pheno.prs[,-c("FID","IID")]) %>%
                summary
        model.r2 <- model$r.squared
        prs.r2 <- model.r2-null.r2
        prs.coef <- model$coeff["SCORE",]
        prs.result %<>% rbind(.,
            data.frame(Threshold=i, R2=prs.r2, 
                        P=as.numeric(prs.coef[4]), 
                        BETA=as.numeric(prs.coef[1]),
                        SE=as.numeric(prs.coef[2])))
    }
    print(prs.result[which.max(prs.result$R2),])
    q() # exit R
    ```

??? note "Which P-value threshold generates the "best-fit" PRS?"
    0.3

??? note "How much phenotypic variation does the "best-fit" PRS explain?"
    0.1638468

