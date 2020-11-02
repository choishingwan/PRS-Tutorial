# Background

PRSice-2 is one of the dedicated PRS programs which automates many of the steps from the previous page that used a sequence of PLINK functions (plus some QC steps). 
On this page you will run a PRS analysis using PRSice-2, which implements the standard C+T method.

## Obtaining PRSice-2
`PRSice-2` can be downloaded from:

| Operating System | Link |
| -----------------|:----:|
| Linux 64-bit | [v2.3.3](https://github.com/choishingwan/PRSice/releases/download/2.3.3/PRSice_linux.zip) |
| OS X 64-bit | [v2.3.3](https://github.com/choishingwan/PRSice/releases/download/2.3.3/PRSice_mac.zip) |

and can be directly used after extracting the file. 

In this tutorial, you will only need `PRSice.R` and `PRSice_XXX` where XXX is the operation system

## Required Data

This analysis assumes that you have the following files (or you can download it from [here](https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view?usp=sharing)): 

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post QC base data file. While PRSice-2 can automatically apply most filtering on the base file, it cannot remove duplicated SNPs|
|**EUR.QC.bed**| This file contains the genotype data that passed the QC steps |
|**EUR.QC.bim**| This file contains the list of SNPs that passed the QC steps |
|**EUR.QC.fam**| This file contains the samples that passed the QC steps |
|**EUR.height**| This file contains the phenotype data of the samples |
|**EUR.cov**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the principal components (PCs) of the samples |

## Running PRS analysis
To run PRSice-2 we need a single covariate file, and therefore our covariate file and PCs file should be combined. This can be done with `R` as follows:

=== "without data.table"    

    ```R
    covariate <- read.table("EUR.cov", header=T)
    pcs <- read.table("EUR.eigenvec", header=F)
    colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
    cov <- merge(covariate, pcs, by=c("FID", "IID"))
    write.table(cov,"EUR.covariate", quote=F, row.names=F)
    q()
    ```

=== "with data.table"

    ```R
    library(data.table)
    covariate <- fread("EUR.cov")
    pcs <- fread("EUR.eigenvec", header=F)
    colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
    cov <- merge(covariate, pcs)
    fwrite(cov,"EUR.covariate", sep="\t")
    q()
    ```

which generates **EUR.cov**.

PRSice-2 can then be run to obtain the PRS results as follows:

=== "Linux"
    ```bash
    Rscript PRSice.R \
        --prsice PRSice_linux \
        --base Height.QC.gz \
        --target EUR.QC \
        --binary-target F \
        --pheno EUR.height \
        --cov EUR.covariate \
        --base-maf MAF:0.01 \
        --base-info INFO:0.8 \
        --stat OR \
        --or \
        --out EUR
    ```

=== "OS X"
    ```bash
    Rscript PRSice.R \
        --prsice PRSice_mac \
        --base Height.QC.gz \
        --target EUR.QC \
        --binary-target F \
        --pheno EUR.height \
        --cov EUR.covariate \
        --base-maf MAF:0.01 \
        --base-info INFO:0.8 \
        --stat OR \
        --or \
        --out EUR
    ```

=== "Windows"
    ```bash
    Rscript PRSice.R ^
        --prsice PRSice_win64.exe ^
        --base Height.QC.gz ^
        --target EUR.QC ^
        --binary-target F ^
        --pheno EUR.height ^
        --cov EUR.covariate ^
        --base-maf MAF:0.05 ^
        --base-info INFO:0.8 ^
        --stat OR ^
        --or ^
        --out EUR
    ```

The meaning of the parameters are as follow:

| Paramter | Value | Description|
|:-:|:-:|:-|
|prsice|PRSice_xxx| Informs `PRSice.R` that the location of the PRSice binary |
|base| Height.QC.gz| Informs `PRSice` that the name of the GWAS summary statistic |
|target| EUR.QC| Informs `PRSice` that the input genotype files should have a prefix of `EUR.QC` |
|binary-target| F| Indicate if the phenotype of interest is a binary trait. F for no |
|pheno| EUR.height| Provide `PRSice` with the phenotype file |
|cov| EUR.covariate| Provide `PRSice` with the covariate file |
|base-maf| MAF:0.05| Filter out SNPs with MAF < 0.05 in the GWAS summary statistics, using information in the `MAF` column|
|base-info| INFO:0.8| Filter out SNPs with INFO < 0.8 in the GWAS summary statistics, using information in the `INFO` column|
|stat| OR| Column name of the column containing the effect size|
|or|-| Inform `PRSice` that the effect size is an Odd Ratio|
|out | EUR | Informs `PRSice` that all output should have a prefix of `EUR`|

This will automatically perform "high-resolution scoring" and generate the "best-fit" PRS (in **EUR.best**), with associated plots of the results. 
Users should read Section 4.6 of our paper to learn more about issues relating to overfitting in PRS analyses.  

??? note "Which P-value threshold generates the "best-fit" PRS?"
    0.13995

??? note "How much phenotypic variation does the "best-fit" PRS explain?"
    0.166117
