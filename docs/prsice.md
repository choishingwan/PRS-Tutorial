Over the following three pages you can run three devoted PRS programs, which automate many of the steps from the previous page that used a sequence of PLINK functions (plus some QC steps from earlier pages). On this page you will run a PRS analysis using PRSice-2, which implements the standard C+T method.

This analysis assumes that you have the following files: 

|File Name | Description|
|:-:|:-:|
|**GIANT.height.gz**| The original base data file. PRSice-2 can apply INFO and MAF filtering to these base summary statistic data directly |
|**EUR.QC.bed**| This file contains the genotype data that passed the QC steps |
|**EUR.QC.bim**| This file contains the list of SNPs that passed the QC steps |
|**EUR.QC.fam**| This file contains the samples that passed the QC steps |
|**EUR.valid.sample**| This file contains the samples that passed the QC steps |
|**EUR.height**| This file contains the phenotype data of the samples |
|**EUR.covariate**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the principal components (PCs) of the samples |

And `PRSice-2`, which can be downloaded from:

| Operating System | Link |
| -----------------|:----:|
| Linux 64-bit | [v2.2.5](https://github.com/choishingwan/PRSice/releases/download/2.2.5/PRSice_linux.zip) |
| OS X 64-bit | [v2.2.5](https://github.com/choishingwan/PRSice/releases/download/2.2.5/PRSice_mac.zip) |
| Windows 32-bit | [v2.2.5](https://github.com/choishingwan/PRSice/releases/download/2.2.5/PRSice_win32.zip) |
| Windows 64-bit | [v2.2.5](https://github.com/choishingwan/PRSice/releases/download/2.2.5/PRSice_win64.zip) |

In this tutorial, you will only need `PRSice.R` and `PRSice_XXX` where XXX is the operation system

# Running PRS analysis
To run PRSice-2 we need a single covariate file, and therefore our covariate file and PCs file should be combined. This can be done with `R` as follows:

```R
covariate <- read.table("EUR.covariate", header=T)
pcs <- read.table("EUR.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,"EUR.cov", quote=F, row.names=F)
```
which generates **EUR.cov**.

PRSice-2 can then be run to obtain the PRS results as follows:

```bash tab="Linux"
Rscript PRSice.R \
    --prsice PRSice_linux \
    --base Height.QC.gz \
    --target EUR.QC \
    --keep EUR.valid.sample \
    --binary-target F \
    --pheno-file EUR.height \
    --cov-file EUR.cov \
    --maf-base MAF,0.05 \
    --info-base INFO,0.8 \
    --out EUR
```


```bash tab="OS X"
Rscript PRSice.R \
    --prsice PRSice_mac \
    --base Height.QC.gz \
    --target EUR.QC \
    --keep EUR.valid.sample \
    --binary-target F \
    --pheno-file EUR.height \
    --cov-file EUR.cov \
    --maf-base MAF,0.05 \
    --info-base INFO,0.8 \
    --out EUR
```

```bash tab="Windows"
Rscript PRSice.R ^
    --prsice PRSice_win64.exe ^
    --base Height.QC.gz ^
    --target EUR.QC ^
    --keep EUR.valid.sample ^
    --binary-target F ^
    --pheno-file EUR.height ^
    --cov-file EUR.cov ^
    --maf-base MAF,0.05 ^
    --info-base INFO,0.8 ^
    --out EUR
```

This will automatically perform "high-resolution scoring" and generate the "best-fit" PRS (in **EUR.best**), with associated plots of the results. Users should read Section 4.6 of our paper to learn more about issues relating to overfitting in PRS analyses.  


