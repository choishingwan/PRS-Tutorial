Over the following three pages you can run three dedicated PRS programs, which automate many of the steps from the previous page that used a sequence of PLINK functions (plus some QC steps). 
On this page you will run a PRS analysis using PRSice-2, which implements the standard C+T method.

This analysis assumes that you have the following files (or you can download it from [here](https://drive.google.com/file/d/1_ujJhCFAAHp_fA2U291pBUPTeF_FQLyu/view?usp=sharing)): 

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post QC base data file. While PRSice-2 can automatically apply most filtering on the base file, it cannot remove duplicated SNPs|
|**EUR.QC.bed**| This file contains the genotype data that passed the QC steps |
|**EUR.QC.bim**| This file contains the list of SNPs that passed the QC steps |
|**EUR.QC.fam**| This file contains the samples that passed the QC steps |
|**EUR.height**| This file contains the phenotype data of the samples |
|**EUR.covariate**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the principal components (PCs) of the samples |

And `PRSice-2`, which can be downloaded from:

| Operating System | Link |
| -----------------|:----:|
| Linux 64-bit | [v2.3.2](https://github.com/choishingwan/PRSice/releases/download/2.3.2/PRSice_linux.zip) |
| OS X 64-bit | [v2.3.2](https://github.com/choishingwan/PRSice/releases/download/2.3.2/PRSice_mac.zip) |

In this tutorial, you will only need `PRSice.R` and `PRSice_XXX` where XXX is the operation system

# Running PRS analysis
To run PRSice-2 we need a single covariate file, and therefore our covariate file and PCs file should be combined. This can be done with `R` as follows:


```R tab="without data.table"    
covariate <- read.table("EUR.cov", header=T)
pcs <- read.table("EUR.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,"EUR.covariate", quote=F, row.names=F)
q()
```

```R tab="with data.table"
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

```bash tab="Linux"
Rscript PRSice.R \
    --prsice PRSice_linux \
    --base Height.QC.gz \
    --target EUR.QC \
    --binary-target F \
    --pheno EUR.height \
    --cov EUR.cov \
    --base-maf MAF:0.01 \
    --base-info INFO:0.8 \
    --stat OR \
    --or \
    --out EUR
```


```bash tab="OS X"
Rscript PRSice.R \
    --prsice PRSice_mac \
    --base Height.QC.gz \
    --target EUR.QC \
    --binary-target F \
    --pheno EUR.height \
    --cov EUR.cov \
    --base-maf MAF:0.01 \
    --base-info INFO:0.8 \
    --stat OR \
    --or \
    --out EUR
```

```bash tab="Windows"
Rscript PRSice.R ^
    --prsice PRSice_win64.exe ^
    --base Height.QC.gz ^
    --target EUR.QC ^
    --binary-target F ^
    --pheno EUR.height ^
    --cov EUR.cov ^
    --base-maf MAF:0.05 ^
    --base-info INFO:0.8 ^
    --stat OR ^
    --or ^
    --out EUR
```

This will automatically perform "high-resolution scoring" and generate the "best-fit" PRS (in **EUR.best**), with associated plots of the results. 
Users should read Section 4.6 of our paper to learn more about issues relating to overfitting in PRS analyses.  

??? note "Which P-value threshold generates the "best-fit" PRS?"
    0.1388

??? note "How much phenotypic variation does the "best-fit" PRS explain?"
    0.174227