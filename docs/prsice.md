An alternative to `plink` is `PRSice-2`, which automates much of the PRS analyses.

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

And `PRSice-2`, which can be downloaded from

| Operating System | Link |
| -----------------|:----:|
| Linux 64-bit | [v2.1.11](https://github.com/choishingwan/PRSice/releases/download/2.1.11/PRSice_linux.zip) |
| OS X 64-bit | [v2.1.11](https://github.com/choishingwan/PRSice/releases/download/2.1.11/PRSice_mac.zip) |
| Windows 32-bit | [v2.1.11](https://github.com/choishingwan/PRSice/releases/download/2.1.11/PRSice_win32.zip) |
| Windows 64-bit | [v2.1.11](https://github.com/choishingwan/PRSice/releases/download/2.1.11/PRSice_win64.zip) |

In this tutorial, you will only need `PRSice.R` and `PRSice_XXX` where XXX is the operation system

# Running PRS analysis
It is simple to run PRSice-2. First, we need a single covariate file. This can be done with `R`:

```R
covariate <- read.table("EUR.covariate", header=T)
pcs <- read.table("EUR.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
cov <- merge(covariate, pcs, by=c("FID", "IID"))
write.table(cov,"EUR.cov", quote=F, row.names=F)
```
which generates **EUR.cov**

PRSice-2 can then be run to obtain the PRS results:

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

This will automatically perform a "high-resolution scoring" and generate the "best" PRS (in **EUR.best**) and relevant graphs    


