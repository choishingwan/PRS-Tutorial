Another popular PRS software is [LDpred](https://github.com/bvilhjal/ldpred), which instead of performing p-value thresholding,
infers the posterior mean effect size of each marker by using a prior on effect sizes and LD information from an external reference panel, 
thus allow for a better $R^2$.

!!! note
    Python 3 and other packages need to be installed before running LDpred. Please refer
    to the github for instructions of installation

!!! note
    Current script is based on version 1.0.0

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


# Running PRS analysis
LDpred does not support in place filtering of sample and SNPs, therefore we need to generate a new QCed genotype file using `plink`

``` bash
# We also add the height phenotype for convenience
plink \
    --bfile EUR.QC \
    --pheno EUR.height \
    --keep EUR.valid.sample \
    --make-bed \
    --out EUR.ldpred
```

LDpred can then performs PRS analysis in three steps

1. Preprocessing the summary statistic file
```bash
# There are 253,288 samples in the Height GWAS
python LDpred.py coord \
    --rs SNP \
    --A1 A1 \
    --A2 A2 \
    --pos BP \
    --chr CHR \
    --pval P \
    --eff OR \
    --ssf-format CUSTOM \
    --N 253288 \
    --ssf Height.QC.gz \
    --out EUR.coord \
    --gf EUR.ldpred
```

2. Adjust the effect size
``` bash
# LDpred recommend radius to be Total number of SNPs in target / 3000
 python LDpred.py gibbs \
    --cf EUR.coord \
    --ldr 183 \
    --ldf EUR.ld \
    --out EUR.weight \
    --N 253288;
```

3. Calculate the PRS
```bash 
python LDpred.py score \
    --gf EUR.ldpred \
    --rf EUR.weight \
    --out EUR.score \
    --pf EUR.height \
    --pf-format LSTANDARD 
```

!!! note
    To obtain the $R^2$ of PRS obtained from different parameters, and / or 
    adjust for covariate, you might need to use `R` (see [here](https://github.com/bvilhjal/ldpred/wiki/Q-and-A#im-having-trouble-with-covariates-can-you-help-me))