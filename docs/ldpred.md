Here we use another PRS program, [LDpred](https://github.com/bvilhjal/ldpred), that uses a Bayesian approach to polygenic risk scoring.

!!! note
    Python 3 and other packages need to be installed before running LDpred. 
    Please refer to their website for installation instructions.
    
    If you have Python installed, you should be able to install LDpred with the following command:
    ```
    pip install ldpred
    ```

!!! note
    The script used here is based on LDpred version 1.0.6

    Due to current update in LDpred, this tutorial is likely outdated. We will try to update it once we've tested the new version.

We assume that you have the following files:

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.covariate**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the PCs of the samples |


# Running PRS analysis
LDpred can then perform PRS analysis in three steps:

1. Preprocessing the base data file:
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
    --gf EUR.QC
```

2. Adjust the effect size estimates:
``` bash
# LDpred recommend radius to be Total number of SNPs in target / 3000
 python LDpred.py gibbs \
    --cf EUR.coord \
    --ldr 183 \
    --ldf EUR.ld \
    --out EUR.weight \
    --N 253288
```

3. Calculate the PRS:
```bash 
python LDpred.py score \
    --gf EUR.ldpred \
    --rf EUR.weight \
    --out EUR.score \
    --pf EUR.height \
    --pf-format LSTANDARD 
```

!!! note
    To obtain the PRS $R^2$ according to the use of different parameters, and / or 
    to adjust for covariates, you may need to use `R` (see [here](https://github.com/bvilhjal/ldpred/wiki/Q-and-A#im-having-trouble-with-covariates-can-you-help-me)).
