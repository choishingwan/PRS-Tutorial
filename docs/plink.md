
# Remove Ambiguous SNPs
If the base and target data were generated using different genotyping chips and the chromosome strand (+/-) for either is unknown, then it is not possible to match ambiguous SNPs (i.e. those with complementary alleles, either C/G or A/T) across the data sets, because it will be unknown whether the base and target data are referring to the same allele or not. 

Ambiguous SNPs can be obtained by examining the bim file:
```bash
awk '!( ($5=="A" && $6=="T") || \
        ($5=="T" && $6=="A") || \
        ($5=="G" && $6=="C") || \
        ($5=="C" && $6=="G")) {print}' \
        EUR.QC.unrel.het.bim > EUR.unambig.snp 
```

The ambiguous SNPs can then be removed by

```bash
plink \
    --bfile EUR.QC.unrel.het \
    --extract EUR.unambig.snp \
    --make-bed \
    --out EUR.QC.unrel.het.NoAmbig
```

??? note "How many ambiguous SNPs were there?"
    There are `17,260` ambiguous SNPs

# Strand Flipping
Alternatively, when there is a non-ambiguous mismatch in allele coding between the data sets, such as A/C in the base
and G/T in the target data, then this can be resolved by ‘flipping’ the alleles in the target data to their complementary alleles. 
This can be done using `R`

# Update Effect Size
When odd ratios (OR) instead of BETA are provided, the PRS might have to calculated using a multiplicative model.
To simplify the calculation, we usually take the natural logarithm of the OR such that an additive model can be applied. 
We can obtain the transformed summary statistics with `R`:

```R
dat <- read.table(gzfile("Height.QC.gz"), header=T)
dat$OR <- log(dat$OR)
write.table(dat, "Height.QC.Transformed", quote=F, row.names=F)
```

!!! warning
    While you can also do the log transformation using `awk`, the resulting transformation will only have precision upto 7th digit. 
    The imprecision can accumulate and lead to slightly less accurate results. 
    Therefore it is best to do the transformation in `R` or allow the PRS software to do the transformation for you. 

# Clumping
Linkage disequilibrium introduce a strong correlation structure across the genome, makes identifying the independent
genetic effects extremely challenging. 
One simple method is to perform clumping, which preferentially selects SNPs most
associated with the trait under study when removing SNPs in LD. 

```bash
plink \
    --bfile EUR.QC.unrel.het \
    --clump-p1 1 \
    --clump-p2 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Height.QC.transformed \
    --out EUR
```

# Generate PRS

# Accounting for Population Stratification

# Finding the "Best" P-value threshold

# Plotting the Results