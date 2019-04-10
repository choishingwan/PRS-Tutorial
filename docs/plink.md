In previous sections, we have generated the following files

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.valid.sample**| This file contains the samples that passed all the QC |

Here, we will try to calculate polygenic risk score using `plink`. 
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
    There are `17,260` ambiguous SNPs

??? note  "Strand Flipping"
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
    --bfile EUR.QC \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Height.QC.transformed \
    --keep EUR.valid.sample \
    --extract EUR.unambig.snp \
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
| keep | EUR.valid.sample | Samples for LD calculation |
| extract | EUR.unambig.snp | Only consider these SNPs for clumping |
| clump-snp-field | SNP | Specify that the column `SNP` contains the SNP IDs |
| clump-field | P | Specify that the column `P` contains the P-value information |

A more detailed document can be found [here](https://www.cog-genomics.org/plink/1.9/postproc#clump)

!!! note
    The $r^2$ values computed by `--clump` are based on maximum likelihood haplotype frequency estimates


This will generate **EUR.clumped**, containing the index SNPs after clumping is performed.
We can extract the index SNP ID by doing

```bash
awk 'NR!=1{print $3}' EUR.clumped >  EUR.valid.snp
```

> `$3` because the third column contains the SNP ID


!!! note
    If your target sample is small (e.g. <500), you can try using the 1000 Genome samples for the LD calculation.
    Make sure you use the population that best represents your sample.

# Generate PRS
`plink` provide a handy function `--score` and `--q-score-range` for calculating polygenic score.

We will need three files

1. The summary statistic file: **Height.QC.Transformed**
2. A file containing SNP ID and their corresponding p-value (`$1` because SNP ID is located at the first column; `$8` because P-value is located at the eigth column)
```bash
awk '{print $1,$8}' Height.QC.Transformed > SNP.pvalue
```
3. A file containing the P-value thresholds we are testing. Here we will only test a few for illustration purposes
```bash
echo "0.001 0 0.001" > range_list
echo "0.05 0 0.05" >> range_list
echo "0.1 0 0.1" >> range_list
echo "0.2 0 0.2" >> range_list
echo "0.3 0 0.3" >> range_list
echo "0.4 0 0.4" >> range_list
echo "0.5 0 0.5" >> range_list
```
The format of the **range_list** file should be as follow

|Name of Threshold|Lower bound| Upper Bound|
|:-:|:-:|:-:|

!!! note
    The boundary are inclusive. For example, for the `0.05` threshold, we include all SNPs with P-value from 
    `0` to `0.05`, **including** any SNPs with P-value equal to `0.05`

We can then calculate the PRS with the following `plink` command:

```bash
plink \
    --bfile EUR.QC \
    --extract EUR.valid.snp \
    --keep EUR.valid.sample \
    --score Height.QC.Transformed 1 4 11 header \
    --q-score-range range_list SNP.pvalue \
    --out EUR
```
Meaning of the new parameters are as follow

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
    (Assuming the effect size of SNP $i$ is $S_i$;  the number of effective allele observed in sample $j$ is $G_j$; the ploidy of the sample is $P$ (It should be 2 for human); and the number of samples and SNPs included in the PRS be $N$ and $M$ respectively)
    $$
    PRS_j = \sum_i^N\frac{S_i*G_{ij}}{P*M}
    $$

    If sample has a missing genotype for SNP $i$, the population minor allele frequency times ploidy ($MAF_i*P$)is then used inplace of $G_{ij}$

# Accounting for Population Stratification

Population structure is the principal source of confounding in GWAS and are usually accounted for by incorporating the principal components (PCs) as a covariate. 
Similarly, we can incorporate PCs in our PRS analysis to account for population stratification.

Again, we can calculate the PCs using `plink` 
```bash
# First, we need to perform prunning
plink \
    --bfile EUR.QC \
    --keep EUR.valid.sample \
    --extract EUR.valid.snp \
    --indep-pairwise 200 50 0.25 \
    --out EUR
# Then we calculate the first 20 PCs
plink \
    --bfile EUR.QC \
    --keep EUR.valid.sample \
    --extract EUR.prune.in \
    --pca 20 \
    --out EUR
```

!!! important
    If the base and target samples are collected from different population (e.g. Caucasian vs African ), the results from PRS analysis will be biased (see [Martin et al](https://www.ncbi.nlm.nih.gov/pubmed/28366442)).


# Finding the "Best" P-value threshold

# Plotting the Results