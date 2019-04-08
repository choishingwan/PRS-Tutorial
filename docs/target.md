Next, we'd like to perform basic quality controls (QC) on the target genotype data. 

In this tutorial, we've simulated some samples using the 1000 genome european genotypes. 
You can download the data [here](https://github.com/choishingwan/PRS-Tutorial/raw/master/resources/EUR.zip). 

Unzip the data as follow:

```bash
unzip EUR.zip
```

??? note "What's the md5sum of the genotype files?"

    |File|md5sum|
    |:-:|:-:|
    |**EUR.bed**           |96ce8f494a57114eaee6ef9741676f58|
    |**EUR.bim**           |852d54c9b6d1159f89d4aa758869e72a|
    |**EUR.covariate**     |afff13f8f9e15815f2237a62b8bec00b|
    |**EUR.fam**           |8c6463c0d8f32f975cdc423b1b80a951|
    |**EUR.height**        |052beb4cae32ac7673f1d6b9e854c85b|

!!! note
    We assume PLINK is installed in your PATH directory, which allow us to use `plink` instead of `./plink`.
    If PLINK is not in your PATH directory, replace all instance of `plink` in the tutorial to `./plink` assuming
    the PLINK executable is located within your working directory

# Basic filterings

```bash
plink --bfile EUR \
    --maf 0.05 \
    --hwe 1e-6 \
    --geno 0.01 \
    --mind 0.01 \
    --make-bed \
    --out EUR.QC
```
??? note "How many SNPs were filtered?"
    A total of `266` SNPs were removed due to Hardy-Weinberg exact test results

# Remove Ambiguous SNPs
If the base and target data were generated using different genotyping chips and the chromosome strand (+/-) for either is unknown, then it is not possible to match ambiguous SNPs (i.e. those with complementary alleles, either C/G or A/T) across the data sets, because it will be unknown whether the base and target data are referring to the same allele or not. 

Ambiguous SNPs can be obtained by examining the bim file:
```bash
awk '!( ($5=="A" && $6=="T") || \
        ($5=="T" && $6=="A") || \
        ($5=="G" && $6=="C") || \
        ($5=="C" && $6=="G")) {print}' EUR.QC.bim > EUR.unambig.snp 
```

The ambiguous SNPs can then be removed by

```bash
plink \
    --bfile EUR.QC \
    --extract EUR.unambig.snp \
    --make-bed \
    --out EUR.QC.NoAmbig
```

??? note "How many ambiguous SNPs were there?"
    There are `17,260` ambiguous SNPs

# Check for mis-matched Sex information

# Filter related samples
```bash
    plink --bfile EUR --indep-pairwise 200 50 0.25 --out EUR
    plink --bfile EUR --extract EUR.prune.in --rel-cutoff 0.125 --out EUR
```