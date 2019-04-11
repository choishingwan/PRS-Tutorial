In previous sections, we have generated the following files

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.valid.sample**| This file contains the samples that passed all the QC |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.covariate**| This file contains the covariates of the samples |

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

# Strand Flipping
Alternatively, when there is a non-ambiguous mismatch in allele coding between the data sets, such as A/C in the base
and G/T in the target data, then this can be resolved by ‘flipping’ the alleles in the target data to their complementary alleles. 
This has to be done with mulitpe steps

1. Get the correct A1 alleles for the bim file
```R
bim <- read.table("EUR.QC.bim", header=F)
colnames(bim) <- c("CHR", "SNP", "CM", "BP", "B.A1", "B.A2")
height <- read.table(gzfile("Height.QC.gz"), header=T)
# Avoid complicated factor problem
height$A1 <- toupper(as.character(height$A1))
height$A2 <- toupper(as.character(height$A2))
bim$B.A1 <- toupper(as.character(bim$B.A1))
bim$B.A2 <- toupper(as.character(bim$B.A2))
info <- merge(bim, height, by=c("SNP", "CHR", "BP"))

# Function for calculating the complementary allele
complement <- function(x){
    switch (x,
        "A" = "T",
        "C" = "G",
        "T" = "A",
        "G" = "C",
        return(NA)
    )
}
# Now get SNPs that has exact match between base and target
info.match <- subset(info, A1==B.A1 & A2==B.A2)
# Check for complementary matchs
info$C.A1 <- sapply(info$B.A1, complement)
info$C.A2 <- sapply(info$B.A2, complement)
info.complement <- subset(info, A1==C.A1 & A2==C.A2)
# Update these allele coding in the bim file 
bim[bim$SNP %in% info.complement$SNP, ]$B.A1 <- sapply(bim[bim$SNP %in% info.complement$SNP, ]$B.A1, complement)
bim[bim$SNP %in% info.complement$SNP, ]$B.A2 <- sapply(bim[bim$SNP %in% info.complement$SNP, ]$B.A2, complement)
# identify SNPs that need flipping 
info.flip <- subset(info, A1==B.A2 & A2==B.A1)
# identify SNPs that need flipping & complement
info.cflip <- subset(info, A1==C.A2 & A2==C.A1)
# Update these allele coding in the bim file 
bim[bim$SNP %in% info.cflip$SNP, ]$B.A1 <- sapply(bim[bim$SNP %in% info.cflip$SNP, ]$B.A1, complement)
bim[bim$SNP %in% info.cflip$SNP, ]$B.A2 <- sapply(bim[bim$SNP %in% info.cflip$SNP, ]$B.A2, complement)
# Get list of SNPs that need to change the A1 encoding
flip <- rbind(info.flip, info.cflip)
flip.snp <- data.frame(SNP=flip$SNP, A1=flip$A1)
write.table(flip.snp, "EUR.update.a1", quote=F, row.names=F)
write.table(bim, "EUR.QC.adj.bim", quote=F, row.names=F, col.names=F)
# And we want to remove any SNPs that do not match with the base data
mismatch <- bim$SNP[!(bim$SNP %in% info.match$SNP | bim$SNP %in% info.complement$SNP | bim$SNP %in% flip$SNP)]
write.table(mismatch, "EUR.mismatch", quote=F, row.names=F, col.names=F)
```

The above script will generate three files: **EUR.QC.adj.bim**, **EUR.update.a1** and **EUR.mismatch**. We want to replace
**EUR.QC.bim** with **EUR.QC.adj.bim**:

```bash
# Make a back up
mv EUR.QC.bim EUR.QC.bim.bk
ln -s EUR.QC.adj.bim EUR.QC.bim
```

We can then generate a new genotype file with the correct genetic encodings
```bash
plink \
    --bfile EUR.QC \
    --a1-allele EUR.update.a1 \
    --make-bed \
    --keep EUR.valid.sample \
    --extract EUR.unambig.snp \
    --exclude EUR.mismatch \
    --out EUR.QC.flipped
```

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
    --bfile EUR.QC.flipped \
    --clump-p1 1 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump Height.QC.transformed \
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
    --bfile EUR.QC.flipped \
    --extract EUR.valid.snp \
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
    (Assuming the effect size of SNP $i$ is $S_i$;  the number of effective allele observed in sample $j$ is $G_{ij}$; the ploidy of the sample is $P$ (It should be 2 for human); the number of samples included in the PRS be $N$; and the number of non-missing SNPs observed in sample $j$ be $M_j$)
    $$
    PRS_j =\frac{ \sum_i^NS_i*G_{ij}}{P*M_j}
    $$

    If sample has a missing genotype for SNP $i$, the population minor allele frequency times ploidy ($MAF_i*P$) is used inplace of $G_{ij}$

# Accounting for Population Stratification

Population structure is the principal source of confounding in GWAS and are usually accounted for by incorporating the principal components (PCs) as a covariate. 
Similarly, we can incorporate PCs in our PRS analysis to account for population stratification.

Again, we can calculate the PCs using `plink` 
```bash
# First, we need to perform prunning
plink \
    --bfile EUR.QC.flipped \
    --extract EUR.valid.snp \
    --indep-pairwise 200 50 0.25 \
    --out EUR
# Then we calculate the first 6 PCs
plink \
    --bfile EUR.QC.flipped \
    --extract EUR.prune.in \
    --pca 6 \
    --out EUR
```

!!! note
    One way to select the appropriate number of PCs is to perform GWAS on the trait of interest with different number of PCs.
    [LDSC](https://github.com/bulik/ldsc) analysis can then be performed on each of the resulted GWAS summary statistics. 
    By observing the estimated intercept, one can select the number of PCs that provide an intercept estimate closer to 1, which might
    suggest a smaller influence of population stratification.

The eigen-vector (PCs) are stored in **EUR.eigenvec** and can be used as a covariate in the regression model to account for population stratification.

!!! important
    If the base and target samples are collected from different population (e.g. Caucasian vs African ), the results from PRS analysis will be biased (see [Martin et al](https://www.ncbi.nlm.nih.gov/pubmed/28366442)).


# Finding the "Best" P-value threshold
The "best" p-value threshold for PRS construction are usually not known. 
To identify the "best" PRS, we can perform a regression between the calculated PRS and the 
sample phenotype and select the PRS that explains most of the phenotypic variation. 
This can be achieved using `R`.

```R tab="detail"
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
covariate <- read.table("EUR.covariate", header=T)
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
```

```R tab="quick"
p.threshold <- c(0.001,0.05,0.1,0.2,0.3,0.4,0.5)
phenotype <- read.table("EUR.height", header=T)
pcs <- read.table("EUR.eigenvec", header=F)
colnames(pcs) <- c("FID", "IID", paste0("PC",1:6)) 
covariate <- read.table("EUR.covariate", header=T)
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
```

??? note "Which p-value threshold generate the "best" PRS?"
    0.2

??? note "How much phenotypic variation does the "best" PRS explains?"
    0.04128065

# Plotting the Results
We can also visualize our results using `R`

!!! note
    We will be using `prs.result` generated in [previous section](#finding-the-best-p-value-threshold)


```R tab="ggplot2"
# ggplot2 is a handy package for plotting
library(ggplot2)
# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) &
                       prs.result$print.p == 0] <-
    format(prs.result$P[!is.na(prs.result$print.p) &
                            prs.result$print.p == 0], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
# Initialize ggplot, requiring the threshold as the x-axis (use factor so that it is uniformly distributed)
ggplot(data = prs.result, aes(x = factor(Threshold), y = R2)) +
    # Specify that we want to print p-value on top of the bars
    geom_text(
        aes(label = paste(print.p)),
        vjust = -1.5,
        hjust = 0,
        angle = 45,
        cex = 4,
        parse = T
    )  +
    # Specify the range of the plot, *1.25 to provide enough space for the p-values
    scale_y_continuous(limits = c(0, max(prs.result$R2) * 1.25)) +
    # Specify the axis labels
    xlab(expression(italic(P) - value ~ threshold ~ (italic(P)[T]))) +
    ylab(expression(paste("PRS model fit:  ", R ^ 2))) +
    # Draw a bar plot
    geom_bar(aes(fill = -log10(P)), stat = "identity") +
    # Specify the colors
    scale_fill_gradient2(
        low = "dodgerblue",
        high = "firebrick",
        mid = "dodgerblue",
        midpoint = 1e-4,
        name = bquote(atop(-log[10] ~ model, italic(P) - value),)
    ) +
    # Some beautification of the plot
    theme_classic() + theme(
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(size = 14),
        legend.title = element_text(face = "bold", size =
                                        18),
        legend.text = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust =
                                       1)
    )
# save the plot
ggsave("EUR.height.bar.png", height = 7, width = 7)
```

```R tab="base plot"
# We strongly recommend the use of ggplot2. Only follow this code if you
# are desperate.
# Specify that we want to generate plot in EUR.height.bar.png
png("EUR.height.bar.png",
      height=10, width=10, res=300, unit="in")
# First, obtain the colorings based on the p-value
col <- suppressWarnings(colorRampPalette(c("dodgerblue", "firebrick")))
# We want the color gradient to match the ranking of p-values
prs.result <- prs.result[order(-log10(prs.result$P)),]
prs.result$color <-  col(nrow(prs.result))
prs.result <- prs.result[order(prs.result$Threshold),]
# generate a pretty format for p-value output
prs.result$print.p <- round(prs.result$P, digits = 3)
prs.result$print.p[!is.na(prs.result$print.p) & prs.result$print.p == 0 ] <-
    format(prs.result$P[!is.na(prs.result$print.p) & prs.result$print.p == 0 ], digits = 2)
prs.result$print.p <- sub("e", "*x*10^", prs.result$print.p)
# Generate the axis labels
xlab <- expression(italic(P) - value ~ threshold ~ (italic(P)[T]))
ylab <- expression(paste("PRS model fit:  ", R ^ 2))
# Setup the drawing area
layout(t(1:2), widths=c(8.8,1.2))
par( cex.lab=1.5, cex.axis=1.25, font.lab=2, 
    oma=c(0,0.5,0,0),
    mar=c(4,6,0.5,0.5))
# Plotting the bars
b<- barplot(height=prs.result$R2, 
            col=prs.result$color, 
            border=NA, 
            ylim=c(0, max(prs.result$R2)*1.25), 
            axes = F, ann=F)
# Plot the axis labels and axis ticks
odd <- seq(0,nrow(prs.result)+1,2)
even <- seq(1,nrow(prs.result),2)
axis(side=1, at=b[odd], labels=prs.result$Threshold[odd], lwd=2)
axis(side=1, at=b[even], labels=prs.result$Threshold[even],lwd=2)
axis(side=1, at=c(0,b[1],2*b[length(b)]-b[length(b)-1]), labels=c("","",""), lwd=2, lwd.tick=0)
# Write the p-value on top of each bar
text( parse(text=paste(
    prs.result$print.p)), 
    x = b+0.1, 
    y =  prs.result$R2+ (max(prs.result$R2)*1.05-max(prs.result$R2)), 
    srt = 45)
# Now plot the axis lines
box(bty='L', lwd=2)
axis(2,las=2, lwd=2)
# Plot the axis titles
title(ylab=ylab, line=4, cex.lab=1.5, font=2 )
title(xlab=xlab, line=2.5, cex.lab=1.5, font=2 )
# Generate plot area for the legend
par(cex.lab=1.5, cex.axis=1.25, font.lab=2, 
      mar=c(20,0,20,4))
prs.result <- prs.result[order(-log10(prs.result$P)),]
image(1, -log10(prs.result$P), t(seq_along(-log10(prs.result$P))), col=prs.result$color, axes=F,ann=F)
axis(4,las=2,xaxs='r',yaxs='r', tck=0.2, col="white")
# plot legend title
title(bquote(atop(-log[10] ~ model, italic(P) - value), ), 
          line=2, cex=1.5, font=2, adj=0)
# write the plot to file
dev.off()
```
