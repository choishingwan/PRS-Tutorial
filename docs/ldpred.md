# Background

[LDpred-2](https://privefl.github.io/bigsnpr/articles/LDpred2.html) is one of the dedicated PRS programs which is an `R` package that uses a Bayesian approach to polygenic risk scoring.

## Installing LDpred-2

!!! note
    The script used here is based on LDpred 2 implemented under bigsnpr version 1.4.7

!!! note
    For more details, please refer to [LDpred 2's homepage](https://privefl.github.io/bigsnpr/articles/LDpred2.html)

You can install `LDpred` and its dependencies in `R` with the following command:

```R
install.packages("remotes")
library(remotes)
remotes::install_github("https://github.com/privefl/bigsnpr.git")
```

!!! note
    For mac users, you might need to follow the guide [here](https://thecoatlessprofessor.com/programming/cpp/openmp-in-r-on-os-x/) to be able to install LDpred2

## Required Data

We assume that you have the following files (or you can download it from [here](https://drive.google.com/file/d/1x_G0Gxk9jFMY-PMqwtg6-vdEyUPp5p5u/view?usp=sharing)):

|File Name | Description|
|:-:|:-:|
|**Height.QC.gz**| The post-QCed summary statistic |
|**EUR.QC.bed**| The genotype file after performing some basic filtering |
|**EUR.QC.bim**| This file contains the SNPs that passed the basic filtering |
|**EUR.QC.fam**| This file contains the samples that passed the basic filtering |
|**EUR.height**| This file contains the phenotype of the samples |
|**EUR.cov**| This file contains the covariates of the samples |
|**EUR.eigenvec**| This file contains the PCs of the samples |


## Running PRS analysis
!!! warning
    While we do provide a rough guide on how to perform LDpred on bed files separated into individual chromosomes, this script is untested and extra caution is required

### 0. Prepare workspace
On some server, you might need to first use the following code in order to run LDpred with multi-thread

=== "prepare workspace and load bigsnpr"

    ```R
    library(bigsnpr)
    options(bigstatsr.check.parallel.blas = FALSE)
    options(default.nproc.blas = NULL)
    ```

### 1. Read in the phenotype and covariate files

=== "read in phenotype and covariates"

    ```R 
    library(data.table)
    library(magrittr)
    phenotype <- fread("EUR.height")
    covariate <- fread("EUR.cov")
    pcs <- fread("EUR.eigenvec")
    # rename columns
    colnames(pcs) <- c("FID","IID", paste0("PC",1:6))
    # generate required table
    pheno <- merge(phenotype, covariate) %>%
        merge(., pcs)
    ``` 
### 2. Obtain HapMap3 SNPs
LDpred2 authors recommend restricting the analysis to only the HapMap3 SNPs

=== "load HapMap3 SNPs"

    ```R
    info <- readRDS(url("https://github.com/privefl/bigsnpr/raw/master/data-raw/hm3_variants.rds"))
    ```
### 3. Load and transform the summary statistic file 

=== "Load summary statistic file"
    
    ``` R
    # Read in the summary statistic file
    sumstats <- bigreadr::fread2("Height.QC.gz") 
    # LDpred 2 require the header to follow the exact naming
    names(sumstats) <-
        c("chr",
        "pos",
        "rsid",
        "a1",
        "a0",
        "n_eff",
        "beta_se",
        "p",
        "OR",
        "INFO",
        "MAF")
    # Transform the OR into log(OR)
    sumstats$beta <- log(sumstats$OR)
    # Filter out hapmap SNPs
    sumstats <- sumstats[sumstats$rsid%in% info$rsid,]
    ```

    !!! Warning
        
        Here, we know the exact ordering of the summary statistics file. 
        However, in many cases, the ordering of the summary statistics differ, 
        thus one must rename the columns according to their actual ordering

### 3. Calculate the LD matrix

=== "Genome Wide bed file"

    ```R
    # Get maximum amount of cores
    NCORES <- nb_cores()
    # Open a temporary file
    tmp <- tempfile(tmpdir = "tmp-data")
    on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
    # Initialize variables for storing the LD score and LD matrix
    corr <- NULL
    ld <- NULL
    # We want to know the ordering of samples in the bed file 
    fam.order <- NULL
    # preprocess the bed file (only need to do once for each data set)
    snp_readBed("EUR.QC.bed")
    # now attach the genotype object
    obj.bigSNP <- snp_attach("EUR.QC.rds")
    # extract the SNP information from the genotype
    map <- obj.bigSNP$map[-3]
    names(map) <- c("chr", "rsid", "pos", "a1", "a0")
    # perform SNP matching
    info_snp <- snp_match(sumstats, map)
    # Assign the genotype to a variable for easier downstream analysis
    genotype <- obj.bigSNP$genotypes
    # Rename the data structures
    CHR <- map$chr
    POS <- map$pos
    # get the CM information from 1000 Genome
    # will download the 1000G file to the current directory (".")
    POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
    # calculate LD
    for (chr in 1:22) {
        # Extract SNPs that are included in the chromosome
        ind.chr <- which(info_snp$chr == chr)
        ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
        # Calculate the LD
        corr0 <- snp_cor(
                genotype,
                ind.col = ind.chr2,
                ncores = NCORES,
                infos.pos = POS2[ind.chr2],
                size = 3 / 1000
            )
        if (chr == 1) {
            ld <- Matrix::colSums(corr0^2)
            corr <- as_SFBM(corr0, tmp)
        } else {
            ld <- c(ld, Matrix::colSums(corr0^2))
            corr$add_columns(corr0, nrow(corr))
        }
    }
    # We assume the fam order is the same across different chromosomes
    fam.order <- as.data.table(obj.bigSNP$fam)
    # Rename fam order
    setnames(fam.order,
            c("family.ID", "sample.ID"),
            c("FID", "IID"))
    ```

=== "Chromosome separated bed files"

    ```R
    # Get maximum amount of cores
    NCORES <- nb_cores()
    # Open a temporary file
    tmp <- tempfile(tmpdir = "tmp-data")
    on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
    # Initialize variables for storing the LD score and LD matrix
    corr <- NULL
    ld <- NULL
    # We want to know the ordering of samples in the bed file 
    info_snp <- NULL
    fam.order <- NULL
    for (chr in 1:22) {
        # preprocess the bed file (only need to do once for each data set)
        # Assuming the file naming is EUR_chr#.bed
        snp_readBed(paste0("EUR_chr",chr,".bed"))
        # now attach the genotype object
        obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,".rds"))
        # extract the SNP information from the genotype
        map <- obj.bigSNP$map[-3]
        names(map) <- c("chr", "rsid", "pos", "a1", "a0")
        # perform SNP matching
        tmp_snp <- snp_match(sumstats[sumstats$chr==chr,], map)
        info_snp <- rbind(info_snp, tmp_snp)
        # Assign the genotype to a variable for easier downstream analysis
        genotype <- obj.bigSNP$genotypes
        # Rename the data structures
        CHR <- map$chr
        POS <- map$pos
        # get the CM information from 1000 Genome
        # will download the 1000G file to the current directory (".")
        POS2 <- snp_asGeneticPos(CHR, POS, dir = ".")
        # calculate LD
        # Extract SNPs that are included in the chromosome
        ind.chr <- which(tmp_snp$chr == chr)
        ind.chr2 <- tmp_snp$`_NUM_ID_`[ind.chr]
        # Calculate the LD
        corr0 <- snp_cor(
                genotype,
                ind.col = ind.chr2,
                ncores = NCORES,
                infos.pos = POS2[ind.chr2],
                size = 3 / 1000
            )
        if (chr == 1) {
            ld <- Matrix::colSums(corr0^2)
            corr <- as_SFBM(corr0, tmp)
        } else {
            ld <- c(ld, Matrix::colSums(corr0^2))
            corr$add_columns(corr0, nrow(corr))
        }
        # We assume the fam order is the same across different chromosomes
        if(is.null(fam.order)){
            fam.order <- as.data.table(obj.bigSNP$fam)
        }
    }
    # Rename fam order
    setnames(fam.order,
            c("family.ID", "sample.ID"),
            c("FID", "IID"))
    ```

### 4. Perform LD score regression

===  "Perform LD score regression"

    ```R
    df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
    ldsc <- snp_ldsc(   ld, 
                        length(ld), 
                        chi2 = (df_beta$beta / df_beta$beta_se)^2,
                        sample_size = df_beta$n_eff, 
                        blocks = NULL)
    h2_est <- ldsc[["h2"]]
    ```

### 5. Calculate the null R2

=== "Calculate the null R2 (quantitative trait)"
    
    ```R
    # Reformat the phenotype file such that y is of the same order as the 
    # sample ordering in the genotype file
    y <- pheno[fam.order, on = c("FID", "IID")]
    # Calculate the null R2
    # use glm for binary trait 
    # (will also need the fmsb package to calculate the pseudo R2)
    null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~Sex+", .) %>%
        as.formula %>%
        lm(., data = y) %>%
        summary
    null.r2 <- null.model$r.squared
    ```

=== "Calculate the null R2 (binary trait)"
    
    ```R
    library(fmsb)
    # Reformat the phenotype file such that y is of the same order as the 
    # sample ordering in the genotype file
    y <- pheno[fam.order, on = c("FID", "IID")]
    # Calculate the null R2
    # use glm for binary trait 
    # (will also need the fmsb package to calculate the pseudo R2)
    null.model <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~Sex+", .) %>%
        as.formula %>%
        glm(., data = y, family=binomial) %>%
        summary
    null.r2 <- fmsb::NagelkerkeR2(null.model)
    ```

!!! important
    Scripts for binary trait analysis only serve as a reference as we have not simulate any binary traits. 
    In addition, Nagelkerke $R^2$ is biased when there are ascertainment of samples. For more information, please refer to [this paper](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.21614)

### 6. Obtain LDpred adjusted beta
=== "infinitesimal model"
    
    ```R
    beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
    ```

=== "grid model"

    ```R
    # Prepare data for grid model
    p_seq <- signif(seq_log(1e-4, 1, length.out = 17), 2)
    h2_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
    grid.param <-
        expand.grid(p = p_seq,
                h2 = h2_seq,
                sparse = c(FALSE, TRUE))
    # Get adjusted beta from grid model
    beta_grid <-
        snp_ldpred2_grid(corr, df_beta, grid.param, ncores = NCORES)
    ```

=== "auto model"

    ```R
    # Get adjusted beta from the auto model
    multi_auto <- snp_ldpred2_auto(
        corr,
        df_beta,
        h2_init = h2_est,
        vec_p_init = seq_log(1e-4, 0.9, length.out = NCORES),
        ncores = NCORES
    )
    beta_auto <- sapply(multi_auto, function(auto)
        auto$beta_est)
    ```

### 7. Obtain model PRS
#### Using Genome wide bed file

=== "infinitesimal model"

    ```R
    if(is.null(obj.bigSNP)){
        obj.bigSNP <- snp_attach("EUR.QC.rds")
    }
    genotype <- obj.bigSNP$genotypes
    # calculate PRS for all samples
    ind.test <- 1:nrow(genotype)
    pred_inf <- big_prodVec(    genotype,
                                beta_inf,
                                ind.row = ind.test,
                                ind.col = info_snp$`_NUM_ID_`)
    ```

=== "grid model"

    ```R
    if(is.null(obj.bigSNP)){
        obj.bigSNP <- snp_attach("EUR.QC.rds")
    }
    genotype <- obj.bigSNP$genotypes
    # calculate PRS for all samples
    ind.test <- 1:nrow(genotype)
    pred_grid <- big_prodMat(   genotype, 
                                beta_grid, 
                                ind.col = info_snp$`_NUM_ID_`)
    ```

=== "auto model"

    ```R
    if(is.null(obj.bigSNP)){
        obj.bigSNP <- snp_attach("EUR.QC.rds")
    }
    genotype <- obj.bigSNP$genotypes
    # calculate PRS for all samples
    ind.test <- 1:nrow(genotype)
    pred_auto <-
        big_prodMat(genotype,
                    beta_auto,
                    ind.row = ind.test,
                    ind.col = info_snp$`_NUM_ID_`)
    # scale the PRS generated from AUTO
    pred_scaled <- apply(pred_auto, 2, sd)
    final_beta_auto <-
        rowMeans(beta_auto[,
                    abs(pred_scaled -
                        median(pred_scaled)) <
                        3 * mad(pred_scaled)])
    pred_auto <-
        big_prodVec(genotype,
            final_beta_auto,
            ind.row = ind.test,
            ind.col = info_snp$`_NUM_ID_`)
    ```

#### Using chromosome separated bed files

=== "infinitesimal model"

    ```R
    pred_inf <- NULL
    for(chr in 1:22){
        obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,".rds"))
        genotype <- obj.bigSNP$genotypes
        # calculate PRS for all samples
        ind.test <- 1:nrow(genotype)
        # Extract SNPs in this chromosome
        chr.idx <- which(info_snp$chr == chr)
        ind.chr <- info_snp$`_NUM_ID_`[chr.idx]
        tmp <- big_prodVec(genotype,
                                beta_inf[chr.idx],
                                ind.row = ind.test,
                                ind.col = ind.chr)
        if(is.null(pred_inf)){
            pred_inf <- tmp
        }else{
            pred_inf <- pred_inf + tmp
        }
    }
    ```

=== "grid model"

    ```R
    pred_grid <- NULL
    for(chr in 1:22){
        obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,"_.rds"))
        genotype <- obj.bigSNP$genotypes
        # calculate PRS for all samples
        ind.test <- 1:nrow(genotype)
        # Extract SNPs in this chromosome
        chr.idx <- which(info_snp$chr == chr)
        ind.chr <- info_snp$`_NUM_ID_`[chr.idx]

        tmp <- big_prodMat( genotype, 
                            beta_grid[chr.idx], 
                            ind.col = ind.chr)

        if(is.null(pred_grid)){
            pred_grid <- tmp
        }else{
            pred_grid <- pred_grid + tmp
        }
    }
    ```

=== "auto model"

    ```R
    pred_auto <- NULL
    for(chr in 1:22){
        obj.bigSNP <- snp_attach(paste0("EUR_chr",chr,"_.rds"))
        genotype <- obj.bigSNP$genotypes
        # calculate PRS for all samples
        ind.test <- 1:nrow(genotype)
        # Extract SNPs in this chromosome
        chr.idx <- which(info_snp$chr == chr)
        ind.chr <- info_snp$`_NUM_ID_`[chr.idx]

        tmp <-
            big_prodMat(genotype,
                        beta_auto[chr.idx],
                        ind.row = ind.test,
                        ind.col = ind.chr)
        # scale the PRS generated from AUTO
        pred_scaled <- apply(tmp, 2, sd)
        final_beta_auto <-
            rowMeans(beta_auto[chr.idx,
                        abs(pred_scaled -
                            median(pred_scaled)) <
                            3 * mad(pred_scaled)])
        tmp <-
            big_prodVec(genotype,
                final_beta_auto,
                ind.row = ind.test,
                ind.col = ind.chr)
        if(is.null(pred_auto)){
            pred_auto <- tmp
        }else{
            pred_auto <- pred_auto + tmp
        }
    }
    ```

### 8. Get the final performance of the LDpred models

=== "infinitesimal model"

    ```R
    reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~PRS+Sex+", .) %>%
        as.formula
    reg.dat <- y
    reg.dat$PRS <- pred_inf
    inf.model <- lm(reg.formula, dat=reg.dat) %>%
        summary
    (result <- data.table(
        infinitesimal = inf.model$r.squared - null.r2,
        null = null.r2
    ))
    ```

=== "grid model"

    ```R
    reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~PRS+Sex+", .) %>%
        as.formula
    reg.dat <- y
    max.r2 <- 0
    for(i in 1:ncol(pred_grid)){
        reg.dat$PRS <- pred_grid[,i]
        grid.model <- lm(reg.formula, dat=reg.dat) %>%
            summary  
        if(max.r2 < grid.model$r.squared){
            max.r2 <- grid.model$r.squared
        }
    }
    (result <- data.table(
        grid = max.r2 - null.r2,
        null = null.r2
    ))
    ```

=== "auto model"

    ```R
    reg.formula <- paste("PC", 1:6, sep = "", collapse = "+") %>%
        paste0("Height~PRS+Sex+", .) %>%
        as.formula
    reg.dat <- y
    reg.dat$PRS <- pred_auto
    auto.model <- lm(reg.formula, dat=reg.dat) %>%
        summary
    (result <- data.table(
        auto = auto.model$r.squared - null.r2,
        null = null.r2
    ))
    ```

??? note "How much phenotypic variation does the PRS from each model explain?"
    Infinitesimal = 0.0100
    
    Grid Model = 0.00180

    Auto Model = 0.171