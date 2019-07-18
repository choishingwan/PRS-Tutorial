# Overview 
The aim of this tutorial is to provide a step by step guide as to how to perform basic polygenic risk score analyses, therefore
allow for a better understanding of the inner mechanism implemented in most polygenic risk score software. 

The tutorial is separated into four main sections

1. [How to perform basic filtering on the summary statistic file (Base)](base.md)
2. [How to perform quanilty controls on the target genotype file](target.md)
3. [The details steps involved in calculating PRS (using `plink`)](plink.md)
4. [Visualizing PRS results](plink_visual.md)

We also provided a brief example on how to perform PRS using the three polygenic risk score software: [PRSice-2](prsice.md), [LDpred](ldpred.md) and [lassosum](lassosum.md)

If you are only interested in how to perform PRS, you can directly skipped to [step 3](plink.md). Links to download the required data are provided under each section.

!!! note

    This tutorial is based on linux and OS X systems. Users who would like to perform polygenic risk score analyses
    on their Windows machine will need to change some of the commands.

!!! note
    Throughout the tutorial, you might see some codes with tab on top:

    ```bash tab="A"
    echo "Tab A"
    ```

    ```bash tab="B"
    echo "Tab B"
    ```

    You can click on the tab to change to relevant codes (e.g. different operation system)

# Requirements
To follow the tutorial, you will need the following programs installed:

1. [R](https://www.r-project.org/) (**version 3.2.3+**)
2. [PLINK 1.9](https://www.cog-genomics.org/plink2)

## Citation
If you find this tutorial helpful, then please cite:

!!! important "Citation"

    A guide to performing Polygenic Risk Score analyses, 
    
    Shing Wan Choi, Timothy Shin Heng Mak, Paul O'Reilly 
    
    bioRxiv 416545 (2018). doi:10.1101/416545
