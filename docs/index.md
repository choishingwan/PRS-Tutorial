# Overview 
This tutorial provides a step-by-step guide to performing basic polygenic risk score (PRS) analyses and accompanies our paper on the topic (https://doi.org/10.1101/416545). Our hope is that this tutorial will allow new users to get started on PRS analyses and provide existing users with a better understanding of the processes and implemention underlying popular PRS software.

The tutorial is separated into four main sections

1. [QC of Base GWAS Summary Data](base.md)
2. [QC of Target Individual-Level Data](target.md)
3. [Running PRS Analyses (using `plink`)](plink.md)
4. [Visualizing PRS Results](plink_visual.md)

We provide brief examples of performing PRS analyses using four software for polygenic risk score analyses: [PLINK](plink.md), [PRSice-2](prsice.md), [LDpred](ldpred.md) and [lassosum](lassosum.md)

If you are only interested in how to perform PRS on previously QC'ed data then you can directly skip to [step 3](plink.md). Links to download the required data are provided under each section.

!!! note

    This tutorial is written for Linux and OS X operating systems. 
    Windows users will need to change some commands accordingly.

!!! note
    Throughout the tutorial you will see tabs above some of the code:

    ```bash tab="A"
    echo "Tab A"
    ```

    ```bash tab="B"
    echo "Tab B"
    ```

    You can click on the tab to change to an alternative code (e.g to a different operation system)

# Requirements
To follow the tutorial, you will need the following programs installed:

1. [R](https://www.r-project.org/) (**version 3.2.3+**)
2. [PLINK 1.9](https://www.cog-genomics.org/plink2)

## Citation
If you find this tutorial helpful for a publication, then please consider citing:

!!! important "Citation"

    Choi SW, Mak TSH, O'Reilly PF.  A guide to performing Polygenic Risk Score analyses. bioRxiv 416545 (2018). https://doi.org/10.1101/416545
