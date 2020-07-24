# Overview 
This tutorial provides a step-by-step guide to performing basic polygenic risk score (PRS) analyses and accompanies our [PRS Guide paper](https://doi.org/10.1101/416545). 
The aim of this tutorial is to provide a simple introduction of PRS analyses to those new to PRS, while equipping existing users with a better understanding of the processes and implementation "underneath the hood" of popular PRS software.

The tutorial is separated into four main sections and reflects the structure of our [guide paper](https://doi.org/10.1101/416545): 
the first two sections on QC correspond to Section 2 of the paper and constitute a 'QC checklist' for PRS analyses, the third section on calculating PRS (here with examples using [PLINK](plink.md), [PRSice-2](prsice.md), [LDpred-2](ldpred.md) and [lassosum](lassosum.md)) corresponds to Section 3 of the paper, while the fourth section, which provides some examples of visualising PRS results, accompanies Section 4 of the paper.

1. [Quality Control (QC) of Base Data](base.md)
2. [Quality Control (QC) of Target Data](target.md)
3. [Calculating and analysing PRS](plink.md)
4. [Visualising PRS Results](plink_visual.md)

We will be referring to our [guide paper](https://doi.org/10.1101/416545) in each section and so you may find it helpful to have the paper open as you go through the tutorial.

If you are only interested in how to perform PRS on previously QC'ed data then you can skip to [Step 3](plink.md). Links to download the required data are provided under each section.
!!! warning

    Data used in this tutorial are simulated and is intended for demonstration purposes only. Result from this tutorial will not reflect the true performance of different software. 

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

    You can click on the tab to change to an alternative code (eg. to a different operation system)
# Datasets
1. [Base data](https://drive.google.com/file/d/1RWjk49QNZj9zvJHc9X_wyZ51fdy6xQjv/view?usp=sharing): Modified summary statistics file from the GIANT consortium study on height
2. [Target data](https://drive.google.com/file/d/1uhJR_3sn7RA8U5iYQbcmTp6vFdQiF4F2/view?usp=sharing): Simulated data based on the 1000 Genomes Project European samples

# Requirements
To follow the tutorial, you will need the following programs installed:

1. [R](https://www.r-project.org/) (**version 3.2.3+**)
2. [PLINK 1.9](https://www.cog-genomics.org/plink2)

## Citation
If you find this tutorial helpful for a publication, then please consider citing:

!!! important "Citation"

    Choi SW, Mak TSH, O'Reilly PF.  A guide to performing Polygenic Risk Score analyses. 
    bioRxiv 416545 (2018). https://doi.org/10.1101/416545
