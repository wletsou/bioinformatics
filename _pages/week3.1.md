---
layout: "single"
title: "Genome-wide association testing (2)"
toc: true
toc_label: "Contents"
sidebar:
  nav: "side"
permalink: /assignments/week3.1/
---

### To turn in: ###

1. Import the [map](https://raw.githubusercontent.com/wletsou/bioinformatics/master/docs/CEU.simulation.chr1.map) and [ped](https://raw.githubusercontent.com/wletsou/bioinformatics/master/docs/CEU.simulation.chr1.ped) files into R.&nbsp; Find the code in [Week 3](https://wletsou.github.io/bioinformatics/assignments/week3) to convert the ped and map file into a gds file using the <kbd>SNPRelate</kbd> package.&nbsp; Paste a summary of the gds output.
2. The sixth column of your ped object is a vector of phenotypes, with 1 indicating control and 2 indicating case.&nbsp; Determine the number of cases and controls in the study.
3. Import the [matrix](https://raw.githubusercontent.com/wletsou/bioinformatics/master/docs/CEU.simulation.disease_snp_matrix.txt) containg the genotypes (i.e., number of alternative alleles) of the subjects at three disease SNPs.&nbsp; Find and adapt the code to compute the allele frequency of each SNP.
4. Estimate the odds ratio\\[\text{OR}=\frac{P\left(\text{Allele} = 1\mid \text{Case}\right)P\left(\text{Allele} = 0\mid \text{Control}\right)}{P\left(\text{Allele} = 0\mid \text{Case}\right)P\left(\text{Allele} = 1\mid \text{Control}\right)}\\]for each SNP by extracting the rows of the genotype matrix corresponding to cases and controls.
5. From the following Manahattan plot<figure class="align-center"><img src="https://github.com/wletsou/bioinformatics/raw/master/docs/Simulation%20Manahattan%20plot.png" alt="Simulation Manhattan plot"><center><i>Simulation Mahattan plot</i></center></figure>discuss which SNPs you think were detected in the study.
