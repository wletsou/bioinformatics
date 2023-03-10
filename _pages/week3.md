---
layout: "single"
title: "Genome-wide association testing"
toc: true
toc_label: "Contents"
permalink: /assignments/week3/
---

In this assignemnt we will combine the data-cleaning steps we learned in [Week 1](https://wletsou.github.io/bioinformatics/assignments/week1) and [Week 2](https://wletsou.github.io/bioinformatics/assignments/week2) to a simulated case-control study of disease.  We will create a simulated dataset from one of the 1KG populations, declare some alleles to be risk alleles, and find the association of each SNP with the simulated phenotype.  First we will give an overview of logistic regression and linear mixed models.

### Statistics ###

#### The odds ratio ####

The logic of any genetic association study is to see if an allele is enriched in subjects affected with disease.  In other words, we want to see if the allele is *associated* with disease.  We do this by collecting may individuals with disease and a similar number of healthy controls from the general popuation.  An individual's risk is the probability \\(p=P\left(\text{Disease}\mid\text{Allele}\right)\\) of developing disease over a lifetime, and the *relative risk* or RR is the ratio of the risk to carriers to the risk to non-carriers of the allele.  However, the RR is difficult to measure because it involves waiting a long time for a potentially smaller number of cases to develop.  If we pre-select cases and controls, we can instead calculate the probability \\(P\left(\text{Allele}\mid\text{Disease}\right)\\) of carrying the allele conditioned on an individual's disease status.  Related to probability is the *odds* \\(\frac{p}{1-p}\\) of the the event, and it turns out that the *odds ratio*\\[\text{OR}=\frac{\frac{P\left(\text{Allele}\mid\text{Cases}\right)}{1-P\left(\text{Allele}\mid\text{Cases}\right)}}{\frac{P\left(\text{Allele}\mid\text{Controls}\right)}{1-P\left(\text{Allele}\mid\text{Controls}\right)}}\tag{1}\\]is invariant to whether we collect subjects prospectively or retrospectively.  Hence we use the OR as a convenient measure the effect of an allele on disease risk in case-control studies. However, the crude measure (1) cannot be adjusted for other factors like sex, age, and race which may also have an effect on disease risk.  For this type of analysis, we need logistic regression.

#### The logistic model of disease risk ####

Normally we test the association between two variables using regression analysis, but doing so requires that both the dependent and independent variables be continuous.  In genetic association studies, neither the outcome (disease) nor the predictor (number of risk alleles) is continuous.  However, an individual's unobserved probability \\(p\\) is a continuous variable that ranges between 0 and 1; furthermore, the individual's *log-odds* \\(\log{\frac{p}{1-p}}\\) of disease is a continuous value that ranges between \\(-\infty\\) and \\(+\infty\\).  Thus if we assume that the log-odds of disease can be represented by the equation \\[\log{\frac{p}{1-p}}=\beta_0+\beta_1X_1,\tag{2}\\] where \\(X_1\\) is number of risk alleles and \\(\beta_1\\) is the *log-odds-ratio*, then we can in principle fit a line and estimate its slope.  This slope \\(\beta_1\\) would then be interpretted as the multiplicative increase in the log-odds of disease.

Now, since we cannot observe \\(p_i\\) for each subject \\(i\\), we cannot actually fit (2) using linear regression.  We can, however, use the concept of *maximum likelihood*.  In statistics, the likelihood of a disease model like (2) is the probability of the data being generated the model, or\\[\mathcal{L}\left(\text{Model}\mid\text{Data}\right)=P\left(\text{Data}\mid\text{Model}\right).\tag{3}\\]Here the model is the set of parameters \\(\beta_0,\beta_1\\) required to predict disease risk, and we can find estimates for the parameters by maximizing (3), i.e., by finding the model most consistent with the data.  Furthermore, if the likelihood function has approximately the shape of a normal distribution, then we can estimate the statistical significance of our estimates by computing their standard error, got from the *curvature* or second derivative of \\(\mathcal{L}\\) near the \\(\beta\\) which maximize it.

The binomial distribution is a good approximation to the normal distribution, so if we model the likelihood of the observed data as\\[\mathcal{L}=\prod_ip_i^{y_i}\left(1-p_i\right)^{1-y_i},\\] where \\(y_i=0,1\\) is an indicator of disease status, then the *log-likelihood* is\\[\begin{align}\ell&=\sum_iy_i\log{p_i}+\left(1-y_i\right)\log{\left(1-p_i\right)}\\\\\\ &=\sum_iy_i\log{\frac{p_i}{1-p_i}}+\log{\left(1-p_i\right)}\end{align}\\]or using (2)\\[\ell=\sum_iy_i\left(\beta_0+\beta_1X_{i1}\right)-\log{\left(1+e^{\beta_0+\beta_1X_{i1}}\right)}.\tag{4}\\]Eq. (4) is a function of the parameter \\(\beta_1\\).  Thus we can find the *maximum-likelihood estimate* \\(\hat{\beta_1}\\) of \\(\beta_1\\) by solving \\(\frac{\partial \ell}{\partial \beta_1}=0\\) and get its standard error \\(\frac{-1}{\frac{\partial^2\ell}{\partial \beta_1^2}\bigr\rvert_{\beta_1=\hat{\beta_1}}}\\)by evaluating the curvature of the log-likelihood at the best estimate of \\(\beta_1\\).  From these we can also get an estimate of the statistical significance.

#### Linear mixed models ####

In genome-wide association studies (GWAS) we'd like to estimate the odds ratio \\(e^{\beta_1}\\) for each SNP to see if any SNPs are associated with disease.  However, there are millions of SNPs and only thousands of subjects, so we cannot fit all the parameters simultaneously.  Instead, one SNP effect \\(\beta_1\\) is fit at a time together with other the covariate effects \\(\beta_j\\) against a background of the composite effect of all the remaining SNPs together, so that our model has two components:\\[Y_i=\log{\frac{p_i}{1-p_i}}=\sum_jX_{ij}\beta_j+\sum_jZ_{ij}u_j+\varepsilon_i.\tag{5}\\]Here, \\(\mathbf{Z}\\) is an \\(n\times m\\) matrix of the (standardized) genotypes of \\(n\\) individuals at \\(m\\) SNPs and \\(u_j\\) is the effect of SNP \\(j\\) on the log-odds for individual \\(i\\).  The mixed-model framework assumes the \\(u\\) come from a normal distribution with mean 0 and standard deviation \\(\sigma\\), so that each variant has but a small effect on disease risk.  

Now, the variance of the log-odds \\(Y\\) about its mean \\(\mathbf{X}\beta\\) becomes\\[\begin{align}\left(Y-\mathbf{X}\beta\right)\left(Y-\mathbf{X}\beta\right)^T&=\mathbf{Z}uu^T\mathbf{Z}^T+\varepsilon\varepsilon^T\\\\\\ &=\left(\mathbf{Z}\mathbf{Z}^T+\mathbf{I}\right)\sigma^2=\mathbf{V}.\end{align}\\]Rearranging and differntiating with respect to \\\(\beta_k\\) obtains (with summation over \\(i\\), \\(j\\), and \\(l\\))\\[\begin{align}V_{li}^{-1}X_{ik}\left(Y_l-X_{lj}\beta_j\right)^T+V_{il}^{-1}\left(Y_l-X_{lj}\beta_j\right)X^T_{ki}&=0\\\\\\ \left(X^T_{ki}V_{il}^{-1}Y_l-X^T_{ki}V_{il}^{-1}X_{lj}\beta_j\right)^T+\left(X^T_{ki}V_{il}^{-1}Y_l-X^T_{ki}V_{il}^{-1}X_{lj}\beta_j\right)&=0,\end{align}\\]since \\(\mathbf{V}\\) and hence \\(\mathbf{V}^{-1}\\) is a symmetric matrix.  And because the \\(k\\)<sup>th</sup> entry of a transposed vector is equal to the \\(k\\)<sup>th</sup> entry of the original vector, we get the maximum-likelihood solution\\[\mathbf{X}^T\left(\mathbf{I}+\mathbf{Z}\mathbf{Z}^T\right)^{-1}\mathbf{X}\hat{\beta}=\mathbf{X}^T\left(\mathbf{I}+\mathbf{Z}\mathbf{Z}^T\right)^{-1}Y,\tag{6}\\]where \\(\frac{1}{m}\mathbf{Z}\mathbf{Z}^T\\) is the genomic relationship matrix (GRM) we used to compute principle components in [Week 1](https://wletsou.github.io/assignments/week1).  Thus we can estimatimate the *fixed effects* \\(\beta\\)&mdash;including the SNP effect \\(\beta_1\\)&mdash;and their standard errors without fitting the *random effects* of every other SNP simultaneously.

### Simulating genotypes and phenotypes ###

Today we will be simulating a case-control study using the package <kbd>sim1000G</kbd> and the logistic model.  Then we will perform a "genome-wide" association study to discover the SNPs that affect disease risk.  We will need the following packages

```
library(sim1000G)
library(SNPRelate)
library(GENESIS)
library(GWASTools)
library(SeqVarTools)
library(data.table)
```

#### Importing data and simulating genotypes ####

First get the file of individuals and their population codes:

```
indivs <- read.table("path/to/file/CHB+YRI+CEU.txt",header = FALSE)
colnames(indivs) <- c("id","pop")
```

Then import the 1KG vcf file for one of the chromosomes using the sim1000G function <kbd>readVCF</kbd>.  Take 100 variants with minor allele frequencies between 0.10 and 0.90 to simulate genotypes comprising commmon variants across the chromosome. 

```
vcf <- readVCF("OneDrive - New York Institute of Technology/Courses/BIOL 350 Spring 2023/CHB+YRI+CEU.chr1.vcf.gz", maxNumberOfVariants = 100 , min_maf = 0.10 , max_maf = 0.90)
```

Make a table of variants as we did [previously](https://wletsou.github.io/bioinformatics/assignments/week2/#making-a-table-of-variants); call it <kbd>variants</kbd>.  We will use this table to make a map file.

Next simulate 2000 individuals from one of the three populations, also as we did [before](https://wletsou.github.io/bioinformatics/assignments/week2/#simulating-individuals).  It is important that we use a large number of subjects so that we have enough statistical power to detect ORs near 1.  At the end of this step you should have two data tables <kbd>dt.gt1.allele</kbd> and <kbd>dt.gt2.allele</kbd> of alleles on each of the two chromosomes of 2000 subjects.  We will use these tables to make a ped file.

You should also have two 2000-by-100 data tables <kbd>dt.gt1</kbd> and <kbd>dt.gt2</kbd> of 0's and 1's, corresponding to the number of copies of the alternative allele on each chromosome.  Turn these into a matrix from which we can compute allele frequencies:

```
genotype.matrix <- as.matrix((dt.gt1 + dt.gt2)[,1:ncol(dt.gt1),,with = FALSE])
colnames(genotype.matrix) <- 1:ncol(genotype.matrix)
rownames(genotype.matrix) <- 1:nrow(genotype.matrix)
```

You can compute the frequency of each minor allele by applying the <kbd>sum</kbd> function to the columns of <kbd>genotype.matrix</kbd> and dividing by *twice* the number of rows.

```
apply(genotype.matrix[,disease.snps],2,sum) / (2 * nrow(genotype.matrix)) # minor allele frequencies
```

These values are called the *minor allele frequencies (MAF)* and will be used to simulate phenotypes.
#### Simulating phenotypes ####

We will simulate a binary phenotype by picking three alleles at random to increase the risk of disease.  We will weight the log-OR \\(\beta\\) by allele frequency according to the formula \\[\beta=\log{\left(10\right)}\times-\log_{10}{\left(\text{MAF}\right)}.\tag{7}\\]Eq. (7) says that the odds of disease increases by a factor of 10 for alleles which occur on only 10% of chromosomes.  This is a much stronger effect than most SNPs show, but implementing it will increase the power of our small study.  Choose and view the effect sizes using:

```
disease.snps <- sample(1:nrow(variants),3,replace = FALSE) # randomly select disease SNPs
beta <- log(10/1) * -log10(apply(genotype.matrix[,disease.snps],2,sum) / (2 * nrow(genotype.matrix))) # weight log-OR by allele frequency
rbind(beta,freq = apply(genotype.matrix[,disease.snps],2,sum) / (2 * nrow(genotype.matrix))) # view frequency and beta values
```

According to Eq. (2), the number of copies of each allele will increase the log-odds of disease for each subject by \\(X_{i1}\beta_1+X_{i2}\beta_2+X_{i3}\beta_3\\).  Carry out this matrix product by

```
logits <- genotype.matrix[,disease.snps] %*% (beta) # genotype matrix multiplied by column of beta values
```

The baseline probability of disease in our study will be 50%, because there will be an approximately equal number of cases and controls, corresponding to a mean log-odds of \\(\log{\left(\frac{0.5}{0.5}\right)}=0\\).  In order that the mean log-odds should be zero, simply subtract the mean of <kbd>logits</kbd> from <kbd>logits</kbd> itself:

```
logits <- logits + (0 - mean(logits)) # log-odds in a balanced case-control study
```

Now we can convert odds into probability by rearrangeing Eq. (2)\\[p_i=\frac{e^{\beta_0+X_i1\beta_1+X_i2\beta_2+X_i3\beta_3}}{1+e^{\beta_0+X_i1\beta_1+X_i2\beta_2+X_i3\beta_3}}.\tag{3}\\]In R, we can do

```
probs <- exp(logits) / (1 + exp(logits)) # disease probability
```

Now we can simulate the binary phenotype by drawing a random number between 0 and 1.  If the number is less than <kbd>probs[i]</kbd> for subject <kbd>i</kbd>, then subject <kbd>i</kbd> is a case; otherwise it is a control.  Simulate the phenotypes by
  
```
pheno <- as.numeric(runif(nrow(probs)) <= probs) # generate disease phenotypes
```
  
Verify that the mean of <kbd>pheno</kbd> is close to 0.5.
  
#### Saving your genotype and phenotypes ####
  
[Last time](https://wletsou.github.io/bioinformatics/assignments/week2/#generating-a-gds-file-from-your-ped-and-map-files) we generated map and ped files from our data.  The map file can be simply got from your <kbd>variants</kbd> table using

```
write.table(data.frame(chr = variants$chr,rsid = variants$rsid,X = 0,pos = variants$pos),col.names = FALSE,row.names = FALSE,quote = FALSE,file = "path/to/file/GWAS.simulation.map") # variant map file
```

To make a ped file with phenotype information, we first need to make a six-column data frame to prepend the genotypes:

```
fam <- data.frame(fid = row.names(genotype.matrix),id = row.names(genotype.matrix), mother = 0,father = 0,sex = sample(c(0,1),nrow(genotype.matrix),replace = TRUE),phenotype = pheno) # variant map file
```

This <kbd>fam</kbd> object is the pedigree of a cohort of unrelated individuals having an approximately equal distributions of cases and controls and men and women.  Now, prepend these six columns to the vector of each subject's alleles on each of its two chromosomes:

```
write.table(cbind(fam,data.frame(cbind(dt.gt1.allele,dt.gt2.allele))[,rep(1:ncol(dt.gt1.allele.ceu),each = 2) + (0:1) * ncol(dt.gt1.allele.ceu)]),col.names = FALSE,row.names = FALSE,quote = FALSE,file = "OneDrive - New York Institute of Technology/Courses/BIOL 350 Spring 2023/CEU.simulation.ped") # genotype .ped file
```

This will be our ped file.  Verify that the object you created has 2000 rows and 206 columns (what do these numbers represent?).

### Association testing ###

In the next part of the analysis, our goal is to use the <kbd>GENESIS</kbd> package to estimate the log-OR associate with each SNP in our dataset.  According to Eq. (6), we need the GRM to do this properly.  GENESIS requires that the GRM be in a specific form, which we will achieve by following the protocol for PC-AiR and PC-Relate.  First of all, we need principal components from the implementation of KING in SNPRelate.

#### Principal components analysis and LD pruning ####
  
