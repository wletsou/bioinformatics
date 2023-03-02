---
layout: "single"
title: "Genome-wide association testing"
toc: true
toc_label: "Contents"
permalink: /assignments/week3/
---

In this assignemnt we will combine the data-cleaning steps we learned in [Week 1](https://wletsou.github.io/assignments/week1) and [Week 2](https://wletsou.github.io/assignments/week2) to a simulated case-control study of disease.  We will create a simulated dataset from one of the 1KG populations, declare some alleles to be risk alleles, and find the association of each SNP with the simulated phenotype.  First we will give an overview of logistic regression and linear mixed models.

### Statistics ###

#### The logistic model of disease risk ####

The logic of any genetic association study is to see if an allele is enriched in subjects affected with disease.  In other words, we want to see if the allele is *associated* with disease.  Normally we test the association between two variables using regression analysis, but doing so requires that both the dependent and independent variables be continuous.  In genetic association studies, neither the outcome (disease) nor the predictor (number of risk alleles) is continuous.  However, an individual's unobserved probability \\(p\\) is a continuous variable that ranges between 0 and 1; furthermore, the individual's *log-odds* \\(\log{\frac{p}{1-p}}\\) of disease is a continuous value that ranges between \\(-\infty\\) and \\(+\infty\\).  Thus if we assume that the log-odds of disease can be represented by the equation \\[\log{\frac{p}{1-p}}=\beta_0+\beta_1X_1,\tag{1}\\] where \\(X_1\\) is number of risk alleles and \\(\beta_1\\) is the *log-odds-ratio*, then we can in principle fit a line and estimate its slope.  This slope \\(\beta_1\\) would then be interpretted as the multiplicative increase in the log-odds of disease.
