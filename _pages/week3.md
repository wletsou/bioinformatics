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

Now, since we cannot observe \\(p_i\\) for each subject \\(i\\), we cannot actually fit (1) using linear regression.  We can, however, use the concept of *maximum likelihood*.  In statistics, the likelihood of a disease model like (1) is the probability of the data being generated the model, or\\[\mathcal{L}\left(\text{Model}\mid\text{Data}\right)=P\left(\text{Data}\mid\text{Model}\right).\tag{2}\\]Here the model is the set of parameters \\(\beta_0,\beta_1\\) required to predict disease risk, and we can find estimates for the parameters by maximizing (2), i.e., by finding the model most consistent with the data.  Furthermore, if the likelihood function has approximately the shape of a normal distribution, then we can estimate the statistical significance of our estimates by computing their standard error, got from the *curvature* or second derivative of \\(\mathcal{L}\\) near the \\(\beta\\) which maximize it.

The binomial distribution is a good approximation to the normal distribution, so if we model the likelihood of the observed data as\\[\mathcal{L}=\prod_ip_i^{Y_i}\left(1-p_i\right)^{1-Y_i},\\] where \\(Y_i=0,1\\) is an indicator of disease status, then the *log-likelihood* is\\[\begin{align}\ell&=\sum_iY_i\log{p_i}+\left(1-Y_i\right)\log{1-p_i}\\ &=\sum_iY_i\log{\frac{p_i}{1-p_i}}+\log{1-p_i}\end{align}\\]or using (1)\\[\ell=\sum_iY_i\left(\beta_0+\beta_1X_{1i}\right)-\log{\left(1+e^{\beta_0+\beta_1X_{1i}}\right)}.\tag{3}\\]Eq. (3) is a function of the parameter \\(\beta_1\\).  Thus we can find the *maximum-likelihood estimate* \\(\hat{\beta_1}\\
) of \\(\beta_1\\) by solving \\(\frac{\partial \ell}{\partial \beta_1}=0\\) and get its standard error \\(\frac{-1}{\frac{\partial^2\ell}{\partial \beta_1^2}\bigr\rvert_{\beta_1=\hat{\beta_1}}}\\)by evaluating the curvature of the log-likelihood at the best estimate of \\(\beta_1\\).  From these we can also get an estimate of the statistical significance.

#### Linear mixed models ####


