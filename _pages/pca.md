---
layout: "single"
title: "Principal components analysis of 1000GP populations"
toc: true
toc_label: "Contents"
permalink: /assignments/pca/
---

### Importing the data

First we'll see how individuals can be separated by genetic ancestry using principal components analysis.  First we'll need to load the libraries <kbd>SNPRelate</kbd> and <kbd>SeqArray</kbd>:

```
library(SNPRelate)
library(SeqArray)
```

Download one of the vcf files from the [Files](https://wletsou.github.io/bioinformatics/files) page of a single chromosome containing just the CEU, YRI, and CHB populations.  Once you have the file, store its name as a variable:

```
vcf <- "path/to/file.vcf.gz"
```

Now we'll convert the vcf format to gds format.  I recommend simply changing the vcf.gz extension to gds.  This may take a minute to complete.  Then we'll import the gds file as a gds object.

```
seqVCF2GDS(vcf.fn = vcf,"path/to/file.gds")
genofile <- seqOpen("path/to/file.gds")
```

You can can see the various fields under <kbd>genofile</kbd> by printing it.  To access the data in one of the fields, do

```
seqGetData(genofile,"sample.id") # view the sample ids
```

where the name of the field is enclosed in quotes.

### Pricipal components analysis

The populations in our dataset can be separated into clusters based on their genotypes.  The inferred groups help control for confoundng due to ancestry and are also more reliable than self-reported race in association studies.  To see how it works, suppose \\(\mathbf{X}\\) is an \\(n\times m\\) (standardized) genotype matrix with individuals down the rows and (standardized) genotypes across the columns.  Principal components analysis (PCA) says we can find an \\(n\times m\\) matrix \\(\mathbf{V}\\), a diagonal \\(n\times n\\) matrix \\(\mathbf{\Sigma}\\), and an \\(n\times n\\) matrix \\(\mathbf{U}\\) satisfying \\[\mathbf{X}=\mathbf{U}\mathbf{\Sigma}\mathbf{V}^T.\\] If we think of \\(\mathbf{V}\\) as the the (standardized) SNP genotypes of an "ideal" person represneting a certain ancestry, then \\[x_{j\cdot}_v_{i\cdot}=u_{ji}\lambda_i\\] represents the amount of idealized person \\(i\\) in actual person \\(i\\), upt to some proportionality constant \\(\lambda_i\\) that depends on the individual.  Then row \\(i\\) of \\(\mathbf{U}\\) are the *ancestries* of person \\(i\\), and column \\(j\\) of \\(\mathbf{U}\\) are the ancestries of each individual on ancestry \\(j\\).  If we rearrange Eq. (1) and use the fact that the columns of \\(\mathbf{U}\\) and \\(\mathbf{V}\\) are *orthonormal*, we can find that \\[\mathbf{X}\mathbf{X}^T\mathbf{U}=\mathbf{\Sigma}^2\mathbf{U}\\] or that the columns of \\(\mathbf{U}\\) are the eigenvectors of the matrix \\[\frac{1}{m}\mathbf{X}\mathbf{X}^T\\]] whose \\(\left(i,j\right)\\) entry is the genetic correlation between individuals \\(i\\) and \\(j\\), sometimes known as the *genomic relationship matrix* or GRM.  PCA works by finding the first several eigenvectors of the GRM and plotting each individual's ancestry along each orthogonal vector in a rectangular grid. Clusters of individuals in this grid represent distinct ancetry groups.
