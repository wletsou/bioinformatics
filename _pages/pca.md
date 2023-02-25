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

### Pricipal components analysis: theory

The populations in our dataset can be separated into clusters based on their genotypes.  The inferred groups help control for confoundng due to ancestry and are also more reliable than self-reported race in association studies.  To see how it works, suppose \\(\mathbf{X}\\) is an \\(n\times m\\) (standardized) genotype matrix with individuals down the rows and (standardized) genotypes across the columns.  Principal components analysis (PCA) says we can find an \\(n\times m\\) matrix \\(\mathbf{V}\\), a diagonal \\(n\times n\\) matrix \\(\mathbf{\Sigma}\\), and an \\(n\times n\\) matrix \\(\mathbf{U}\\) satisfying \\[\mathbf{X}=\mathbf{U}\mathbf{\Sigma}\mathbf{V}^T.\tag{1}\\] If we think of \\(\mathbf{V}\\) as the the (standardized) SNP genotypes of an "ideal" person represneting a certain ancestry, then \\[x_{j\cdot} v_{i\cdot}=u_{ji}\lambda_i\\] represents the amount of idealized person \\(i\\) in actual person \\(i\\), upt to some proportionality constant \\(\lambda_i\\) that depends on the individual.  Then row \\(i\\) of \\(\mathbf{U}\\) are the *ancestries* of person \\(i\\), and column \\(j\\) of \\(\mathbf{U}\\) are the ancestries of each individual on ancestry \\(j\\).  If we rearrange Eq. (1) and use the fact that the columns of \\(\mathbf{U}\\) and \\(\mathbf{V}\\) are *orthonormal*, we can find that \\[\mathbf{X}\mathbf{X}^T\mathbf{U}=\mathbf{\Sigma}^2\mathbf{U}\\] or that the columns of \\(\mathbf{U}\\) are the eigenvectors of the matrix \\[\frac{1}{m}\mathbf{X}\mathbf{X}^T\tag{2}\\] whose \\(\left(i,j\right)\\) entry is the genetic correlation between individuals \\(i\\) and \\(j\\), sometimes known as the *genomic relationship matrix* or GRM.  PCA works by finding the first several eigenvectors of the GRM and plotting each individual's ancestry along each orthogonal vector in a rectangular grid. Clusters of individuals in this grid represent distinct ancetry groups.

### Principal components analysis: practice

To run PCA in R, simply do

```
pca <- snpgds(genofile) # runs PCA
```

to create an objects with 32 eignevectors of the GRM.  Make a data frame of the first several eigenvectors along with subject ids:

```
df.pca <- data.frame(sample = pca$sample.id,EV1 = pca$eigenvect[,1],EV2 = pca$eigenvect[,3],EV1 = pca$eigenvect[,3],stringsAsFactors = FALSE)
```

We'll plot individuals along EV1, EV2, and EV3 in several two-dimensional projections.  But we'll want to see how the clustering done by PCA corresponds to individuals' self-reported race; for that we'll need another column in our data frame.

#### Getting population labels ####

The [file](https://raw.githubusercontent.com/wletsou/bioinformatics/master/docs/CHB%2BYRI%2BCEU.txt) contains each subject id along with its 1KG population group. Let's import the file now

```
indivs <- read.table("path/to/CHB+YRI+CEU.txt",header = TRUE)
```

The second field of this table is <kbd>pop</kbd>, an assignment to each <kbd>id</kbd> of one of three 1KG population groups.  We want to match the right id in <kbd>indivs</kbd> to the right id in <kbd>df.pca</kbd> so that we can color our PCA plots by population.  If the tables are in the same order, matching will be easy, but it not, we have to use the <kbd>match(x,y)</kbd> function, which finds the positions in <kbd>y</kbd> corresponding to the same items in <kbd>x</kbd>.  Thus we can make a new column <kbd>pop</kbd> in <kbd>df.pca</kbd> with the corresponding <kbd>pop</kbd> values from <kbd>indivs</kbd> by

```
df.pca$pop[match(indivs$id,df.pca$sample.id)] <- indivs$pop # find the population group of each individual in df.pca
```
