---
layout: "single"
title: "Principal components analysis of 1000GP populations"
permalink: /assignments/pca/
---

### Importing the data

First we'll see how individuals can be separated by genetic ancestry using principal components analysis.  First we'll need to load the libraries <kbd>SNPRelate</kbd> and <kbd>SeqArray</kbd>:

```
library(SNPRelate)
library(SeqArray)

Download one of the vcf files from the [Files](https://wletsou.github.io/bioinformatics/files) page of a single chromosome containing just the CEU, YRI, and CHB populations.  Once you have the file, store the file name as a variable:

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
