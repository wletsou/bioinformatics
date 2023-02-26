---
type: "single"
title: "Simulating genotypes and kinship analysis"
toc: true
toc_label: "Contents"
permalink: /assignments/week2/
---

For this assignment we need to load the following packages from Bioconductor. See the [Assignments page](https://wletsou.github.io/bioinformatics/assignments) for more details on installing packages for the first time.

```
library(sim1000G)
library(SNPRelate)
library(GENESIS)
library(GWASTools)
```

We'll also want the package <kbd>data.table</kbd> from CRAN.  <kbd>data.table</kbd> is a more flexible structure than <kbd>data.frame</kbd>:

```
install.packages("data.table") # your first installation
library(data.table) # all subsequent loadings
```

We'll also need the list of individuals again, available [here](https://raw.githubusercontent.com/wletsou/bioinformatics/master/docs/CHB%2BYRI%2BCEU.txt):

```
indivs <- read.table("path/to/file.txt",header = FALSE) # subset of 1KG individuals with population codes
colnames(indivs) <- c("id","pop")
```

### Simulating haplotypes from a phased vcf file ###

For the association study we'll be doing next week, we'll need genotypes of unrelated individuals with and without the disease phenotype of interest.  To protect subjects' privacy, real sequencing data are not readily accessible without data-sharing agreements.  Instead, we'll have to simulate haplotypes from the 1KG vcf files.  The simulation program <kbd>sim1000G</kbd> works by drawing \\(m\\)-SNP haplotypes from an \\(m\\)-dimensional multivariate normal distirbution such that the *allele frquency* of each variant and *linkage disequilibrium* between each pair of variants is preserved.  <kbd>sim1000G</kbd> uses the observed correlation between variants to know what should be the pairwise correlation between any two SNPs in the simulated data.  In order for the correlation to be measured, the input data need to be *phased*, so that we know not only each subject's genotypes, but also which allele appears on which of the subject's two chromosomes.  For example, a line of a phased vcf file might look like

```
1	14930	rs75454623	A	G	100	PASS	AC=284;AF=0.482228;AN=620;NS=2504;DP=42231;EAS_AF=0.4137;AMR_AF=0.5231;AFR_AF=0.4811;EUR_AF=0.5209;SAS_AF=0.4857;AA=a|||;VT=SNP	GT	1|0	1|0	1|0	1|0	0|1	1|1	0|1	0|1	1|1	1|0	0|1
```

with 0 indicating the *reference allele* and 1 the *alternative allele*.  The vertical bars indicate the data have been computationally phased, so that alleles on the left side of the bar all reside on one chromosome and the alleles on the right on the other.  The correlation between alleles \\(i\\) and \\(j\\) is then measured by \\[r_{ij} = \frac{p_{i,j}-p_ip_j}{\sqrt{p_i\left(1-p_i\right)p_j\left(1-p_j\right)}},\tag{1}\\] where \\(p_i\\) is the sample allele frequency.  The measured correlation matrix is then used to generate haplotypes from an \\(m\\)-dimensional multivariate normal distirbution.

#### Running sim1000G ####

To use <kbd>sim1000G</kbd> we first need to read in a vcf file.  vcf files can be very large (hundreds of MB to GB), so we should ony read a small portion in.  <kbd>sim1000G</kbd> lets us choose the number of variants we'd like to simulate and the minimum and maximum allele frequency.  Download one of the 1KG files for one of the chromomsomes [here](https://wletsou.github.io/bioinformatics/files) and read it into R using the <kbd>readVCF</kbd> function:

```
vcf <- readVCF("path/to/file/CHB+YRI+CEU.chr1.vcf.gz", maxNumberOfVariants = 2000 , min_maf = 0.05 , max_maf = 0.95) # read 2000 variants from chr1 with MAF between 5% and 95%
```

This step may take a few minutes to complete depending on the size of the file, so be patient.
