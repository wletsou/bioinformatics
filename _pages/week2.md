---
layout: "single"
title: "Genotype simulation and kinship analysis"
toc: true
toc_label: "Contents"
sidebar:
  nav: "side"
permalink: /assignments/week2/
---

For this assignment we need to load the following packages from Bioconductor.&nbsp;  See the [Assignments page](https://wletsou.github.io/bioinformatics/assignments) for more details on installing Bioconductor packages for the first time.

```
library(sim1000G)
library(SNPRelate)
library(GENESIS)
library(GWASTools)
```

We'll also want the packages <kbd>data.table</kbd> and <kbd>fields</kbd> from CRAN.&nbsp; <kbd>data.table</kbd> is a more flexible structure than <kbd>data.frame</kbd>, and <kbd>fields</kbd> is for plotting images.

```
install.packages("data.table") # your first installation
library(data.table) # all subsequent loadings

install.packages("fields")
library(fields)
```

We'll also need the list of individuals again, available [here](https://raw.githubusercontent.com/wletsou/bioinformatics/master/docs/CHB%2BYRI%2BCEU.txt):

```
indivs <- read.table("path/to/file.txt",header = FALSE) # subset of 1KG individuals with population codes
colnames(indivs) <- c("id","pop")
```

### Simulating haplotypes from a phased vcf file ###

For the association study we'll be doing next week, we'll need genotypes of unrelated individuals with and without the disease phenotype of interest.&nbsp; To protect subjects' privacy, real sequencing data are not readily accessible.&nbsp; Instead, we'll have to simulate haplotypes from the 1KG vcf files.&nbsp; The simulation program [<kbd>sim1000G</kbd>](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2611-1) works by drawing \\(m\\)-SNP haplotypes from an \\(m\\)-dimensional multivariate normal distirbution such that the *allele frquency* of each variant and *linkage disequilibrium* between each pair of variants is preserved.&nbsp; <kbd>sim1000G</kbd> uses the observed correlation between variants to infer what should be the pairwise correlation between any two SNPs in the simulated data.&nbsp; In order for the correlation to be measured, the input data need to be *phased*, so that we know not only each subject's genotypes, but also which allele appears on which of the subject's two chromosomes.  For example, a line of a phased vcf file might look like

```
1	14930	rs75454623	A	G	100	PASS	AC=284;AF=0.482228;AN=620;NS=2504;DP=42231;EAS_AF=0.4137;AMR_AF=0.5231;AFR_AF=0.4811;EUR_AF=0.5209;SAS_AF=0.4857;AA=a|||;VT=SNP	GT	1|0	1|0	1|0	1|0	0|1	1|1	0|1	0|1	1|1	1|0	0|1
```

with 0 indicating the *reference allele* and 1 the *alternative allele*.&nbsp; The vertical bars indicate the data have been computationally phased, so that alleles on the left side of the bar all reside on one chromosome and the alleles on the right on the other.&nbsp; The correlation between alleles \\(i\\) and \\(j\\) is then measured by \\[r_{ij} = \frac{p_{i,j}-p_ip_j}{\sqrt{p_i\left(1-p_i\right)p_j\left(1-p_j\right)}},\tag{1}\\] where \\(p_i\\) is the sample allele frequency and \\(p_{i,j}\\) is the frequency with which two appear together on the same chromosome.&nbsp;  The measured correlation matrix is then used to generate haplotypes from an \\(m\\)-dimensional multivariate normal distirbution.

#### Reading in a vcf file ####

To use <kbd>sim1000G</kbd> we first need to read in a vcf file.  vcf files can be very large (hundreds of MB to GB), so we should ony read a small portion in.&nbsp; <kbd>sim1000G</kbd> lets us choose the number of variants we'd like to simulate and the minimum and maximum allele frequency.&nbsp; Download the  the 1KG chr1 vcf files [here](https://raw.githubusercontent.com/wletsou/bioinformatics/master/docs/CHB+YRI+CEU.chr1.vcf.gz) and read it into R using the <kbd>readVCF</kbd> function:

```
vcf <- readVCF("path/to/file/CHB+YRI+CEU.chr1.vcf.gz", maxNumberOfVariants = 2000 , min_maf = 0.05 , max_maf = 0.95) # read 2000 variants from chr1 with MAF between 5% and 95%
```

This step may take a few minutes to complete depending on the size of the file, so be patient.

#### Making a table of variants ####

Before we generate haplotypes from this vcf file, we'll need to make PLINK *map* and *ped* files, explained below, to be converted into gds objects. 

A map file is simply a manifest of the SNPs, their locations, and their alleles.&nbsp;  You can access the sampled variants in the <kbd>$varid</kbd> field of your <kbd>vcf</kbd> object.  We're going to split up these lines into "chr", "pos", "rsid", "ref",and "alt", and then <kbd>rbind</kbd> the rows together.  Do

```
variants <- data.frame(Reduce(rbind,lapply(vcf$varid,function(X) unlist(strsplit(X,split = " "))))) # variant info, chr, pos, rsid, ref, alt
rownames(variants) <- 1:nrow(variants)
colnames(variants) <- c("chr","pos","rsid","ref","alt") # rename the columns
```

Check that the first few lines agree with the first few entries in <kbd>vcf$varid</kbd>.&nbsp; Your new table of variants should have exactly 2000 records, one for each SNP we sampled from the vcf file.&nbsp; For later use, save the reference and alternate alleles in a variable:

```
ref <- variants[,4]
alt <- variants[,5]
```

Now we'll convert this table to map format.&nbsp; Per the [PLINK file format reference page](https://www.cog-genomics.org/plink/1.9/formats#map), your table should contain

> 1. Chromosome code. PLINK 1.9 also permits contig names here, but most older programs do not.
> 2. Variant identifier
> 3. Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
> 4. Base-pair coordinate

We can achieve this format from our <kbd>variants</kbd> object using

```
df.variants <- data.frame(chr = variants$chr,rsid = variants$rsid,X = rep(0,nrow(variants)),pos = variants$pos) # four columns of the PLINK map file
write.table(df.variants,col.names = FALSE,row.names = FALSE,quote = FALSE,file = "path/to/file/CHB+YRI+CEU.simulation.chr1.map") # save the file with a name you'll use throughout
```

#### Simulating individuals ####

Now we're going to simulate 100 individuals from each of the three populations.&nbsp; We need to do the simulation step separately for each population and then combine the genotypes in a single PLINK ped file.&nbsp; A line from a ped file has the following format:

```
1 1 0 0 1 1 C C A C G G G G A A A A T T C C G G G G C C T T G G G T C C C C T T C C C C G G G G A G G G G G C C G G C C A A G G G A G G C C T T A G C C A G A A
```

where the first six fields are numbers and the seventh through last fields are alleles.&nbsp; The numeric fields should be

1. Family id
2. Individual id
3. Paternal id
4. Maternal id
5. Sex
6. Phenotype

Since we're simulating unrelated individuals, fields 3 and 4 can both be blank (i.e., the individual has no parents in the dataset), and fields 1 and 2 can be identical (i.e., each individual belongs to its own family).&nbsp; Sex is coded as 1 for male and 2 for female, and phenotype is coded as 1 for unaffected and 2 for affected.&nbsp; So the recod above for individual 1 corresponds to an unaffected male.&nbsp;  The 4000 alleles for each of the 2000 SNPs come afterward in pairs.&nbsp; There is no explicit phase information in ped files, but we can assume that the odd-numbered alleles are paternal and the evens maternal.

The challenge will be to get our simulated data into ped format.

To begin our simulation, run

```
N <- 100 # simulate 100 individuals
SIM$reset()
startSimulation(vcf,totalNumberOfIndividuals = N,subset = indivs[indivs$pop == "CEU",]$id) # only take the CEU individuals from the vcf
ids <- generateUnrelatedIndividuals(N) # run the simulation
```

You have now created a <kbd>SIM</kbd> object with fields <kbd>$gt1</kbd> and <kbd>$gt2</kbd> containing 100-by-2000 matrices containing, respectively, the maternal and paternal alleles of the one hundred simulated individuals, coded as 0 for REF and 1 for ALT.&nbsp; To turn these matrices into a ped file, we (1) have to interleave the two matrices and (2) replace the 0's and 1's with A, T, G, and C.

First let's convert our matrices into data tables to make them easier to work with:

```
dt.gt1 <- data.table(SIM$gt1) # first chromosomes alleles
dt.gt2 <- data.table(SIM$gt2) # second chromosome alleles
```

Now, remember those <kbd>ref</kbd> and <kbd>alt</kbd> variables we created earlier?&nbsp;  We'll use them them to write a function that can be applied to each column of the data table:

```
dt.gt1.allele.ceu <- dt.gt1[,mapply(function(X,Y,Z) ifelse(X == 0,Y,Z),.SD,ref,alt),.SDcols = 1:ncol(dt.gt1),by = .I] # convert 0/1 to ATCG
dt.gt2.allele.ceu <- dt.gt2[,mapply(function(X,Y,Z) ifelse(X == 0,Y,Z),.SD,ref,alt),.SDcols = 1:ncol(dt.gt2),by = .I] # convert 0/1 to ATCG
```

This function converts the \\(i\\)<sup>th</sup> column from a number (0 or 1) to the \\(i\\)<sup>th</sup> entry of <kbd>ref</kbd> or <kbd>alt</kbd> depending on its value. &nbsp; The arguments <kbd>.SDcols</kbd> and <kbd>by = .I</kbd> tell <kbd>data.table</kbd> to apply the function rowwise at each column individually.

You can assign the correct variant names to the columns by running

```
colnames(dt.gt1.allele.ceu) <- variants$rsid
colnames(dt.gt2.allele.ceu) <- variants$rsid
```

Before we interleave the columns of these tables to make a ped file, **repeat the above procedure for the YRI and CHB populations** and generate appropriately named data tables.&&nbsp; Be sure to give your tables unique names, or else you will have to run the simulation all over again.

#### Generating a gds file from your ped and map files ####

Once you have the six data tables <kbd>dt.gt1.allele.ceu</kbd> to <kbd>dt.gt2.allele.chb</kbd>, you're ready to make a ped file and convert it into a gds object suitable for <kbd>SNPRelate</kbd> and <kbd>GENESIS</kbd>.&nbsp; First let's create the six numeric columns of the file for our 300 individuals.&nbsp;  All we need is to number individuals/families from 1 to 300 (<kbd>=nrow(dt.gt1.allele.chb) + nrow(dt.gt1.allele.yri) + nrow(dt.gt1.allele.ceu)</kbd>) and ensure they have no parents in the data.&nbsp;  We'll randomly select the subject's sex and for now assume that each individual is unaffected.

```
fam <- data.frame(fid = 1:(nrow(dt.gt1.allele.chb) + nrow(dt.gt1.allele.yri) + nrow(dt.gt1.allele.ceu)),id = 1:(nrow(dt.gt1.allele.chb) + nrow(dt.gt1.allele.yri) + nrow(dt.gt1.allele.ceu)), mother = 0,father = 0,sex = sample(c(0,1),nrow(dt.gt1.allele.chb) + nrow(dt.gt1.allele.yri) + nrow(dt.gt1.allele.ceu),replace = TRUE),phenotype = 1) # numeric columns of ped file
```

This should look something like

```
head(fam)
  fid id mother father sex phenotype
1   1  1      0      0   1         1
2   2  2      0      0   0         1
3   3  3      0      0   1         1
4   4  4      0      0   1         1
5   5  5      0      0   0         1
6   6  6      0      0   1         1
```

Now, to interleave the columns of our <kbd>dt.gt1.allele</kbd> and <kbd>dt.gt2.allele</kbd> tables for each population, we need to first glue the tables together side-by-side and take, in pairs, columns 1 and 2001, 2 and 2002, 3 and 2003, etc.  We can do this by repeating each number in 1 to 2000 twice (i.e., 1 1 2 2 3 3 ...) and then adding 2000 to every other value:

```
cbind(dt.gt1.allele.chb,dt.gt2.allele.chb)[,rep(1:ncol(dt.gt1.allele.chb),each = 2) + (0:1) * ncol(dt.gt1.allele.chb)]
```

Then we'll group these tables (one for each population) in an R <kbd>list</kbd> and <kbd>Reduce</kbd> the <kbd>rbind</kbd> function over the three entries of the list to make a table with 300 rows, which we'll finally join to our <kbd>fam</kbd> object above.&nbsp;  Save this table as file using a similar naming scheme you used for the map file.

```
write.table(cbind(fam,Reduce(rbind,list(data.frame(cbind(dt.gt1.allele.chb,dt.gt2.allele.chb))[,rep(1:ncol(dt.gt1.allele.chb),each = 2) + (0:1) * ncol(dt.gt1.allele.chb)],data.frame(cbind(dt.gt1.allele.yri,dt.gt2.allele.yri))[,rep(1:ncol(dt.gt1.allele.yri),each = 2) + (0:1) * ncol(dt.gt1.allele.yri)],data.frame(cbind(dt.gt1.allele.ceu,dt.gt2.allele.ceu))[,rep(1:ncol(dt.gt1.allele.ceu),each = 2) + (0:1) * ncol(dt.gt1.allele.ceu)]))),col.names = FALSE,row.names = FALSE,quote = FALSE,file = "path/to/file/CHB+YRI+CEU.simulation.chr1.ped") # save .ped file
```

Finally, let's convert the map and ped files to gds format an create a gds object.  Use the <kbd>snpgdsSummary</kbd> function to verify that you have the right number of samples and variants.

```
snpgdsPED2GDS("path/to/file/CHB+YRI+CEU.simulation.chr1.ped","path/to/file/CHB+YRI+CEU.simulation.chr1.map","path/to/file/CHB+YRI+CEU.simulation.chr1.gds")
snpgdsSummary("path/to/file/CHB+YRI+CEU.simulation.chr1.gds")
genofile <- SNPRelate::snpgdsOpen("path/to/file/CHB+YRI+CEU.simulation.chr1.gds")
```

If you have any trouble with this step, you can close all open gds objects with

```
showfile.gds(closeall=TRUE)
```

### To turn in (1): ###

On your slides, print out:

1. The first ten rows of your <kbd>variants</kbd> table
2. A ten-by-ten sample of your <kbd>dt.gt1.allele.ceu</kbd> file
3. A ten-by-ten sample of the corresponding <kbd>dt.gt2.allele.ceu</kbd> file
4. A ten-by-sixteen sample of the object which became your ped file
5. The <kbd>genofile</kbd> summary

Explain how these files are related.

### Kinship analysis ###

The kinship analysis program [KING](https://academic.oup.com/bioinformatics/article/26/22/2867/228512?login=false) is used to estimate whether any individuals are related to each other.&nbsp;  Cryptic relatedness can bias the results of association studies, so we should correct for it.  KING works by comparing the alleles of each SNP for a pair of individuals.&nbsp;  Before we get into computing the genetic relationship coefficient, we first need to reduce our SNP set to a list of approximately independent alleles.&nbsp;  This first step is known as *LD-pruning*.

#### LD-pruning ####

The <kbd>SNPRelate</kbd> package computes the correlation coefficient (Eq. (1)) for each pair of SNPs and removes SNPs which are correlated by more than a certain threshold.&nbsp; To get a subset of quasi-independent SNPs, run

```
snpset <- snpgdsLDpruning(genofile,method = "corr", slide.max.bp = 10e6,ld.threshold = sqrt(0.1), verbose = FALSE) # remove correlated SNPs
pruned <- unlist(snpset, use.names=FALSE) # get the pruned list
```

This function computes the correlation coefficient \\(r_{i,j}\\) between each pair of SNPs and removes those with \\(r < \sqrt{0.1}\\) or LD \\(r^2 < 0.1\\) with an index SNPs.&nbsp;  The computation is restricted to SNPs in ten-million-bp windows.&nbsp;  To see that your SNPs are indepdents, make and plot a matrix of the \\(r^2\\) values for each pair of SNPs using

```
LD.mat <- snpgdsLDMat(genofile, snp.id = pruned, slide = 0,method = "corr") # creates LD matrix of r values
fields::image.plot(1:nrow(LD.mat$LD),1:ncol(LD.mat$LD),LD.mat$LD ^ 2,main = "LD r2",xlab = "SNPs",ylab = "SNPs",asp = 1,frame.plot = FALSE) # plots r^2 values
```

Be sure to square each entry of <kbd>LD.mat</kbd> to get \\(r^2\\) instead of \\(r\\).  How many SNPs are left in your pruned set?

#### Kinship: theory ####

Kinship can be defined as the expected fraction of alleles that two individuals got from the same ancestor(s).&nbsp; We say that two individuals share an allele of a SNP *identical-by-descent* or *IBD* if they inherited the same copy of the allele from a common ancestor.&nbsp; IBD-sharing is different from simply carrying the same allele of a gene (known as *identical-by-state* or *IBS*-sharing), which unrelated individuals may do if the allele is common enough in the population.&nbsp; The *degree* \\(R\\) of relationship may be defined as the effective number of meioses separating the relatives through the equation \\[\frac{1}{2^R}=\frac{1}{2^{R_1}}+\frac{1}{2^{R_2}},\tag{1}\\]where \\(R_i\\) is the number of meioses separating the relatives through the first relative's \\(i\\)<sup>th</sup> parent.&nbsp; For example, sibs are connected by two meioses through two parents, while a parent and child are connected by one meiosis through one parent: both relationships are degree-1.  

The probability that two relatives share an allele IBD is \\(\frac{1}{2^R}\\), as there is a \\(\frac{1}{2}\\) chance that an allele is passed on in any meiosis, and \\(R\\) is the effective number of meioses or steps between the relatives.&nbsp; If \\(2\times\frac{1}{2^R}\\) is the expected number of alleles shared IBD at any given locus, then the fraction of the genome shared by any two relatives is\\[r=\frac{2\times\\frac{1}{2^R}}{2}=\frac{1}{2^R}.\tag{2}\\]

However, genomic sharing can be realized in different ways depending on the probabilities \\(\pi_0\\), \\(\pi_1\\), and \\(\pi_2\\) that individuals share zero, one, or two copies IBD at a locus.&nbsp; The probability that the relatives inherit both both copies IBD\\[\pi_2=P\left(\text{share 2 IBD}\right)=2^2\frac{1}{2^{R_1}2^{R_2}}\tag{3a}\\] is simply the product of the probabilities of sharing through both parents.&nbsp; The probability of sharing exactly one allele IBD\\[\pi_1=P\left(\text{share 1 IBD}\right)=2\left(\frac{1}{2^{R_1}}+\frac{1}{2^{R_2}}\right)-2^3\frac{1}{2^{R_1}2^{R_2}}\tag{3b}\\] is the got by finding the probability \\(2\left(\frac{1}{2^{R_1}}+\frac{1}{2^{R_2}}\right)-2^2\frac{1}{2^{R_1}2^{R_2}}\\) of sharing at least one allele IBD less the probability \\(\pi_2\\) of sharing two.&nbsp; Finally, the probability of sharing at zero alleles IBD\\[\pi_0=P\left(\text{share 0 IBD}\right)=1-2\frac{1}{2^{R_1}}-2\frac{1}{2^{R_2}}+2^2\frac{1}{2^{R_1}2^{R_2}}\tag{3c}\\]is got by subtracting the probability \\(\pi_1+\pi_2\\) of sharing at least one allele IBD from 1.&nbsp; The coefficients account for the fact that there are \\(2\\) alleles at each locus and \\(2^2\\) that can be shared.

From (3a)&ndash;(3c), the fraction of the genome shared IBD is\\[r=\frac{2\pi_2+1\pi_1}{2}=\frac{1}{2^{R_1}}+\frac{1}{2^{R_2}}=\frac{1}{2^R}.\tag{4}\\]But the same value of \\(r\\) can obtain from different values of \\(\pi_1\\) and \\(\pi_2\\).&nbsp; For example, full sibs have a 25% probability of sharing two alleles and a 50% chance of sharing one allele at a locus, for a total fraction \\(r=0.5\\) shared.&nbsp; But a parent and child have 0% chance of sharing two alleles and a 100% chance of sharing one, also giving \\(r=0.5\\).&nbsp; Thus, your parent does not equal your sibling, despite your sharing equal amounts of your genome with each of them.&nbsp; Put another way, parents cannot pass on their genotypes to their offspring.

#### KING ####

KING computes both the probability \\(\pi_0\\) that two relatives share 0 alleles IBD as well as the coefficient of relatedness \\(\phi=\frac{r}{2}\\), defined as the probability that two alleles taken one from each relative are IBD at a locus (the maximum probability is \\(\frac{1}{2}\\) because there is a 50% chance that the alleles chosen come from different parents).&nbsp; The idea is to compare the counts \\(X\\) and \\(Y\\) of the alternative alleles that two individuals each have at a genetic locus.&nbsp; If the two individuals are from the same population, the expected values and variances of the allele counts are \\(\mathbb{E}\left(X\right)=\mathbb{E}\left(Y\right)=2p\\) and \\(\sigma^2_X=\mathbb{E}\left(X^2\right)-\mathbb{E}\left(X\right)^2=\\(\mathbb{E}\left(Y^2\right)-\mathbb{E}\left(Y\right)^2=2p\left(1-p\right)\\).&nbsp; Thus the expected value of the difference \\(\mathbb{E}\left(X^2-Y^2\right)=\mathbb{E}\left(X^2\right)+\mathbb{E}\left(Y^2\right)-2\mathbb{E}\left(XY\right)\\) is\\[\frac{\mathbb{E}\left(X^2-Y^2\right)}{\sigma_X^2+\sigma_Y^2}=1-\frac{\sigma_{XY}}{\sigma_X\sigma_Y}=1-r,\tag{5}\\]where \\(\sigma_{XY}=\mathbb{E}\left(XY\right)-\mathbb{E}\left(X\right)\mathbb{E}\left(Y\right)\\) is the covariance of the genotype counts and \\(r=2\phi\\) is the genetic correlation between two individuals, also interprettted as the amount of the genome shared IBD.

KING estimates \\(\phi\\) from Eq. (5) by counting the number \\(N\\) of loci at which two individuals are heterozygous \\(Aa,Aa\\) or opposite homozygous \\(AA,aa\\), as well as the total number of alleles at which each individual is heterozygous \\(Aa\\), using the equation:\\[\hat{\phi_{ij}}=\frac{N_{Aa,Aa}-2N_{AA,aa}}{N_{Aa}^{\left(i\right)}+N_{Aa}^{\left(j\right)}}.\tag{6}\\]From (6) it can be seen that shared heterozygous sites increase the estimated relatedness, and that unshared homozygous sites decrease relatedness.&nbsp;  Eq. (6) is called a "robust" estimator because it measures relatedness in a purely pairwise fashion; it does not rely on population estimates of allele frequencies.&nbsp; However, if the individuals are not of the same genetic background, the allele frequency \\(p\\) is not well-defined and Eq. (5) does not hold, leading to negative estimates of \\(\phi\\); this feature is not necessarily a problem, as it helps us to distinguish different ancestries within a single population.

To run KING we need only a gds object and a set of SNPs.&nbsp;  We will use the LD-pruned set <kbd>pruned</kbd> we computed above and the <kbd>genofile</kbd> containing simulated haplotypes from CHB, YRI, and CEU individuals.&nbsp;  Running

```
ibd <- snpgdsIBDKING(genofile,snp.id = pruned) # run KING
colnames(ibd$kinship) <- ibd$sample.id # label rows and columns of kinship matrix
rownames(ibd$kinship) <- ibd$sample.id
```

will gives us an object <kbd>ibd</kbd> that contains two matrics <kbd>$IBS0</kbd> and <kbd>$kinship</kbd> that contain the probabilities of sharing zero alleles identical by state (IBS, not IBD) and the estimated kinship coefficients.&nbsp; Make a plot of <kbd>IBS0</kbd> vs. <kbd>kinship</kbd> to see if subjects who share few alleles IBS are unrelated:

```
plot(ibd$IBS0,ibd$kinship,ylab = "Kinship coeffecient",xlab = "IBS0 proportion",main = "KING relatedness estimation")
```

You should see that most values are negative, indicating that individuals come from different populations.&nbsp;  Recalling that individuals 1 to 100 are CHB, 101 to 200 are YRI, and 201 to 300 are CEU, make plots for the populations separately to see if you get fewer negative values.

We are done using the <kbd>genofile</kbd> object, so close it now:

```
closefn.gds(genofile)
```

#### PC-AiR and PC-Relate ####

Your kinship plots probably look very messy.&nbsp; This is okay since we are using simulated data at only a few markers.&nbsp; To improve our estimates, we will use two related packages: <kbd>pcair</kbd> and <kbd>pcrelate</kbd>.&nbsp; These packages improve our estimates of kinship by first running PCA on a subset of ancestry-representative individuals (PC-Air) and then using the principal components to separate sharing due to common ancestry from sharing due to kinship (PC-Relate).

To use these programs we also need the package <kbd>GWASTools</kbd> which creates an object that can be used for association testing in [Week 3](https://wletsou.github.io/bioinformatics/assignments/week3).  First we read in our gds file and create a <kbd>GenotypeData</kbd> object:

```
geno <- GdsGenotypeReader("path/to/file/CHB+YRI+CEU.simulation.chr1.gds")
genoData <- GenotypeData(geno)
```

Now we'll run PC-AiR to identify and compute PCs of a subset of unrelated individuals who represent each ancestry group.&nbsp; We need the initial results of our KING analysis above as well as the list of LD-pruned SNPs.&nbsp; Then using the <kbd>genoData</kbd> object, run

```
mypcair <- pcair(genoData,kinobj = ibd$kinship,divobj = ibd$kinship,snp.include = pruned) # genotype principal components based on a subset of unrelated individuals
```

You can quickly generate plots of the principal components using:

```
plot(mypcair,vx = 1, vy = 2) # plot of PC2 vs PC1
```

The black dots represent the "ancestry-representative subset" for whom PCs have been computed and the blue dots individuals whose PCs have been estimated based on their similarity to the reprentatives.&nbsp; Do the PCs look similar to your output from last week?

The important part of <kbd>mypcair</kbd> is the field <kbd>$vectors</kbd> containing new PCs for this unrelated subset of individuals.&nbsp; Next we will create a <kbd>GenotypeBlockIterator</kbd> needed for PC-Relate.&nbsp; This object is important for running genome-wide analyses in parallel for millions of SNPs; we only have about 1000 SNPs, so parallelization is unnecessary.

```
genoData.iterator <- GenotypeBlockIterator(genoData,snpInclude = pruned)
```

Now we'll use PC-Relate to update the kinship coefficients for the entire sample based on the fact that the GRM is biased when genotypes are standardized to allele frequencies measured among related individuals.&nbsp;  To get a corrected GRM, we supply PC-Relate with the PC results from PC-AiR and the list of unrelated individuals for training the model:

```
mypcrel <- pcrelate(genoData.iterator,pcs = mypcair$vectors[,1:3],training.set = mypcair$unrels) # kinship based on unrelated individuals
myGRM <- pcrelateToMatrix(mypcrel,scaleKin = 2) # genotype relatedness matrix 
```

The <kbd>mypcrel</kbd> object contagins fields <kbd>$kinBtwn\$k0</kbd> and <kbd>$kinBtwn\$kin</kbd> for the recomputed IBS0 proportions and kinship coefficients.  Make a new plot using

```
plot(mypcrel$kinBtwn$k0,mypcrel$kinBtwn$kin,xlab = "IBD0 proportion",ylab = "Kinship coefficient",main = "PC-relate relatedness estimation") # kinship vs. IBD0
```

Also plot subsets of the data frame <kbd>$kinBtwn</kbd> for values of <kbd>ID1</kbd> and <kbd>ID2</kbd> restricted to each of the three populations as we did above, for example

```
plot(mypcrel$kinBtwn[mypcrel$kinBtwn$ID1 %in% 1:100 & mypcrel$kinBtwn$ID2 %in% 1:100,]$k0,mypcrel$kinBtwn[mypcrel$kinBtwn$ID1 %in% 1:100 & mypcrel$kinBtwn$ID2 %in% 1:100,]$kin,xlab = "IBD0 proportion",ylab = "Kinship coefficient",main = "PC-relate relatedness estimation (CHB)")
```

Does the plot look cleaner?&nbsp;  You can get the unique individuals who have kinship closer than \\(R = -log_2{\left(2\phi\right)}\\) by doing

```
unique(c(subset(mypcrel$kinBtwn[mypcrel$kinBtwn$ID1 %in% 1:100 & mypcrel$kinBtwn$ID2 %in% 1:100,],kin > phi)$ID1,subset(mypcrel$kinBtwn[mypcrel$kinBtwn$ID1 %in% 1:100 & mypcrel$kinBtwn$ID2 %in% 1:100,],kin > phi)$ID2)) # unique individuals having kinship coefficient greater than phi
```

Your new GRM <kbd>myGRM</kbd> contains estimates of the genetic relatedness between study subjects which we will use next week for association testing.

Since we are done for now, you may close your genoData object:
  
```
close(genoData)
```

### To turn in (2): ###

Generate PC-Relate plots for each of the CHB, YRI, and CEU populations.&nbsp; For each group, answer:
1. What is the approximate genetic relationship of the most closely related pair of individuals in each population?
2. How many individuals have a third-degree relative or closer in the population?
