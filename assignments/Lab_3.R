#This lab introduces concepts for analyzing case-control whole-genome association studies.  Similar to the previous lab, we will be implementing more than a single software tool in this lab, to better understand both the tools and the application to whole-genome association studies.  This data set includes 43 pediatric patients that are considered controls for the outcome under study.  We are comparing these control children to a much larger cohort of 400 control children.  We would expect few differences between these 43 ‘cases’ and 400 controls.  

#Like the beginning of any SNP analysis, we need to QC the data and conduct data cleansing, similar to what we have done previously.  Then we will construct a vector to help control for population stratification using MDS and compare the output between case-control associations controlling for population stratification and not controlling for this.

#1.)	First, get the case-control_SNP_files.zip file from the website and unzip it.  We will be working with the ped and map files from this zip. Isolate the nr_no_443.ped pedigree file and the snps_1809.map map file.
ped <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/case-control_SNP_files/nr_no_443.ped')
map <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/case-control_SNP_files/snps_1809.map')

#2.)	Run analysis to filter out SNPs based on the following criteria
#a.	<5% missingness rate per SNP
#b.	MAF > 5%
#c.	<10% missingness rate per subject
#d.	HWE significance at p<.001
#Make sure to create a new pedigree and map file from this analysis.  How many SNPs remain?  How many subjects remain?
./plink --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/case-control_SNP_files/nr_no_443.ped --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/case-control_SNP_files/snps_1809.map --out filtered_out --geno 0.05 --maf 0.05 --mind 0.1 --hwe 0.001 --recode
#PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
#(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to filtered_out.log.
#Options in effect:
#  --geno 0.05
#  --hwe 0.001
#  --maf 0.05
#  --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/case-control_SNP_files/snps_1809.map
#  --mind 0.1
#  --out filtered_out
#  --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/case-control_SNP_files/nr_no_443.ped
#  --recode

#16384 MB RAM detected; reserving 8192 MB for main workspace.
#.ped scan complete (for binary autoconversion).
#Performing single-pass .bed write (1809 variants, 443 people).
#--file: filtered_out-temporary.bed + filtered_out-temporary.bim +
#filtered_out-temporary.fam written.
#1809 variants loaded from .bim file.
#443 people (443 males, 0 females) loaded from .fam.
#443 phenotype values loaded from .fam.
#2 people removed due to missing genotype data (--mind).
#IDs written to filtered_out.irem .
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 441 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate in remaining samples is 0.995947.
#29 variants removed due to missing genotype data (--geno).
#--hwe: 15 variants removed due to Hardy-Weinberg exact test.
#93 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#1672 variants and 441 people pass filters and QC.
#Among remaining phenotypes, 41 are cases and 400 are controls.
#--recode ped to filtered_out.ped + filtered_out.map ... done.

#After these filtering steps, just 1672 variants and 441 people (subjects) remain.

#3.)	Using this newly filtered ped and map file, calculate a simple association test using the Fisher’s exact test method.  Plot the results in R in a -log10(p-value) plot similar to how you have before to visualize the major significant SNPs that show association.  How many SNPs are significant at p<.01 and p<.05?

./plink --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.ped --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.map --assoc fisher #fisher's exact test
#PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
#(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to plink.log.
#Options in effect:
#  --assoc fisher
#  --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.map
#  --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.ped

#16384 MB RAM detected; reserving 8192 MB for main workspace.
#.ped scan complete (for binary autoconversion).
#Performing single-pass .bed write (1672 variants, 441 people).
#--file: plink-temporary.bed + plink-temporary.bim + plink-temporary.fam written.
#1672 variants loaded from .bim file.
#441 people (441 males, 0 females) loaded from .fam.
#441 phenotype values loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 441 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.997326.
#1672 variants and 441 people pass filters and QC.
#Among remaining phenotypes, 41 are cases and 400 are controls.
#Writing C/C --assoc report to plink.assoc.fisher ... done.

# reading in association test data
fisher <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/plink.assoc.fisher', header=TRUE)
# getting counts of SNPs at sig levels 0.05 and 0.01
sig0.05 <- fisher[fisher$P < 0.05,]; nrow(sig0.05)
#[1] 54
sig0.01 <- fisher[fisher$P < 0.01,]; nrow(sig0.01)
#[1] 11

# plotting p-values histogram
fisher$logP <- -1 * log10(fisher$P)
par(mar=c(6,4,4,2))
plot(fisher$BP, fisher$logP, type='h', main="Distribution of p-values from Fisher's exact test across chromosome 6", axes=FALSE, xlab='Position', ylab='-log10(p-value)', cex=0.5, lwd=0.8)
options(scipen=999)
axis(1, las=2, at=seq(min(fisher$BP), max(fisher$BP), length.out=10), cex.axis=0.5)
axis(2, cex.axis=0.9)
abline(h=2, col='red', lty=2)

#There are just 11 variants that are significant at p<0.01. On the other hand, there are 54 variants significant at p<0.05.

#4.)	Now calculate MDS on this pedigree file in Plink and plot the first 2 eigenvectors against each other in R.  Make sure to color the 2 groups different colors and add a legend.

# MDS vectors for ‘case’ and control data
./plink --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.ped --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.map --genome
#PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
#(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to plink.log.
#Options in effect:
#  --genome
#  --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.map
#  --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.ped

#16384 MB RAM detected; reserving 8192 MB for main workspace.
#.ped scan complete (for binary autoconversion).
#Performing single-pass .bed write (1672 variants, 441 people).
#--file: plink-temporary.bed + plink-temporary.bim + plink-temporary.fam written.
#1672 variants loaded from .bim file.
#441 people (441 males, 0 females) loaded from .fam.
#441 phenotype values loaded from .fam.
#Using up to 8 threads (change this with --threads).
#Before main variant filters, 441 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.997326.
#1672 variants and 441 people pass filters and QC.
#Among remaining phenotypes, 41 are cases and 400 are controls.
#IBD calculations complete.  
#Finished writing plink.genome .

./plink --noweb --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.ped --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.map --read-genome plink.genome --cluster --mds-plot 2
#PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
#(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to plink.log.
#Options in effect:
#  --cluster
#  --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.map
#  --mds-plot 2
#  --noweb
#  --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.ped
#  --read-genome plink.genome

#Note: --noweb has no effect since no web check is implemented yet.
#16384 MB RAM detected; reserving 8192 MB for main workspace.
#.ped scan complete (for binary autoconversion).
#Performing single-pass .bed write (1672 variants, 441 people).
#--file: plink-temporary.bed + plink-temporary.bim + plink-temporary.fam written.
#1672 variants loaded from .bim file.
#441 people (441 males, 0 females) loaded from .fam.
#441 phenotype values loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 441 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.997326.
#1672 variants and 441 people pass filters and QC.
#Among remaining phenotypes, 41 are cases and 400 are controls.
#Clustering... done.                        
#Cluster solution written to plink.cluster1 , plink.cluster2 , and
#plink.cluster3 .
#Performing multidimensional scaling analysis (SVD algorithm, 2
#dimensions)... done.
#MDS solution written to plink.mds .

#plot first 2 vectors from mds for rejector/control data
mds <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/plink.mds', header=T)
sam <- c(rep(1,41), rep(2,400))
plot(mds[,4:5], col=sam, pch=16, xlab='p1', ylab='p2', main='MDS plot of "case" and control samples')
legend('topright', c('"case"', 'control'), col=c(1,2), pch=16, cex=0.9)

#5.)	Write these first 2 eigenvectors out in a covariate file using the correct Plink format (Family ID, Individual ID, covariate 1, covariate 2) and run the linear regression association test in Plink.  In R, plot the same -log10 plot as before.  Make sure to only extract that ADD column since this indicates the coefficient in the model that we are testing for.  How many SNPs are significant at p<.01 and p<.05?

write.table(data.frame(mds$FID, mds$IID, mds$C1, mds$C2), file='/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/mycov.txt', sep='\t', col.names=FALSE, row.names=FALSE)

# plink code for linear (logistic) regression association test
./plink --noweb --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.ped --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.map --logistic --covar /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/mycov.txt
#PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
#(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to plink.log.
#Options in effect:
#  --covar /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/mycov.txt
#  --logistic
#  --map /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.map
#  --noweb
#  --ped /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/filtered_out.ped

#Note: --noweb has no effect since no web check is implemented yet.
#16384 MB RAM detected; reserving 8192 MB for main workspace.
#.ped scan complete (for binary autoconversion).
#Performing single-pass .bed write (1672 variants, 441 people).
#--file: plink-temporary.bed + plink-temporary.bim + plink-temporary.fam written.
#1672 variants loaded from .bim file.
#441 people (441 males, 0 females) loaded from .fam.
#441 phenotype values loaded from .fam.
#Using 1 thread (no multithreaded calculations invoked).
#--covar: 2 covariates loaded.
#Before main variant filters, 441 founders and 0 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.997326.
#1672 variants and 441 people pass filters and QC.
#Among remaining phenotypes, 41 are cases and 400 are controls.
#Writing logistic model association results to plink.assoc.logistic ... done.

# reading in association test data
logreg <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/plink.assoc.logistic', header=TRUE)
logreg_add <- logreg[logreg$TEST == 'ADD',]

# getting counts of SNPs at sig levels 0.05 and 0.01
nrow(logreg_add[logreg_add$P < 0.05,])
#[1] 48
nrow(logreg_add[logreg_add$P < 0.01,])
#[1] 9

#There are 48 SNPs that are significant at p<0.05. Only 9 SNPs are significant at p<0.01.

# plotting p-values histogram
logreg_add$logP <- -1 * log10(logreg_add$P)
par(mar=c(6,4,4,2))
plot(logreg_add$BP, logreg_add$logP, type='h', main="Distribution of p-values from logistic regression association test across chromosome 6", axes=FALSE, xlab='Position', ylab='-log10(p-value)', cex=0.5, lwd=0.8)
options(scipen=999)
axis(1, las=2, at=seq(min(logreg_add$BP), max(logreg_add$BP), length.out=10), cex.axis=0.5)
axis(2, cex.axis=0.9)
abline(h=2, col='red', lty=2)
