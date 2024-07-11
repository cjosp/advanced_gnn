#This lab introduces concepts for analyzing family-based SNP data.  We will be implementing multiple software tools in this lab, to better understand both the tools and the application to family-based SNP data.  First, we will begin with some important data cleansing concepts such as SNP filtering based on HWE criterion, missingness, and MAFs.  Then we will run an IBD analysis for verification of parental genotypes.  This is important because sometimes the parental information is not from a biological parent, so this subject needs to be removed from the analysis since this person did not transmit the genotype information to the offspring under study.  We will then calculate a TDT and plot the results.

#1.)	First, get the family_SNP_files.zip file from the website and unzip it.  We will be working with the ped and map files from this zip.
ped <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/family_SNP_files/ped_final_re_MHC_mod.ped', header=FALSE)
map <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/family_SNP_files/ped_final_re_MHC_mod.map', header=FALSE)

#2.)	Run analysis to filter out SNPs based on the following criteria
#a.	<5% missingness rate per SNP
#b.	MAF > 10%
#c.	<20% missingness rate per subject
#d.	HWE significance at p<.001
#Make sure to create a new pedigree and map file from this analysis.  How many SNPs remain?  How many subjects remain?
./plink --noweb --file /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/family_SNP_files/ped_final_re_MHC_mod --out filtered_ped_final_re_MHC_mod --geno 0.05 --maf 0.1 --mind 0.2 --hwe 0.001 --recode
#PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
#(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to filtered_ped_final_re_MHC_mod.log.
#Options in effect:
#  --file /Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/family_SNP_files/ped_final_re_MHC_mod
#  --geno 0.05
#  --hwe 0.001
#  --maf 0.1
#  --mind 0.2
#  --noweb
#  --out filtered_ped_final_re_MHC_mod
#  --recode

#Note: --noweb has no effect since no web check is implemented yet.
#16384 MB RAM detected; reserving 8192 MB for main workspace.
#.ped scan complete (for binary autoconversion).
#Performing single-pass .bed write (1813 variants, 201 people).
#--file: filtered_ped_final_re_MHC_mod-temporary.bed +
#filtered_ped_final_re_MHC_mod-temporary.bim +
#filtered_ped_final_re_MHC_mod-temporary.fam written.
#1813 variants loaded from .bim file.
#201 people (134 males, 67 females) loaded from .fam.
#201 phenotype values loaded from .fam.
#7 people removed due to missing genotype data (--mind).
#IDs written to filtered_ped_final_re_MHC_mod.irem .
#Using 1 thread (no multithreaded calculations invoked).
#Before main variant filters, 132 founders and 62 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate in remaining samples is 0.99132.
#46 variants removed due to missing genotype data (--geno).
#--hwe: 1 variant removed due to Hardy-Weinberg exact test.
#344 variants removed due to minor allele threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#1422 variants and 194 people pass filters and QC.
#Among remaining phenotypes, 30 are cases and 164 are controls.
#--recode ped to filtered_ped_final_re_MHC_mod.ped +
#filtered_ped_final_re_MHC_mod.map ... done.

#1422 SNPs (variants) and 194 subjects (people) remain after filtering for this criteria using PLINK.

#3.)	Using this newly filtered ped and map file, run IBD to identify any inconsistencies in parental genotypes using a p<.01 threshold.  How many paternal or maternal subjects have a p<.01 when compared to their offspring?  How many paternal or maternal subjects have p<.01 compared to each other?  We would expect paternal and maternal subjects to not be ‘genetically related’ and thus have a significant p-value, but for many paternal and maternal subjects, this is not the case, why do you think this is not occurring using this data set?
./plink --file ./filtered_ped_final_re_MHC_mod --genome rel-check
#PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
#(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to plink.log.
#Options in effect:
#  --file ./filtered_ped_final_re_MHC_mod
#  --genome rel-check

#16384 MB RAM detected; reserving 8192 MB for main workspace.
#.ped scan complete (for binary autoconversion).
#Performing single-pass .bed write (1422 variants, 194 people).
#--file: plink-temporary.bed + plink-temporary.bim + plink-temporary.fam written.
#1422 variants loaded from .bim file.
#194 people (129 males, 65 females) loaded from .fam.
#194 phenotype values loaded from .fam.
#Using up to 8 threads (change this with --threads).
#Before main variant filters, 132 founders and 62 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.99267.
#1422 variants and 194 people pass filters and QC.
#Among remaining phenotypes, 30 are cases and 164 are controls.
#IBD calculations complete.  
#Finished writing plink.genome .

geno <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/plink.genome', header=TRUE)
po <- geno[geno$RT=='PO',]; ot <- geno[geno$RT=='OT',]
po_sig <- po[po$PPC<0.01,]; ot_sig <- ot[ot$PPC<0.01,]
nrow(po_sig); nrow(ot_sig)
#[1] 0
#[1] 4
#There are 0 paternal or maternal subjects that have a p<0.01 when compared to their offspring, but there are 4 paternal or maternal subjects that have p<0.01 when compared to each other. 

#I think a possible reason for this finding is that, in some cultures/communities, relatives with some considerable degree of genetic relation (e.g., second or third cousins) get married. This leads to increased consanguinity between the spouses than is expected with a marriage “at random”. The samples of this data set might reflect a population(s) where marriage between related individuals is not uncommon.

#4.)	For those parents that differ significantly from their child (from above, if any), remove them (using the –-remove flag and instructions on either the plink website or the lecture notes) and run the TDT with output of 95% confidence intervals.
# no parents appeared to differ significantly from their child (i.e., nrow(po_sig)) was 0)
# so no PO rows to remove

# running TDT with output of 95% CI
./plink --file filtered_ped_final_re_MHC_mod --tdt --model perm --ci 0.95
#PLINK v1.90b7.2 64-bit (11 Dec 2023)           www.cog-genomics.org/plink/1.9/
#(C) 2005-2023 Shaun Purcell, Christopher Chang   GNU General Public License v3
#Logging to plink.log.
#Options in effect:
#  --ci 0.95
#  --file filtered_ped_final_re_MHC_mod
#  --model perm
#  --tdt

#16384 MB RAM detected; reserving 8192 MB for main workspace.
#.ped scan complete (for binary autoconversion).
#Performing single-pass .bed write (1422 variants, 194 people).
#--file: plink-temporary.bed + plink-temporary.bim + plink-temporary.fam written.
#1422 variants loaded from .bim file.
#194 people (129 males, 65 females) loaded from .fam.
#194 phenotype values loaded from .fam.
#Using up to 8 threads (change this with --threads).
#Before main variant filters, 132 founders and 62 nonfounders present.
#Calculating allele frequencies... done.
#Total genotyping rate is 0.99267.
#1422 variants and 194 people pass filters and QC.
#Among remaining phenotypes, 30 are cases and 164 are controls.
#Writing --model report to plink.model ... done.                    
#3718 (adaptive) permutations complete.
#Permutation test report written to plink.model.best.perm .
#--tdt: Report written to plink.tdt .

#5.)	Now we want to visualize the regions along chromosome 6 where there are significant p-values representing significant transmission of alleles from parent to offspring.  Read in the TDT results and create the plot below applying your knowledge of the plot function in R.  Be sure to add the axes after you plot the lines and plot structure.  You may also have to use a for loop to add the lines to the plot.  How many SNPs have p<.05?  Load the qvalue R package and compute false discovery rate adjusted p-values using fdr.level=0.05.  How many SNPs have FDR adjusted p-values of p<.05?  What does this latter result indicate to you?
# reading in TDT data
tdt <- read.table('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/plink/plink.tdt', header=TRUE)
tdt$logP <- -1 * log10(tdt$P)

# plotting p-values histogram
par(mar=c(6,4,4,2))
plot(tdt$BP, tdt$logP, type='h', main='Distribution of p-values from TDT across chromosome 6', axes=FALSE, xlab='Position', ylab='-log10(p-value)', cex=0.5, lwd=0.8)
options(scipen=999)
axis(1, las=2, at=seq(min(tdt$BP), max(tdt$BP), length.out=10), cex.axis=0.5)
axis(2, cex.axis=0.9)
abline(h=2, col='red', lty=2)

sum(tdt$P<0.05)
#[1] 120

#There appear to be 120 SNPs with p-value<0.05. (121 SNPs have p-values with -log10(p-value)<0.05)

# computing q-values at FDR level = 0.05
library(qvalue)
qvals <- qvalue(tdt$P, fdr.level=0.05)
sum(qvals$qvalues<0.05)
#[1] 0

#There are no (0) SNPs with FDR-adjusted p-value<0.05 at the FDR level of 0.05. 

#This suggests that at the FDR level of 0.05, the 120 SNPs previously found to be significant according to their raw p-values with the alpha level of 0.05, are not significant at the same alpha level when considering their FDR-adjusted p-values. 

#In other words, those 120 SNPs might not robust enough to withstand the adjustment for multiple comparisons at this FDR level. This doesn't definitively prove they are all false positives, as some potentially true positives with weaker evidence may have been mistakenly excluded by this FDR correction. 
