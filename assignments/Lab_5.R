#In this lab, we will be analyzing expression quantitative trait loci (eQTL) data. The data matrix is taken from a study where the investigators are profiling both whole genome transcripts and genome-wide SNPs on normal healthy controls and subjects diagnosed with bipolar disorder. The specimen taken for measurement is post mortem brain tissue from the subjects.
#We will be focusing on the issues of data manipulation and analysis that go into computing eQTLs between gene expression microarray and SNP array data. There are 83 subjects in each data matrix (i.e. expression matrix and pedigree/map files). The differential expression analysis has already been conducted for genes that discriminate between bipolar and control subjects. We have identified the gene DICER1 as a top candidate for being differentially expressed between the two conditions. Now we want to identify any SNPs, both trans (limited to chromosome 14) and cis, that may be correlated with this gene.

#1.) Get the eqtl_files.zip file from the course website and extract all files. Read in all files but the demographic file into R, since we will not be working with this one for this lab. For the pedigree file, take note that it is space delimited and not tab delimited. There are no headers on the ped or map files.
setwd('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN')
ped <- read.table("eqtl_files/bp_ct_83_subjects.ped",header=F,sep=" ")
map <- read.table("eqtl_files/bp_ct_83_subjects.map",header=F,sep="\t")
exp <- read.table("eqtl_files/bp_ct_u133A.txt",sep="\t",header=T)

#2.) The Affymetrix probe ID for DICER1 is below, as are the SNPs within a 100 kb distance of the DICER1 transcript. Read them into R and find the positions of these SNPs in the map file. You will need to find the actual SNPs in the pedigree file. To do this, use the formula:

#Ped position2 = (Map position*2)+6
#Ped position1 = Ped position2 – 1
#DICER1 gene and SNPs cis to DICER1
#probe <- "216280_s_at"
#snps <- c("SNP_A-2240938","SNP_A-2249234","SNP_A-2120212","SNP_A-1864785","SNP_A-1944204", "SNP_A-2101022","SNP_A-2192181","SNP_A-4296300","SNP_A-1961028","SNP_A-2215249","SNP_A-2172180","SNP_A-1945117","SNP_A-1941491")

# DICER1 gene and SNPs cis to APOE
probe <- "216280_s_at"
snps <- c("SNP_A-2240938","SNP_A-2249234","SNP_A-2120212","SNP_A-1864785","SNP_A-1944204",
"SNP_A-2101022","SNP_A-2192181","SNP_A-4296300","SNP_A-1961028","SNP_A-2215249","SNP_A-2172180","SNP_A-1945117","SNP_A-1941491")
#find positions of SNPs cis to DICER1 in the map file
cis.snps <- snps
map.cis.pos <- match(cis.snps, map$V2)
#calculate the positions of those SNPs in the ped file
ped.cis.pos2 <- (map.cis.pos*2) + 6
ped.cis.pos1 <- ped.cis.pos2 - 1
map.cis.pos; ped.cis.pos2; ped.cis.pos1
# [1] 2337 2338 2339 2340 2341 2342 2343 2344 2345 2346 2347 2348 2362
# [1] 4680 4682 4684 4686 4688 4690 4692 4694 4696 4698 4700 4702 4730
# [1] 4679 4681 4683 4685 4687 4689 4691 4693 4695 4697 4699 4701 4729

#3.) Now, subset the pedigree file by these positions of the 13 cis-acting SNPs and apply the following function below to convert the SNPs from a 2 allele genotype code to a single 1 number genotype code (i.e. 11=0, 12/21=1, 22=2, and 00=NA). Hint: use the apply statement on the pedigree file that you just subset.
#recode <- function(x) {
  #x <- as.numeric(x)
  #x1 <- seq(1,length(x),by=2)
  #x2 <- seq(2,length(x),by=2)
  #geno <- paste(x[x1],x[x2],sep="")
  #geno[geno=="00"] <- NA
  #geno[geno=="11"] <- 0
  #geno[geno=="12" | geno=="21"] <- 1
  #geno[geno=="22"] <- 2
  #geno
#}
#create ped subset for the 13 cis-acting SNPs
ped.cis.subset <- ped[,sort(c(ped.cis.pos1, ped.cis.pos2))]
#reorder cis SNPs to make sure they agree with sorting by positions above
cis.snps <- snps[order(pos2)]
recode <- function(x) {
  x <- as.numeric(x)
  x1 <- seq(1,length(x),by=2)
  x2 <- seq(2,length(x),by=2)
  geno <- paste(x[x1],x[x2],sep="")
  geno[geno=="00"] <- NA
  geno[geno=="11"] <- 0
  geno[geno=="12" | geno=="21"] <- 1
  geno[geno=="22"] <- 2
  geno
}
# convert the SNPs from a 2 allele genotype code to a single 1 number genotype code
ped.cis.converted <- apply(ped.cis.subset, 1, recode)

#4.) Calculate a linear regression model between the gene expression values for DICER1 (using the probe variable and the gene expression matrix) and the single numeric coded genotypes for the 13 cis-acting SNPs. Hint: use the lin.mod() function below with an apply statement.
#lin.mod <- function(x,y) {
  #x <- as.numeric(x)
  #dat <- data.frame(y,x)
  #out <- lm(y~x,data=dat)
  #outx <- summary(out)
  #return(outx$coefficients["x",4])
#}
lin.mod <- function(x,y) {
  x <- as.numeric(x)
  dat <- data.frame(y,x)
  out <- lm(y~x,data=dat)
  outx <- summary(out)
  return(outx$coefficients["x",4])
}
cis.pvals <- apply(ped.cis.converted, 1, lin.mod, y=as.numeric(exp[probe,]))
cis.pvals
#The eQTL p-values for the 13 cis-acting SNPs are as follows:
#> cis.pvals
#[1] 0.29890757 0.66356446 0.37822627 0.72387428 0.42564052 0.41066079 0.37232864 0.24248884
#[9] 0.71146674 0.58259739 0.37334937 0.53435622 0.01159989

#5.) Which of these 13 SNPs has the lowest p-value? What is the p-value?
#SNP_A-1941491 has the lowest p-value; the p-value is 0.01159989.

#6.) Now extract all of the remaining SNPs only on chromosome 14. These are the trans-acting SNPs. Go back to questions #2-4 and recalculate the eQTL p-values, but this time, you will be using 303 trans-acting SNPs, instead of the 13 cis-acting SNPs to DICER1. Hint: use setdiff() to get the remaining SNPs on chromosome 14, different from the 13 SNPs, after you have subset the map file. Then use these SNP names with the match() function on the original map file.
# get the remaining SNPs on chromosome 14 (trans-acting SNPs)
map14 <- subset(map, V1 == "14")
trans.snps <- setdiff(map14$V2, cis.snps)
# get the SNPs’ positions in the original map file
map.trans.pos <- match(trans.snps, map$V2)
#calculate the positions of these SNPs in the ped file
ped.trans.pos2 <- (map.trans.pos*2) + 6
ped.trans.pos1 <- ped.trans.pos2 - 1
head(map.trans.pos); head(ped.trans.pos2); head(ped.trans.pos1)
# [1] 156 157 158 159 160 161
# [1] 318 320 322 324 326 328
# [1] 317 319 321 323 325 327
#create ped subset for the 303 trans-acting SNPs
ped.trans.subset <- ped[,sort(c(ped.trans.pos1, ped.trans.pos2))]
#reorder trans SNPs to make sure they agree with sorting by positions above
trans.snps <- trans.snps[order(pos2)]
# convert the SNPs from a 2 allele genotype code to a single 1 number genotype code
ped.trans.converted <- apply(ped.trans.subset, 1, recode)
#calculate a linear regression model on the trans-acting SNPs and get eQTL p-values
trans.pvals <- apply(ped.trans.converted, 1, lin.mod, y=as.numeric(exp[probe,]))
head(trans.snps); head(trans.pvals)
# [1] "SNP_A-1882555" "SNP_A-1956152" "SNP_A-1782659" "SNP_A-1846216" "SNP_A-2253929"
# [6] "SNP_A-1903236"
# [1] 0.5366836 0.8634538 0.4652868 0.8634538 0.2908562 0.6085535

#7.) Use match() on the original map file and identify the physical distances for the 13 cis-acting and 303 trans-acting SNPs. Then plot 2 plots (use the par(mfrow=c(2,1)) function): one with the –log10(p-value) vs. physical distance for the cis-acting SNPs and the other for the trans-acting SNPs.
# get physical distances for the obtained cis-acting and trans-acting SNPs
cis.pd <- map[match(cis.snps, map$V2), 'V4']
trans.pd <- map[match(obtained.trans.snps, map$V2), 'V4']
# plotting for cis- and trans-acting SNPs of gene DICER1
par(mfrow=c(2,1))
plot(cis.pd, -log10(cis.pvals), xlab="Physical distance ",
ylab="-log10(p-value)", main="Cis-acting SNPs of DICER1")
plot(trans.pd, -log10(trans.pvals), xlab="Physical distance",
ylab="-log10(p-value)", main="Trans-acting SNPs of DICER1")
