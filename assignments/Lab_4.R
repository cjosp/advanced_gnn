#In this lab, we will be analyzing copy number alterations from a set of 9 arrays. These arrays were taken from a prostate cancer cell line study, where each array represents a different cell line. I have already processed the raw CEL files using the APT software and provided you the formatted output data for all 23 chromosomes for these 9 cell lines, to save you time. You will be implementing some tools in various R packages to visualize CN alterations and understand the similarity between these 9 cell lines. Then you will use ABSOLUTE to calculate the ploidy and tumor purity for one of the cell lines in the dataset.

#1.) Get the copy number file from the website called CNV_processed_files.zip. This contains both CN state values and log2 ratio values for the 9 samples.
setwd('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN')
cn_states <- read.table('CNV_processed_files/cn_states.txt')
log_ratios <- read.table('CNV_processed_files/log2ratios.txt')

#2.) Subset the data by chromosome 13, run smoothing and segmentation with the DNAcopy library. Then plot the curves for all 9 samples for chromosome 13.
cn_states <- cn_states[cn_states$Chromosome == '13',]
log_ratios <- log_ratios[log_ratios$Chromosome == '13',]
library(DNAcopy)
#create a CNA object
CNA_object <- CNA(cbind(log_ratios$Log2Ratio, log_ratios$m.Log2Ratio, log_ratios$m.Log2Ratio.1, log_ratios$m.Log2Ratio.2, log_ratios$m.Log2Ratio.3, log_ratios$m.Log2Ratio.4, log_ratios$m.Log2Ratio.5, log_ratios$m.Log2Ratio.6, log_ratios$m.Log2Ratio.7), log_ratios$Chromosome, log_ratios$Position, data.type='logratio', sampleid=c('Sample 1','Sample 2','Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9'))
#run smoothing
dat_smoothed <- smooth.CNA(CNA_object)
#run circular binary segmentation with default parameters
dat_segment <- segment(dat_smoothed, verbose=1)
#plot curves for the 9 samples for chromosome 13
pdf("CNV_plot.pdf")
plot(dat_segment, plot.type="chrombysample", pt.cex=0.5, lwd=0.5)
dev.off()

#3.) Now do the same type of plots using the plot() function for chromosome 13 using the CN state values for the 9 samples. Make sure to title and label things appropriately.
#create another CNA object
CNA_object2 <- CNA(cbind(cn_states$CNState, cn_states$m.CNState, cn_states$m.CNState.1, cn_states$m.CNState.2, cn_states$m.CNState.3, cn_states$m.CNState.4, cn_states$m.CNState.5, cn_states$m.CNState.6, cn_states$m.CNState.7), cn_states$Chromosome, cn_states$Position, data.type='logratio', sampleid=c('Sample 1','Sample 2','Sample 3', 'Sample 4', 'Sample 5', 'Sample 6', 'Sample 7', 'Sample 8', 'Sample 9'))
#run smoothing
dat_smoothed2 <- smooth.CNA(CNA_object2)
#run circular binary segmentation with default parameters
dat_segment2 <- segment(dat_smoothed2, verbose=1)
#plot curves for the 9 samples for chromosome 13
pdf("CNV_plot2.pdf")
plot(dat_segment2, plot.type="chrombysample", pt.cex=0.5, lwd=0.5)
dev.off()

#4.) Using the CNseg(), getRS(), and rs() functions in the CNTools package on the output from the segment() function, create a data matrix, run correlation between the 9 samples, and plot this correlation heatmap using the image() function. Be sure to label the x and y-axes appropriately and title the plot.
library(CNTools)
#create a data matrix of the combined CN vectors
seg <- CNSeg(dat_segment$output)
rs.region <- getRS(seg, by="region", imput=FALSE, XY=FALSE, what='mean')
mat <- rs(rs.region)
mat_vals <- data.matrix(mat[,4:12], rownames.force=NA)
mat.cor <- cor(mat_vals, method='pearson')
#plot correlation heatmap
library(gplots)
layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 5, 2, byrow=TRUE))
cx = rev(colorpanel(25,"red","white","blue"))
leg = seq(min(mat.cor,na.rm=T), max(mat.cor,na.rm=T), length=10)
par(oma=c(5,5,1,1))
image(mat.cor, main='Correlation Heatmap for 9 Cell Lines; Chromosome 13', axes=F, col=cx)
axis(1, at=seq(0,1,length=ncol(mat.cor)), label=dimnames(mat.cor)[[2]], cex.axis=0.9, las=2)
axis(2, at=seq(0,1,length=ncol(mat.cor)),label=dimnames(mat.cor)[[2]], cex.axis=0.9, las=2)
par(mar=rep(2, 4))
image(as.matrix(leg), col=cx, axes=F)
tmp <- round(leg,2)
axis (1, at=seq(0,1,length=length(leg)), labels=tmp, cex.axis=1)

#5.) Now go back to the log2 ratio file and extract the first 5 columns, which includes just the first cell line using awk, or some other linux command. Then read this 5 column file into R and load the ABSOLUTE and numDeriv packages. The ABSOLUTE R package is available on the course website.
#on CLI…. awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' log2ratios.txt > log2ratios_5.txt
#in R:
log_ratios5 <- read.csv('CNV_processed_files/log2ratios_5.txt', header=TRUE, sep='\t', col.names=c('ProbeSetName', 'Chromosome', 'Position', 'Log2Ratio'))
library(numDeriv)
library(ABSOLUTE)

#6.) Run CNA(), smoothing, and segmentation with the DNAcopy library across the genome for this single cell line.
#create another CNA object
CNA_object3 <- CNA(cbind(log_ratios5[,4]), log_ratios5$Chromosome, log_ratios5$Position, data.type='logratio', sampleid='Sample 1')
#run smoothing
dat_smoothed3 <- smooth.CNA(CNA_object3)
#run circular binary segmentation with default parameters
dat_segment3 <- segment(dat_smoothed3, verbose=1)

#7.) Use the output of the segmentation to run ABSOLUTE with the following parameters (we are coding this as a sequencing technology, though it’s really an array technology):
#sigma.p <- 0
#max.sigma.h <- 0.02
#min.ploidy <- 0.95
#max.ploidy <- 10
#max.as.seg.count <- 1500
#max.non.clonal <- 0
#max.neg.genome <- 0
#genome <- "hg19"
#platform <- "Illumina_WES"
#This step will take approximately a few hours to run, so plan accordingly. Now looking just at the pdf output generated from the RunAbsolute() function, for the optimal solution determined, what is tumor purity calculated for this cell line as well as approximate ploidy state?
#...
# output from DNAcopy segment()
dat <- dat_segment3$output
names(dat) <-c("ID", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")
dat.fn <-"test.abs.dat"
write.table(file=dat.fn, dat, quote=F, col.names=T, row.names=F, sep="\t")
sigma.p <- 0
max.sigma.h <- 0.02
min.ploidy <- 0.95
max.ploidy <- 10
max.as.seg.count <- 1500
max.non.clonal <- 0
max.neg.genome <- 0
genome <- "hg19"
platform <- "Illumina_WES"
sample.name <- "Sample 1"
primary.disease <- "cancer"
RunAbsolute(dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform, sample.name, ".", max.as.seg.count, max.non.clonal, max.neg.genome, "total", verbose=TRUE)
#The optimal solution determined suggests that the tumor purity calculated for this cell line (sample 1) is 0.26 and the approximate ploidy state is 2.18.

#External reference: https://www.genepattern.org/analyzing-absolute-data#gsc.tab=0
