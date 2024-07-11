#Analysis and Interpretation of splice variants

#This lab introduces concepts for analyzing exon array data. The dataset that we are using has both exon- and gene-level summaries, separated into 2 distinct data matrices. The data was taken from a study that was designed to investigate alternative splicing between pediatric patients that have experienced liver rejection after an organ transplant and pediatric patients that have not experienced such rejection. These are blood samples taken from these patients, so the effect sizes are rather modest. As a result, we will be a little more liberal in our exon probe filtering methods to retain as much information as possible.

#1.) Obtain and load in the exon-rma-sketch.summary.txt, gene-rma-sketch.summary.txt, dabg.summary.txt, and HuEx-1_0-st-v2.na24.hg18.probeset_abbr.csv files into R.
unzip('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/exon-rma-sketch.summary.zip', list=TRUE)
# Name Length Date
#1 exon-rma-sketch.summary.txt 62216797 2009-08-05 10:16:00
unzip('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/gene-rma-sketch.summary(1).zip', list=TRUE)
# Name Length Date
#1 gene-rma-sketch.summary.txt 4769999 2009-08-06 12:05:00
unzip('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/dabg.summary.zip', 'dabg.summary.txt', list=TRUE)
# Name Length Date
#1 dabg.summary.txt 62068469 2009-08-05 16:32:00
unzip('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/HuEx-1_0-st-v2.na24.hg18.probeset_abbr (1).csv.zip', 'HuEx-1_0-st-v2.na24.hg18.probeset_abbr.csv', list=TRUE)
# Name Length Date
#1 HuEx-1_0-st-v2.na24.hg18.probeset_abbr (1).csv 115592023 2022-01-31 11:02:00
e <- read.table(unzip('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/exon-rma-sketch.summary.zip', 'exon-rma-sketch.summary.txt'), header=TRUE, row.names=1, skip=70, sep='\t')
g <- read.table(unzip('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/gene-rma-sketch.summary(1).zip', 'gene-rma-sketch.summary.txt'), header=TRUE, row.names=1, skip=70, sep='\t')
p <- read.table(unzip('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/dabg.summary.zip', 'dabg.summary.txt'), header=TRUE, row.names=1, skip=70, sep='\t')

map <- read.csv(unzip('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN/HuEx-1_0-st-v2.na24.hg18.probeset_abbr (1).csv.zip', 'HuEx-1_0-st-v2.na24.hg18.probeset_abbr (1).csv'))

#2.) Continue to follow code in the lecture, up to the line of: you should filter probes based on guidelines…
#define class membership
r <- c(1,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1)
#intersect matching probes between annotation file and data matrix (filtering for core probesets)
x <- intersect(dimnames(e)[[1]], dimnames(map)[[1]])
#subset the rows in the exon, p-value matrix, and annotation file to the intersecting probes
e <- e[x,]
p <- p[x,]
map <- map[x,]
#get unique transcript cluster IDs (gene IDs) from annotation file
u <- unique(as.character(map$transcript_cluster_id))
u <- intersect(u, dimnames(g)[[1]])

#3.) Filter probes in the exon matrix using the detection above background (DAGB) p-values. If >50% of the subjects in BOTH groups have a p-value > .05, remove this exon probe. Note that this is a little more liberal than the thresholds that we explained in the lecture. The effects are small for this dataset, so this is why we are being a little less conservative with our filtering. How many exon probes remain after this filtering? Subset the exon matrix by those probes that pass the filter.
#define a function to count how many samples have a detection P > 0.05 and apply to each group separately
count.undet <- function(x){length (which(x>0.05))}
group1.undet <- apply(p[as.logical(r)], 1, count.undet)
group2.undet <- apply(p[!as.logical(r)], 1, count.undet)
#remove probesets if >50% of the subjects in BOTH groups have P > 0.05
joint_undet <- intersect(which(group1.undet>0.5*length(p[as.logical(r)])), which(group2.undet>0.5*length(p[!as.logical(r)])))
e.rm <- e[joint_undet,]
e.fil <- e[!(rownames(e) %in% rownames(e.rm)), ]
e <- e.fil
nrow(e)
#[1] 163972
#163,972 exon probes remain after filtering out those that are “undetected” at this threshold level (i.e. where the DABG p-value > 0.05 in >50% of BOTH groups).

#4.) Go back to the annotation table (i.e. map) and subset it by the shorter set of exon probes that you will now be using. Get the unique transcript cluster IDs from this new annotation table and intersect them with the gene matrix. You should now have a shorter u variable. How many unique transcript cluster IDs do you have?
map <- map[rownames(e),]
u <- unique(as.character(map$transcript_cluster_id))
u <- intersect(u, dimnames(g)[[1]])
length(u)
#[1] 16659
#16,659 unique transcript cluster IDs are identified after this additional filtering step.

#5.) Now, using the code in the lecture, write a loop that loops through all values of the variable u (unique transcript cluster ID). In each loop iteration, you need to:
#a. Get the correct mapping between the ith transcript cluster ID and the exon probe IDs
#b. Subset the gene and exon matrices by the appropriate transcript cluster ID or exon probe IDs
#c. Call the exon.ni function
#d. Store the output p-values, test statistics, and CIs from this function for each transcript cluster ID (and exon probe IDs) 
#For this loop, you also need to add an if statement to only conduct this calculation when you have at least 2 rows in the d.exon matrix (i.e. at least 2 exon probe IDs for a single transcript cluster ID). This loop will take up to an hour to run (depending on your computer speed), so plan accordingly.
#exon ni value calculation and t - test function
#takes transcript cluster values, exon probes for this transcript cluster, and classification vector
exon.ni <- function(genex, exonx, rx) {
  ni <- t(t(exonx) - genex)
  ttest <- t(apply(ni, 1, t.two, sam=as.logical(r), v=F))
  return(ttest)
}
#two-sample t-test function (called by an apply() function)
t.two <- function(x, sam, v=F) {
  x <- as.numeric(x)
  out <- t.test(as.numeric(x[sam]), as.numeric(x[!sam]), alternative="two.sided", var.equal=v)
  o <- as.numeric(c(out$statistic, out$p.value, out$conf.int[1], out$conf.int[2]))
  names(o) <- c("test_statistic", "pv", "lower_ci", "upper_ci")
  return(o)
}
ni.df <- data.frame()
for (i in u) {
  ex <- dimnames(map[map$transcript_cluster_id==i,])[[1]]
  d.exon <- e[ex,]
  d.gene <- g[i,]
  if (nrow(d.exon) >= 2) {
    ni.out <- exon.ni(genex=as.numeric(d.gene), exonx=d.exon, rx=r)
    ni.df <- rbind(ni.df, ni.out)
  }
}
dim(ni.df)
#[1] 157492 4

#6.) Sort the matrix where you stored all of the summary statistics by p-value and get the transcript cluster ID with the lowest p-value.
#sorting the summary stats matrix by p-value
ni.df.sorted <- ni.df[order(ni.df$pv),]
#getting the transcript cluster ID with the lowest p-value
ni.df.sorted[1,]
ti <- as.character(map[rownames(ni.df.sorted[1,]),]$transcript_cluster_id)
ti
#[1] "2426951"
#The transcript cluster ID with the lowest p-value is 2426951.

#7.) Run the exon plot function with this transcript cluster ID and provide the plot in your solutions. You will need to redefine d.exon and d.gene to do this.
plot.exons <- function(exonx,genex,rx,ti) {
  rr <- rx
  rx <- rep(rx,nrow(exonx))
  rx[rx==1] <- "A"
  rx[rx==0] <- "B"
  rx <- as.factor(rx)
  ni <- t(t(exonx)-genex)
  exonx <- as.data.frame(t(ni))
  ex.stack <- stack(exonx)
  d <- data.frame(ex.stack,rx)
  names(d) <- c("exon_values","exon_id","class")
  d$exon_id <- as.factor(d$exon_id)
  d$class <- as.factor(d$class)
  genex.title <- as.character(map[match(ti,as.character(map$transcript_cluster_id)),"gene_assignment"])
  plot(c(.5,(ncol(exonx)+.5)),range(d[,1]),type="n",axes=F,xlab="",ylab="")
  boxplot(exon_values~exon_id,add=T,subset=d$class=="A",d,col="salmon",border='red',cex.axis=.75,las=2,xlab='Exon probe ID',ylab='Log2 normalized intensity',main=paste("Gene ID:",ti,"\n",genex.title),boxwex=0.4)
  boxplot(exon_values~exon_id,subset=d$class=="B",d,add=T,col="green",border='darkgreen',axes=F,boxwex=0.4, at=c(1:ncol(exonx))+0.1)
  legend("topleft", fill=c('salmon', 'green'), legend=c("Class A/1", "Class B/0"), cex=0.7, bty='n')
}
ex <- dimnames(map[map$transcript_cluster_id==ti,])[[1]]
d.exon <- e[ex,]
d.gene <- g[ti,]
plot.exons(exonx=d.exon, genex=as.numeric(d.gene), rx=r, ti=ti)
