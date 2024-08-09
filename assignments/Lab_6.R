#In this lab, we will be working with a ChIP-Seq dataset to identify TF binding sites and chromatin modifications from the authors of PeakSeq. However, we will be using MACS, rather than PeakSeq to analyze it. As bioinformaticians, we often find that when working with multiple software tools, we have to format and reformat data to accommodate the requirements of the tool being used. You will be using some simple shell commands and R (or another language of your choice) to format this data into a format that MACS will accept. We will then run MACS to identify the peaks in the dataset, and plot the results in R.

#1.) Get the file called FC305JN_s_1_eland_result.txt.gz, unzip it and extract only lines with “chr6.fa” since we will only be analyzing the sequence information from chromosome 6. How many lines are in the new file that you created?
# in R:
setwd('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN')
# read in file and extract only lines with ‘chr6.fa’
txt <- readLines("FC305JN_s_1_eland_result.txt")
txt_chr6 <- txt[grep("chr6.fa", txt)]
# count number of lines in subset
length(txt_chr6)
#[1] 166056
# write to a new file
write(txt_chr6, file='chipseq_chr6.txt')
#There are 166,056 lines in the new file that contains only lines with “chr6.fa”.

#2.) Now since we don’t need all of the columns, extract columns 1,2,7,8,9 from this chromosome 6 file that you created in #1. These columns correspond to the read name, read sequence, chr6.fa, locus, and strand. One way to do this is with the awk command, although there are other ways to do this as well.
# on command line:
#awk '{print $1,$2,$7,$8,$9}' FS='\t' OFS='\t' chipseq_chr6.txt > chipseq_chr6_select.txt

#3.) In R, or a language of your choice, read in the 5 column file you created in #2 and count up the number of unique locus values for each strand separately. Then output a MACS formatted file (see lecture slides) that contains the appropriate fields for each strand. You will need to change F/R to +/- in addition to providing the frequencies for each locus (read). You may also need to cast a character field or two to numeric with something like as.numeric(as.matrix(x)) to appropriately order the locus values in ascending order. Be sure to check your tag counts and coordinates in R to make sure that they are not accidentally cast to character.
# reading in new 5 column file
txt_chr6_select <- read.table('chipseq_chr6_select.txt', sep='\t')
colnames(txt_chr6_select) <- c('read_name', 'read_sequence', 'chr6.fa', 'locus', 'strand')
# a) counts of unique locus values for strands separately
library(dplyr)
txt_chr6_select %>%
  group_by(strand, locus) %>%
  summarise(count = n(),
    seq_length = nchar(first(read_sequence)),
    .groups = "drop") %>%
  ungroup() %>%
  arrange(locus) -> n_loci
nrow(n_loci[n_loci$strand=='F',]); nrow(n_loci[n_loci$strand=='R',])
#[1] 68183
#[1] 68328
# b) counts of unique loci values for strands separately (considering read seq uniqueness)
txt_chr6_select %>%
  group_by(strand, read_sequence, locus) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  arrange(locus) -> n_loci2
nrow(n_loci2[n_loci2$strand=='F',]); nrow(n_loci2[n_loci2$strand=='R',])
#[1] 74452
#[1] 74405
#Not considering the read sequence values, there are 68,183 unique loci values for the forward strand and 68,328 unique loci values for the reverse strand. However, when also taking into consideration uniqueness of the read sequences, there are 74,452 and 74,405 unique loci values for the forward and reverse strands, respectively.
# creating the MACS input file
macs_input <- matrix(nrow=nrow(n_loci), ncol=6, dimnames=list(NULL, c('Chr', 'start', 'end', '??', 'tags', 'sense')))
macs_input[,1] <- 'chr6'
macs_input[,2] <- n_loci$locus
macs_input[,3] <- n_loci$locus + n_loci$seq_length
macs_input[,4] <- 0
macs_input[,5] <- n_loci$count
macs_input[,6] <- ifelse(n_loci$strand=='F', '+', '-')
write.table(macs_input, 'macs_input.txt', sep='\t', quote=FALSE, row.names=FALSE, col.names=FALSE)

#4.) You are now ready to run MACS on the file that you created in #3. Use the –t argument and the –name arguments and run MACS to identify significant regions. We will use the default p-value threshold of 1e-5. Print out the table in the *.xls file of the significant regions (this table should have <400 rows). You can also run the R script that is output with MACS to see the peak distributions.
# on command line: create a virtual environment to run MACS3
conda create -n myenv python=2.7
conda activate myenv
cd MACS-macs_v1
python setup.py install
# run MACS on the new file
macs --treatment ../macs_input.txt --name ../macs_output --pvalue 1e-5
# MACS running updates
INFO @ Thu, 14 Mar 2024 15:57:38:
# ARGUMENTS LIST:
# name = ../macs_output
# format = AUTO
# ChIP-seq file = ../macs_input.txt
# control file = None
# effective genome size = 2.70e+09
# band width = 300
# model fold = 10,30
# pvalue cutoff = 1.00e-05
# Large dataset will be scaled towards smaller dataset.
# Range for calculating regional lambda is: 10000 bps

#INFO @ Thu, 14 Mar 2024 15:57:38: #1 read tag files...
#INFO @ Thu, 14 Mar 2024 15:57:38: #1 read treatment tags...
#INFO @ Thu, 14 Mar 2024 15:57:38: Detected format is: BED
#INFO @ Thu, 14 Mar 2024 15:57:39: #1 tag size is determined as 28 bps
#INFO @ Thu, 14 Mar 2024 15:57:39: #1 tag size = 28
#INFO @ Thu, 14 Mar 2024 15:57:39: #1 total tags in treatment: 136511
#INFO @ Thu, 14 Mar 2024 15:57:39: #1 user defined the maximum tags...
#INFO @ Thu, 14 Mar 2024 15:57:39: #1 filter out redundant tags at the same location and the same strand by allowing at most 1 tag(s)
#INFO @ Thu, 14 Mar 2024 15:57:39: #1 tags after filtering in treatment: 136511
#INFO @ Thu, 14 Mar 2024 15:57:39: #1 Redundant rate of treatment: 0.00
#INFO @ Thu, 14 Mar 2024 15:57:39: #1 finished!
#INFO @ Thu, 14 Mar 2024 15:57:39: #2 Build Peak Model...
#INFO @ Thu, 14 Mar 2024 15:57:39: #2 number of paired peaks: 0
#WARNING @ Thu, 14 Mar 2024 15:57:39: Too few paired peaks (0) so I can not build the model! Broader your MFOLD range parameter may erase this error. If it still can't build the model, please use --nomodel and --shiftsize 100 instead.
#WARNING @ Thu, 14 Mar 2024 15:57:39: Process for pairing-model is terminated!
#WARNING @ Thu, 14 Mar 2024 15:57:39: #2 Skipped...
#WARNING @ Thu, 14 Mar 2024 15:57:39: #2 Use 100 as shiftsize, 200 as fragment length
#INFO @ Thu, 14 Mar 2024 15:57:39: #3 Call peaks...
#INFO @ Thu, 14 Mar 2024 15:57:39: #3 shift treatment data
#INFO @ Thu, 14 Mar 2024 15:57:39: #3 merge +/- strand of treatment data
#INFO @ Thu, 14 Mar 2024 15:57:39: #3 call peak candidates
#INFO @ Thu, 14 Mar 2024 15:57:40: #3 use self to calculate local lambda and filter peak candidates...
#INFO @ Thu, 14 Mar 2024 15:57:41: #3 Finally, 349 peaks are called!
#INFO @ Thu, 14 Mar 2024 15:57:41: #4 Write output xls file... ./../macs_output_peaks.xls
#INFO @ Thu, 14 Mar 2024 15:57:41: #4 Write peak bed file... ./../macs_output_peaks.bed
#INFO @ Thu, 14 Mar 2024 15:57:41: #4 Write summits bed file... ./../macs_output_summits.bed
#INFO @ Thu, 14 Mar 2024 15:57:41: #5 Done! Check the output files!

# checking *.xls file output
#cat ../macs_output_peaks.xls
#The *.xls table has 349 rows. No R script from the output was locatable.

#5.) In R, make a plot similar to that below where the blue peaks are the read counts based on the unique locus calculations you conducted in #4 and the red peaks are those found using MACS. Use the polygon function for those regions identified by MACS since these span multiple physical positions. Also only make your y-axis up to 100 tags, so we can see the patterns better.
# read in MACS peaks output data
macs_output <- read.table('macs_output_peaks.xls', skip=19, header=TRUE)
# plotting peaks
plot(n_loci$locus, n_loci$count, type="h", col='blue', ylim=c(0, 100), main='Chr 6 - Peaks identified by MACS', axes=FALSE, xlab='Position', ylab='No. of tags')
for (i in 1:length(macs_output$start)) {
  polygon(c(macs_output$start[i], macs_output$start[i], macs_output$end[i], macs_output$end[i]),
          c(0, macs_output$tags[i], macs_output$tags[i], 0), density=1, col='red')
}
options(scipen=999)
axis(1, las=2, at=seq(min(n_loci$locus), max(n_loci$locus), length.out=40), cex.axis=0.5)
axis(2, cex.axis=0.9)

#6.) Why is it likely that the 4 largest blue peaks were not identified as significant regions?
#The 4 largest blue peaks may not have been identified as significant regions by the MACS algorithm possibly because those they might not have a high enough fold enrichment relative to the background, especially since there appears to be elevated read density overall across the chromosome. Since MACS uses statistical tests to determine significance, it also depends on the p-value cutoff setting that we use. Additionally, according to the HBC Training Program’s page on peak calling with MACS, the algorithm relies on detecting bimodal enrichment patterns of the tags and may shift the tags to improve the accuracy of peak calling (and better representing the most likely protein-DNA interaction sites).

#7.) Why do you think there is a gap at the physical position ~60,000,000 for this chromosome?
#Given that chromosome 6 is approximately 170 Mbs long and submetacentric, my guess is that its centromere is located at around the 60 Mbs region. In that case, the 60 Mbs region may be heterochromatic and rich with repetitive DNA elements, so it cannot to be captured well by ChIP-seq. The Biostars thread below clarifies that this region would appear as a gap where no reads map to when using certain sequencing techniques, such as ChIP-seq.

#Reference:
#https://www.biostars.org/p/11133/#32826
