#In this lab, we will be analyzing copy number alterations from a set of 9 arrays. These arrays were taken from a prostate cancer cell line study, where each array represents a different cell line. I have already processed the raw CEL files using the APT software and provided you the formatted output data for all 23 chromosomes for these 9 cell lines, to save you time. You will be implementing some tools in various R packages to visualize CN alterations and understand the similarity between these 9 cell lines. Then you will use ABSOLUTE to calculate the ploidy and tumor purity for one of the cell lines in the dataset.

#1.) Get the copy number file from the website called CNV_processed_files.zip. This contains both CN state values and log2 ratio values for the 9 samples.
setwd('/Users/josephchristina/Desktop/JHU_Y1/AdvancedGGN')
cn_states <- read.table('CNV_processed_files/cn_states.txt')
log_ratios <- read.table('CNV_processed_files/log2ratios.txt')

#2.) Subset the data by chromosome 13, run smoothing and segmentation with the DNAcopy library. Then plot the curves for all 9 samples for chromosome 13.
cn_states <- cn_states[cn_states$Chromosome == '13',]
