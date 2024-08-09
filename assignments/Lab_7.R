#In this lab, we will be working with an E.coli genome reference sequence. We will be mapping a set of 1 million paired end reads that I have simulated from this reference sequence (using wgsim from Samtools) in fastq format. Obviously, since these reads were not generated from an instrument, the quality values will all be the same, but we will still work with FASTQ files instead of FASTA files, so you have an idea of the code that you will have to use. Each read is 70mer and the insert size between the 2 mates is 500 bps. You will be using Bowtie to perform this mapping and Samtools to produce various summaries of the information.
#Then we will work with human whole exome data truncated to only chromosome 22. We will utilize VarScan and a web annotation engine to call somatic variants between a liver tumor and normal healthy liver tissue specimen.
#1.) We first need to index the reference sequence, so using the NC_008253.fa file, create the necessary index files for bowtie. Also get the 2 paired end FASTQ files called out1.fq and out2.fq.
# create the necessary index files from the reference sequence
gunzip NC_008253.fa.gz
bowtie-build -f NC_008253.fa e_coli
# get the 2 paired end .fq files
unzip out1.fq.zip
gunzip out1.fq.gz
unzip out2.fq.zip
gunzip out2.fq.gz

