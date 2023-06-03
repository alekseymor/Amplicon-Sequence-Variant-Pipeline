Amplicon Sequence Variant (ASV) processing pipeline

Tools:
MultiQC
DADA2
phyloseq

## Load necessary packages
library(dada2)
library(ggplot2)
library(phyloseq)

# Set filepath
path <- 'Path/to/files'
list.files(path) # verify all of the files of interest are present.

# Sort files into forward and reverse reads.
F_reads <- sort(list.files(path, pattern='_R1_001.fastq'))
R_reads <- sort(list.files(path, pattern='_R2_001.fastq'))

# Obtain sample names from the filenames
sample.names <- sapply(strsplit(F_reads, '_'), `[`, 1)

# Set full path for read files
F_reads <- file.path(path, F_reads)
R_reads <- file.path(path, R_reads)

# Create directory for filtered and trimmed read files
filtered_path <- file.path(path, 'filtered')
# Rename the newly filtered and trimmed files
F_filtered <- file.path(filtered_path, paste0(sample.names, '_F_filtered.fastq.gz'))
R_filtered <- file.path(filtered_path, paste0(sample.names, '_R_filtered.fastq.gz'))

## Filtering and Trimming
output <- filterAndTrim(F_reads, F_filtered, R_reads, R_filtered, 
                        truncLen=c(290,290), # truncLen should be set based on fastQC results
                        maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE)
head(output)

## Estimating Error Rate
errorRate_F <- learnErrors(F_filtered, multithread=TRUE)
errorRate_R <- learnErrors(R_filtered, multithread=TRUE)
# Visualizing the error rates
plotErrors(errorRate_F, nominalQ=TRUE)
plotErrors(errorRate_R, nominalQ=TRUE)


## Dereplicate reads to improve computation
derep_F <- derepFastq(F_filtered, verbose=TRUE)
derep_R <- derepFastq(R_filtered, verbose=TRUE)
names(derep_F) <- sample.names
names(derep_R) <- sample.names


## Sample Inference
F_inference <- dada(derep_F, err=errorRate_F, multithread=TRUE)
R_inference <- dada(derep_R, err=errorRate_R, multithread=TRUE)


## Merging Reads (exact overlap)
merged <- mergePairs(F_inference, derep_F, R_inference, derep_F, verbose=TRUE)


## Create Sequence Table
seq_table <- makeSequenceTable(merged)
table(nchar(getSequences(seq_table))) # Provides count of sequence variants of various lengths

