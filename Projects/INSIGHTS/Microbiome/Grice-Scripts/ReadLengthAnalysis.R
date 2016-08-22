# File: ReadLengthAnalysis.R
# Written by: loesche
#
# Takes in the read length summary & size cutoffs and creates a summary of the
# reads, the filtering, and makes a density plot in the working directory.
#
# Arguments: sum_file, min, max

# Read in command line arguments
arguments <- commandArgs(T)
sum_file <- arguments[1]
min <- as.integer(arguments[2])
max <- as.integer(arguments[3])

# Get project name
proj <- gsub(sum_file, pattern = '_specific_raw_SizeStats.txt', replacement = '')

# Read in the output summary
reads <- read.delim(sum_file)

# Attach data
attach(reads)

## Make density plot
# Turn summary into a long format so that density can "count" the data points
read_length_counts <- rep(Size, Num.Seqs)
# Open PDF
pdf(paste(proj,'_read_length_density_plots.pdf', sep = ''))
# Make density plot showing all reads
plot(density(read_length_counts), main = 'Density of Read Lengths', xlab = 'Read Lengths')
# Make density plot zoomed in to the selected for region
plot(density(read_length_counts), main = 'Density of Read Lengths Zoomed In', xlim = c(min-30,max+30), xlab = 'Read Lengths')
# Make lines showing the min and max
abline(v = min, col = 'red'); abline(v=max, col='red')
dev.off() # Closes PDF

# Make summary stats
total_reads <- length(read_length_counts)
bad_reads <- sum(read_length_counts < min | read_length_counts > max)
good_reads <- sum(read_length_counts >= min & read_length_counts <= max)
prop_bad <- bad_reads/total_reads
prop_good <- good_reads/total_reads

# Print out the results
cat('Size Filtering Analysis:\n==============================\nTotal Reads:',
    total_reads,'\nReads Kept:',good_reads,'\nReads Discarded:',bad_reads,
    '\n\nProportion of reads kept:',prop_good,'\nProportion of reads discarded:',
    prop_bad,'\n\nMade density plot\n\n',file = paste('./',proj, '_reads_analyses.txt', sep = ''))

