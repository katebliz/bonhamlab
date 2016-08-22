# R - Script: SampleQualityandRarefaction.R
# Written by: Michael Loesche
#
# Description: These are commands used to generate stats and graphs to visualize output
#              files generated during the psoriasis baseline analysis QIIME run.
#
# Arguments: map_file, biom_summary, factors...

# Load required libraries or install them if missing
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE, repos="http://cran.rstudio.com/")
  require(p, character.only = TRUE)
}
usePackage('qiimer')
usePackage('ggplot2')
usePackage('scales')
usePackage('getopt')
usePackage('dplyr')

# Make matrix specifying the available flags
spec <-matrix(c('help','h',0,'logical',
                'project_name','p',1,'character',
                'map','m',1,'character',
                'biom_summary','b',1,'character',
                'collated_alpha','a',1,'character',
                'factors','f',1,'character'), ncol=4,byrow=T)
# Read in command line arguments
opt <- getopt(spec)

# Very crude help menu
if ( !is.null(opt$help) ) {
  cat(getopt(spec, usage=TRUE))
  q(save = 'no',status=1)
}
# Quit if missing required flags
if(is.null(opt$project_name) | is.null(opt$map) | is.null(opt$biom_summary) | is.null(opt$collated_alpha)){
  cat('\n**A required option is missing!**\nCheck help file for more info -h\n\n');cat(getopt(spec, usage=TRUE));q('no',1)
}

# Store arguments into variables
proj <- opt$project_name
map_file <- opt$map
biom_summary <- opt$biom_summary
collated_path <- opt$collated_alpha
# Split input factors by comma. If no factors were entered create a dummy factor
# that will be used to refer to whole data set.
if(!is.null(opt$factors)) {
  factors <- factor(tmp<-unlist(strsplit(opt$factors,',',fixed=T)))
} else { factors <- as.factor('Data') }

# Read in mapping file
map <- read.table(map_file,sep='\t',header=T,colClasses='factor',comment.char='')
names(map)[1] <- 'SampleID'
# Create dummy factor column with one level refering to whole data set
if(is.null(opt$factors)) map$Data <- as.factor(rep('All Data',nrow(map)))

#############################################################################
# Sample Quality Analysis
#############################################################################
## Sink output
sink(paste('./',proj,'_reads_analyses.txt',sep=''), append = T)

# Read in the biom summary table for sequences per sample.
reads <- read.table(file = biom_summary, skip = 16, col.names = c('SampleID', 'reads'), sep = ':', header = F)

# Make vector of samples missing from BIOM table but present in map.
diff <- setdiff(map$SampleID, reads$SampleID)
cat('Sample Quality Analysis:\n==============================\nSamples in map file missing from biom summary:', diff,'\n\n')

# Merge mapping file with BIOM summary
map_reads <- merge(map, reads, by = 'SampleID')
rm(reads)
attach(map_reads)

# Get reads/sample by percentiles
summary(reads)
percentiles <- c(0.01, 0.05, 0.10, 0.15, 0.2, 0.25)
cat('Reads/Sample Percentiles:\n')
print(quantile(reads, percentiles))

# Get percentiles of reads/sample
cat('\nProportion of samples above reads/sample cutoff:\n')
cutoff_lengths <- c(500, 1000, 1500, 2000, 2500, 5000, 10000, 15000)
print(sapply(cutoff_lengths, FUN = function(x){
  # Gets percentile, rounds to one decimal and then adds a '%' to each number.
  percent <- round(sum(reads >= x)/length(reads),3)
  names(percent) <- x
  return(percent)
}))

## Plots

pdf(paste('./',proj,'_reads_per_sample_and_collectors_curve_plots.pdf',sep=''))
# Make density plot of reads per sample
plot(density(log10(reads)), xlim=c(0,max(log10(reads))), main = 'Density Plot of Reads per Sample', xlab = 'log10 of Reads per Sample')
cat('\nMade density plot\n')

# Make plots of reads per sample by Sequencing Run
for(f in factors) {
  print(ggplot(map_reads, aes_string(x = f, y = 'reads', fill = f)) + 
          geom_boxplot(notch = F, outlier.size = 1.3) + 
          scale_y_continuous(trans=log10_trans(), name='Reads/Sample') + 
          ggtitle('Boxplot of Reads per Sample'))
}
cat('Made boxplot(s) by factor\n')

# Return output to console
sink()

# Detach data frame
detach(map_reads)

#############################################################################
# Rarefaction plots
#############################################################################

## Read in collated alpha tables and format data

# Get metrics used from collated file names
metrics <- gsub(list.files(collated_path), pattern = '.txt', replacement = '')
names(metrics) <- metrics

# Read in rarefaction data
rare <- lapply(metrics, FUN = function(m) {
  # Read in rarefaction data
  tmp <- read_qiime_rarefaction(paste(collated_path, m, '.txt',sep=''))
  # Fix SampleID name problem: remvoes leading "X"
  levels(tmp$SampleID) <- gsub(levels(tmp$SampleID), pattern = 'X', replacement = '')
  # Collapse multiple rarefactions into means
  rarefaction_stats(tmp)
})

# Merge map with rarefaction table
map_rare <- lapply(rare, merge, map, by = 'SampleID')
# Sort data frame so that ggplot paths go in order
map_rare <- lapply(map_rare, arrange, SampleID, sequences_per_sample)

## Make rarefaction plots

# Loop through the different metrics
tmp <-lapply(names(map_rare), FUN = function(m) {
  mr <- map_rare[[m]]
  # Loop through the factors of interest
  lapply(factors, FUN = function(f) {
    # Turn f into character
    f <- as.character(f)
    # Make a group_by object
    grouped <- group_by_(mr, .dots = list(as.symbol('sequences_per_sample'), as.symbol(f)))
    # Get group means and SE's
    df <- summarise(grouped, diversity = mean(diversity.mean), se = sd(diversity.mean)/sqrt(n()))
    # Make rarefaction plot
    print(ggplot(df, aes_string(x='sequences_per_sample',y='diversity',group=f,color=f))+
      geom_path(size=1)+geom_errorbar(aes(ymin=diversity-se,ymax=diversity+se))+
        ggtitle(paste('Collector\'s Curve:',m)) + ylab(m))
  })
})
cat('Made collector\'s curve(s)\n\n',file = paste('./',proj, '_reads_analyses.txt', sep = ''),append=T)

# Turn off PDF
dev.off()

