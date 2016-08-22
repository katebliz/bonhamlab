# Script to read in OTU table

# Read in Basic Mapping Info ------------------------------------------------------------------

# read in metadata
source('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/dfu_metadata.R')

library(biom)
library(qiimer)
library(loescher)
library(plyr)
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(ape)
library(picante)
library(GUniFrac)

# Set wd
setwd('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/')


# Read in OTU Table ------------------------------------------------------------------------------

# Read OTU tables
b <- read_biom('./uclust_no_singletons.biom')
b <- read_biom('./dfu_otu_table.biom')

# Extract OTU by sample counts
otu <- as.data.frame(as.matrix(biom_data(b)))
# Put samples in order of mapping file
otu <- otu[,as.character(dfu$SampleID[dfu$has_raw_sample])]

# Remove samples or OTUs as necessary here
otu <- otu[rowSums(otu) != 0,]

# Extract OTU taxonomy assignments
#otu_taxa <- pretty_taxa(biom_taxonomy(b))
#otu_taxa <- otu_taxa[row.names(otu),]

# Read in the blast assignments subset those in OTU table  
blast <- read.csv('./blast_assignments_gg_no_singletons.csv')
blast <- read.csv('./blast_assignments_gg.csv')
# Make otu taxonomy translation data frame
otu_taxa <- gsub(pattern = '.__', replacement = '', x = nameify(blast$Taxa,blast$OTU)) %>%
  strsplit(split = '; ') %>%
  lapply(function(x) x[x != '']) %>%
  pretty_taxa(add_rank = T)
# Subset otu_taxa to OTUs in otu table
otu_taxa <- filter(otu_taxa, OTU %in% row.names(otu))
# Remove OTUs that could not be classified by blast
otu <- otu[otu_taxa$OTU,]



# Taxonomic Table -----------------------------------------------------------------------------

taxa <- aggregate(as.matrix(otu), by = list(otu_taxa$Genus, otu_taxa$Species), sum)

otu_taxa <- otu_taxa %>% mutate(unique = paste(Genus,Species)) %>% 
  select(-OTU) %>%
  group_by(unique) %>% 
  summarise_each(funs(first))
otu_taxa$OTU <- paste('tax',1:nrow(otu_taxa))

taxa <- taxa %>% mutate(unique = paste(Group.1,Group.2)) %>%
  arrange(unique) %>%
  select(-unique,-Group.1,-Group.2)
otu <- taxa
row.names(otu) <- paste('tax',1:nrow(otu))

# Sum together duplicates ---------------------------------------------------------------------

# Make list of duplicates
dup <- filter(dfu, !is.na(Notes), SubjectID != 'Control')$SampleID

# Merge samples by summation
otu$merge178_1 <- otu[,dup[1]] + otu[,dup[2]]
otu$merge178_6 <- otu[,dup[3]] + otu[,dup[4]]

# Remove original samples from OTU table
otu <- otu[,!colnames(otu) %in% dup]

# Combine samples metadata by removing the second and replacing the sampleid
dfu <- filter(dfu, is.na(Notes) | Notes != 'Second')
dfu$SampleID[dfu$SubjectID == '178' & dfu$Visit == 1] <- 'merge178_1'
dfu$SampleID[dfu$SubjectID == '178' & dfu$Visit == 6] <- 'merge178_6'

# rearrange otu to be in order of mapping file
otu <- otu[,dropNA(dfu$SampleID)]

# Add total sample reads to 
dfu$RawReads <- NA
dfu$RawReads[dfu$has_raw_sample] <- colSums(otu)


# Look at controls ----------------------------------------------------------------------------

swab <- otu[,c('115087','115088')] %>% as.data.frame %>% 
  add_rownames('OTU') %>% 
  left_join(otu_taxa, by = 'OTU') %>%
  rename(swab1 = `115087`, swab2 = `115088`) %>%
  group_by(Genus,Species) %>%
  summarise(swab1 = sum(swab1), swab2 = sum(swab2)) %>%
  ungroup %>% 
  filter(swab1 > 0 | swab2 > 0) %>%
  arrange(desc(swab1 + swab2))

# Plot relative abundances of swabs - filter at 0.01
swab %>% mutate(swab1 = swab1/sum(swab1), swab2 = swab2/sum(swab2)) %>% 
  mutate(cs1 = cumsum(swab1), cs2 = cumsum(swab2)) %>%
  filter(swab1 >= 0.01, swab2 >= 0.01) %>% 
  gather('Sample','Proportion',swab1:swab2) %>% 
  ggplot(aes(x = Sample, y = Proportion, fill = Genus)) + 
  geom_bar(stat = 'identity')

# How much do those taxa represent other samples?
cbind(Taxa = otu_taxa$Genus, as.data.frame(otu)) %>% 
  filter(Taxa %in% c('g__Geobacillus','g__Bacillus','g__Lactococcus')) %>%
  group_by(Taxa) %>%
  summarise_each(funs(sum)) %>%
  gather('SampleID','Reads', -Taxa) %>%
  group_by(Taxa) %>%
  summarise(Reads = sum(Reads)) %>%
  mutate(Proportion = Reads/sum(otu))

rm(swab)

# Remove Contaminants -------------------------------------------------------------------------

# Remove lactococcus from table
otu <- otu[otu_taxa$Genus != 'g__Lactococcus',]
otu_taxa <- otu_taxa[otu_taxa$Genus != 'g__Lactococcus',]

# Remove Geobacillus OTU from table
#otu <- otu[!otu_taxa$Genus %in% c('f__Bacillaceae','g__Geobacillus','g__Bacillus'),]
#otu_taxa <- otu_taxa[!otu_taxa$Genus %in% c('f__Bacillaceae','g__Geobacillus','g__Bacillus'),]
otu <- otu[!otu_taxa$Species %in% c('g__Geobacillus','g__Bacillus'),]
otu_taxa <- otu_taxa[!otu_taxa$Species %in% c('g__Geobacillus','g__Bacillus'),]

# Add new Reads column to metadata
dfu$Reads <- NA
dfu$Reads[dfu$has_raw_sample] <- colSums(otu)
summary(dfu$Reads/dfu$RawReads)


# Finalize Data -----------------------------------------------------------------------------

# Make column indicating having a sample above the cutoff
dfu$has_sample <- this2that(dfu$Reads >= 1200, NA, FALSE)

# Remove "failed" reactions
otu <- otu[,(colSums(otu) >= 1200)]
# Remove OTUs with no reads after filtering
otu <- otu[rowSums(otu) != 0,]
# Remove taxa that were filtered
otu_taxa <- otu_taxa[otu_taxa$OTU %in% row.names(otu),]

# Subsample table
set.seed(1234)
rare <- t(rrarefy(t(otu), 1200))
rare <- rare[rowSums(rare) > 0,]

# Make phylum, genus, and species level tables
taxa <- lapply(nameify(c('Phylum','Family','Genus','Species')), function(rank) {
  tax <- aggregate(otu, by = otu_taxa[,rank], sum)
  row.names(tax) <- tax$Group.1
  tax$Group.1 <- NULL
  return(tax)
})


# Add microbiome fields to metadata -----------------------------------------------------------

# Define genera of interest
want <- c('g__Staphylococcus','g__Streptococcus','g__Corynebacterium','g__Anaerococcus',
          'p__Proteobacteria','p__Firmicutes','s__aureus','s__pettenkoferi')
rank <- c(rep('Genus',4),rep('Phylum',2),rep('Species',2))

# Extract proportions of wanted OTUs/genera and add to DFU table
dfu <- mapply(function(w,r) {
  norm <- apply(otu, 2, function(x) x/sum(x))  
  otus <- otu_taxa$OTU[otu_taxa[[r]] == w]
    if(length(otus) == 1) {
      norm[otus,]
    } else  {
      colSums(norm[otus,])
    }
  }, want, rank, SIMPLIFY = T) %>% 
  as.data.frame %>% 
  add_rownames('SampleID') %>%
  left_join(x = dfu, y = ., by = 'SampleID')
# remove underscores from taxa names
names(dfu) <- gsub(names(dfu), pattern = '.__', replacement = '')


# Anaerobe information ------------------------------------------------------------------------

# read in oxygen tolerance information
anaerobes <- read.csv('./oxygen_tolerance.csv') %>%
  # define strict anaerobes and lenient anaerobes
  mutate(Anaerobe = Oxygen.Tolerance == 'Anaerobe',
         Anaerobe.lenient = Oxygen.Tolerance %in% c('Anaerobe','Anaerobe/Facultative')) %>%
  # remove others
  filter(Anaerobe) %>%
  # add OTU identifier 
  inner_join(otu_taxa, by = c('Genus','Species')) %>% 
  select(Anaerobe,OTU)

# Create a column in dfu for the relative abundance of anaerobes in each sample
tmp <- colSums(otu[anaerobes$OTU,])/colSums(otu) %>% 
  data.frame(Anaerobes = .)
dfu <- tmp %>% add_rownames('SampleID') %>%
  left_join(x = dfu, y = ., by = 'SampleID')


# DMN Clusters --------------------------------------------------------------------------------

dfu <- read.csv('./dfu_dmn_assignments.csv', colClasses = c('character','factor')) %>% 
  left_join(x = dfu, y = ., by = 'SampleID')


# Alpha diversity -----------------------------------------------------------------------------

# Read in phylogenetic tree
# phylo <- read.tree('./dfu_rep_set_no_singletons.tre')
# phylo <- root(phylo, 1, resolve.root = T)

# Calculate phylogenetic distance and observed species
# alpha <- pd(samp = t(rare), tree = phylo, include.root = F) %>% 
#   add_rownames('SampleID') %>%
#   rename(Richness = SR)

# Calculate shannon diversity
alpha <- apply(rare, 2, function(x) sum(x > 0)) %>% data.frame(SampleID = names(.), Richness = .)
alpha$Shannon <- diversity(rare, MARGIN = 2, index = 'shannon')

# Merge with dfu data frame
dfu <- left_join(dfu, alpha, by = 'SampleID')

# Beta Diversity ------------------------------------------------------------------------------

# Not sure why but need to remove otus to match phylogeny
#tmp <- otu[row.names(otu) %in% phylo$tip.label,]
tmp <- otu

# Make bray-curtis distance matrix
bc <- vegdist(x = as.matrix(t(tmp)), method = 'bray') %>% as.matrix
jc <- vegdist(x = as.matrix(t(tmp)), method = 'jaccard') %>% as.matrix
#uf <- GUniFrac(otu.tab = t(tmp), tree = phylo)
#dist_mat <- list(bc = bc, jc = jc, wuf = uf$unifracs[,,'d_1'], uuf = uf$unifracs[,,'d_UW'])
dist_mat <- list(bc = bc, jc = jc)
rm(bc,jc,uf)

# Make NMDS object
beta <- ldply(dist_mat, function(d) {
  metaMDS(d,k=2)$points %>% as.data.frame %>% add_rownames('SampleID')
}, .id = 'Beta') %>% left_join(x = dfu, y = ., by = 'SampleID') 


# Relative Abundance -------------------------------------------------------------------------------

# Normalize table
#abund <- apply(taxa$Genus, MARGIN = 2, function(x) x/sum(x))
row.names(taxa$Genus) <- taxa$Genus$Genus
abund <- apply(taxa$Genus[,-1], MARGIN = 2, function(x) x/sum(x))
# Add relative abundances to mapping file and make into long format
abund <- merge(dfu, t(abund), by.x = 'SampleID', by.y = 'row.names') %>%
  gather('Taxa','Proportion',ncol(dfu)+(1:nrow(abund))) %>%
  filter(SubjectID != 'Control')

# Make table of 20 most present genera
taxa_freq <- abund %>% group_by(Taxa) %>% 
  mutate(PropNA = this2that(Proportion, 0, NA)) %>%
  summarise(num.samples = sum(Proportion > 0, na.rm = T),
            mean.all = mean(Proportion, na.rm = T),
            mean.present = mean(PropNA, na.rm = T),
            sparsity = 1 - (mean.all/mean.present),
            var.all = var(Proportion, na.rm = T),
            var.present = var(PropNA, na.rm = T)) %>%
  arrange(desc(sparsity))

# Make list of wanted taxa: >= 1% overall or 5% when present
want <- taxa_freq %>% filter(mean.present >= 0.05 | mean.all >= 0.01, num.samples >= 10)

# Make dataframe with taxa comprising more than 0.5% of all
top <- abund %>% mutate(SampleID = as.factor(SampleID)) %>%
  filter(Taxa %in% want$Taxa) %>%
  take_top_taxa('Taxa', 0.001, order_taxa = T, order_samples = T, add_other = F)


# Temporal Variation --------------------------------------------------------------------------

# Calculate subject pairwise distnaces in chronological order
# Assumes samp and id are in chronological order
sequential_dist <- function(samp, dist, id = NULL) {
  # check for empty factor levels
  if(any(table(id) == 0)) stop('Empty factor levels present in id.')
  # define comparison vectors of sample ids
  from <- head(samp,-1)
  to <- tail(samp,-1)
  # Extract pairwise distance for all sequential samples
  distances <- mapply(function(from, to) dist[from,to], from, to)
  # make results data frame
  results <- data.frame(from, to, distances, row.names = NULL)
  # if id grouping variable, put NAs at the beginning of each group
  if(!is.null(id)) {
    # make sure they are the same length
    if(length(id) != length(samp)) stop('The variables samp and id of different lengths')
    # find indices of first, remove the last element because of the n-1 issue for distances
    replace <- head(tapply(1:length(id), id, tail, 1), -1)
    results$distances[replace] <- NA
  }
  # return results
  return(results)
}

# Calculate temporal variation for each distance metric and merge with dfu
## loop through metrics, take in subsetted data frame as argument
dfu <- ldply(dist_mat, function(d, df) {
  sequential_dist(samp = df$SampleID, id = df$SubjectID, dist = d)
  ## hand subsetted data frame to ldply function
}, filter(dfu, has_sample, SubjectID != 'Control') %>% droplevels %>% arrange(SubjectID, Visit)) %>% 
  # spread out the data frame so that each distance metric gets one column
  spread('.id','distances') %>%
  # remove the from column
  select(-contains("from")) %>%
  # merge into full data frame
  left_join(x = dfu, y = ., by = c('SampleID' = 'to'))


# clean up
rm(alpha,want,rank,b,dup,anaerobes,tmp,blast)
detach("package:dplyr", unload=TRUE)
detach("package:loescher", unload=TRUE)
detach("package:GUniFrac", unload=TRUE)
detach("package:picante", unload=TRUE)
detach("package:ape", unload=TRUE)
detach("package:vegan", unload=TRUE)
detach("package:biom", unload=TRUE)
detach("package:qiimer", unload=TRUE)
detach("package:tidyr", unload=TRUE)
detach("package:ggplot2", unload=TRUE)
