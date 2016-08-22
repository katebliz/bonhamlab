# Explore contaminants ------------------------------------------------------------------------

# Numer of total reads and total otus
sum(otu) #10,095,728 reads
nrow(otu) #11,321

## Geobacillus vs Bacillaceae
mutate(otu_taxa, Reads = rowSums(otu)) %>%
  filter(Family == 'f__Bacillaceae') %>% 
  group_by(Species) %>% 
  summarise(nOTU = length(Species), Genus = first(Genus), Reads = sum(Reads)) %>%
  select(Genus, Species, nOTU, Reads)

# Distribution of Bacillaceae OTUs
tmp <- otu[otu_taxa$Genus == 'f__Bacillaceae',] %>% 
  rowSums %>% 
  sort(decreasing = T) %>% 
  data.frame(OTU = names(.), Reads = ., row.names = NULL, Index = 1:length(.)) 
head(tmp)
# Linear plot
g <- ggplot(tmp, aes(x = Index, y = Reads)) + geom_point() + 
  geom_line() + ggtitle('Genus: f__Bacillaceae')
# Logarithmic plot
g2 <- ggplot(tmp, aes(x = Index, y = Reads)) + geom_point() + 
  geom_line() + ggtitle('Genus: f__Bacillaceae') +
  scale_y_log10(breaks = 10^(1:5))
multiplot(g,g2, cols = 2)
# Distribution of Geobacillus OTUs
tmp <- otu[otu_taxa$Genus == 'g__Geobacillus',] %>% 
  rowSums %>% 
  sort(decreasing = T) %>% 
  data.frame(OTU = names(.), Reads = ., row.names = NULL, Index = 1:length(.)) 
head(tmp)
# Linear plot
g <- ggplot(tmp, aes(x = Index, y = Reads)) + geom_point() + 
  geom_line() + ggtitle('Genus: g__Geobacillus')
# Logarithmic plot
g2 <- ggplot(tmp, aes(x = Index, y = Reads)) + geom_point() + 
  geom_line() + ggtitle('Genus: g__Geobacillus') +
  scale_y_log10(breaks = 10^(1:5))
multiplot(g,g2, cols = 2)


## Lactococcus vs Streptococcacaeae
mutate(otu_taxa, Reads = rowSums(otu)) %>%
  filter(Family == 'f__Streptococcaceae') %>% 
  group_by(Species) %>% 
  summarise(nOTU = length(Species), Genus = first(Genus), Reads = sum(Reads)) %>%
  select(Genus, Species, nOTU, Reads)


## Correlations of contaminants and total reads
tmp <- mutate(otu, Genus = otu_taxa$Genus, Species = otu_taxa$Species) %>%
  filter(Genus %in% c('f__Bacillaceae','g__Geobacillus','g__Streptococcus','g__Lactococcus')) %>%
  group_by(Species) %>% 
  select(-Genus) %>%
  summarise_each(funs(sum))
row.names(tmp) <- tmp$Species
tmp$Species <- NULL
tmp <- as.data.frame(t(tmp))
tmp$Total <- rowSums(tmp)
tmp <- tmp %>% add_rownames('SampleID') %>% gather('Taxa','Reads',-SampleID, -Total)
# linear scale plot
g <- ggplot(tmp, aes(x = Total, y = Reads, color = Taxa, group = Taxa)) + 
  geom_point() +
  scale_color_manual(values = pick_colors(20)) +
  ylab('Taxa Reads') + xlab('Total Reads') +
  guides(col = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(legend.position = 'none')
# log scale plot
g2 <- ggplot(tmp, aes(x = Total, y = Reads, color = Taxa, group = Taxa)) + 
  geom_point() +
  scale_y_log10() + scale_x_log10() +
  scale_color_manual(values = pick_colors(20)) +
  ylab('Taxa Reads (log10)') + xlab('Total Reads (log10)') +
  guides(col = guide_legend(nrow = 3, byrow = TRUE)) +
  theme(legend.position = 'bottom')
multiplot(g,g2, cols = 1)

# Negative Controls ---------------------------------------------------------------------------

# Isolate the control swabs
neg <- otu[,dfu$SampleID[dfu$SubjectID == 'Control']]
neg <- neg[rowSums(neg) > 0,]

# get sample depth
colSums(neg)

# Add taxa information
neg <- merge(neg, otu_taxa[,c('OTU','Genus','Species')], by.x = 'row.names', by.y = 'OTU') %>%
  rename(OTU = Row.names)

# sort by otu counts
arrange(neg, desc(rowSums(neg[,2:3]))) %>% head(20)


# Cutoff Thresholds ---------------------------------------------------------------------------

## Assess distribution of reads per sample
# plot showing number of reads per sample with lines depicting potential cutoffs
tmp <- data.frame(Reads = sort(colSums(otu)), Samples = (1:ncol(otu))/ncol(otu))
cutoffs <- c(500,1000,1500,2000,2500,3000,3500,4000,4500,5000)
g <- ggplot(tmp, aes(x = Samples, y = Reads)) + 
  geom_line() +
  scale_x_continuous(breaks = seq(0,1,0.05)) +
  scale_y_sqrt()
# zoomed in on potential cutoff range
g2 <- ggplot(tmp[1:100,], aes(x = Samples, y = Reads)) + 
  geom_hline(yintercept = cutoffs, color = 'red', linetype = 'dashed') +
  geom_text(data = data.frame(cutoffs), aes(label = cutoffs, y = cutoffs), x = .01) +
  geom_line() +
  scale_x_continuous(breaks = seq(0,1,0.05)) +
  scale_y_sqrt()
multiplot(g,g2,cols = 2)

# Make plot to see where we would lose samples at different thresholds
ggplot(subset(dfu, has_raw_sample), aes(x = Visit, y = SubjectID)) + 
  geom_tile(fill = 'red', color = 'white') +
  geom_tile(data = subset(dfu, Reads >= 1500), color = 'white')

# Rarefaction Analysis ------------------------------------------------------------------------

# define cutoffs
cutoffs <- c(seq(100,499,100), seq(500,1499,250), seq(1500,4000,500))

# make a bunch of subsets of varying depths and calculate alpha diversity
rare <- lapply(1:5, function(i) {
  lapply(cutoffs, function(x){
    # Remove samples below cutoff
    sub <- otu[,(colSums(otu) >= x)]
    # Subsample table
    rare <- t(rrarefy(x = t(sub), sample = x))
    # get richness
    richness <- specnumber(rare,MARGIN = 2)
    # get shannon
    shannon <- diversity(rare, MARGIN = 2, index = 'shannon')
    # return results
    data.frame(richness,shannon) %>% 
      gather('Index','Diversity',richness:shannon) %>%
      mutate(Depth = as.numeric(x))
  }) %>% do.call(what = 'rbind')
}) %>% do.call(what = 'rbind') %>%
  group_by(Index,Depth) %>%
  summarise(se = se(Diversity,na.rm=T), Diversity = mean(Diversity,na.rm=T))

# plot results
ggplot(rare, aes(x = Depth, y = Diversity, group = Index)) +
  geom_ribbon(aes(ymin = Diversity - se, ymax = Diversity + se, fill = Index), alpha = 0.6) +
  geom_line(aes(color = Index)) +
  geom_point(aes(color = Index)) +
  facet_wrap(~Index, scales = 'free_y') +
  scale_color_manual(values = pick_colors(20)) +
  scale_fill_manual(values = pick_colors(20)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')



