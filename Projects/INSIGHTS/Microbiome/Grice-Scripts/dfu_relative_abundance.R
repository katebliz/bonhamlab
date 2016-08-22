library(loescher)
library(ggplot2)
library(dplyr)
library(tidyr)

load('./blast.RData')

# Read in dmn data
dmn <- read.csv('./dfu_dmn_assignments.csv') %>% select(SampleID, dmn) %>% mutate(dmn = as.factor(dmn))
dfu <- dfu %>% select(-dmn) %>% left_join(dmn, by = 'SampleID')

# Format the taxa names for the plot
# Format Taxa names data frame to show last known rank
tmp <- otu_taxa %>% select(Genus, Species) %>%
  mutate_each(funs(gsub(pattern = '.__', replacement = '', x = .))) %>%
  mutate(Genus_species = paste(Genus, Species))
tmp$Genus_species[tmp$Genus == tmp$Species] <- tmp$Genus[tmp$Genus == tmp$Species]
otu_taxa$Genus_species <- tmp$Genus_species
rm(tmp)

# aggregate by new assignemnts
gs <- aggregate(otu, by = list(Taxa = otu_taxa[['Genus_species']]), sum)

normalize <- function(x) x/sum(x)
make_rowname <- function(x, col, as.matrix = F, drop = T) {
  row.names(x) <- x[,col]
  if(drop) x[,col] <- NULL
  if(as.matrix) x <- as.matrix(x)
  return(x)
}

# Normalize the taxa and subset to samples being used
gs <- gs %>% make_rowname('Taxa') %>% 
  apply(2, function(x) x/sum(x)) %>% t %>%
  as.data.frame %>% add_rownames('SampleID') %>%
  filter(SampleID %in% dfu$SampleID)

# make a list of taxa that contribute greater than 1%
keep <- gs %>% select(-SampleID) %>%
  summarise_each(funs(mean)) %>% 
  gather('Taxa','Mean') %>%
  filter(Mean >= 0.01)

# get species means
gs %>% gather('Taxa','Mean', -SampleID) %>%
  right_join(dfu, by = 'SampleID') %>%
  filter(!is.na(dmn)) %>%
  mutate(dmn = paste0('CT',dmn)) %>%
  group_by(dmn, Taxa) %>%
  summarise(Proportion = mean(Mean, na.rm = T)) %>% 
  filter(Taxa %in% keep$Taxa) %>% 
  ggplot(aes(x = dmn, y = Proportion, fill = Taxa)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = pick_colors()) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(ncol = 2)) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  xlab('')
  
ggsave('./paper_plots/dmn_means_plot.pdf', width = 3.5, height = 8)
