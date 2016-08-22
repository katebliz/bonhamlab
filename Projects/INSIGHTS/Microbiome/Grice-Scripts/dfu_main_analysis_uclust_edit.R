# Project: DFU Analysis - New
# Written by: Michael Loesche

# Set Up --------------------------------------------------------------------------------------

#source('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/dfu_load_data.R')
#load('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/data.RData')
load('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/blast.RData')
# Read in file
setwd('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/')

# Load required packages
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(loescher)
library(grid) # makes unit object work in ggplot2

# Set default ggplot theme
theme_set(theme_light())


# Alpha Diversity -----------------------------------------------------------------------------

# Make list of alpha metrics
alpha_metrics <- c('Shannon','Richness','PD')

# Baseline Samples Split by Outcome
## heals 12 weeks
g <- dfu %>% filter(Visit == 0) %>%
  gather_('Metric','Diversity',alpha_metrics) %>%
  ggplot(aes(x = heals.12wks, y = Diversity, fill = heals.12wks)) + 
  geom_boxplot(notch = T) +
  facet_wrap(~Metric, scales = 'free_y') +
  scale_fill_manual(values = c('white','gray60')) +
  theme(legend.position = 'none')
## develops complications
g2 <- dfu %>% filter(Visit == 0) %>%
  gather_('Metric','Diversity',alpha_metrics) %>%
  ggplot(aes(x = comp, y = Diversity, fill = comp)) + 
  geom_boxplot(notch = T) +
  facet_wrap(~Metric, scales = 'free_y') +
  scale_fill_manual(values = c('white','gray60')) +
  theme(legend.position = 'none')
multiplot(g,g2)
## stats
dfu %>% 
  filter(Visit == 0) %>%
  gather_('Metric','Diversity',alpha_metrics) %>%
  group_by(Metric) %>%
  summarise(heals.12wks = test.dplyr(Diversity,heals.12wks,test = 'wilcox'), 
            comp = test.dplyr(Diversity,comp,test = 'wilcox'))


# Look at baseline samples split by various metadata
want <- c('Sex','race_cat','EOSReas','ulcerloc_cat','Plate','Offloading','v.abx')
g <- lapply(want, function(x) {
      dfu %>% filter(Visit == 0) %>%
        gather_('Metric','Diversity',alpha_metrics) %>%
        ggplot(aes_string(x = x, y = 'Diversity', fill = x)) + 
          geom_boxplot() +
          facet_wrap(~Metric, scales = 'free_y') +
          scale_fill_manual(values = pick_colors()) +
          theme(legend.position = 'none',
                axis.text.x = element_text(angle = 45, hjust = 1))
     })
# multiplot 2x2 grids
lapply(1:ceiling(length(g)/4), function(x) {
  start <- 4*(x-1) + 1
  end <- x * 4
  if(x == ceiling(length(g)/4) & end > length(g)) end <- length(g)
  multiplot(plotlist = g[start:end], cols = 2)
})
## stats
sapply(alpha_metrics, function(a) {
  sapply(want, function(x) {
    df <- filter(dfu, Visit == 0, has_sample)
    kruskal.test(x = df[[a]], g = df[[x]])$p.value
  })
}) %>% as.data.frame %>% 
  add_rownames('Variable') %>% 
  gather('Alpha','pval',Shannon:PD) %>% 
  filter(pval <= 0.05)

# do the same for continuous variables
want <- c('DiabDur','UlcerDur','Age','TBPI','ABI','hgba1c','crp','glucose','esr','wbc','meantisox','surarea','depth','log_burden','Reads','RawReads')
g <- lapply(want, function(x) {
  dfu %>% filter(Visit == 0) %>%
    gather_('Metric','Diversity',alpha_metrics) %>%
    ggplot(aes_string(x = x, y = 'Diversity')) + 
    geom_point() +
    facet_wrap(~Metric, scales = 'free_y') +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45, hjust = 1))
})
# multiplot 2x2 grids
lapply(1:ceiling(length(g)/4), function(x) {
  start <- 4*(x-1) + 1
  end <- x * 4
  if(x == ceiling(length(g)/4) & end > length(g)) end <- length(g)
  multiplot(plotlist = g[start:end], cols = 2)
})
## stats
lapply(alpha_metrics, function(a) {
  sapply(want, function(x) {
    df <- filter(dfu, Visit == 0, has_sample)
    cor.test(x = df[[a]], y = df[[x]], method = 'spearman')$p.value
  }) %>% data.frame(var = names(.), pval = ., row.names = NULL, alpha = a)
}) %>% do.call(what = 'rbind') %>% filter(pval <= 0.05)

# All samples Alpha by visit
dfu %>% filter(SubjectID != 'Control', has_sample) %>% droplevels %>%
  ggplot(aes(x = Visit, y = Shannon)) + 
  geom_line() +
  geom_point(aes(color = v.comp, shape = v.abx)) +
  facet_wrap(~EOSReas+SubjectID) +
  scale_color_manual(values = pick_colors()) + 
  scale_shape_manual(values = c(16,21)) 

# Longitudinal outcomes
lapply(c('v.b4.healing','v.b4.comp','v.comp','v.b4.abx','v.abx'), function(x) {
  dfu %>% gather('Alpha','Diversity',Shannon,PD,Richness) %>%
    ggplot(aes_string(x = x, y = 'Diversity', fill = x)) +
    geom_boxplot() +
    facet_wrap(~Alpha, scales='free_y') +
    scale_fill_manual(values = pick_colors(), guide = 'none')
})
## stats
dfu %>% gather('Alpha','Diversity',Shannon,PD,Richness) %>%
  group_by(Alpha) %>%
  summarise(v.b4.healing = test.dplyr(Diversity, v.b4.healing, test = 'wilcox'),
            v.b4.comp = test.dplyr(Diversity, v.b4.comp, test = 'wilcox'),
            v.comp = test.dplyr(Diversity, v.comp, test = 'wilcox'),
            v.b4.abx = test.dplyr(Diversity, v.b4.abx, test = 'wilcox'),
            v.abx = test.dplyr(Diversity, v.abx, test = 'wilcox'))


# Beta Diversity ------------------------------------------------------------------------------

# plots
beta %>% filter(has_sample) %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = v.abx, shape = v.comp), size = 3) +
  scale_color_manual(values = pick_colors()) +
  scale_shape_manual(values = c(21,16)) +
  facet_wrap(~Beta)

# stats
sapply(dist_mat, function(d) {
  no_neg <- subset_beta(as.matrix(d), subset = SubjectID != 'Control', data = dfu)
  vc <- adonis(no_neg ~ v.comp, data = filter(dfu, has_sample, SubjectID != 'Control'))
  va <- adonis(no_neg ~ v.abx, data = filter(dfu, has_sample, SubjectID != 'Control'))
  return(c(v.comp = vc[[1]]$`Pr(>F)`[1], v.abx = va[[1]]$`Pr(>F)`[1]))
})

# colored with taxon cutoffos
beta %>% filter(has_sample) %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(alpha = 0.2, size = 2.5, color = 'black') +
  geom_point(aes(alpha = Firmicutes > 0.5, color = 'Firmicutes'), size = 2.5) +
  geom_point(aes(alpha = Streptococcus > 0.5, color = 'Streptococcus'), size = 2.5) +
  geom_point(aes(alpha = Staphylococcus > 0.5, color = 'Staphylococcus'), size = 2.5) +
  geom_point(aes(alpha = Corynebacterium > 0.5, color = 'Corynebacterium'), size = 2.5) +
  geom_point(aes(alpha = Anaerococcus > 0.5, color = 'Anaerococcus'), size = 2.5) +
  geom_point(aes(alpha = Proteobacteria > 0.5, color = 'Proteobacteria'), size = 2.5) +
  scale_alpha_manual(values = c(0,1), guide = 'none') +
  facet_wrap(~Beta) +
  scale_color_manual(name='Taxa > 0.5', values = c('Firmicutes' = 'purple4', 
    'Streptococcus' = 'blue', 'Staphylococcus' = 'red', 'Corynebacterium' = 'purple', 
    'Anaerococcus' = 'orange', 'Proteobacteria' = 'darkgreen'))


# look at metadata
pdf('./beta_metadata.pdf')
want <- c('Sex','race_cat','EOSReas','ulcerloc_cat','Plate','Offloading','v.abx','DiabDur','UlcerDur','Age','TBPI','ABI','hgba1c','crp','glucose','esr','wbc','meantisox','surarea','depth','log_burden','Reads','RawReads')
lapply(want, function(var) {
  beta %>% filter(has_sample) %>%
    ggplot(aes(x = MDS1, y = MDS2)) + 
    geom_point(aes_string(color = var), size = 2) +
    facet_wrap(~Beta)
})
dev.off()

# Anaerobes
ggplot(filter(beta, has_sample), aes(x = MDS1, y = MDS2, color = sqrt(Anaerobes))) + 
  geom_point() + 
  facet_wrap(~Beta) + 
  scale_color_gradient(low = 'red', high = 'green')


# Relative Abundance --------------------------------------------------------------------------

# plot showing different cutoff thresholds for top taxa
taxa_freq %>% filter(mean.all >= 0.005) %>%
  ggplot(aes(x = (1 - sparsity), y = mean.present)) + 
  geom_line(data=data.frame(x = seq(0,1,0.01), y = 0.005/seq(0,1,0.01)), 
            aes(x = x, y = y), color = 'gray45', linetype = 'dashed') +
  geom_line(data=data.frame(x = seq(0,1,0.01), y = 0.01/seq(0,1,0.01)), 
            aes(x = x, y = y), color = 'gray10', linetype = 'dashed') +
  geom_errorbar(aes(ymin = mean.present - var.present, 
                    ymax = mean.present + var.present), alpha = .25) +
  geom_text(aes(label = Taxa, color = mean.all >= 0.01)) +
  scale_color_manual(values = c('red3','forestgreen'), guide = F) +
  scale_y_continuous(limits = c(0,.31))

# Open the PDF
# Make a plot showing all samples
ggplot(top, aes(x = SampleID, y = Proportion, fill = Taxa, color = Taxa)) +  
  geom_bar(stat='identity') + 
  facet_grid(~ EOSReas, space = 'free_x', scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90, size = 10)) + 
  scale_x_discrete(name="") +
  scale_color_manual(values = pick_colors()) + 
  scale_fill_manual(values = pick_colors()) + 
  theme(legend.key.height=unit(2,"line"), axis.ticks = element_blank(), 
        axis.text.x = element_blank(), legend.text = element_text(size = 13))

# Take means by EOSReas at baseline
top %>% filter(Visit == 0, !is.na(Outcome)) %>%
  group_by(Taxa, EOSReas) %>% 
  summarise(Proportion = mean(Proportion,na.rm=T)) %>%
  ggplot(aes(x = EOSReas, y = Proportion, fill = Taxa)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = pick_colors()) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
## stats
top %>% group_by(Taxa) %>% 
  filter(SubjectID != 'Control', EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  summarise(pval = test.dplyr(Proportion, EOSReas, test = 'kruskal')) %>%
  mutate(padj = p.adjust(pval, method = 'fdr')) %>%
  filter(pval <= 0.05) %>% arrange(pval)
  
# Take means by EOSReas and split by baseline or not
top %>% group_by(Taxa, EOSReas, Visit != 0) %>% 
  filter(SubjectID != 'Control', EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  summarise(Proportion = mean(Proportion,na.rm=T)) %>%
  ggplot(aes(x = `Visit != 0`, y = Proportion, fill = Taxa)) + 
  geom_bar(stat = 'identity') + 
  facet_wrap(~EOSReas) +
  xlab('Baseline Visit On Left') +
  scale_fill_manual(values = pick_colors()) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
# stats
top %>% group_by(Visit == 0, Taxa) %>% 
  filter(SubjectID != 'Control', EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  summarise(pval = test.dplyr(Proportion, EOSReas, test = 'kruskal')) %>%
  mutate(padj = p.adjust(pval, method = 'fdr')) %>%
  filter(padj <= 0.05) %>% arrange(`Visit == 0`, pval)

# Make plot showing the healed breakdown
top %>% filter(Visit == 0, healed) %>%
  group_by(Outcome, Taxa) %>%
  summarise(Proportion = mean(Proportion,na.rm=T)) %>%
  ggplot(aes(x = Outcome, y = Proportion, fill = Taxa)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = pick_colors()) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
top %>% filter(Visit == 0, healed) %>%
  group_by(Taxa) %>%
  summarise(kw = test.dplyr(Proportion, Outcome, 'kruskal')) %>%
  mutate(fdr = p.adjust(kw, 'fdr')) %>%
  filter(fdr <= 0.1)

top %>% group_by(Taxa,v.comp) %>% 
  summarise(Proportion = mean(Proportion,na.rm=T)) %>%
  ggplot(aes(x = v.comp, y = Proportion, fill = Taxa)) +
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values = pick_colors())
top %>% group_by(Taxa) %>%
  summarise(wilcox = test.dplyr(Proportion, v.comp, 'wilcox')) %>%
  filter(wilcox <= 0.05)

## Make plot of major taxa changing over time
spec_taxa <- c('Streptococcus','Staphylococcus','Proteobacteria','Anaerococcus',
               'Corynebacterium','Anaerobes')
g <- lapply(spec_taxa, function(tax) {
  dfu %>% filter(has_sample, !is.na(Outcome)) %>%
    ggplot(aes_string(x = 'Visit', y = tax, group = 'SubjOrder')) +
      geom_line(size = 0.75) + 
      geom_point(size = 3, shape = 21, fill = 'white') + 
      geom_point(aes(alpha = v.comp, color = 'Comp'), size = 3) + 
      geom_point(aes(alpha = v.abx, color = 'Abx'), size = 3) + 
      geom_point(aes(alpha = v.comp & v.abx, color = 'Both'), size = 3) + 
      geom_point(aes(alpha = v.b4.healing | last.visit, color = 'Last Visit'), size = 3) + 
      scale_alpha_manual(values = c(0,.75), guide = 'none') +
      scale_color_manual(name = 'Event', values = c('Comp'='red','Abx'='blue','Both'='purple',
                                                    'Last Visit'='green')) +
      facet_grid(Outcome~., scales = 'free_x')
})
multiplot(plotlist = g, cols = 2)

# heatmap timeline
top %>% filter(has_sample, EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  arrange(idx.last.visit,SubjectID,Visit) %>%
  mutate(SampOrder = orderFactor(paste(SubjectID,Visit))) %>%
  ggplot(aes(x = SampOrder, y = Taxa)) + 
  geom_tile(aes(fill = sqrt(Proportion))) +
  facet_grid(~EOSReas, scales = 'free_x', space = 'free_x') +
  theme(panel.margin = unit(.5, 'lines'),
        axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position='none')

# Mega timeline
# Make annotation data frame
annot <- mutate(dfu,Taxa=NA)
# Make plot of each patient's timeline
g <- ggplot(top, aes(x = Visit, y = Proportion, fill = Taxa)) + 
  geom_bar(stat = 'identity') + 
  geom_point(data = annot, aes(x = Visit, y = 1, alpha = v.comp, shape = v.wd), size = 2) +
  geom_point(data = annot, aes(x = Visit, y = 1.15, alpha = v.abx),shape=18,color='red',size=2) +
  facet_grid(EOSReas + SubjOrder~., scales = 'free') +
  scale_fill_manual(values = pick_colors()) +
  scale_alpha_manual(values = c(0,1)) +
  scale_shape_manual(values = c(21,16)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1.25)) +
  guides(fill = guide_legend(reverse = TRUE))

ggsave(g, filename = './mega_timeline.pdf', width = 9, height = 120, limitsize = F)



# Temporal Stability Metrics ------------------------------------------------------------------

# Make heat map plot showing the change over time
dfu %>% filter(EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  filter(!(healed & last.visit), Visit != 0) %>%
  gather('Beta','Distance',bc:wuf) %>%
  ggplot(aes(x = Visit, y = SubjOrder)) + 
  geom_tile(aes(fill = Distance), color = 'white') +
  geom_point(aes(alpha = v.comp | v.abx, color = v.abx)) +
  facet_grid(EOSReas~Beta, scales = "free_y", space = 'free_y') +
  scale_alpha_manual(values = c(0,1)) +
  scale_color_manual(values = c('green','red'))

# Make a line plot showing the change over time
dfu %>% filter(EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  filter(!(healed & last.visit), Visit != 0, has_sample) %>%
  gather('Beta','Distance',bc:wuf) %>%
  ggplot(aes(x = Visit, y = Distance, group = SubjOrder)) + 
  geom_line() +
  geom_point(shape = 21) + 
  geom_point(aes(alpha = v.comp | v.abx, color = v.abx)) +
  facet_grid(EOSReas~Beta, scales = "free_y", space = 'free_y') +
  scale_alpha_manual(values = c(0,1)) +
  scale_color_manual(values = c('green','red'))

# Boxplot of EOSReas distances
dfu %>% filter(EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  gather('Beta','Distance',bc:wuf) %>%
  ggplot(aes(x = EOSReas, y = Distance)) + 
  geom_boxplot(notch = T) +
  facet_wrap(~Beta, scales = 'free_y')
#
dfu %>% filter(EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  ungroup %>%
  mutate(Visit = as.factor(Visit)) %>%
  ggplot(aes(x = Visit, y = wuf, fill = EOSReas)) +
  geom_boxplot(notch = F) +
  facet_wrap(~EOSReas)
#
dfu %>% filter(healed, has_sample) %>%
  ggplot(aes(x = Visit, y = wuf, group = SubjectID)) + 
  geom_line() +
  geom_point(aes(size = sqrt(surarea), color = v.abx, shape = v.comp), fill = 'white') +
  scale_shape_manual(values = c(21, 16)) +
  facet_wrap(~idx.last.visit)

# Tests
dfu %>% gather('Beta','Distance',bc:wuf) %>%
  group_by(Beta) %>%
  mutate(VisitType = paste(v.comp, v.abx, sep = ':')) %>%
  summarise(healed = test.dplyr(Distance, healed, 'wilcox'),
            comp = test.dplyr(Distance, comp, 'wilcox'),
            abx = test.dplyr(Distance, antibiotics, 'wilcox'),
            v.abx = test.dplyr(Distance, v.abx, 'wilcox'),
            v.comp = test.dplyr(Distance, v.comp, 'wilcox'),
            VisitType = test.dplyr(Distance, VisitType, 'kruskal'))
#
dfu %>% filter(!is.na(v.comp), !is.na(v.abx)) %>%
  mutate(VisitType = paste(v.comp, v.abx, sep = ':')) %>%
  with(pairwise.wilcox.test(x = wuf, g = VisitType))
# plot findings
dfu %>% gather('Beta','Distance',bc:wuf) %>%
  group_by(Beta) %>%
  mutate(VisitType = paste(v.comp, v.abx, sep = ':')) %>%
  ggplot(aes(x = VisitType, y = Distance, fill = VisitType)) +
  geom_boxplot() +
  facet_wrap(~Beta, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  xlab('Antibiotics:Complication')


