# DFU Correlations and Associations
library(loescher)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

load('./blast.rdata')


# Create species columns ----------------------------------------------------------------------

# aggregate OTUs by taxa
norm <- taxa$Species %>% apply(2,function(x) x/sum(x))
want <- norm[c('s__aureus','s__pettenkoferi','g__Staphylococcus'),] %>% t %>% 
  as.data.frame %>% add_rownames()
names(want) <- c('SampleID','S.aureus','S.pettenkoferi','S.other')
dfu <- left_join(dfu, want, by = 'SampleID')


clinical <- c('hgba1c','meantisox','crp','UlcerDur','surarea','wbc','depth')
microbiome <- c('Staphylococcus','S.aureus','S.pettenkoferi','S.other','Streptococcus','Anaerobes','Anaerococcus','Corynebacterium','Proteobacteria','Firmicutes','Shannon','PD','Richness')

correlation <- filter(dfu,Visit == 0)[,c(clinical,microbiome)] %>%
  multi.cor.test(method = 'spearman', plot = F)

# Make correlation heat map of baseline samples
corr <- correlation[microbiome,clinical] %>% as.data.frame %>% 
  add_rownames('Microbiome') %>% 
  gather('Clinical','Correlation',-Microbiome)
sig <- t(correlation[clinical,microbiome]) %>% as.data.frame %>% 
  add_rownames('Microbiome') %>% 
  gather('Clinical','p.val',-Microbiome)
comb <- full_join(corr,sig, by = c('Microbiome','Clinical')) %>%
  # address p-values of zero
  mutate(p.val = this2that(p.val, 0, 0.0001))
#pdf('./correlations/heatmap_all_samples.pdf')
ggplot(comb, aes(x = Clinical, y = Microbiome, fill = Correlation)) +
  geom_tile() +
  geom_point(aes(size = -log10(p.val), alpha = p.val <= 0.05), shape = 42, color = 'white') +
  scale_fill_gradient2(low = 'red', mid = 'black', high = 'green', midpoint = 0) +
  scale_alpha_manual(values = c(0,1)) +
  scale_size_continuous(range = c(3,10))
#dev.off()

# Make combined table for ease of display
comb <- full_join(corr,sig, by = c('Microbiome','Clinical')) %>%
  rename(rho = Correlation)

# get the heatmap hierarchical clustering
tmp <- comb %>% select(-p.val) %>%
  spread(Clinical,rho)
comb %>% select(-rho) %>%
  spread(Clinical,p.val) %>% View
row.names(tmp) <- tmp$Microbiome
tmp$Microbiome <- NULL
pdf('./paper_plots/correlations_heatmap.pdf')
heatmap(as.matrix(tmp), margins = c(9,9), col = colorRampPalette(c("red", "black","green"))(100))
dev.off()


# Longitudinal Correlations -------------------------------------------------------------------

library(piecewiseSEM)
library(lme4)

# Make function to estimate p-values for parameters
get_lmer_pval <- function(obj) {
  # Extract the t-statistic for the parameter
  tStat <- data.frame(coef(summary(obj)))$t.value
  # Get a p-value assuming the normal approximation, assuming two-tailed
  pVal <- 2 * (1 - pnorm(abs(tStat)))
  # return result
  return(pVal)
}

# Make data frame of associations
long <- ldply(microbiome, function(m) {
  ldply(clinical, function(c) {
    form <- as.formula(paste(c,'~',m,'+ (',m,'|SubjectID)'))
    fit <- lmer(form, data = dfu, REML = F)
    coef <- coef(summary(fit))[2,1]
    pval <- get_lmer_pval(fit)[2]
    r2 <- sem.model.fits(fit)[c('Marginal','Conditional')]
    data.frame(clinical = c, microbiome = m, 
               coef = coef, marginal = r2[1], conditional = r2[2], pval = pval)
  })
})

# Make modified marginal to include sign of coef
long <- long %>% mutate(Marginal.Sign = ifelse(coef > 0, Marginal, Marginal * -1),
                        Conditional.Sign = ifelse(coef > 0, Conditional, Conditional * -1))

# make heatmaps
ggplot(long, aes(x = clinical, y = microbiome)) + 
  geom_tile(aes(fill = Marginal.Sign)) +
  geom_tile(aes(fill = Conditional.Sign)) +
  geom_point(aes(size = -log10(pval), alpha = pval <= 0.05), shape = 42, color = 'white') +
  scale_fill_gradient2(mid = 'black', high = 'green', low = 'red') +
  scale_alpha_manual(values = c(0,1)) +
  scale_size_continuous(range = c(3,8)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Categorical Variables -----------------------------------------------------------------------

sapply(microbiome, function(m) {
  form <- as.formula(paste(m,'~ ulcerloc_cat'))
  kruskal.test(form, data = subset(dfu, Visit == 0))$p.value
})

sapply(clinical, function(c) {
  form <- as.formula(paste(c,'~ dmn'))
  plot(form, dfu)
  #kruskal.test(form, data = subset(dfu, Visit == 0))$p.value
  kruskal.test(form, data = dfu)$p.value
})

