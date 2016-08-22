# pca analysis

library(vegan)
library(ggplot2)

# normalize table
norm <- apply(taxa$Species, 2, function(x) x/sum(x))

pca <- rda((norm), scale = T)

biplot(pca, scaling = -1)


filter(beta, !is.na(Beta)) %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = dmn), size = 2.5) +
  scale_color_manual(values = pick_colors()) +
  facet_wrap(~Beta)


# Beta Diversity ------------------------------------------------------------------------------

# Not sure why but need to remove otus to match phylogeny
#tmp <- otu[row.names(otu) %in% phylo$tip.label,]
tmp <- taxa$Species

# Make bray-curtis distance matrix
bc <- vegdist(x = as.matrix(t(tmp)), method = 'bray') %>% as.matrix
jc <- vegdist(x = as.matrix(t(tmp)), method = 'jaccard') %>% as.matrix
dist_mat <- list(bc = bc, jc = jc)
rm(bc,jc)

# Make NMDS object
beta.species <- ldply(dist_mat, function(d) {
  metaMDS(d,k=2)$points %>% as.data.frame %>% add_rownames('SampleID')
}, .id = 'Beta') %>% left_join(x = dfu, y = ., by = 'SampleID') 

#
beta.species$Staphylococcus <- beta.species$Staphylococcus + 0.0001
beta.species$Streptococcus <- beta.species$Streptococcus + 0.0001
beta.species$aureus <- beta.species$aureus + 0.0001
beta.species$pettenkoferi <- beta.species$pettenkoferi + 0.0001

# Make plot
filter(beta.species, !is.na(Beta), Beta == 'bc') %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = dmn), size = 2.5) +
  scale_color_manual(values = pick_colors()) +
  facet_wrap(~Beta)

filter(beta.species, !is.na(Beta), Beta == 'bc') %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = log2(Staphylococcus/Streptococcus), 
                 size = Staphylococcus+Streptococcus)) +
  facet_wrap(~Beta) +
  scale_color_gradient2(low = 'deepskyblue2', high = 'deeppink2', mid = 'black', midpoint = 0)

filter(beta.species, !is.na(Beta), Beta == 'bc') %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = log2(aureus/pettenkoferi), 
                 size = Staphylococcus)) +
  facet_wrap(~Beta) +
  scale_color_gradient2(low = 'orange', high = 'purple', mid = 'black', midpoint = 0)
  
filter(beta.species, !is.na(Beta), Beta == 'bc') %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = aureus/Staphylococcus, size = Staphylococcus)) +
  facet_wrap(~Beta) +
  scale_color_gradient(low = 'purple', high = 'orange')

filter(beta.species, !is.na(Beta), Beta == 'bc') %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = pettenkoferi/Staphylococcus, size = Staphylococcus)) +
  facet_wrap(~Beta) +
  scale_color_gradient(low = 'purple', high = 'lightseagreen')

filter(beta.species, !is.na(Beta), Beta == 'bc') %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = (Staphylococcus-pettenkoferi-aureus+0.0002)/Staphylococcus, 
                 size = Staphylococcus)) +
  facet_wrap(~Beta) +
  scale_color_gradient('Unclassified Staph',low = 'purple', high = 'limegreen')

filter(beta.species, !is.na(Beta), Beta == 'bc') %>%
  ggplot(aes(x = MDS1, y = MDS2)) + 
  geom_point(aes(color = Anaerobes, size = Anaerobes)) +
  facet_wrap(~Beta) +
  scale_color_gradient(low = 'black', high = 'red')


filter(beta.species, !is.na(Beta), Beta == 'bc', !is.na(Outcome), has_sample, !is.na(wuf),
       idx.last.visit > 2) %>%
  arrange(SubjectID,Visit) %>%
  ggplot(aes(x = MDS1, y = MDS2, group = SubjectID)) + 
  geom_path(aes(color = Visit)) +
  geom_point(aes(fill = paste(v.abx,v.comp)), shape = 21, size = 2.5) +
  scale_color_gradient(low = 'red', high = 'green') +
  facet_wrap(~SubjectID + EOSReas)


filter(beta.species, !is.na(Beta), Beta == 'bc', has_sample, !is.na(wuf), !is.na(abx)) %>%
  left_join(abx, by = c('SubjectID','Visit')) %>%
  ggplot(aes(x = MDS1, y = MDS2, group = SubjectID)) + 
  geom_point(aes(color = abx), shape = 16, size = 2.5) +
  facet_wrap(~antibiotics)
