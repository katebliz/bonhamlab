library(loescher)
library(dplyr)
library(tidyr)

# Subsample table
set.seed(1234)
taxa.rare <- t(rrarefy(t(taxa$Species), 1200))
taxa.rare <- taxa.rare[rowSums(taxa.rare) > 0,]

# Calculate shannon diversity
alpha <- apply(taxa.rare, 2, function(x) sum(x > 0)) %>% data.frame(SampleID = names(.), Richness = .)
alpha$taxaShannon <- diversity(taxa.rare, MARGIN = 2, index = 'shannon')
alpha <- rename(alpha, taxRichness = Richness)

# Merge with dfu data frame
dfu <- left_join(dfu, alpha, by = 'SampleID')


# Make plots
dfu %>% gather('Alpha','Diversity',taxaRichness,taxaShannon) %>%
  ggplot(aes(x = complications, y = Diversity)) +
    geom_boxplot(aes(fill = (Anaerobes > 0.4))) +
    facet_wrap(~Alpha, scales = 'free_y')

dfu %>% gather('Alpha','Diversity',taxaRichness,taxaShannon) %>%
  ggplot(aes(x = abx, y = Diversity)) +
  geom_boxplot() +
  facet_wrap(~Alpha, scales = 'free_y')