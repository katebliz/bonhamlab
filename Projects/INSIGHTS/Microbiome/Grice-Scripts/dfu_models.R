# Set Up for Modeling -------------------------------------------------------------------------
library(loescher)
library(dplyr)
library(lme4)
library(Gmisc)

# # Create a column for staph aureus specifically
# dfu <- taxa$Species %>% apply(2, function(x) x/sum(x)) %>% 
#   .['s__aureus',] %>% 
#   data.frame(SampleID = names(.), S.aureus = .) %>% 
#   right_join(dfu, by = 'SampleID')
# # s. pettenkoferi
# dfu <- taxa$Species %>% apply(2, function(x) x/sum(x)) %>% 
#   .['s__pettenkoferi',] %>% 
#   data.frame(SampleID = names(.), S.pettenkoferi = .) %>% 
#   right_join(dfu, by = 'SampleID')


# Open PDF
pdf('./outcome_models/forest_plots.pdf', width = 8, height = 3.5)

# Make vector of variable factors for analysis
parameters <- c('Anaerobes','Proteobacteria','Staphylococcus','Streptococcus','S.aureus',
                'Corynebacterium','Shannon','Richness','PD')

# Make vector of variable names that you don't want in your odds tables
remove <- c('(Intercept)','UlcerDur','depth','meantisox','ulcerloc_catMidfoot','ulcerloc_catHindfoot','ulcerloc_catForefoot','wbc','hgba1c','crp','Visit')


# Baseline: healing ---------------------------------------------------
# Baseline univariate analysis: heals.12wks
m.unadj <- 
  lapply(parameters, FUN = function(x) {
    form <- as.formula(paste('heals.12wks ~',x))
    glm(form, family = 'binomial', data = dfu, subset = Visit == 0 & has_sample)
  })
# Make odds ratio table
or.unadj <- make_odds_table(m.unadj)

# Baseline multivariate analysis: heals.12wks
m.adj <- 
  lapply(parameters, FUN = function(x) {
    glm(substitute(heals.12wks ~ i + UlcerDur + depth + meantisox + 
                     ulcerloc_cat + wbc + hgba1c + crp, list(i=as.name(x))), family = 'binomial',
        data = dfu, subset = Visit == 0)
  })
# Make odds-ratio table
or.adj <- make_odds_table(m.adj, remove)

# Plot odds
or.unadj$adj <- 'Unadjusted'; or.adj$adj <- 'Adjusted'
or.comb <- rbind(or.unadj, or.adj)
or.comb$adj <- factor(or.comb$adj, levels = c('Unadjusted','Adjusted'))
plot_odds(or.comb,facet = '~adj', title = 'Healed in 12 Weeks')
# Write table of odds
or.comb$outcome <- 'Healed by 12 Weeks'
#write.table(or.comb,'./outcome_models/OR_table.txt',sep='\t',append=F,row.names=F,quote=F,col.names=T)


# Baseline: complications -------------------------------------------------
# Baseline single variable analysis: overall_comp
m.unadj <- 
  lapply(parameters, FUN = function(x) {
    glm(substitute(comp ~ i, list(i=as.name(x))), family = 'binomial',
        data = dfu, subset = Visit == 0 & has_sample)
  })
# Make odds ratio table
or.unadj <- make_odds_table(m.unadj)

# Baseline multivariate analysis: overall_comp
m.adj <- 
  lapply(parameters, FUN = function(x) {
    glm(substitute(comp ~ i + UlcerDur + depth + meantisox + ulcerloc_cat + wbc + 
                     hgba1c + crp, list(i=as.name(x))), family = 'binomial', data = dfu, 
        subset = Visit == 0 & has_sample)
  })
# Make odds ratio table
or.adj <- make_odds_table(m.adj, remove)

# Plot odds
or.unadj$adj <- 'Unadjusted'; or.adj$adj <- 'Adjusted'
or.comb <- rbind(or.unadj, or.adj)
or.comb$adj <- factor(or.comb$adj, levels = c('Unadjusted','Adjusted'))
plot_odds(or.comb,facet = '~adj', title = 'Complications Overall')
# Write table of odds
or.comb$outcome <- 'Overall Complications'
write.table(or.comb,'./outcome_models/OR_table.txt',sep='\t',append=T,row.names=F,quote=F,col.names=F)


# Longitudinal: Healing ---------------------------------------------------

# Add mean.dist to remove group
parameters <- c(parameters, 'wuf')
#remove <- c(remove, 'mean.dist')

# Longitudinal single variable analysis: v.b4.healing ~ x + Visit + (SubjectID)
m.unadj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.b4.healing ~ i + Visit + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.unadj <- make_odds_table(m.unadj, remove)

# Longitudinal multivariate analysis: v.b4.healing ~ x + ... + Visit
m.adj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.b4.healing ~ i + UlcerDur + depth + meantisox  + 
                       ulcerloc_cat + wbc + hgba1c + crp + Visit + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.adj <- make_odds_table(m.adj, remove)

# Plot odds
or.unadj$adj <- 'Unadjusted'; or.adj$adj <- 'Adjusted'
or.comb <- rbind(or.unadj, or.adj)
or.comb$adj <- factor(or.comb$adj, levels = c('Unadjusted','Adjusted'))
plot_odds(or.comb,facet = '~adj', title = 'Visit Before Healing')
# Write table of odds
or.comb$outcome <- 'Visit Before Healing'
write.table(or.comb,'./outcome_models/OR_table.txt',sep='\t',append=T,row.names=F,quote=F,col.names=F)


# Longitudinal: Complication Exclusive with NA's --------------------------------

# Longitudinal multivariate analysis: v.b4.comp ~ Visit + x + (SubjectID)
m.unadj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.b4.comp.na ~ Visit + i + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.unadj <- make_odds_table(m.unadj, remove)

# Longitudinal multivariate analysis: v.b4.comp ~ Visit + x + ... + (SubjectID)
m.adj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.b4.comp.na ~ i + UlcerDur + depth + meantisox + ulcerloc_cat + 
                       wbc + hgba1c + crp + Visit + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.adj <- make_odds_table(m.adj, remove)

# Plot odds
or.unadj$adj <- 'Unadjusted'; or.adj$adj <- 'Adjusted'
or.comb <- rbind(or.unadj, or.adj)
or.comb$adj <- factor(or.comb$adj, levels = c('Unadjusted','Adjusted'))
plot_odds(or.comb,facet = '~adj', title = 'Visit Before Complication: Exclusive with NAs')
# Write table of odds
or.comb$outcome <- 'Visit Before Complication: Exclusive with NAs'
write.table(or.comb,'./outcome_models/OR_table.txt',sep='\t',append=T,row.names=F,quote=F,col.names=F)


# Longitudinal: Wound Deterioration ---------------------------------------

# Longitudinal multivariate analysis: v.wd ~ Visit + x + (SubjectID)
m.unadj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.wd ~ Visit + i + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.unadj <- make_odds_table(m.unadj, remove)

# Longitudinal multivariate analysis: v.wd ~ Visit + x + ... + (SubjectID)
m.adj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.wd ~ i + UlcerDur + depth + meantisox + ulcerloc_cat + 
                       wbc + hgba1c + crp + Visit + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.adj <- make_odds_table(m.adj, remove)

# Plot odds
or.unadj$adj <- 'Unadjusted'; or.adj$adj <- 'Adjusted'
or.comb <- rbind(or.unadj, or.adj)
or.comb$adj <- factor(or.comb$adj, levels = c('Unadjusted','Adjusted'))
plot_odds(or.comb,facet = '~adj', title = 'Visit of Wound Deterioration')
# Write table of odds
or.comb$outcome <- 'Visit of Wound Deterioration'
write.table(or.comb,'./outcome_models/OR_table.txt',sep='\t',append=T,row.names=F,quote=F,col.names=F)


# Longitudinal: Antibiotics -------------------------------------------------------------------

# Longitudinal multivariate analysis: v.abx ~ Visit + x + (SubjectID)
m.unadj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.abx ~ Visit + i + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.unadj <- make_odds_table(m.unadj, remove)

# Longitudinal multivariate analysis: v.abx ~ Visit + x + ... + (SubjectID)
m.adj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.abx ~ i + UlcerDur + depth + meantisox + ulcerloc_cat + 
                       wbc + hgba1c + crp + Visit + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.adj <- make_odds_table(m.adj, remove)

# Plot odds
or.unadj$adj <- 'Unadjusted'; or.adj$adj <- 'Adjusted'
or.comb <- rbind(or.unadj, or.adj)
or.comb$adj <- factor(or.comb$adj, levels = c('Unadjusted','Adjusted'))
plot_odds(or.comb,facet = '~adj', title = 'Visit of Antibiotics')
# Write table of odds
or.comb$outcome <- 'Visit of Wound Deterioration'
write.table(or.comb,'./outcome_models/OR_table.txt',sep='\t',append=T,row.names=F,quote=F,col.names=F)


# Longitudinal: Wound Deterioration and Antibiotics ---------------------------------------

dfu <- dfu %>% mutate(v.wd.abx = v.wd | v.abx)

# Longitudinal multivariate analysis: v.wd ~ Visit + x + (SubjectID)
m.unadj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.wd.abx ~ Visit + i + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.unadj <- make_odds_table(m.unadj, remove)

# Longitudinal multivariate analysis: v.wd ~ Visit + x + ... + (SubjectID)
m.adj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.wd.abx ~ i + UlcerDur + depth + meantisox + ulcerloc_cat + 
                       wbc + hgba1c + crp + Visit + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.adj <- make_odds_table(m.adj, remove)

# Plot odds
or.unadj$adj <- 'Unadjusted'; or.adj$adj <- 'Adjusted'
or.comb <- rbind(or.unadj, or.adj)
or.comb$adj <- factor(or.comb$adj, levels = c('Unadjusted','Adjusted'))
plot_odds(or.comb,facet = '~adj', title = 'Visit of WD or Abx')
# Write table of odds
or.comb$outcome <- 'Visit of Wound Deterioration'
write.table(or.comb,'./outcome_models/OR_table.txt',sep='\t',append=T,row.names=F,quote=F,col.names=F)


# Longitudinal: Visit before Wound Deterioration ---------------------------------------

# Longitudinal multivariate analysis: v.wd ~ Visit + x + (SubjectID)
m.unadj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.b4.wd ~ Visit + i + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.unadj <- make_odds_table(m.unadj, remove)

# Longitudinal multivariate analysis: v.wd ~ Visit + x + ... + (SubjectID)
m.adj <- 
  lapply(parameters, FUN = function(x) {
    glmer(substitute(v.b4.wd ~ i + UlcerDur + depth + meantisox + ulcerloc_cat + 
                       wbc + hgba1c + crp + Visit + (1|SubjectID), list(i=as.name(x))), 
          family = binomial, data = dfu)
  })
# Make odds-ratio table
or.adj <- make_odds_table(m.adj, remove)

# Plot odds
or.unadj$adj <- 'Unadjusted'; or.adj$adj <- 'Adjusted'
or.comb <- rbind(or.unadj, or.adj)
or.comb$adj <- factor(or.comb$adj, levels = c('Unadjusted','Adjusted'))
plot_odds(or.comb,facet = '~adj', title = 'Visit Before Wound Deterioration')
# Write table of odds
or.comb$outcome <- 'Visit Before Wound Deterioration'
write.table(or.comb,'./outcome_models/OR_table.txt',sep='\t',append=T,row.names=F,quote=F,col.names=F)

# Close PDF
dev.off()