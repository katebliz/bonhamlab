# Set Up --------------------------------------------------------------------------------------

# Read in file
#source('~/Coding/Club_Grice/scripts/loesche/dfu_uclust/dfu_load_data.R')
load('./blast.RData')
# set ggplot theme
library(ggplot2)
library(dplyr)
library(tidyr)
library(loescher)
theme_set(theme_light())


# Perform Clustering --------------------------------------------------------------------------

library(DirichletMultinomial)
library(lattice) #visualization
library(parallel) #to use multiple cores with cross-validation

# # Subset to 100 most abundant OTUs
# sub <- otu[order(rowSums(otu), decreasing = T)[1:100],] %>% as.matrix %>% t
# sub <- sub[rownames(sub) %in% dfu$SampleID,]

# Collapse at species level
sub <- aggregate(otu, by = list(otu_taxa$Species), sum)
row.names(sub) <- sub$Group.1
sub$Group.1 <- NULL
# remove empty taxa if any
sub <- sub[rowSums(sub) > 0,]
# remove taxa not present in at least 10% of samples
sub <- sub[rowSums(sub > 0) >= floor(nrow(sub) * 0.1),]

# get idea of taxa distribution
cnts <- log10(rowSums(sub))
densityplot(cnts, xlim=range(cnts),xlab="Taxon representation (log 10 count)")

# Do Dirichlet-Multinomial clustering with parallelized lapply
fit <- mclapply(c(1:6), dmn, count = t(sub), verbose = TRUE, seed = 1234)
#fit2 <- mclapply(c(7:8), dmn, count = t(sub), verbose = TRUE, seed = 1234)
#fit <- c(fit,fit2)
# Extract fit measure from models
lplc <- sapply(fit, laplace)
pdf('./paper_plots/laplace_plot.pdf')
plot(lplc, type = 'b', xlab = 'Dirichlet Components',ylab='Laplace Estimate')
dev.off()
# Pick best model
best <- fit[[which.min(lplc)]]
# View pi and theta (homogeneity) values
mixturewt(best)
# View Dirichlet component estimates
#mixture(best)
# Visualize taxa contributions to each component
#splom(log(fitted(best)))
# Compare 1 component model to best fit
p0 <- fitted(fit[[1]], scale = TRUE)
pbest <- fitted(best, scale = TRUE)
colnames(pbest) <- paste('m',1:ncol(pbest),sep='')
(meandiff <- colSums(abs(pbest - as.vector(p0))))
#sum(meandiff)
# Look at top 10 discriminating taxa from the mean
diff <- rowSums(abs(pbest - as.vector(p0)))
o <- order(diff, decreasing=TRUE)
cdiff <- cumsum(diff[o]) / sum(diff)
(df <- head(cbind(Mean=p0[o], pbest[o,], diff=diff[o], cdiff), 10))

# Add to DFU data frame
cluster <- sapply(1:nrow(best@group), FUN=function(x) which.max(best@group[x,]))
dfu$dmn <- NA
dfu$dmn[match(colnames(sub),dfu$SampleID)] <- cluster
dfu$dmn <- as.factor(dfu$dmn)



# Visualize Sample Means ----------------------------------------------------------------------

top$dmn <- NULL

top <- ungroup(dfu) %>%
  select(SampleID, dmn) %>% 
  right_join(ungroup(top), by = 'SampleID')

top %>% filter(has_sample) %>%
  group_by(dmn, Taxa) %>%
  summarise(Proportion = mean(Proportion,na.rm=T)) %>%
  ggplot(aes(x = dmn, y = Proportion, fill = Taxa)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values = pick_colors())


# Visualize Timeline ---------------------------------------------------------------------------

# Make plot
dir.create('./dmn/',F)
pdf('./dmn/dmn_over_time.pdf', height = 10, width = 8)
# Make plot of 
dfu %>% filter(!(healed & idx.last.visit == Visit)) %>%
  filter((healed & idx.last.visit > 1) | (!healed & idx.last.visit > 0)) %>%
  filter(EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  ggplot(aes(x = Visit, y = SubjOrder, group = SubjOrder)) +
  geom_line(linetype='dashed',color = 'darkgray') +
  geom_point(aes(color = dmn, shape = has_sample), size = 3, fill = 'white') +
  facet_grid(EOSReas~.,scales = 'free',space='free') +
  scale_shape_manual(values = c(21,16)) +
  xlab('Visit Number') + ylab('Subject Number') +
  ggtitle('Dirichlet Cluster Over Time')
# As tiles
dfu %>% filter(!(healed & idx.last.visit == Visit)) %>%
  filter((healed & idx.last.visit > 1) | (!healed & idx.last.visit > 0)) %>%
  filter(EOSReas %in% c('Healed','Unhealed','Amputation')) %>%
  ggplot(aes(x = Visit, y = SubjOrder, group = SubjOrder)) +
  geom_vline(xintercept = c(5,10), linetype='dashed', color = 'darkgray') +
  geom_tile(aes(fill = dmn), color = 'white') +
  geom_point(aes(alpha = v.abx), shape = 18, size = 3) +
  facet_grid(EOSReas~.,scales = 'free',space='free') +
  scale_alpha_manual(values = c(0,1)) +
  scale_x_continuous(expand = c(0, 0)) +
  xlab('Visit Number') + ylab('Subject Number') +
  ggtitle('Dirichlet Cluster Over Time')
dev.off()


# DMN Run Length Analysis -------------------------------------------------------------------------

dfu.raw <- dfu
dfu <- filter(dfu.raw, SampleID %in% colnames(sub), SubjectID != 'Control') %>% droplevels

# Make list of dmn runs by patient id
num.changes <- 
  sapply(unique(dfu$SubjectID)[-18], simplify = F, function(x) {
    # Subset data frame
    df <- subset(dfu, SubjectID == x)
    
    # Make new dmn column filled with NAs
    d <- rep(NA, df$idx.last.visit[1]+1)
    
    # Replace NAs with dmn column from visits with samples
    have <- which(0:df$idx.last.visit[1] %in% df$Visit)
    d[have] <- df$dmn
    
    # Initialize variables
    last <- df$idx.last.visit[1] + 1
    runs <- numeric()
    run <- 0
    i <- 1
    
    ## Count runs of consecutive dmn assignments
    # Iterate over all visits
    for(i in 1:last) {
      # If NA append NA to runs and reset
      if(is.na(d[i])) {
        runs <- c(runs, NA)
        run <- 0
      } else { # If not NA then...
        # Add 1 to the run count
        run <- run + 1
        # If next visit is NA...
        if(is.na(d[i+1])){
          # Break runs, resets run above
          runs <- c(runs,run)
          # If not NA, and two visits not same break runs and reset
        } else if(d[i] != d[i+1]) {
          runs <- c(runs,run)
          run <- 0
        }
        # If not NA, and visits same, no action required. Adds 1 to count when loops.
      }
    }
    # Return results
    return(runs)
    #Close lapply loop
  })

# Take overall mean and median
mean(unlist(num.changes),na.rm = T)
median(unlist(num.changes),na.rm = T)

# Take means by EOSReas
split.ids <- tapply(dfu$SubjectID, dfu$EOSReas, unique)
sapply(split.ids, function(x) mean(unlist(num.changes[x]),na.rm=T)) %>% as.data.frame
sapply(split.ids, function(x) median(unlist(num.changes[x]),na.rm=T))

# Take means by healing in 12 weeks or not
split.ids <- tapply(dfu$SubjectID, dfu$heals.12wks, unique)
sapply(split.ids, function(x) mean(unlist(num.changes[x]),na.rm=T))
sapply(split.ids, function(x) median(unlist(num.changes[x]),na.rm=T))
wilcox.test(x = unlist(num.changes[split.ids[[1]]]), y = unlist(num.changes[split.ids[[2]]]),na.action = na.remove)

# Take means by overall complication or not
split.ids <- tapply(dfu$SubjectID, dfu$comp, unique)
sapply(split.ids, function(x) mean(unlist(num.changes[x]),na.rm=T))
sapply(split.ids, function(x) median(unlist(num.changes[x]),na.rm=T))
wilcox.test(x = unlist(num.changes[split.ids[[1]]]), y = unlist(num.changes[split.ids[[2]]]),na.action = na.remove)

# Take means by overall complication or not
split.ids <- tapply(dfu$SubjectID, dfu$antibiotics, unique)
sapply(split.ids, function(x) mean(unlist(num.changes[x]),na.rm=T))
sapply(split.ids, function(x) median(unlist(num.changes[x]),na.rm=T))
wilcox.test(x = unlist(num.changes[split.ids[[1]]]), y = unlist(num.changes[split.ids[[2]]]),na.action = na.remove)
  

# Transitions by Visit ------------------------------------------------------------------------

# Make a wide data frame of the dmns for each visit
dmn <- dfu %>% filter(has_sample, !is.na(Outcome)) %>%
#dmn <- dfu %>% filter(has_sample, !is.na(Outcome), healed, Visit <= 7) %>%
  # Make sure data frame is in order.
  arrange(SubjectID, Visit) %>%
  # Columns can't be named after numbers
  mutate(Visit = paste0('v',formatC(Visit, width=2, flag="0"))) %>%
  select(SubjectID, Visit, dmn, EOSReas) %>%
  spread(key = Visit, value = dmn)

# Make a plot
mapply(function(a,b) a == b, select(dmn, v00:v12), select(dmn, v01:v13)) %>% 
  cbind(select(dmn, EOSReas)) %>%
  group_by(EOSReas) %>%
  summarise_each(funs(mean(.,na.rm=T))) %>%
  gather("Visit", "TransFreq", v00:v12) %>%
  mutate(TransFreq = 1-TransFreq) %>%
  ggplot(aes(x = Visit, y = TransFreq, group = EOSReas)) + 
  geom_line(position = position_dodge(width = .5)) + 
  geom_point(aes(color = EOSReas), size = 3, position = position_dodge(width = 0.5))

# Proportion of visits that are 
dfu %>% group_by(EOSReas,dmn) %>% 
  filter(!is.na(dmn), !is.na(Outcome)) %>% 
  summarise(Count = length(dmn)) %>% 
  mutate(norm = Count/sum(Count)) %>% 
  ggplot(aes(x = dmn, y = EOSReas, fill = norm)) + 
  geom_tile()


# Associations with outcomes ------------------------------------------------------------------

dfu %>% filter(!is.na(Outcome)) %>% 
  droplevels %>% 
  with(table(dmn, EOSReas)) %>% 
  fisher.test(simulate.p.value = T, B = 1000000)

dfu %>% filter(!is.na(Outcome)) %>% 
  droplevels %>% 
  with(table(dmn, EOSReas)) %>% 
  as.data.frame %>% 
  group_by(EOSReas) %>%
  mutate(Proportion = Freq/sum(Freq)) %>%
  ggplot(aes(x = dmn, y = EOSReas, fill = Proportion)) +
  geom_tile() +
  scale_fill_gradient(low = 'black', high = 'skyblue')
  #scale_fill_gradient2(low = 'red', mid = 'black', high = 'green', midpoint = 0.25)



# Species medians -----------------------------------------------------------------------------

# this is showing the median
taxa$Species %>% apply(2, function(x) x/sum(x)) %>% t %>%
  as.data.frame %>% add_rownames('SampleID') %>% 
  right_join(select(dfu,SampleID,dmn), by = 'SampleID') %>%
  filter(!is.na(dmn)) %>%
  select(-SubjectID) %>%
  gather('Taxa','Proportion',-dmn,-SampleID) %>%
  mutate(dmn = paste0('ct',dmn)) %>%
  group_by(dmn, Taxa) %>%
  summarise(Proportion = mean(Proportion, na.rm = T)) %>% 
  spread(dmn, Proportion) %>% 
  mutate(sum = ct1 + ct2 + ct3 + ct4) %>%
  arrange(desc(sum)) %>%
  write.csv('./cluster_means.csv', row.names = F)

# this is showing the medians for genera strep and staph
taxa$Genus %>% apply(2, function(x) x/sum(x)) %>% t %>%
  as.data.frame %>% add_rownames('SampleID') %>% 
  right_join(select(dfu,SampleID,dmn), by = 'SampleID') %>%
  filter(!is.na(dmn)) %>%
  select(-SubjectID) %>%
  gather('Taxa','Proportion',-dmn,-SampleID) %>%
  mutate(dmn = paste0('ct',dmn)) %>%
  group_by(dmn, Taxa) %>%
  summarise(Proportion = mean(Proportion, na.rm = T)) %>% 
  spread(dmn, Proportion) %>% 
  mutate(sum = ct1 + ct2 + ct3 + ct4) %>%
  filter(Taxa %in% c('g__Staphylococcus','g__Streptococcus','g__Corynebacterium'))

# this is showing the median
df <- taxa$Species %>% apply(2, function(x) x/sum(x)) %>% t %>%
  as.data.frame %>% add_rownames('SampleID') %>% 
  right_join(select(dfu,SampleID,dmn), by = 'SampleID') %>%
  filter(!is.na(dmn)) %>%
  select(-SubjectID) %>%
  gather('Taxa','Proportion',-dmn,-SampleID) %>%
  mutate(dmn = paste0('ct',dmn)) %>%
  group_by(dmn, Taxa) %>%
  summarise(Proportion = mean(Proportion, na.rm = T)) %>% 
  spread(dmn, Proportion) 
row.names(df) <- df$Taxa
df$Taxa <- NULL
df <- as.matrix(df)
df <- apply(df + 0.00001, 2, logit) %>% apply(1, function(x) x - mean(x)) %>% t

# make distance matrix
df.dist <- dist(df, method = 'canberra')
# perform hclust
df.hc.row <- hclust(df.dist, method = 'complete')
# define sample colors
heatmap(df, col = redgreen(50), Rowv = as.dendrogram(df.hc.row))
