library(loescher)
library(plyr)
library(dplyr)
library(tidyr)

load('./blast.RData')

# Read in dmn data
dmn <- read.csv('./dfu_dmn_assignments.csv') %>% select(SampleID, dmn) %>% mutate(dmn = as.factor(dmn))
dfu <- dfu %>% select(-dmn) %>% left_join(dmn, by = 'SampleID')

# Split data frame into patients that had complications and those that did not
df <- split.data.frame(dfu, dfu$heals.12wks)
names(df) <- c('no.comp','comp')
#df <- list(dfu)

# Function that will count transitions in a factor vector
count_transitions <- function(x, count_bookends = T, keep_na = T) {
  # make sure x is a factor
  if(class(x) != 'factor') x <- as.factor(x)
  # create a named vector with possible transitions
  result <- rep(0, times = (nlevels(x)+1)^2)
  names(result) <-  paste(rep(c(levels(x),NA), each = nlevels(x) + 1),
                          rep(c(levels(x),NA), times = nlevels(x) + 1),
                          sep = '->')
  # get last index
  end <- length(x)
  # initialize transition vector
  trans <- paste(NA, x[1], sep = '->')
  # loop through time points
  for(i in 1:end) {
    trans <- c(trans, paste(x[i],x[i+1],sep='->'))
  }
  # remove bookends if count_bookends = F
  if(!count_bookends) trans <- trans[2:(length(trans)-1)]
  # count transitions
  count <- table(trans)
  # add to result data frame
  result[names(count)] <- count
  # return result data frame
  if(keep_na) {
    return(result)
  } else {
    return(result[grep('NA', names(result),invert = T)])
  }
}
# Count transitions
trans <- lapply(df, function(d) {
           sapply(unique(d$SubjectID), function(n) {
             filter(d, SubjectID == n) %>% .$dmn %>% 
               count_transitions(count_bookends = T, keep_na = F)
           }) %>% rowSums %>% data.frame(freq = .)
          })
# Split label
trans <- lapply(trans, function(tr) {
           tr[,c('current','next')] <- row.names(tr) %>% 
             as.character %>% strsplit(split = '->') %>% do.call('rbind',.)
           return(tr)
         })
# Spread data into wide format: rows are current, columns are next
trans <- lapply(trans, function(tr) {
           tr <- tr %>% spread(`next`, freq)
           row.names(tr) <- tr$current
           tr <- tr[,-1]
         })
# Normalize counts into proportions by row
trans.n <- lapply(trans, function(tr) t(apply(tr, 1, function(x) x/sum(x))))

# Calculate 
t_persist <- sapply(trans.n, function(tr) -1/log(diag(tr)))

# test for differences
ldply(trans, function(x) {
  x %>% as.data.frame %>% 
    add_rownames('Start') %>% 
    gather('End','Count',2:5)
}) %>% spread(.id, Count) %>% select(no.comp,comp) %>% 
  fisher.test(simulate.p.value = T, B = 1000000)
  #fisher.test(simulate.p.value = F, workspace = 10000000)


# Make plot for no comp -----------------------------------------------------------------------

# Libraries
library(markovchain)
library(igraph)
library(loescher)

# Set current data frame
d <- df$no.comp
tr <- trans$no.comp
tr.n <- trans.n$no.comp
# d <- df[[1]]
# tr <- trans[[1]]
# tr.n <- trans.n[[1]]

# Create markov chain plot object
mc <- new("markovchain", states=row.names(tr.n),
              transitionMatrix = tr.n, name="DFU DMN")
netMC <- markovchain:::.getNet(mc, round = TRUE)

# Set arrow weights to proportions
wts <- E(netMC)$weight/100
# Create variables containing transition types
edgel <- get.edgelist(netMC)
elcat <- paste(edgel[,1], edgel[,2])
elrev <- paste(edgel[,2], edgel[,1])

# Make graphing parameter to figure out if arrow edges should be curved
# to show both directions clearly
edge.curved <- sapply(elcat, function(x) x %in% elrev)
# Make color scale for vertices based on proportion of 
vert.clrs <- pick_colors(4)
# Make vertex size based on frequency of cluster
vert.size <- rowSums(tr)/sum(tr) * 150
# Make loop edge angles
edge.loop.angle <- c(3.925,0,0,0,0,2.35,0,0,0,0,.785,0,0,0,-.785)
# Make layout matrix
layout.matrix <- matrix(c(-1,1,-1,-1,1,-1,1,1), byrow = T, ncol = 2)*.5
# Make edge labels
edge.label <- as.character(round(E(netMC)$weight/100,2))
edge.label[-c(1,6,11,15)] <- NA
# Make plot
pdf('./markovplot_split_heals12wks.pdf', height = 6, width = 6)
par(mar=rep(1,4))#,oma=rep(1,4))
plot.igraph(netMC, 
            edge.arrow.size = 0.8,
            edge.arrow.width = 2,
            edge.width = wts*15,
            edge.color = 'gray47',
            edge.curved=edge.curved,
            edge.loop.angle = edge.loop.angle,
            edge.label = edge.label,
            edge.label.color = 'gray25',
            edge.label.font = 2,
            edge.label.cex = .9,
            edge.label.family = 'Courier',
            edge.label.y = c(1,0,0,0,0,-.95,0,0,0,0,-.95,0,0,0,.90),
            vertex.color=vert.clrs,
            vertex.size=vert.size,
            vertex.label.font = 2, 
            vertex.label.family = 'Courier',
            vertex.label.cex = c(0.9,0.9,0.85,0.9),
            vertex.label.color = 'white', 
            vertex.frame.color = NA, 
            layout=layout.matrix,
            frame = F,
            rescale = F,
            main = 'No healing within 12 weeks')


# Make plot for comp -----------------------------------------------------------------------

# Set current data frame
d <- df$comp
tr <- trans$comp
tr.n <- trans.n$comp

# Create markov chain plot object
mc <- new("markovchain", states=row.names(tr.n),
          transitionMatrix = tr.n, name="DFU DMN")
netMC <- markovchain:::.getNet(mc, round = TRUE)

# Set arrow weights to proportions
wts <- E(netMC)$weight/100
# Create variables containing transition types
edgel <- get.edgelist(netMC)
elcat <- paste(edgel[,1], edgel[,2])
elrev <- paste(edgel[,2], edgel[,1])

# Make graphing parameter to figure out if arrow edges should be curved
# to show both directions clearly
edge.curved <- sapply(elcat, function(x) x %in% elrev)
# Make color scale for vertices based on proportion of 
vert.clrs <- pick_colors(4)
# Make vertex size based on frequency of cluster
vert.size <- rowSums(tr)/sum(tr)*150
# Make loop edge angles
edge.loop.angle <- c(3.925,0,0,0,0,2.35,0,0,0,0,.785,0,0,-.785)
# Make layout matrix
layout.matrix <- matrix(c(-1,1,-1,-1,1,-1,1,1), byrow = T, ncol = 2)*.5
# Make edge labels
edge.label <- as.character(round(E(netMC)$weight/100,2))
edge.label[-c(1,6,11,14)] <- NA
# Make plot
par(mar=rep(1,4))#,oma=rep(1,4))
plot.igraph(netMC, 
            edge.arrow.size = 0.8,
            edge.arrow.width = 2,
            edge.width = wts*20,
            edge.color = 'gray47',
            edge.curved=edge.curved,
            edge.loop.angle = edge.loop.angle,
            edge.label = edge.label,
            edge.label.color = 'gray25',
            edge.label.font = 2,
            edge.label.cex = 0.9,
            edge.label.family = 'Courier',
            edge.label.y = c(1.05,0,0,0,0,-.95,0,0,0,0,-.92,0,0,.9),
            vertex.color=vert.clrs,
            vertex.size=vert.size,
            vertex.label.font = 2, 
            vertex.label.family = 'Courier',
            vertex.label.cex = 0.9,
            vertex.label.color = 'white', 
            vertex.frame.color = NA, 
            layout=layout.matrix,
            frame = F,
            rescale = F,
            main = 'Heals within 12 weeks')
dev.off()



# Markov Chain Analysis -----------------------------------------------------------------------

mc.heals <- new("markovchain", states=row.names(trans.n$comp),
                transitionMatrix = trans.n$comp, name="Heals < 12 weeks")
mc.noheal <- new("markovchain", states=row.names(trans.n$no.comp),
                transitionMatrix = trans.n$no.comp, name="Heals > 12 weeks")
# 1) the stationary distribution (long term effect) of the four CT categories in fast and slow
# healing patients, which could be used to predict the relative effect of each CT during DFU healing
## calculate the steady state distribution
steadyStates(mc.heals)
steadyStates(mc.noheal)

# 2). the mean expected recurrence time of each CT, which could be used to guide potential
# antibiotic treatment intervals.
2/steadyStates(mc.heals) # number of weeks
2/steadyStates(mc.noheal) # number of weeks
