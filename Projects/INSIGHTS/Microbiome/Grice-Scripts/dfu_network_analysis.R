# Load in data
#load(file = '~/Coding/Club_Grice/scripts/loesche/dfu_uclust/data.RData')
load(file = '~/Coding/Club_Grice/scripts/loesche/dfu_uclust/blast.RData')

# Set up
library(loescher)
library(vegan)
library(dplyr)
library(tidyr)

# Rarify by taxa
taxa.rare <- rrarefy(x = t(taxa$Genus), sample = 1200)
#taxa.rare <- rrarefy(x = t(otu), sample = 1200)

# Make subset dfu table with just the samples
dfu.sub <- droplevels(subset(dfu, has_sample & Visit == 0))
dfu.sub <- droplevels(subset(dfu, has_sample & healed))
dfu.sub <- droplevels(subset(dfu, has_sample))

# Subset taxa table to samples
#tx <- taxa.rare[dfu.sub$SampleID,]
tx <- taxa$Genus[,dfu.sub$SampleID] %>% t

# Subset taxa to most present taxa - greater than 1 sample
tx <- tx[,colSums(tx > 0) > 1]
# Subset taxa to most present taxa - greater than 10%
#tx <- tx[,colSums(tx > 0) >= floor(nrow(tx)*0.1)]


# Spiec-Easi: MB ------------------------------------------------------------------------------

# Make spiec.easi objects
library(SpiecEasi)
se.mb <- spiec.easi(tx, method='mb', lambda.min.ratio=1e-3, #nlambda=20,
                    icov.select.params=list(rep.num=25, ncores = 4, verbose = T))

# Create igraph objects
ig.mb <- adj2igraph(se.mb$refit, rmEmptyNodes = F)
#ig.mb <- graph.adjacency(se.mb$refit, mode='undirected')

beta <- pick_edge_association(obj = se.mb, class = 'mb')$x.mean
# Make data frame containing edge information

# # Remove vertices with no connections
# del <- which(degree(ig.mb) < 1)
# ig.mb <- delete.vertices(ig.mb, v = del)
# node.tax <- colnames(tx)[-del]

# Change layout
am.coord <- layout.fruchterman.reingold(ig.mb)

# Color for edges - green/red is positive/negative
edge.color <- ifelse(beta > 0, 'forest green','red')
# Edge thickness
edge.width <- log2(abs(beta))/2 + 6
# Set size of vertex proportional to clr-mean
vsize <- rowMeans(clr(tx, 1)) + 6
#vsize <- 3 * (rowMeans(clr(tx, base = 2, mar = 1)) + 2)

# Set color of vertices by oxygen tolerance
anaerobes <- read.csv('./oxygen_tolerance.csv') %>% 
  mutate(Oxygen.Tolerance = this2that(Oxygen.Tolerance, "", NA)) %>% 
  droplevels
#want <- otu_taxa$Species[match(colnames(tx), otu_taxa$OTU)]
#anaerobes <- anaerobes$Oxygen.Tolerance[match(want, anaerobes$Species)]
want <- match(colnames(tx), anaerobes$Genus)
anaerobes <- anaerobes$Oxygen.Tolerance[want]

vcolor <- c("blue","green","red","orange","gold","brown","purple")[as.integer(anaerobes)]
vcolor <- this2that(vcolor, NA, "gray")

# Plot out network
par(mar=c(0,0,2,5))
plot(ig.mb,
     #layout = am.coord,
     vertex.size = vsize, 
     vertex.color = vcolor,
     vertex.label.color = 'white',
     vertex.label.cex = 0.65,
     edge.color = edge.color, 
     edge.width = edge.width,
     main = "Spiec-Easi Network: MB")
#legend("topright", levels(anaerobes$Oxygen.Tolerance), 
#       col=c("blue","red","purple"), pch=19, cex = 0.5)



# Spiec-Easi: Glasso --------------------------------------------------------------------------

# Make Spiec-Easi Glsso object
se.gl <- spiec.easi(tx, method='glasso', lambda.min.ratio=1e-2,
                    nlambda=20, icov.select.params=list(rep.num=50))

# Create igraph objects
ig.gl <- graph.adjacency(se.gl$refit, mode='undirected', diag=FALSE)

# Make data frame containing edge information
beta <- pick_edge_association(se.gl, class = 'glasso')$x.mean

# Color for edges - green/red is positive/negative
## Note the inversion because inverse covariance matrix used in the case of glasso
edge.color <- ifelse(beta*-1 > 0, 'forest green','red')
# Edge thickness
edge.width <- log2(abs(beta))/2 + 6

# Make plot of glasso network
plot(ig.gl,
     layout = am.coord,
     vertex.size = vsize, 
     vertex.color = vcolor,
     vertex.label.color = 'white',
     vertex.label.cex = 0.65,
     edge.color = edge.color, 
     edge.width = edge.width,
     main = "Spiec-Easi Network: Glasso")
legend("topright", levels(phyla), col=pick_colors(nlevels(phyla)), pch=19, cex = 0.5)



# Spies-Easi: SPARCC --------------------------------------------------------------------------

# Make SPARCC correlation matrix
sp <- sparcc(tx)
# Define arbitrary threshold for SparCC correlation matrix for the graph
sp.graph <- abs(sp$Cor) >= 0.3

# Create igraph objects
ig.sp <- graph.adjacency(sp.graph, mode='undirected', diag=FALSE)

# Make data frame containing edge information
beta <- pick_edge_association(sp, class = 'sparcc', cutoff = 0.3)$x.mean

# Remove vertices with no connections
#del <- which(degree(ig.sp) < 1)
#ig.sp <- delete.vertices(ig.sp, v = del)
#node.tax <- colnames(tx)[-del]

# Vertex size
vsize <- rowMeans(clr(tx, 1)) + 6
# Color for edges - green/red is positive/negative
edge.color <- ifelse(beta > 0, 'forest green','red')
# Edge thickness
edge.width <- abs(beta) * 6

# Change layout
am.coord <- layout.fruchterman.reingold(ig.sp)

# Make plot of glasso network
plot(ig.sp,
     layout = am.coord,
     vertex.size = vsize, 
     vertex.color = vcolor,
     vertex.label.color = 'white',
     vertex.label.cex = 0.65,
     edge.color = edge.color, 
     edge.width = edge.width,
     main = "Spiec-Easi Network: SPARCC")
 


# Degree Analysis -----------------------------------------------------------------------------

#Lets look at the degree statistics from the networks inferred by each method.

par(mar=c(4,4,4,4))
dd.gl <- degree.distribution(ig.gl)
dd.mb <- degree.distribution(ig.mb)
dd.sp <- degree.distribution(ig.sp)

plot(0:(length(dd.sp)-1), dd.sp, type='b', 
     ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso", "sparcc"), 
       col=c("forestgreen", "red", "black"), pch=1, lty=1)

  #rm(vsize,edge.width,edge.color,beta,del,phyla,vcolor,ig.mb,ig.gl,sp.graph,ig.sp,dd.gl,dd.mb,dd.sp)
  #rm(am.coord,se.mb,se.gl,sp)
  