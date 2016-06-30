##load the ggplot2 package. Don't worry if you get an error message about the package being built under a different version of R.
library(ggplot2)

##load the dataset - we'll be using a default dataset called Indometh and put it in a dataframe called "indo"
data(Indometh)
indo <- Indometh

##check out the documentation for the dataset so we know what we're looking at
?Indometh

##look at the beginning of the data so we can see how it's structured
head(indo)

##define the data mapping for our visualization
a <- ggplot(data = indo, aes(x = time, y = conc))

##create a scatterplot
b <- a + geom_point()

##color each point by subject
a <- ggplot(data = indo, aes(x = time, y = conc, col = Subject))

##create the colored scatterplot
b <- a + geom_point()

##instead of colors, use a different shape for each subject
a <- ggplot(data = indo, aes(x = time, y = conc, shape = Subject))

##create the scatterplot with different shapes for each subject
b <- a + geom_point()

##create a line graph with the data
a <- ggplot(data = indo, aes(x = time, y = conc, col = Subject))
c <- a + geom_line()

##resize the lines
c <- a + geom_line(size = 1)

##relabel axes
d <- c + xlab("Time, hours") + ylab("Plasma concentration, mcg/ml") + ggtitle("Pharmacokinetics of Indomethacin")

##add points on the line chart
e <- d + geom_point(size = 3)
e

##reorder items in the legend
f <- e + scale_colour_hue(breaks = c("1", "2", "3", "4", "5", "6"))
f

##change the location of the legend
g <- f + theme(legend.position = "top")

##place the legend within the graph itself
g <- f + theme(legend.position=c(.9,.8))

##create a tick mark for each time point measured in our dataset
h <- g + scale_x_continuous(breaks = c(0.25, 0.5, .75, 1, 1.25, 2, 3, 4, 5, 6, 8))

##change the tick marks to make them evenly spaced and add more tick marks on the y-axis
h <- g + scale_x_continuous(breaks = c(0.5:8)) + scale_y_continuous(breaks=seq(0, 3, 0.25))

##switch to black and white color scheme
i <- h + theme_bw()
i

##change colors!
i <- h + theme(
  panel.background = element_rect(fill = "lavenderblush1", colour = "thistle4"),
  legend.background = element_rect(fill = "white", colour = "black"),
  legend.key = element_rect(fill = "white"),
  panel.grid.major = element_line(colour = "thistle4"))

##return to the black and white theme and add some customization of text
j <- h + theme_bw() + theme(
  title = element_text(family="serif", colour = "grey38"),
  text = element_text(family="serif", colour = "grey38"),
  plot.title = element_text(face = "bold", size = 24),
 axis.title = element_text(face = "italic", size = 12)
 )
j
