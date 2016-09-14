##load the diamonds data
data(diamonds)
diam <- diamonds

##create a quick histogram of diamond prices
hist <- qplot(price, data = diam, geom="histogram")
hist

##create a scatterplot of price to carat
a <- ggplot(diam, aes(x = price, y = carat)) + geom_point()

##add a facet for cut
b <- a + facet_wrap(~cut, nrow = 1)
b

##update our original scatterplot to add faceting by color
a <- ggplot(diam, aes(x = price, y = carat, col = color)) + geom_point()

##add cut faceting to the color-faceted scatterplot
b <- a + facet_wrap(~cut, nrow = 1)
b

##add additional faceting to include clarity as well
c <- a + facet_grid(clarity~cut)
c

##create a boxplot of cut versus price, with outliers colored in red
d <- ggplot(data = diam, aes(x = cut, y = price)) + geom_boxplot(outlier.colour = "red")
d

##add color as a dimension in our boxplot
e <- ggplot(data = diam, aes(x = cut, y = price)) + geom_boxplot(outlier.colour = "red", aes(fill=color))
e
