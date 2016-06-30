##Exercise 2: Using the indomethacin dataset, create individual plots for each subject using the facet_wrap function.

##read in the data
data(Indometh)
indo <- Indometh

##create the plot
indo$Subject <- factor(indo$Subject, levels = c("1", "2", "3", "4", "5", "6"))
p <- ggplot(data = indo, aes(x = time, y = conc, col = Subject)) + geom_line(size = 1) + facet_wrap(~Subject, nrow=1)
p


##Exercise 3:  Use the geom_bar() function instead of geom_line() to turn the above plot into a bar chart.  Hint: geom_bar() will need the argument stat = "identity".

##create the plot
indo$Subject <- factor(indo$Subject, levels = c("1", "2", "3", "4", "5", "6"))
p <- ggplot(data = indo, aes(x = time, y = conc, fill = Subject)) + geom_bar(stat = "identity") + facet_wrap(~Subject, nrow=1)
p
