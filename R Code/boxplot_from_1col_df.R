## boxplot from a single column of df
male = data.frame(c(127,44,28,83,0,6,78,6,5,213,73,20,214,28,11)) # data from page 66
require(ggplot2)
ggplot(data = male, aes(x = "", y = male)) + 
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 150)) # I set the y axis scale so the plot looks better.
