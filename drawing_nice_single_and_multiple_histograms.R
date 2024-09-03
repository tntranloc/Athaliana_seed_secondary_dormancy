## Functions to draw multiple or single histogram
# single histogram
plot_histogram = function(df, feature) {
  plt = ggplot(df, aes(x=eval(parse(text=feature)))) +
    geom_histogram(aes(y = ..density..), alpha=0.7, fill="darkorange3", color="black") +
    geom_density(alpha=0.3, fill="purple4") +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density") + 
    theme_classic() +
    theme(plot.title=element_text(family = "Lucida Grande", face = "bold", size=13,hjust=0.5,vjust=1), 
          legend.position = "right", 
          legend.title = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          legend.text = element_text(colour = "black", size = 11, family = "Lucida Grande"),
          axis.title.y = element_text(colour = 'black',  size = 13, family = "Lucida Grande", angle=90), 
          axis.title.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          axis.text.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          axis.text.y = element_text(colour = "black", size = 13, family = "Lucida Grande"))
  print(plt)
}

plot_histogram(data, "variable")

## overlaying histogram - simple

plot_histogram = function(df, feature) {
  plt = ggplot(df, aes(x=eval(parse(text=feature)))) +
    geom_histogram(aes(y = ..density..), alpha=0.7, fill="orange3", color="black") +
    geom_density(alpha=0.3, fill="skyblue") +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density") + 
    theme_classic() +
    theme(plot.title=element_text(family = "Lucida Grande", face = "bold", size=13,hjust=0.5,vjust=1), 
          legend.position = "right", 
          legend.title = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          legend.text = element_text(colour = "black", size = 11, family = "Lucida Grande"),
          axis.title.y = element_text(colour = 'black',  size = 13, family = "Lucida Grande", angle=90), 
          axis.title.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          axis.text.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          axis.text.y = element_text(colour = "black", size = 13, family = "Lucida Grande"))
  print(plt)
}


#### overlaying histograms - add colour palette

plot_multi_histogram = function(df, feature, label_column) {
  plt = ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.5, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.5) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=label_column)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(plot.title=element_text(family = "Lucida Grande", face = "bold", size=13,hjust=0.5,vjust=1), 
          legend.position = "right", 
          legend.title = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          legend.text = element_text(colour = "black", size = 11, family = "Lucida Grande"),
          axis.title.y = element_text(colour = 'black',  size = 13, family = "Lucida Grande", angle=90), 
          axis.title.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          axis.text.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          axis.text.y = element_text(colour = "black", size = 13, family = "Lucida Grande"))
}

#options(repr.plot.width = 20, repr.plot.height = 8)
plot_multi_histogram(compare2, "variable", "group")


#####add mean for each set
plot_multi_histogram = function(df, feature, label_column, means) {
  plt = ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    geom_vline(xintercept=means, color="black", linetype="dashed", size=1)
  labs(x=feature, y = "Density")
  plt + guides(fill=guide_legend(title=label_column)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_classic() +
    theme(plot.title=element_text(family = "Lucida Grande", face = "bold", size=13,hjust=0.5,vjust=1), 
          legend.position = "right", 
          legend.title = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          legend.text = element_text(colour = "black", size = 11, family = "Lucida Grande"),
          axis.title.y = element_text(colour = 'black',  size = 13, family = "Lucida Grande", angle=90), 
          axis.title.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          axis.text.x = element_text(colour = "black", size = 13, family = "Lucida Grande"),
          axis.text.y = element_text(colour = "black", size = 13, family = "Lucida Grande"))
}

plot_multi_histogram(data, "variable", "group", c(1,2,3))
