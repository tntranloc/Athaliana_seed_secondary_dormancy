setwd("/your/path/")

library(tidyverse)
library(dplyr)
library(ggplot2)

###### 1. load your pca output from plink ######
pca = read.table("yourfilename.eigenvec", header = F)
eigenval = scan("yourfilename.eigenvec.eigenval")


head(pca) # you have around 20 columns (PCs)

#we need to do some cleaning jobs
pca = pca[,-1] # the first one is redundant with the second column
# name the first column properly, this is your sample IDs (here I call "ind" for short individual)
# all columns after that are different PCs. 
names(pca)[1] = "ind" 
names(pca)[2:ncol(pca)] = paste0("PC", 1:(ncol(pca) -1))

#check again before proceeding
head(pca) 

###### 2. load your sample info ######
samples = read.csv("/path/to/your/file.csv", header = T)
head(samples)
#here I ordered it to make the ID columns in same order between pca and sample data
  # feel free to order different ways but always make sure the sample order is the same
samples$ID = as.character(samples$ID)
pca$ind = as.character(pca$ind) # very important to make them same class before re-ordering
samples = samples[order(samples$ID),]
pca = pca[order(pca$ind),]

  #check by
identical(sdorm$ID, pca$ind) # should return >TRUE

# check how many populations you have, here my column is named "pop"
summary(samples$pop) 

###### 3. combine your data and pca result for plotting ######
#attach the population info to PCA
pca = as_tibble(data.frame(pca, samples$pop))
head(pca)
#set that population info column with proper name
colnames(pca)[22] = "Pop"
pca$Pop = as.factor(pca$Pop)
#check if things are good # how many populations you have and how many individuals each
summary(pca$Pop)
#calculate the Proportion of Variancce Explained
pve =data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
head(pve)


###### 4. To see how many percentage of variance is explained by the PCs ######
pve_plot = ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
  ylab("Percentage of variance explained") + theme_light() 
pve_plot 
#so I see PC1 and PC2 explain 40% of my data, I would just present my data structure with PC1 and PC2

###### 5. plotting ######

p1 = ggplot(pca, aes(PC1, PC2, col = Pop)) + geom_point() +
  #I set colours manually, I seven populations hence 7 colours
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) +
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (",signif(pve$pve[1],3), "%)")) +
  ylab(paste0("PC2 (",signif(pve$pve[2],3), "%)"))
p1

#here are name code of R colour palettes: https://r-graph-gallery.com/38-rcolorbrewers-palettes.html


p2 = ggplot(pca, aes(PC1, PC2, col = Pop)) + geom_point() +
  #you can also just give a palette # limited to 9-12 colours 
  #if you have more populations, just set colours manually
  scale_color_brewer(palette = "Dark2") + 
  coord_equal() + theme_light() +
  xlab(paste0("PC1 (",signif(pve$pve[1],3), "%)")) +
  ylab(paste0("PC2 (",signif(pve$pve[2],3), "%)"))
p2



################## Resources ####################
#here are R colour codes https://r-graph-gallery.com/42-colors-names.html

#here are name code of R colour palettes: https://r-graph-gallery.com/38-rcolorbrewers-palettes.html
#or just display it in R
library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()










