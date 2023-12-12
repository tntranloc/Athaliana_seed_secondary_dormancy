getwd()
setwd("/home/alle/Documents/NhuTran/dormancy02/20230511_remap_1001genome_Frcol/vcfstats")
maf = read.table("20230610_remap_1001_FrCol_mpileup_allvar_filtered_maf.frq", header = T)

maf = read.table("20230611_remap_1001_FrCol_mpileup_allvar_fullsamples_rmindels_filtered_maf_missing_minQ_renamed_MAFcheck.frq", header = T)

head(maf)
h = hist(maf$MAF)
plot(h, freq = FALSE)
library(tidyverse)
pca = read.table("230623_myPCA_remap_samples_sdorm02.eigenvec", header = F)
eigenval = scan("230623_myPCA_remap_samples_sdorm02.eigenval")

setwd("/home/alle/Documents/NhuTran/dormancy02")
pca = read.table("230627_myPCA_remap_samples_sdorm_ld08.eigenvec", header = F)
pve = read.table("230627_myPCA_remap_samples_sdorm_ld08.eigenval")

head(pca)
pca = pca[,-1]
names(pca)[1] = "ind"
names(pca)[2:ncol(pca)] = paste0("PC", 1:(ncol(pca) -1))

head(sdorm)
sdorm = sdorm[order(sdorm$ID),]
rownames(sdorm) = NULL
View(sdorm)
sdorm = sdorm[-c(337,338),]
sdorm = sdorm[c(337:363,1:336),]
sdorm = sdorm[-c(261,357),]
sdorm$Admixture_as_database[1:27] = "france"
sdorm$Admixture_as_database[1:27] = "western_europe"

library(tidyverse)
pca = as_tibble(data.frame(pca, sdorm$Admixture_as_database))
head(pca)
colnames(pca)[22] = "Pop"
pca$Pop = as.factor(pca$Pop)
summary(pca$Pop)
pve =data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
head(pve)

library(ggplot2)
a = ggplot(pca, aes(PC1, PC2, col = Pop)) + geom_point() +
  #scale_color_brewer(palette = "Dark2") + 
  coord_equal() + theme_light() +
  #scale_color_manual(values = c("plum1","lightsalmon","grey0","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) +
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) +
  xlab(paste0("PC1 (",signif(pve$pve[1],3), "%)")) +
  ylab(paste0("PC2 (",signif(pve$pve[2],3), "%)"))
a

d = ggplot(pca, aes(PC2, PC3, col = Pop)) + geom_point() +
  #scale_color_brewer(palette = "Dark2") + 
  coord_equal() + theme_light() +
  scale_color_manual(values = c("plum1","lightsalmon","grey0","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) +
  xlab(paste0("PC2 (",signif(pve$pve[2],3), "%)")) +
  ylab(paste0("PC3 (",signif(pve$pve[3],3), "%)"))
d

e = ggplot(pca, aes(PC1, PC3, col = Pop)) + geom_point() +
  #scale_color_brewer(palette = "Dark2") + 
  coord_equal() + theme_light() +
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) +  
  xlab(paste0("PC1 (",signif(pve$pve[1],3), "%)")) +
  ylab(paste0("PC3 (",signif(pve$pve[3],3), "%)"))
e

b=ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") +
  ylab("Percentage of variance explained") + theme_light() 
b

dim(pca)
dim(sdorm)



library(tidyverse)
setwd("/home/alle/Documents/NhuTran/dormancy02/20230511_remap_1001genome_Frcol")

var_freq= read_delim("20230610_remap_1001_FrCol_mpileup_allvar_subset2_allele_stats.frq", delim = "\t",
                     skip = 1, col_names = c("chr","pos", "nalleles", "nchr", "a1", "a2"))
var_freq= read_delim("20230610_remap_1001_FrCol_mpileup_allvar_subset2_allele_stats.frq", delim = "\t",
                     skip = 1)
head(var_freq)
var_freq = read_tsv("20230610_remap_1001_FrCol_mpileup_allvar_subset2_allele_stats.frq")
var_freq$maf = var_freq %>% select(FREQ1,FREQ2) %>% apply(1, function(z) min(z))

d = ggplot(var_freq, aes(maf)) + geom_density(fill="dodgerblue1", colour = "black",
                                               alpha = 0.3) + theme_light()
d

var_qual= read_delim("20230610_remap_1001_FrCol_mpileup_allvar_subset2_site_qual.lqual", delim = "\t",
                     col_names = c("chr","pos","qual"), skip =1)

a = ggplot(var_qual, aes(qual)) + geom_density(fill="dodgerblue1", colour = "black",
                                               alpha = 0.3) + theme_light() + xlim(50, 80000)
a

var_depth = read_delim("20230610_remap_1001_FrCol_mpileup_allvar_subset2_site_depth.ldepth.mean", delim = "\t", 
                       skip =1, col_names=c("chr","pos","mean_depth","var_depth"))

b = ggplot(var_depth, aes(mean_depth)) + geom_density(fill="dodgerblue1", colour = "black",
                                               alpha = 0.3) + theme_light()
b
summary(var_depth$mean_depth)

var_miss = read_delim("20230610_remap_1001_FrCol_mpileup_allvar_subset2_missing_site.lmiss",  delim = "\t", 
                      skip =1, col_names=c("chr","pos","nchr","nfiltered", "nmiss", "fmiss"))
c = ggplot(var_miss, aes(fmiss)) + geom_density(fill="dodgerblue1", colour = "black",
                                                      alpha = 0.3) + theme_light() 
c

summary(var_miss$fmiss)




#################################################
###gwas run function
GWAS_run <- function(output_gemma, threshold_pvalue="0", highlighted_SNP=""){
  
  # Highlighted_SNP allows to display in green the SNP of interested on the Manahattan plot
  # It can be 1 SNP (e.g. highlighted_SNP="4:10420088") or several SNPs, passed as a vector
  # (e.g. highlighted_SNP=c("4:10420088","5:112000"). No SNP highlighted by default
  
  # Import GEMMA output file
  gwas.results <- read.delim(path.file, sep="\t")
  
  # Plot QQ plot (need to precise the package as lattice has a similar function
  qqman::qq(gwas.results$P, main=file.name)
  
  # One can select SNPs above the Bonferroni corrected p-value threshold
  # by using the argument "bonferroni"
  if(threshold_pvalue == "bonferroni"){
    # Calculate Bonferroni threshold with risk 5%
    ## Get total number of SNPs
    nb_snps <- dim(gwas.results)[[1]]
    
    ## Calculate Bonferroni corrected P-value threshold
    bonferroni_threshold <- 0.05/nb_snps
    
    threshold_pvalue <- bonferroni_threshold
  } else {
    # In case the variable was entered as string and is not "bonferroni"
    # convert to numeric. Set to 0 by default if user does not want any threshold
    threshold_pvalue <-  as.numeric(threshold_pvalue)
  }
  
  # Get positions of the chromosome with SNPs having a -log(P) > 5
  gwas_significant <- subset(gwas.results, P < threshold_pvalue)
  
  # Default p-value threshold line commonly used in GWAS -> -log10(5e-8) => red line. 
  # Set genomewideline to False has it makes little sense for Arabidopsis genome
  
  # suggestive line = Bonferroni corrected P-value threshold => blue line
  
  # Plot manhattan plot
  manhattan(gwas.results, highlight=highlighted_SNP, main=file.name, suggestiveline = -log10(threshold_pvalue), genomewideline = FALSE)
  
  #Check if dataframe is not empty (no SNPs above threshold value
  if(dim(gwas_significant)[[1]] != 0){ 
    return(gwas_significant)
  }
}


##loading gwas res
setwd("/home/alle/Documents/NhuTran/dormancy02/20230511_remap_1001genome_Frcol")
gwa_bio9 = read.table("2306_gwas_remap_bio9.assoc.txt", header = T, sep = "\t")
head(gwa_bio9)

gwa_bio9_2 = read.table("mygwas_20230630_remap_1001_FrCol_mpileup_allvar_fullsamples_filtered_SDORM02_bio9.assoc.txt", header = T, sep = "\t")
gwa_bio9 = read.table("mygwas_20230611_remap_1001_FrCol_mpileup_allvar_fullsamples_filtered_DORM01_bio9.assoc.txt", header = T, sep = "\t")
gwa_bio9_3 = read.table("mygwas_20230701_remap_1001_FrCol_mpileup_allvar_fullsamples_filtered_SDORM02_bio9.assoc.txt", header = T, sep = "\t")

gwa_sdorm = read.table("2306_gwas_remap_sdorm.assoc.txt", header = T, sep = "\t")
gwa_sdorm_2 = read.table("mygwas_20230630_remap_1001_FrCol_mpileup_allvar_fullsamples_filtered_SDORM02_dorm.assoc.txt", header = T, sep = "\t")
gwa_sdorm_3 = read.table("mygwas_20230701_remap_1001_FrCol_mpileup_allvar_fullsamples_filtered_SDORM02_sdorm.assoc.txt", header = T, sep = "\t")

gwa_bio9_3$rs = paste(gwa_bio9_3$chr, gwa_bio9_3$ps, sep= "_")
head(gwa_bio9_3)
gwa_sdorm_3$rs = paste(gwa_sdorm_3$chr, gwa_sdorm_3$ps, sep= "_")
head(gwa_sdorm1)


gwa_bio9_chr1to5 = gwa_bio9
gwa_bio9_chr1to5$chr = as.numeric(gwa_bio9_chr1to5$chr)
gwa_bio9_chr1to5 = na.omit(gwa_bio9_chr1to5)

gwa_sdorm_chr1to5 = gwa_sdorm
gwa_sdorm_chr1to5$chr = as.numeric(gwa_sdorm_chr1to5$chr)
#gwa_sdorm_chr1to5 = na.omit(gwa_sdorm_chr1to5)

threshold_pval_sdorm = 0.05/dim(gwa_sdorm_chr1to5)[[1]]
threshold_pval_sdorm
-log10(threshold_pval_sdorm)

threshold_pval_sdorm = 0.05/dim(gwa_sdorm_3)[[1]]


threshold_pval_bio9 = 0.05/dim(gwa_bio9_chr1to5)[[1]]
threshold_pval_bio9 = 0.05/dim(gwa_bio9)[[1]]
threshold_pval_bio9

library(qqman)

manhattan(gwa_bio9_chr1to5, chr="chr", bp="ps", p="p_score", snp = "rs", suggestiveline=-log10(threshold_pval_bio9), genomewideline = FALSE)
manhattan(gwa_sdorm_3, chr="chr", bp="ps", p="p_score", snp = "rs", suggestiveline = 7, genomewideline = FALSE, annotatePval = 5)
manhattan(gwa_sdorm1, chr="chr", bp="ps", p="p_score", snp = "rs", suggestiveline = 7.3, genomewideline = FALSE, annotatePval = 7.3)

manhattan(gwa_bio9_1, chr="chr", bp="ps", p="p_score", snp = "rs", suggestiveline=7.3, genomewideline = FALSE, annotatePval = 7.3)


qq(gwa_sdorm_chr1to5$p_score)
qq(gwa_bio9_chr1to5$p_score)


plot(myGwas_south,
     plotType = "manhattan", trait = "rate",
     col = c("#88CCEE", "#44AA99", "#117733", "#999933","#DDCC77"))

library(dplyr)
don = gwa_bio9_1 %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(ps)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwa_bio9_1, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, ps) %>%
  mutate( pscum=ps+tot) %>%
  # Add highlight and annotation information
  mutate( is_annotate=ifelse(-log10(p_score)>7.3, "yes", "no")) 

axisdf = don %>% group_by(chr) %>% summarize(center=( max(pscum) + min(pscum) ) / 2 )

library(ggplot2)
a = ggplot(don, aes(x=pscum, y=-log10(p_score))) +
  # Show all points
  geom_point( aes(color=as.factor(chr)), alpha=0.8, size=1.3) +
  scale_color_manual(values = c("#88CCEE", "#44AA99", "#117733", "#999933","#DDCC77")) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  geom_point(data=subset(don, is_annotate=="yes"), size=2, color = "red") +
  ylim(0,9) + 
  ylab("-log10(p)") +
  xlab("Chromosome") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
a
setwd("/home/alle/Documents/NhuTran/dormancy01")
gwa_bio9_1 = read.table("mygwas_20230703_remap_1001_FrCol_mpileup_allvar_fullsamples_filtered_DORM01_bio9.assoc.txt", header = T, sep = "\t")

gwa_sdorm1 = read.table("mygwas_20230701_remap_1001_FrCol_mpileup_allvar_fullsamples_filtered_SDORM02_as_in_DORM01_sdorm.assoc.txt", header = T, sep = "\t")


gwa_bio9_1$rs = paste(gwa_bio9_1$chr, gwa_sdorm1$ps, sep= "_")
gwa_sdorm1$chr = as.numeric(gwa_sdorm1$chr)

threshold_pval_sdorm1 = 0.05/dim(gwa_sdorm1)[[1]]
-log10(threshold_pval_sdorm1)
manhattan(gwa_bio9_1, chr="chr", bp="ps", p="p_score", snp = "rs", suggestiveline = -log10(threshold_pval_sdorm1), genomewideline = FALSE, annotatePval = -log10(threshold_pval_sdorm1))

gwa_bio9_1 = read.table("mygwas_20230701_remap_1001_FrCol_mpileup_allvar_fullsamples_filtered_SDORM02_as_in_DORM01_bio9.assoc.txt", header = T, sep = "\t")
gwa_bio9_1$rs = paste(gwa_bio9_1$chr, gwa_bio9_1$ps, sep= "_")

trr = read.csv("/home/alle/Documents/NhuTran/dormancy02/externals/TRR_acessions_385.csv", header = T)
head(trr)
trr$No = as.character(trr$No)

sdorm$ID = as.character(sdorm$ID)
sdorm$Internal_id = as.character(sdorm$Internal_id)
colnames(trr)[1] = "Internal_id"
df1 = full_join(sdorm, trr, by = "Internal_id")
head(df1)
View(df1)

dim(df1) 
sum(is.na(df1$ID))
write.table(df1, "TRR_384_accessions_withID.csv", row.names = F, col.names = T)





write.table(b1$ID, "2023037_sdorm02_id_only_in_dorm01.txt", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(b1$rate, "2023037_sdorm02_only_in_dorm01_pheno_sdorm.txt", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(b1$bio9, "2023037_sdorm02_only_in_dorm01_pheno_bio9.txt", row.names = F, col.names = F, sep = "\t", quote = F)




sd = subset(a, Treatment == "Secondary_dormancy")
View(sd)
head(sd)
rownames(sd) = NULL
sd = sd[order(sd$ID),]
sd = sd[-c(267,268),]
sd = sd[-c(263,191),]
sd = sd[c(265:291,1:264),]

write.table(sd$ID, "20230629_dorm01_IDs.txt", row.names = F, col.names = F, sep = "\t", quote = F)
write.table(sd1$Germination_rate, "20230629_dorm01_pheno_sdorm_matchingvcf.txt", row.names = F, col.names = F, quote = F)
write.table(sd1$bio9, "20230629_dorm01_pheno_bio9_matchingvcf.txt", row.names = F, col.names = F, quote = F)
sd$bio9 = sd$bio9/10

sd1 = sd[-c(68,137,203,236,242),]

###################################################################################

setwd("/home/alle/Documents/NhuTran/dormancy01/")
a = read.csv("20230629_dorm01_IDs.txt", header = F)
head(a)

dorm = read.csv("seed_id_eu_filtered.for.viability.csv", header = T)
head(a1)
a1 = subset(dorm, Treatment == "Secondary_dormancy")
a1 = subset(a2, a1$ID %in% b1$ID)
a2 = subset(a2, a2$ID %in% b1$ID)
a2 = subset(dorm, Treatment == "primary_dormancy")

a1 = a1[order(a1$ID),]
a2 = a2[order(a2$ID),]


rownames(a2) = NULL
View(a1)
a1 = a1[c(260:285,1:259),]
a2 = a2[c(260:285,1:259),]

b = subset(sdorm, sdorm$ID %in% a$V1)
head(b)
b1 = b[,c("ID","rate","bio9")]
head(b1)
sum(is.na(b1))
View(b1)
b1= b1[order(b1$ID),]
rownames(b1) = NULL
b1 = b1[c(260:285,1:259),]

head(b1)
dim(b1)
b2 = subset(sdorm, sdorm$ID %in% b1$ID)
View(b2)
b2 = b2[order(b2$ID),]
rownames(b2) = NULL
b2 = b2[c(260:285, 1:259),]

identical(a1$ID, b2$ID)
plot(a1$bio9, b2$bio9)
plot(a1$Germination_rate, b2$rate)

d = merge(a1, b1, by = "ID")
head(d)

ggplot(d, aes(x=Germination_rate, y = rate, col = Admixture)) + geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) 


head(a1)
View(d)

d1 = subset(d, Admixture == "spain")
d2 = subset(d, Admixture == "sweden")

e1 = merge(a1,a2,by = "ID")
head(e1)
e2 =subset(e1, Admixture.x =="spain")
e3 =subset(e1, Admixture.x =="sweden")

ggplot(d, aes(x=Germination_rate, y = rate, col = Admixture)) + geom_point() +
  theme_minimal() +
  xlab("First exp") +
  ylab("Second exp") +
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) 

par(mfrow=c(2,2))
ggplot(a1, aes(y=Germination_rate, col = Admixture)) + geom_boxplot() +
  theme_minimal() +
  ggtitle("First exp - secondary dormancy") +
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) 
ggplot(a2, aes(y=Germination_rate, col = Admixture)) + geom_boxplot() +
  theme_minimal() +
  ggtitle("First exp - primary dormancy") +
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) 

ggplot(b2, aes(y=rate, col = Admixture_as_database)) + geom_boxplot() +
  theme_minimal() +
  ggtitle("Second exp - secondary dormancy") +
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) 
names(b2)



#primary vs secondary dormancy in first exp
ggplot(e2, aes(x=Germination_rate.x, y = Germination_rate.y, col = Admixture.x)) + geom_point() +
  theme_minimal() +
  xlab("Secondary dormancy") +
  ylab("Primary dormancy") +
  scale_color_manual(values = c("plum1","lightsalmon","firebrick2", "mediumpurple1", "burlywood3","cornflowerblue","palegreen1")) 


##########################################################################
library(readr)
library(ggplot2)
setwd("/home/alle/Documents/NhuTran/dormancy01")
ld = read.table("20230611_remap_1001_FrCol_mpileup_allvar_fullsamples_rmindels_filtered_maf_missing_minQ_renamed_SDORM02_LDCalc.stat.gz", header = F)
head(ld)
colnames(ld) = c("Dist","Mean_r2","Mean_D","Sum_r2","Sum_D", "Numberpairs")

ggplot(ld, aes(x=Dist,y=Mean_r2)) + geom_point(size = 1, col = "darkgreen", alpha = 0.5) + geom_line() +
  theme_minimal() +
  xlab("Genomic distance") + ylab("r²") 
  

###############################################
##############ERRAND DATA ARRANGING############
a = read.table("/home/alle/Downloads/sdorm02_id_list_randomised.txt", header = T)
head(a)

random = a[sample(nrow(a)),]
head(random)
random = as.data.frame(random)
dim(random)

write.table(random,"/home/alle/Downloads/sdorm02_id_list_randomised.txt", row.names = F, quote = F)
#############################################

setwd("/home/alle/Downloads")
a1 = read.csv("arrange1.csv", header = T)
a2 = read.csv("arrange2.csv", header = T)
class(a2$name)
a1 = a1[order(a1$name),]
a2 = a2[order(a2$name),]
head(a2)
rownames(a2) = NULL
View(a1)
View(a2)
a1$name[1] = "MAR2-3"

library(dplyr)
a3 = full_join(a1, a2, by = "name")

a2 = a2[!duplicated(a2$name),]
a3 = a3[!duplicated(a3$name),]

a3 = a3[order(a3$No.),]
View(a3)

write.csv(a3, "20230801_Athaliana_database_internal_ids.csv", row.names = F)






