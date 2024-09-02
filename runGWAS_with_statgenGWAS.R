####If preferred running gwas in R ####
### statgenGWAS allows adding one more covariate beside kinship (emma libraries only allow adding kinship as covariate)
### see complete_mapping_from_fastq_to_vcf_workflow.sh script 
  ### for running GWAS in Bash using Gemma -- high recommended for efficiency #####



###PHENOTYPE data
# "pheno" dataframe with columns of sample and traits
# column 1 is sample ID; other columns are your traits of interest

# prepare covariate matrix
# for example, getting ID and Origin columns our of your pheno df (considerng Origin is a covariate here)
covar = pheno[,c('ID','Origin']
library(textshape)
covar = column_to_rownames(covar, "ID")


###KINSHIP data
load("K_382.rds") #for example, this is your kinship matrix
#make sure you have the same number of sample ID in your data and your kinship matrix here
a = unique(pheno$ID)
K.sub = subset(K, rownames(K) %in% a, colnames(K) %in% a)
dim(K.sub)
b = unique(rownames(K.sub))

pheno.sub = subset(pheno, pheno$ID %in% b)
covar.sub = subset(covar, rownames(covar) %in% b)

rm(a,b)

#and the same order as well
identical(pheno.sub$ID, rownames(K.sub)
#if TRUE, moving on. if FALSE, order your data
#check that names are all of same type, f.e., character or numeric, then order
K.sub = K.sub[order(rownames(K.sub), order(rownames(K.sub)]

colnames(pheno.sub) = c("genotype","FT16") #statgenGWAS requires "genotype" as colname!

##GENOTYPE DATA
# "genotype" matrix
#row names are sample ID, column names are SNP ID
#check if names are all the same and in same order also
identical(rownames(genotype),pheno.sub$genotype)


###MAP data
#marker map with 2 columns - chromosome and position; rownames are SNP marker names


snps.info = data.frame(cbind(colnames(genotype),
                               matrix(nrow=ncol(genotype.sub),ncol=2,
                                      data=unlist(strsplit(colnames(genotype.sub),split='- ')),byrow=T)))
colnames(snps.info)<-c('SNP','chr','pos') #name the col names like this as statgenGWAS requires

snps.info[,2] = as.numeric(snps.info[,2])
snps.info[,3] = as.numeric(snps.info[,3])
snps.info = column_to_rownames(snps.info, "SNP")




#CREATE gDATA with statgenGWAS
library(statgenGWAS)
mygData = createGData(geno = genotype.sub, map = snps.info, pheno = pheno.sub, 
                          kin = K.sub, covar = covar.sub)
summary(mygData) #this can take quite a while and a bit unneccesary  

plot(mygData)

#run single trait GWAS
mygwas = runSingleTraitGwas(gData = gDataSeed,
                               traits = c("sdorm","pdorm","FT16"), #here you can put several traits of interest, they run seperately anyways
                               remlAlgo = "EMMA",
                               GLSMethod="single",
                               useMAF = TRUE,
                               MAF=0.05,
                               thrType = "bonf",
                               sizeInclRegion = 0)

head(mygwas$GWAResult$pheno.sub) 
summary(mygwas)

#print significant snps
print(GWASDor$signSnp$pheno.sub, row.names = FALSE)
#manhattan plot
plot(mygwas,
     plotType = "manhattan", trait = "FT16",
     col = c("#88CCEE", "#44AA99", "#117733", "#999933","#DDCC77")) # here I have only 5 chromosomes to colour


################################### DONE ###################################


