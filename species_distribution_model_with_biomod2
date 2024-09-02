####################### BIOMOD2 Species Distribution Model ########################

######################### Required packages #########################
library(biomod2)
library(raster)
library(rasterVis)
library(terra)
library(sp)
library(geodata)
library(tidyterra)
library(ggtext)

##########################################################################
###### Input: two ecotypes to predict for, and their coordinates (longitude and latitude) ######

datsp1 = ecotype1_df[,c("Longitude","Latitude","Athaliana")] # set up data species df with lon, lat and species presence(/absence)
colnames(datsp1)[3] = "Athaliana.eco1.sdorm"  # f.e. ecotype1
#head(datsp1)

datsp2 = ecotype2_df[,c("Longitude","Latitude","Athaliana")] #set up data species df with lon, lat and species presence(/absence)
colnames(datsp2)[3] = "Athaliana.eco2.sdorm" 
#head(datsp2)

###### Get bioclimatic variables ######

### Current ###
biocur = worldclim_global(var = "bio", res = 10, "/Users/nhutran/Documents/PhD/dormancy_redo/external")

#renaming is necessary!
bio = sprintf("bio%d", 1:19)
names(biocur) = bio
eu_ext = c(-10, 50, 35, 70) # extracting whole europe range, covering Scandinavia and Caspian

#cropping the climatic data 
biocur = crop(biocur, eu_ext) 

#rasterise bio variables # DON'T SKIP!
biocur = raster::stack(biocur)

### Future ###
biofut245 = cmip6_world(model = "CanESM5", ssp = "245", 
                        time = "2041-2060" , var = "bioc", res = 10, 
                        path = "/where/to/download/to")

#biofut370 # different scenarios
#biofut585 # different scenarios

names(biofut245) = bio
biofut245 = crop(biofut245, eu_ext) 
biofut245 = raster::stack(biofut245)


### Past - Last glacial maximum ###
biopas_files = list.files(path= "/Users/nhutran/Documents/PhD/dormancy_redo/external/chelsa_LGM_v1_2B_r10m/10min",  pattern = ".tif$", full.names = TRUE)

biopas = stack(biopas_files)
biopas = crop(biopas, eu_ext) 
biopas = raster::stack(biopas)

#because the order is different, make sure to rename them correctly!
n1 = sprintf("bio%d", 10:19)
n2 = sprintf("bio%d",2:9)
names(biopas) = c("bio1", n1, n2)

#check if everything is alright

summary(biofut245)
summary(biocur)
summary(biopas)

### Some values got multiple by *10 or *100 depending on the database 
## check here for the change factor, f.e. divide bio9, bio8, by 10 for the real value
## http://www.paleoclim.org/methods/

### IMPORTANT: check if coordinates are matching ###
print(crs(biopas))
print(crs(biocur))
print(crs(biofut245))
#or 
identical(crs(biocur), crs(biofut245))
identical(crs(biocur), crs(biopas))

# IF coordinate not matching
#  chelsa one is always matched, ccglmb one is not for example
# then do this for reprojection
#biopas_reproj = projectRaster(biopas, crs = crs(biocur))

#################### DATA for the MODEL #######################
#Define these variables for biomod formating data below

myExpl = stack(biocur$bio3, biocur$bio8, biocur$bio9, biocur$bio15, biocur$bio19)

myRespName1 = "Athaliana.eco1.sdorm"
myResp1 = as.numeric(datsp1[,"Athaliana.eco1.sdorm"])
myRespXY1 = datsp1[,c("Longitude","Latitude")]

myRespName2 = "Athaliana.eco2.sdorm"
myResp2 = as.numeric(datsp2[,"Athaliana.eco2.sdorm"])
myRespXY2 = datsp2[,c("Longitude","Latitude")]



#################### MUST: frmat dataset with BIOMOD_FormatingData #######################

#number of pseudo absence should be half or up to roughly the same as your presence
dim(datsp1) 
myBiomodData1 = BIOMOD_FormatingData(resp.var = myResp1,
                                    expl.var = myExpl,
                                    resp.xy = myRespXY1,
                                    resp.name = myRespName1,
                                    PA.nb.rep = 3,
                                    PA.nb.absences = 100,
                                    PA.strategy = "sre")

myBiomodData1
plot(myBiomodData1)

dim(datsp2)
myBiomodData2 = BIOMOD_FormatingData(resp.var = myResp2,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY2,
                                     resp.name = myRespName2,
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 70,
                                     PA.strategy = "sre")

plot(myBiomodData2)




######################## RUN THE MODEL ##################
## STEP1: RUN INDIVIDUAL MODELS
# data split 80 means 70% used for train data and 30% for test data
myBiomodModelOut1 = BIOMOD_Modeling(
  myBiomodData1,
  models = c("GBM", "GLM", "MAXNET", "RF"),
  #bm.options = myBiomodOption,
  CV.nb.rep=3,
  CV.perc=0.8,
  var.import=3,
  metric.eval = c('TSS','ROC','KAPPA'),
  scale.models  = TRUE,
  modeling.id = paste(myRespName1,"mod1",sep=""))

myBiomodModelOut2 = BIOMOD_Modeling(
  myBiomodData2,
  models = c("GBM", "GLM", "MAXNET", "RF"),
  #bm.options = myBiomodOption,
  CV.nb.rep=3,
  CV.perc=0.8,
  var.import=3,
  metric.eval = c('TSS','ROC','KAPPA'),
  #save.output = TRUE,
  scale.models  = TRUE,
  modeling.id = paste(myRespName2,"mod1",sep=""))

# get all models evaluation
myBiomodModelEval1 = get_evaluations(myBiomodModelOut1)
myBiomodModelEval2 = get_evaluations(myBiomodModelOut2)

## plot models evaluation scores
bm_PlotEvalMean(bm.out = myBiomodModelOut1, 
  metric.eval = c("ROC","TSS"))
## you can plot it with different metrics
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut1, 
                   group.by = c('algo','run'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut1, 
                   group.by = c('PA','run'))
bm_PlotEvalBoxplot(bm.out = myBiomodModelOut1, 
                   group.by = c('PA','algo'))


#or check the scores of specific models, TSS scores of Random Forest
RF= myBiomodModelEval1[myBiomodModelEval1$metric.eval=="TSS" & myBiomodModelEval1$algo == "RF",]
RF #model sensitivity and specificity

######## STEP 1.5: Get Variable Importance ###########

# print variable importance

models_var_import1= get_variables_importance(myBiomodModelOut1)
head(models_var_import1)
v1= tapply(models_var_import1$var.imp, models_var_import1$expl.var, mean)

models_var_import2= get_variables_importance(myBiomodModelOut2)
head(models_var_import2)
v2= tapply(models_var_import2$var.imp, models_var_import2$expl.var, mean)

## A bit of manual work 
v1= as.data.frame(v1)
v2= as.data.frame(v2)
v1$Var = rownames(v1)
v2$Var = rownames(v2)
colnames(v1)[1] = "Var_import"
colnames(v2)[1] = "Var_import"
v12 = rbind(v1, v2)
rownames(v12) = NULL
v12$group = NA
v12[1:5,3] = "Eco1"
v12[6:10,3] = "Eco2"
summary(v12)
v12$group = as.factor(v12$group)
v12$Var = as.factor(v12$Var)

v12$var_import_perc = NA
v12[1,4] = v12[1,1]/sum(v12[1:5,1])
v12[2,4] = v12[2,1]/sum(v12[1:5,1])
v12[3,4] = v12[3,1]/sum(v12[1:5,1])
v12[4,4] = v12[4,1]/sum(v12[1:5,1])
v12[5,4] = v12[5,1]/sum(v12[1:5,1])


v12[6,4] = v12[6,1]/sum(v12[6:10,1])
v12[7,4] = v12[7,1]/sum(v12[6:10,1])
v12[8,4] = v12[8,1]/sum(v12[6:10,1])
v12[9,4] = v12[9,1]/sum(v12[6:10,1])
v12[10,4] = v12[10,1]/sum(v12[6:10,1])


v12 %>% group_by(group) %>% summarise(sum= sum(var_import_perc))
library(ggplot2)
ggplot(v12, aes(fill = Var, y = var_import_perc*100, x = group)) + 
  geom_bar(position = "stack", stat = "identity") +
  ggtitle("Percentage of variable importance in SDM") +
  theme_minimal() +
  #scale_fill_brewer(palette = "Pastel2", name = "Variable") +
  scale_fill_manual(values = c("#B3CDE3","#DECBE4", "#CCEBC5","#FED9A6", "#FBB4AE"), name = "Variable") +
  #facet_wrap(~group) +
  xlab("Ecotype") +
  ylab("Percentage of importance of variables") +
  theme(plot.title=element_text(family = "Lucida Grande", face = "bold", size=13,hjust=0.5,vjust=1), 
        legend.position = "right", 
        strip.text = element_text(colour = "black", size = 11, family = "Lucida Grande"),
        legend.title = element_text(colour = "black", size = 11, family = "Lucida Grande"),
        legend.text = element_text(colour = "black", size = 11, family = "Lucida Grande"),
        axis.title.y = element_text(colour = 'black',  size = 11,vjust = 2, family = "Lucida Grande", angle=90), 
        axis.title.x = element_text(colour = "black", size = 10, vjust = -1,family = "Lucida Grande"),
        axis.text.x = element_text(colour = "black", size = 11, family = "Lucida Grande"),
        axis.text.y = element_text(colour = "black", size = 11, family = "Lucida Grande"))


#Individual model response plot

bm_PlotResponseCurves(bm.out = myBiomodModelOut1,
                      models.chosen = 'all',
                      new.env = get_formal_data(myBiomodModelOut1, "expl.var"),
                      show.variables = get_formal_data(myBiomodModelOut1, "expl.var.names"),
                      fixed.var = "median")

  # or only one algorithm
bm_PlotResponseCurves(bm.out = myBiomodModelOut1,
                      models.chosen = get_built_models(myBiomodModelOut1, algo = "RF"),
                      new.env = get_formal_data(myBiomodModelOut1, "expl.var"),
                      show.variables = get_formal_data(myBiomodModelOut1, "expl.var.names"),
                      fixed.var = "mean")



####### STEP 2: RUN ENSEMBLE models #######
# BIOMOD_EnsembleModeling combines individual models to build some kind of meta-model.
# here decide to exclude models with TSS(AUC) less than 0.8
myBiomodEM1 = BIOMOD_EnsembleModeling(
  myBiomodModelOut1,
  models.chosen = 'all',
  em.by='all',
  metric.select = c("TSS","ROC"),
  metric.eval = c("TSS", "ROC"),
  metric.select.thresh = c(0.8,0.8),
  em.algo = c('EMmean', 'EMcv', 'EMmedian', 'EMca', 'EMwmean'),
  var.import = 0)

myBiomodEM2 = BIOMOD_EnsembleModeling(
  myBiomodModelOut2,
  models.chosen = 'all',
  em.by='all',
  metric.select = c("TSS","ROC"),
  metric.eval = c("TSS", "ROC"),
  metric.select.thresh = c(0.8,0.8),
  em.algo = c('EMmean', 'EMcv', 'EMmedian', 'EMca', 'EMwmean'),
  var.import = 0)

#get evaluation scores
myBiomodEMEval1 = get_evaluations(myBiomodEM1)


bm_PlotEvalMean(bm.out = myBiomodEM1, 
                metric.eval = c("ROC","TSS"))

bm_PlotEvalBoxplot(bm.out = myBiomodEM1, 
                   group.by = c('algo','merged.by.PA'))
bm_PlotEvalBoxplot(bm.out = myBiomodEM1, 
                   group.by = c('algo','merged.by.run'))



#### STEP 3: MODEL PROJECTION - CURRENT CLIMATE ######


###MODEL PROJECTION
###prediction under current conditions
myBiomodProj1 = BIOMOD_Projection(
  myBiomodModelOut1,
  new.env = myExpl,
  new.env.xy = myRespXY1,
  proj.name = 'current',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

myBiomodProj2 = BIOMOD_Projection(
  myBiomodModelOut2,
  new.env = myExpl,
  new.env.xy = myRespXY2,
  proj.name = 'current',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

myBiomodEnsembleProj1 = BIOMOD_EnsembleForecasting(
  myBiomodEM1,
  bm.proj = myBiomodProj1,
  proj.name = 'current',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

myBiomodEnsembleProj2 = BIOMOD_EnsembleForecasting(
  myBiomodEM2,
  bm.proj = myBiomodProj2,
  proj.name = 'current',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')


# only plotting the committee averaging and mean weight here

biomod2::plot(myBiomodProj1, str.grep = 'EMwmean')
plot(myBiomodEnsembleProj1,str.grep = 'EMwmean', main = "Eco1")
plot(myBiomodEnsembleProj2,str.grep = 'EMWmean', main = "Eco2")

#### do not worry, the final plot is yet to come ####


####### STEP 4: PREDICTION UNDER FUTURE CLIMATE #######

###prediction under future conditions
#get future data
myExplFut = stack(biofut245$bio3 ,biofut245$bio8, biofut245$bio9, biofut245$bio15, biofut245$bio19)

myBiomodProjFuture1 = BIOMOD_Projection(
  myBiomodModelOut1,
  new.env = myExplFut,
  new.env.xy = myRespXY1,
  proj.name = 'future',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

myBiomodProjFuture2 = BIOMOD_Projection(
  myBiomodModelOut2,
  new.env = myExplFut,
  new.env.xy = myRespXY2,
  proj.name = 'future',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

myBiomodEMProjFuture1 = BIOMOD_EnsembleForecasting(
  myBiomodEM1,
  bm.proj = myBiomodProjFuture1,
  proj.name = 'future',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = T,
  build.clamping.mask = T,
  output.format = '.img')

myBiomodEMProjFuture2 = BIOMOD_EnsembleForecasting(
  myBiomodEM2,
  bm.proj = myBiomodProjFuture2,
  proj.name = 'future',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = T,
  build.clamping.mask = T,
  output.format = '.img')

## check how projections looks like
biomod2::plot(myBiomodEMProjFuture1, str.grep = 'EMwmean')
biomod2::plot(myBiomodEMProjFuture2, str.grep = 'EMwmean')

### keep moving ####


####### STEP 4: PROJECTION UNDER PAST CLIMATE ######

###prediction under last glacial maximum
# get past data
# remember to check if you have to divide the unit!!
myExplPas = stack(biopas$bio3, biopas$bio8/10, biopas$bio9/10, biopas$bio15, biopas$bio19)
myExplPas = crop(myExplPas, eu_ext)
myExplPas = raster::stack(myExplPas)

## OPTIONAL: see how things are in the past
par(mfrow=c(1,3))

plot(myExplPas[[1]]/10, main = "Last Glacial Maximum")
plot(myExpl[[1]], main = "Current")
plot(myExplFut[[1]], main = "Future")

## or plot it with proper asethetics
library(rasterVis) 

# Define a color palette
myPalette = colorRampPalette(c("green4", "white", "orange3"))

levelplot(myExplPas[[1]],
          main = "Last Glacial Maximum",
          col.regions = myPalette,
          margin = FALSE,  # Remove extra margins
          colorkey = list(space = "bottom", # Position the legend at the bottom
                          labels = list(at = seq(-40, 40, by = 10)),
                          fontfamily = "Lucida Grande"), # legend
          #scales = list(draw = TRUE), # Draw axes
          scales = list(draw = TRUE, # Draw axes
                        tck = c(1,0), # Specify the tick length
                        col = "black", # Specify the color of axis ticks
                        fontfamily = "Lucida Grande"), # Specify the font family of axis ticks
         
          #xlab = "Longitude", ylab = "Latitude" # Label axes
          xlab = list(label = "Latitude", cex = 1.0, col = "black", fontfamily = "Lucida Grande"), # Format x-axis label
          ylab = list(label = "Longitude", cex = 1.0, col = "black", fontfamily = "Lucida Grande"),
          par.settings = list(background = list(col = "#F6F6E9")

))

#################################

myBiomodProjPast1 = BIOMOD_Projection(
  myBiomodModelOut1,
  new.env = myExplPas,
  new.env.xy = myRespXY1,
  proj.name = 'past',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

myBiomodProjPast2 = BIOMOD_Projection(
  myBiomodModelOut2,
  new.env = myExplPas,
  new.env.xy = myRespXY2,
  proj.name = 'past',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = 'xz',
  build.clamping.mask = F,
  output.format = '.img')

myBiomodEMProjPast1 = BIOMOD_EnsembleForecasting(
  myBiomodEM1,
  bm.proj = myBiomodProjPast1,
  proj.name = 'past',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = T,
  build.clamping.mask = T,
  output.format = '.img')

myBiomodEMProjPast2 = BIOMOD_EnsembleForecasting(
  myBiomodEM2,
  bm.proj = myBiomodProjPast2,
  proj.name = 'past',
  models.chosen = 'all',
  metric.binary = 'TSS',
  compress = T,
  build.clamping.mask = T,
  output.format = '.img')

## check how projections looks like
biomod2::plot(myBiomodProjPast1, str.grep = 'EMca')
biomod2::plot(myBiomodEMProjPast1, str.grep = 'EMwmean')
biomod2::plot(myBiomodEMProjPast2, str.grep = 'EMwmean')


############# STEP 5: SPECIES RANGE CHANGE calculation ##########

CurrentProj1 = get_predictions(myBiomodEnsembleProj1,
                              metric.binary = "TSS",
                              model.as.col = TRUE)
FutureProj1 = get_predictions(myBiomodEMProjFuture1,
                             metric.binary = "TSS",
                             model.as.col = TRUE)
PastProj1 = get_predictions(myBiomodEMProjPast1,
                              metric.binary = "TSS",
                              model.as.col = TRUE)

CurrentProj2 = get_predictions(myBiomodEnsembleProj2,
                               metric.binary = "TSS",
                               model.as.col = TRUE)
FutureProj2 = get_predictions(myBiomodEMProjFuture2,
                              metric.binary = "TSS",
                              model.as.col = TRUE)
PastProj2 = get_predictions(myBiomodEMProjPast2,
                            metric.binary = "TSS",
                            model.as.col = TRUE)

#when dimensions are not matching
#PastProj1_resampled = raster::resample(PastProj1, CurrentProj1)
#PastProj2_resampled = resample(PastProj2, CurrentProj2)

myBiomodRangeSize1.1 =
  BIOMOD_RangeSize(proj.current = PastProj1, 
                   proj.future = CurrentProj1)

myBiomodRangeSize1.2 =
  BIOMOD_RangeSize(proj.current = CurrentProj1, 
                   proj.future = FutureProj1)

myBiomodRangeSize2.1=
  BIOMOD_RangeSize(proj.current = PastProj2, 
                   proj.future = CurrentProj2 )

myBiomodRangeSize2.2=
  BIOMOD_RangeSize(proj.current = CurrentProj2, 
                   proj.future = FutureProj2)


plot(myBiomodRangeSize2.1$Diff.By.Pixel, col = c("salmon3","skyblue","lightgrey","#33A02C"))
plot(myBiomodRangeSize1.1$Diff.By.Pixel, col = c("salmon3","skyblue","lightgrey","#33A02C"))
plot(myBiomodRangeSize2.2$Diff.By.Pixel,  col = c("salmon3","skyblue","lightgrey","#33A02C"))
plot(myBiomodRangeSize1.2$Diff.By.Pixel, col = c("salmon3","skyblue","lightgrey","#33A02C"))

write.csv(myBiomodRangeSize1$Compt.By.Models, "sdm_src_perc_strongsdorm.csv")
write.csv(myBiomodRangeSize2$Compt.By.Models, "sdm_src_perc_weaksdorm.csv")

#-2: Represents areas where the species has completely disappeared (loss of habitat).
#-1: Represents areas where the species has lost some of its range (partial loss).
# 0: Represents areas where the species' range has remained stable (no change).
# 1: Represents areas where the species has gained new habitat (expansion).

############## STEP 6: MAKE the PIE CHART of loss gain stable #######

### Load necessary library
library(dplyr)
library(ggplot2)

### prepare the data - a bit manual work

# selecting EMwmean by ROC

a = as.data.frame(myBiomodRangeSize1.1$Compt.By.Models[8,])
b = as.data.frame(myBiomodRangeSize1.2$Compt.By.Models[8,])
c = as.data.frame(myBiomodRangeSize2.1$Compt.By.Models[8,])
d = as.data.frame(myBiomodRangeSize2.2$Compt.By.Models[8,])

abcd = cbind(a,b,c,d)
abcd$Metrics = rownames(abcd)

rownames(abcd) = NULL
names(abcd) = c("smalldiff_lgm", "smalldiff_fut", "largediff_lgm", "largediff_fut", "Metrics")
rm(a,b,c,d)

# reshape data 
library(tidyr)
library(dplyr)
library(stringr)

# Reshape the data
abcd_new = pivot_longer(abcd, cols = -Metrics, 
                        names_to = "ecotype", 
                        values_to = "value")

abcd_new = as.data.frame(abcd_new)

# check and change to factors when necessary
summary(abcd_new)
levels(abcd_new$Metrics)
abcd_new$Metrics = as.factor(abcd_new$Metrics)
abcd_new$ecotype = as.factor(abcd_new$ecotype)

abcd_new$Metrics = factor(abcd_new$Metrics, levels = c("Loss","Gain","Stable1", "Stable0"))
abcd_new = na.omit(abcd_new)
abcd_new1 = abcd_new %>%
  mutate(abcd_new = ifelse(grepl("^smalldiff", ecotype), "Strong dormancy", "Weak dormancy"))
names(abcd_new1)[4] = "Ecotype"
names(abcd_new1)[2] = "group"
abcd_new1$Ecotype = as.factor(abcd_new1$Ecotype)
abcd_new1 = abcd_new1 %>%
  mutate(Period = str_extract(group, "(lgm|fut)"))
abcd_new1$Period = factor(abcd_new1$Period,
                          levels = c("lgm", "fut"),
                          labels = c("LGM-Current", "Current-Future"))


# Calculate the percentages within each ecotype-period combination
abcd_new1$Metrics = gsub("Stable1", "Stable", abcd_new1$Metrics)
abcd_new1$Metrics = gsub("Stable0", "Unoccupied", abcd_new1$Metrics)

###### Aggregating the data ######

aggregated_data = abcd_new1 %>%
  group_by(Ecotype, Period) %>%
  mutate(Total = sum(value)) %>%
  mutate(Percentage = (value / Total) * 100)

aggregated_data$Metrics = factor(aggregated_data$Metrics, levels = c("Unoccupied","Stable","Loss","Gain"))
levels(aggregated_data$Metrics)

### Option2: ignore unoccupied one
aggregated_data2 = abcd_new2 %>%
  group_by(Ecotype, Period) %>%
  mutate(Total = sum(value)) %>%
  mutate(Percentage = (value / Total) * 100)
aggregated_data2$Metrics = factor(aggregated_data2$Metrics, levels = c("Stable","Loss","Gain"))

# PLOTTING
# with unoccupied area
ggplot(aggregated_data, aes(x = "", y = Percentage, fill = Metrics)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_grid(Ecotype ~ Period) +
  theme_minimal() +
  labs(title = "Pie Chart of Loss, Gain, and Stable by Period and Ecotype",
       fill = "Metrics") +
  scale_fill_manual(values =c("lightgrey","skyblue", "salmon3", "#33A02C")) +  # Adjust as desired
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),   # Remove axis text
    panel.grid = element_blank(),  # Remove gridlines
    strip.text = element_text(colour = "black", size = 8.5, family = "Lucida Grande"),
    plot.title = element_text(family = "Lucida Grande", face = "bold", size = 11, hjust = 0.5, vjust = 1),
    legend.position = "bottom"
  )


# without unoccupied area
ggplot(aggregated_data2, aes(x = "", y = Percentage, fill = Metrics)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  facet_grid(Ecotype ~ Period) +
  theme_minimal() +
  labs(title = "Pie Chart of Loss, Gain, and Stable by Period and Ecotype",
       fill = "Metrics") +
  scale_fill_manual(values =c("skyblue", "salmon3", "#33A02C")) +  # Adjust as desired
  theme(
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),   # Remove axis text
    panel.grid = element_blank(),  # Remove gridlines
    strip.text = element_text(colour = "black", size = 8.5, family = "Lucida Grande"),
    plot.title = element_text(family = "Lucida Grande", face = "bold", size = 11, hjust = 0.5, vjust = 1),
    legend.position = "bottom"
  )

  
#see and choose colour from R palette
library(RColorBrewer)
display.brewer.pal(n = 8, name = "Paired") # personal favourite and colourblind friendly
#heximal code
brewer.pal(n = 8, name = "Paired")


########### Before you go: BIOMDO2 syntax changed very often, don't panic if any command suddenly doesn't work. Keep calm and check the latest manual #########

######################################################### GOOD LUCK!! #########################################################
