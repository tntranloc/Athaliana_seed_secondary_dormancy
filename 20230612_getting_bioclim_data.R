library(sp)
library(raster)
library(geodata)
library(terra)

###### 1. load your sample info ######
samples = read.csv("/path/to/your/file.csv", header = T)
head(samples)

#extract your coordinate data
xy = samples[,c("Longitude", "Latitude")] 
#create spatial point dataframe
spdf = SpatialPointsDataFrame(xy, samples, 
                              proj4string = CRS("+proj=longlat +datum=WGS84 +ELLPS=WGS84 +towgs84=0,0,0"))

#check - you should have an extra column as first column now named "coordinates", 
#with longitude and latitude values
head(spdf)

#Use this spatial point object to extract climatic variables
biodata = worldclim_global(var = "bio", res = 10, "/path/to/where/these/data/are/downloaded/to")
biodata # inspect the data
#cleaning up the names of 19 variables 
bionames = sprintf("bio%d", 1:19)

biodata_extract = extract(biodata[[1:19]], xy, df = T)
summary(biodata_extract) #you should have bio1 - bio19 columns now
#attach it to the original df
samples_bio = cbind(samples, biodata_extract)
#work with samples_bio now

######################## done ###################################


