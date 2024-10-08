library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)


#make an sf object for geom_sf()
eu = ne_countries(scale ="medium", continent = "europe", type = "countries",returnclass="sf")


library(maptools)
library(ggplot2)
library(raster)
library(RColorBrewer)
library(ggspatial)

#The dataframe for geom_point(data=...) should have columns Longitude and Latitude. The rest doesn't matter

theme_set(theme_bw())
g = ggplot(data = eu)+
  geom_sf() +
  geom_point(data=yourdata,aes(x=Longitude,y=Latitude),pch=20,
             fill="forestgreen",col="forestgreen",size=2)+
  ##name axes
  labs(title="Sampling sets of ... accessions")+
  xlab("Longitude") + ylab("Latitude")+
  ##set font size
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10))+
  ##Insert compass and scale
  annotation_scale(location = "bl", width_hint = 0.5) + #set the map scale
  annotation_north_arrow( 
    location = "br", which_north = "true", 
    pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm"),
    style = north_arrow_fancy_orienteering)+. #set the compass, locations can be br, bl, tr, or tl
  #set map limits
  xlim(-10,40) + #longitude range
  ylim(34,70) + #latitude range
  theme(panel.grid.major = element_line(color = gray(.5), 
                                        linetype = "dashed", size = 0.5), 
        panel.background = element_rect(fill = "aliceblue")) + #preference for the map background
  theme(plot.title=element_text(family = "Optima", face = "bold", size=13,hjust=0.5,vjust=1), 
        axis.title.y = element_text(colour = 'black',  size = 11, family = "Optima", face= "bold", angle=90), 
        axis.title.x = element_text(colour = "black", size = 11, family = "Optima", face = "bold"),
        axis.text.x = element_text(colour = "black", size = 11, family = "Optima"),
        axis.text.y = element_text(colour = "black", size = 11, family = "Optima")) #for the font preference


g
