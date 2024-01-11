
library(raster)
library(dplyr)
library(pastclim)
library(rnaturalearth)


#### FIGURE 3 - CLIMATE OF EA MSA ####

setwd("/Volumes/Lucy/Habitability_modelling")

points <- c(30,55,-9,20) # eastern Africa

setwd("/Volumes/Lucy/climate_data/final_msa_files")
sites <- read.csv("CB2_MSA.csv")
sites <- sites %>% distinct(MIN_AGE, MAX_AGE, .keep_all=TRUE)
mid_ages <- sites$MID_AGE #extract mid ages of sites

# RAW DATA
bio01_rasterbrick <- pastclim:::region_series(time_bp = seq(-320000, -21000, 1000), 
                                              bio_variables = "bio01", 
                                              dataset = 'Krapp2021', ext = points)
bio12_rasterbrick <- pastclim:::region_series(time_bp = seq(-320000, -21000, 1000), 
                                              bio_variables = "bio12", 
                                              dataset = 'Krapp2021', ext = points)
clims <- list(bio01 =  stack(bio01_rasterbrick$bio01),
              bio12 = stack(bio12_rasterbrick$bio12))
clims1 <- lapply(clims, function(x) crop(x, points))
clims2 <- lapply(clims1, function(x) stack(x))


clims3 <- list()
for(i in 1:length(clims2)){
  dataset <-clims2[[i]]
  names(dataset@layers) <- 320:21
  dataset@layers <- dataset@layers[order(as.numeric(names(dataset@layers)))]
  crs(dataset) <-  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  # dataset <- projectRaster(dataset, crs ="+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
  clims3[[i]] <- dataset
}

names(clims3) <- names(clims)

raster_list <- list()
for(i in 1:length(clims3)){
  var <- clims3[[i]]
  raster_stack <- stack(var)
  raster_list[[i]] <- raster_stack
}
names(raster_list) <- names(clims)


npp_stack <- pastclim:::region_series(time_bp = seq(-320000, -21000, 1000), 
                                      bio_variables = "npp", 
                                      dataset = 'Krapp2021', ext = points)

npp_stack <- stack(npp_stack$npp)

names(npp_stack@layers) <- 320:21
npp_stack@layers <- npp_stack@layers[order(as.numeric(names(npp_stack@layers)))]

crs(npp_stack) <-  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#npp_stack <- projectRaster(npp_stack, crs ="+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

npp_stack <- resample(npp_stack, clims3[[1]])

res(npp_stack) == res(clims3[[1]]) # check same resolution


# Altitude 

relief_rast <- pastclim:::download_etopo_subset(bio01_rasterbrick[[1]])
crs(relief_rast) <-  "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#relief_rast <- projectRaster(raster(relief_rast), crs ="+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")


# Water
setwd("/Volumes/Lucy/climate_data/gis_data/final_PhD_GIS_files")
library(rgdal)
library(rgeos)
af_lakes <- readOGR(dsn = "lakes", layer = "lakes")
af_lakes <-spTransform(af_lakes, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
af_lakes2 <- as(af_lakes, "SpatialLines") 

rivers10 <- readOGR(dsn=("10m-rivers-lake-centerlines"),
                    layer="ne_10m_rivers_lake_centerlines") 
rivers10 <- crop(rivers10, points)
rivers10 <-spTransform(rivers10, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

water_bodies <- gUnion(rivers10 , rivers10)
water_bodies <- gUnion(water_bodies, af_lakes2)

mean_bio12 <- calc(raster_list[['bio12']], fun=mean) #mean precipitation
mean_bio01 <- calc(raster_list[['bio01']], fun=mean) #mean temperature
mean_npp <- calc(npp_stack, fun=mean) #mean NPP

require("RColorBrewer")
require("rasterVis")
require("gridExtra")
require("scales")
require("latticeExtra")

sites_sp <- SpatialPointsDataFrame(sites[,c("E", "N")], sites, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
#sites_sp<- spTransform(sites_sp, CRSobj = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
sites$ID <- seq(1, nrow(sites))
  
colr <- colorRampPalette(brewer.pal(7, 'RdYlGn'))
colr2 <- colorRampPalette(rev(c(heat.colors(20))))
colr3 <- colorRampPalette(brewer.pal(7, 'RdBu'))
colr4 <- colorRampPalette(c(terrain.colors(20)))


p1 <- levelplot(mean_bio12, margin=FALSE,                       # suppress marginal graphics
                colorkey=list(
                  space='bottom',                   # plot legend at bottom
                  labels=list(font=4),
                  title = "mm"),      # legend ticks and labels     
                par.settings=list(
                  axis.line=list(col='transparent')), # suppress axes and legend outline
                scales=list(xlab=NULL,ylab=NULL,raw=FALSE),            # suppress axis labels
                col.regions=colr3 ,  # colour ramp 
                at = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300), # colour ramp breaks  
                main = "Total annual precipitation", # title
                maxpixels=ncell(mean_bio12)) +
  layer(sp.polygons(af_lakes, col = "white", fill = "white")) +
  layer(sp.polygons(water_bodies, col = "white", fill = "white"))+
layer(sp.points(sites_sp, pch = 20, col = "black"))


p2 <-levelplot(mean_bio01, margin=FALSE,                       # suppress marginal graphics
               colorkey=list(
                 space='bottom',                   # plot legend at bottom
                 labels=list(font=4),
                 title = "°C"),      # legend ticks and labels     
               par.settings=list(
                 axis.line=list(col='transparent')), # suppress axes and legend outline
               scales=list(xlab=NULL,ylab=NULL,raw=FALSE),            # suppress axis labels
               col.regions=colr2 ,  # colour ramp 
               at = c(5, 6,7, 8,9, 10, 11,12,13, 14, 15, 16,17,  18,19, 20,21, 22, 23, 24, 25, 26, 27), # colour ramp breaks  
               main = "Mean annual temperature", # title
               maxpixels=ncell(mean_bio01))   +
  layer(sp.polygons(af_lakes, col = "white", fill = "white")) +
  layer(sp.polygons(water_bodies, col = "white", fill = "white"))+
  layer(sp.points(sites_sp, pch = 20, col = "black"))


p3 <-levelplot(mean_npp, margin=FALSE,                       # suppress marginal graphics
               colorkey=list(
                 space='bottom',                   # plot legend at bottom
                 labels=list(font=4),
                 title = "gCm-2year-1"),      # legend ticks and labels     
               par.settings=list(
                 axis.line=list(col='transparent')), # suppress axes and legend outline
               scales=list(xlab=NULL,ylab=NULL,raw=FALSE),            # suppress axis labels
               col.regions=colr ,  # colour ramp 
               at = c(0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300), # colour ramp breaks  
               main = "Net primary productivity", # title
               maxpixels=ncell(mean_bio01))   +
  layer(sp.polygons(af_lakes, col = "white", fill = "white")) +
  layer(sp.polygons(water_bodies, col = "white", fill = "white"))+
  layer(sp.points(sites_sp, pch = 20, col = "black"))



p4 <- levelplot(relief_rast, margin=FALSE,                       # suppress marginal graphics
                colorkey=list(
                  space='bottom',                   # plot legend at bottom
                  labels=list(font=4),
                  title = "masl"),      # legend ticks and labels     
                par.settings=list(
                  axis.line=list(col='transparent')), # suppress axes and legend outline
                scales=list(xlab=NULL,ylab=NULL,raw=FALSE),            # suppress axis labels
                col.regions=colr4 ,  # colour ramp 
                at = c(-50,0,200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400), # colour ramp breaks  
                main = "Altitude", # title
                maxpixels=ncell(relief_rast)) +
layer(sp.polygons(af_lakes, col = "white", fill = "white")) +
  layer(sp.polygons(water_bodies, col = "white", fill = "white")) +
  layer(sp.points(sites_sp, pch = 20, col = "black"))

gridExtra:: grid.arrange(p1, p2, p3, p4, nrow = 1)


#### FIGURE 1 - MAP OF SPECIFIC VS GENERIC MSA ####
setwd("/Users/lucytimbrell/Documents/LT_documents/publications/BIEA_paper")
africa <- rnaturalearth::ne_countries(continent = 'africa') #africa shape
africa <- africa[-38,] #removes madagascar
africa <- terra::vect(africa) #to spatvector
africa <- terra::aggregate(africa, dissolve=T) #disolved internal divisions

alt  <- pastclim:::region_series(time_bp = 0,  bio_variables = "altitude",  dataset = 'Krapp2021', crop = africa)
raster(alt$altitude)
plot(alt$altitude)

data <- read.csv("MSA sites Figure 1.csv")

require("RColorBrewer")
require("rasterVis")
require("gridExtra")
require("scales")
require("latticeExtra")

sites_sp <- SpatialPointsDataFrame(data[,c("x..long.", "y..lat.")], data, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
#sites_sp<- spTransform(sites_sp, CRSobj = "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")


colr4 <- colorRampPalette(c(terrain.colors(10)))

ptsCols <- c("green", 'blue', "purple", "black","darkorange1",  "#E7298A", "red") 
ptspch <- c(20, 20, 20, 1, 20, 20, 20)

levelplot(raster(alt$altitude), margin = FALSE,                      # suppress marginal graphics
                colorkey= NULL,      # legend ticks and labels     
                par.settings=list(
                axis.line=list(col='transparent')),# suppress axes and legend outline
                scales=list(xlab=NULL,ylab=NULL,raw=FALSE),            # suppress axis labels
                col.regions= colr4  ,  # colour ramp 
                main = NULL, # title
                maxpixels=ncell(alt$altitude),
                key = list(space = 'right', text = list(levels(as.factor(data$X))), 
                     points = list(col = ptsCols, pch = ptspch, cex = 2))) + 
         layer(sp.points(sites_sp,col = 
                           ifelse(sites_sp$X == "Still Bay", "red",
                                  ifelse(sites_sp$X == "Howiesons Poort", "blue", 
                                         ifelse(sites_sp$X == "Pietersburg", "darkorange1",
                                                ifelse(sites_sp$X == "Sangoan", "#E7298A", 
                                                       ifelse(sites_sp$X == "MSA", "black", 
                                                              ifelse(sites_sp$X == "Lupemban", "purple", 
                                                                     ifelse(sites_sp$X == "Aterian", "green", "red"))))))),
                         pch = ifelse(sites_sp$X == "MSA", 1,20), cex = ifelse(sites_sp$X == "MSA", 1.5,3)))





#### FIGURE 5 - BIOMES ####

####load data####
points <- c(30,55,-9,20) # eastern Africa

africa <- rnaturalearth::ne_countries(continent = 'africa') #africa shape
africa <- africa[-38,] #removes madagascar
africa <- terra::vect(africa) #to spatvector
africa <- terra::aggregate(africa, dissolve=T) #disolved internal divisions

africa <- crop(africa, points)
biomes_rasterbrick <- pastclim:::region_series(time_bp = seq(-320000, -21000, 1000), 
                                               bio_variables = "biome", 
                                               dataset = 'Krapp2021', ext = points, crop = africa)

biomes <- data.frame(ID=1:28, biome <- c("Tropical evergreen forest",
                                         "Tropical semi-deciduous forest",
                                         "Tropical deciduous forest/woodland",
                                         "Temperate deciduous forest",
                                         "Temperate conifer forest",
                                         "Warm mixed forest",
                                         "Cool mixed forest",
                                         "Cool conifer forest",
                                         "Cold mixed forest",
                                         "Evergreen taiga/montane forest",
                                         "Deciduous taiga/montane forest",
                                         "Tropical savanna",
                                         "Tropical xerophytic shrubland",
                                         "Temperate xerophytic shrubland",
                                         "Temperate sclerophyll woodland",
                                         "Temperate broadleaved savanna",
                                         "Open conifer woodland",
                                         "Boreal parkland",
                                         "Tropical grassland",
                                         "Temperate grassland",
                                         "Desert",
                                         "Steppe tundra",
                                         "Shrub tundra",
                                         "Dwarf shrub tundra",
                                         "Prostrate shrub tundra",
                                         "Cushion forb lichen moss tundra",
                                         "Barren",
                                         "Land ice")) #biome names

biome_stack <- list()
for(i in 1:300){
  biomes_x <-biomes_rasterbrick$biome[[i]]
  biomes_x <- raster::ratify(raster(biomes_x))
  levels(biomes_x)<-biomes
  biome_stack[[i]] <- biomes_x
} #different data management required to attribute biome values


af_biome <- stack(biome_stack) #reorder
names(af_biome @layers) <- 320:21 #reorder
af_biome@layers <- af_biome @layers[order(as.numeric(names(af_biome@layers)))] #reorder from present to past

####modal biome####

af_biome_us <- unstack(af_biome)
df1 <- matrix(nrow = length(af_biome_us[[1]]), ncol=length(af_biome_us))

for(i in 1:length(af_biome_us)){df1[,i] <- values(af_biome_us[[i]])}
df2 <- data.frame(df1)

levels(as.factor(unlist(df2)))

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

df3 <- apply(df2, 1, function(x)getmode(x)) #unique values per row excluding NaNs
af_biome3 <- af_biome_us[[1]]
values(af_biome3) <- df3

col14 <- c("#1B5E20", #trop evergeen)
                    "#4CAF50", #trop decid
                    "#B2DFDB", #temp conifer teal100
                    "#FFD54F", #tropop xero shrub
                    "#CDDC39", #temp xero shrub lime 500
                    "#F0F4C3", #temp grass lime100
                    "#FFEB02", #trop grasslanf orange200
                    "red") #steppe tundra) #barren
                    
biome_names <- c('Tropical evergreen forest',
                 'Temperate conifer forest',
                 'Warm mixed forest',
                 'Tropical xerophytic shrubland',
                 'Temperate sclerophyll woodland',
                 'Open conifer woodland',
                 'Desert',
                 "Steppe Tundra")

modal_biome_95 <- af_biome3
modal_biome_95_2 <- modal_biome_95
mb95 <- values(modal_biome_95_2)
table(values(modal_biome_95_2)) # find ones to delete

mb95[is.nan(mb95)] <- NA
values(modal_biome_95_2) <- as.numeric(as.factor(mb95))

table(values(modal_biome_95_2)) # find ones to delete

plot(modal_biome_95_2, col=col14, legend=F, axes = F)
plot(africa, add=T, col="transparent", border = "grey50")
legend("bottomright", legend = biome_names, fill=col14, xpd=T)


p1 <- levelplot(modal_biome_95_2, margin=FALSE,                       # suppress marginal graphics
          colorkey=NULL,      # legend ticks and labels     
          par.settings=list(
          axis.line=list(col='transparent')), # suppress axes and legend outline
          scales=list(xlab=NULL,ylab=NULL,raw=FALSE),            # suppress axis labels
          col.regions=col14 ,  # colour ramp 
          main = "Modal biome", # title
          maxpixels=ncell(modal_biome_95_2),
          key = list(space = 'right', text = list(biome_names, cex = 1), 
                     points = list(fill = col14, pch = 22)))

#### biome count ####
bluescale <- colorspace::diverge_hcl(7, c=100, l=c(50,90), power=1)
grn_purp <- rev(colorspace::diverge_hcl(7, "Purple-Green"))

af_biome_us <- unstack(af_biome)
df1 <- matrix(nrow = length(af_biome_us[[1]]), ncol=length(af_biome_us))
for(i in 1:length(af_biome_us)){
  df1[,i] <- getValues(af_biome_us[[i]])
}
df2 <- data.frame(df1)
df3 <- apply(df2, 1, function(x)length(unique(x[!is.na(x)]))) #unique values per row excluding NaNs
af_biome2 <- af_biome[[1]]
values(af_biome2) <- df3

count_biome_95 <- af_biome2
values(count_biome_95)[values(count_biome_95)==0]<- NA
plot(count_biome_95, col=grn_purp,  legend = F, axes = FALSE)
plot(africa, add=T, col="transparent", border = "grey50")
legend("bottomright", legend = c(1:6),
       fill=grn_purp,
       #inset = c(-0.15,1),
       xpd=T)

p2 <- levelplot(count_biome_95, margin=FALSE,                       # suppress marginal graphics
                colorkey=NULL,      # legend ticks and labels     
                par.settings=list(
                axis.line=list(col='transparent')), # suppress axes and legend outline
                scales=list(xlab="",ylab="",raw=FALSE),            # suppress axis labels
                col.regions=grn_purp ,  # colour ramp 
                main = "Number of biomes over time", # title
                maxpixels=ncell(count_biome_95),
                key = list(space = 'right', text = list(c("1","2","3","4","5","6", "7"), cex = 1.5), 
                           points = list(fill = grn_purp, pch = 22, cex = 1.5)))



#### rate of biome change through time ####

grn_or <- colorspace::diverging_hcl(5, "Blue-Yellow2", power=0.75)

af_biome_chg <- list()
for(i in 1:(length(af_biome@layers)-1)){
  af_biome_chg[[i]] <- af_biome[[i]]-af_biome[[(i+1)]] # differences between consecutive timeslices
  af_biome_chg[[i]][af_biome_chg[[i]]>0] <- 1 #to binary for change or not
  af_biome_chg[[i]][af_biome_chg[[i]]<0] <- 1 #to binary for change or not
}
af_biome_chg_total <- sum(stack(af_biome_chg)) # adds all changes together, i.e. count of total changes between consecutive timeslices
af_biome_chg_prop <- (af_biome_chg_total/(length(af_biome@layers)-1))*100 #there are 55 changes between 56 datasets

change_biome_95 <- af_biome_chg_prop
plot(change_biome_95, col=grn_or, breaks = c(-0.01, 0.01, 6, 12, 24, 48), legend = FALSE, axes = FALSE)
plot(africa, add=T, col="transparent", border = "grey50")
legend("bottomright", legend = rev(c("No Change", "60-20kyrs", "20-10kyrs", "10-5kyrs", "<5kyrs")),
       fill=rev(grn_or),
       #inset = c(-0.15,1),
       xpd=T)

p3 <- levelplot(change_biome_95 , margin=FALSE,                       # suppress marginal graphics
                colorkey=NULL,      # legend ticks and labels     
                par.settings=list(
                  axis.line=list(col='transparent')), # suppress axes and legend outline
                scales=list(xlab=NULL,ylab=NULL,raw=FALSE),            # suppress axis labels
                col.regions=rev(grn_or),  # colour ramp 
                at = c(-0.01, 0.01, 6, 12, 24, 48),
                main = "Rate of biome flux", # title
                maxpixels=ncell(change_biome_95),
                key = list(space = 'right', text = list(rev(c("No Change", "60-20kyrs", "20-10kyrs", "10-5kyrs", "<5kyrs")), cex = 1.5), 
                           points = list(fill = rev(grn_or), pch = 22, cex = 1.5)))

####ecotonality####

af_biome6 <- af_biome5
values(af_biome6)[values(af_biome6)==0]<- NA

af_biome2 <- af_biome5
for(j in 1:af_biome6@data@nlayers){ #this takes ages; read from saved rasters
  for(i in 1:length(af_biome6[[j]])){
    cells_ad <- adjacent(af_biome6[[j]], i, directions=8, pairs=F)
    queen <- c(extract(af_biome6[[j]],cells_ad), extract(af_biome6[[j]],i))
    af_biome2[[j]][i] <- length(unique(na.omit(queen)))
    print(i)
  }}

setwd("/Users/lucytimbrell/Desktop/Refugia")
af_biome4 <- as.list(af_biome2)
lapply(seq_along(af_biome4), function(x) {writeRaster(af_biome4[[x]], paste0("df", x, ".tif"), format='GTiff', overwrite=T)})

#read from saved rasters; set path to chosen directory
rastlist <- list.files(path = "/Users/lucytimbrell/Desktop/Refugia", pattern='.tif', 
                       all.files=TRUE, full.names=FALSE)
af_biome2 <- list()
af_biome2 <-  lapply(rastlist, raster)
af_biome2 <- stack(af_biome2)

af_biome3 <- sum(af_biome2)/length(af_biome2@layers)
af_biome3 <- crop(af_biome3, af_map)#crop to africa
af_biome3 <- mask(af_biome3, af_map)#crop to africa

af_biome_95 <- af_biome3
values(af_biome_95)[values(af_biome_95)==0]<- NA

plot(af_biome_95, breaks = c(1,1.001, 1.5, 1.999, 2), col=col4, legend = F, axes = F)
plot(africa, add=T, col="transparent", border = "grey50")
legend("bottomright", legend = c("Homogeneuos", "Low Variability", "Impersistent Ecotone", "Persistent Ecotone"),
       fill=col4,
       #inset = c(-0.15,1),
       xpd=T)

p4 <- levelplot(af_biome_95 , margin=FALSE,                       # suppress marginal graphics
                colorkey=NULL,      # legend ticks and labels     
                par.settings=list(
                axis.line=list(col='transparent')), # suppress axes and legend outline
                scales=list(xlab=NULL,ylab=NULL,raw=FALSE),            # suppress axis labels
                col.regions=col4,  # colour ramp 
                at = c(1,1.001, 1.5, 1.999, 2),
                main = "Ecotones", # title
                maxpixels=ncell(af_biome_95),
                key = list(space = 'right', text = list(c("Homogeneuos", "Low Variability", "Impersistent Ecotone", "Persistent Ecotone"), cex = 1.5), 
                           points = list(col = col4, pch = 20, cex = 1.5)))


p1
p2
