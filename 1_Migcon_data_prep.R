
# 0. House keeping ####


setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")
devtools::install_github("SMBC-NZP/MigConnectivity")

pacman::p_load("tidyverse", "sp", "rgeos", "plyr", "rgdal", "tidylog",
                "sf", "kableExtra", "ggthemes", "rasterVis", "viridis") 
               
# source functions
source("scripts/migcon_AUX.R")

# 1. Load and explore data ####

## representative colonies ####
boot_out_summary <- read.csv("embarcadero_SDMs/data/colony_bootstrap_summary.csv",
                             row.names = NULL)

representative_colonies <- boot_out_summary %>%
  dplyr::filter(SampleSize > 4 & Representativity > 70) %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) 

  

all_complete <- read.csv("data/complete_fenology.csv", row.names = NULL)

all_complete %<>%
  dplyr::filter(unique_ID_Bird %!in% c("L38927_first", "6195227_2014",
                                "6143080_2012")) %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) %>% 
  dplyr::filter(Colony %in% representative_colonies$Colony) %>% 
  dplyr::select(Species, Colony, unique_ID_Bird, Date_Time, Lat, Lon, fen) 

## projections ####
load("data/projections.Rdata")

# worldmap
world_simple <- rworldmap::getMap(resolution = "high")
world_simple <- gBuffer(world_simple, byid = F, width = 0)
world_simplekM <- spTransform(world_simple, CRS(projections$EquidistConic))
world_simple_sf <- st_as_sf(world_simple)

## breeding colonies data ####
breeding_pops <- read.csv("data/colonies_info_100km_for_ArcGis.csv", 
                          row.names = NULL)

breeding_info <- breeding_pops %>% 
  filter(!is.na(Pop.Min)) %>% # remove colonies without pop.size
  ddply(~Radius100, summarize, Tot.pairs = sum(Pop.Min)) %>% # sum pop.size of all colonies inside the radius ("population")
  full_join(breeding_pops, by = "Radius100") %>%  # merge the sum with the original dataset
  filter(!is.na(Sample.Size)) %>% 
  rename_at(c("Radius100", "SampledColony"), ~c("Population", "Colony")) %>% 
  mutate(Tot.pop = Tot.pairs*2)  %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) 


all_complete %<>% 
  left_join(breeding_info, by = c("Colony", "Species")) %>% 
  droplevels()

## wintering polygons ####
win_areas <- read_sf(dsn = "data/shapefiles", 
                     layer = "wintering_areas_1")

## breeding by sp ####
sp_breeding <- read_sf(dsn = "data/shapefiles", layer = "sp_breeding")

## breeding populations ####
bre_pops <- read_sf(dsn = "data/shapefiles", layer =
                      "locs_CopyFeatures_Buffer1")

names(bre_pops)[1:4] <- c("Population", "Longitude", "Latitude", "Species")
bre_pops$Area <- c("IberianCoast", "Alboran", "Azores", "Azores", 
                   "Macaronesia", "Macaronesia", "IberianCoast", "Alboran", 
                   "Macaronesia", "Azores", "WestMediterranean", 
                   "WestMediterranean", "CentralMediterranean",
                   "WestMediterranean",
                   "CentralMediterranean", "WestMediterranean", 
                   "EastMediterranean", "EastMediterranean",
                   "EastMediterranean",
                   "CapeVerdeIslands", "CapeVerdeIslands")

bre_colonies <- read_sf(dsn = "data/shapefiles", layer = "Colonies_buffers")
 
#' We are going to use the package MigConnectivity and follow 
#' 
#' >Cohen EB, Hostetler JA, Hallworth MT, Rushing CS, Sillett TS, Marra PP. 2017 Quantifying the strength of migratory connectivity. Methods Ecol. Evol. , 1-12. (doi: https://doi.org/10.1111/2041-210X.12916)
#' 
#' First, we create a list and use it to store all the elements needed for the estimation of migratory connectivity, and then we calculate it. We are going to calculate it for each species at a Colony, Population and Area level, and then for the three species together at a species level (so we'll obtain a measure of spatial segregation between species). We can see a summary of colonies, populations, areas and species here: 

breeding_info %>%
  kable(digits = 2) %>%
  kable_styling()

ggplot() + 
  geom_sf(data = world_simple_sf, col = "darkgray", fill = "lightgray") + 
  geom_sf(data = bre_colonies, 
          aes(col = Sp, fill = Sp)) + 
  scale_fill_tableau() + 
  scale_color_tableau() + 
  coord_sf(xlim = st_bbox(bre_colonies)[c(1,3)], 
           ylim = st_bbox(bre_colonies)[c(2,4)], 
           expand = T) + 
  theme_bw()

# 2.  Create lists and store data necessary for the analysis ####

## For _Calonectris diomedea_ ####
 
### Distance matrix between breeding locations ####
#' 
#' **WARNING: coordinates MUST BE in lat/lon** 
#' 

CALDIO <- list()

#### For colonies ####
brCD <- breeding_info %>% 
  filter(Species == "CALDIO") %>% 
  filter(Colony != "Chafarinas") %>% 
  filter(Colony %in% representative_colonies$Colony) 

locsCD <- brCD %>% 
  dplyr::select(Colony, Longitude, Latitude) %>% 
  dplyr::rename_at(c("Colony", "Longitude", "Latitude"), 
                   ~c("name", "lon", "lat")) 

distBR_CD_col <- GeoDistanceInMetresMatrix(locsCD)

levelplot((distBR_CD_col)/1000, col.regions = viridis(100), 
          contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

 
#### For populations  ####
bre_pops_CD <- as_Spatial(bre_pops[bre_pops$Species == "CALDIO",])
bre_pops_CD <- bre_pops_CD[bre_pops_CD$Population != "NWAfrica",]

centroidsPop <- data.frame(name = bre_pops_CD$Population,
                           lon = rgeos::gCentroid(bre_pops_CD, 
                                                  byid = T)@coords[,1],
                           lat = rgeos::gCentroid(bre_pops_CD, 
                                                  byid = T)@coords[,2])


distBR_CD_pop <- GeoDistanceInMetresMatrix(centroidsPop)
levelplot((distBR_CD_pop)/1000, col.regions = viridis(100), 
          contour = F, margin = F, 
          scales = list(x = list(rot = 45)))


 
#### For areas ####

bre_pops_sf <- st_as_sf(bre_pops_CD)

br_areas_sf_merged <- bre_pops_sf %>% 
  group_by(Area) %>%
  summarise_all(first) %>% 
  dplyr::select("Area")

bre_pops_AREA <- as(br_areas_sf_merged, "Spatial")

centroidsAr <- data.frame(name = bre_pops_AREA$Area, 
                          lon = rgeos::gCentroid(bre_pops_AREA,
                                                 byid = T)@coords[,1],
                          lat = rgeos::gCentroid(bre_pops_AREA, 
                                                 byid = T)@coords[,2])

distBR_CD_Ar <- GeoDistanceInMetresMatrix(centroidsAr)

levelplot((distBR_CD_Ar)/1000, col.regions = viridis(100), 
          contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

distBR_CD <- list(distBR_CD_col, distBR_CD_pop, distBR_CD_Ar)
names(distBR_CD) <- c("Colony", "Population", "Area")
CALDIO[1] <- list(distBR_CD)
names(CALDIO) <- c("originDist")

### Relative abundance of the breeding sites ####

#' This will be a numeric vector of length equal to number of breeding sites, and that adds up to 1
#' 
#### For colonies ####

colAB_CD <- brCD$Tot.pop/sum(brCD$Tot.pop)

names(colAB_CD) <- brCD$Colony

relAB <- list()
relAB[1] <- list(colAB_CD)

#### For populations ####

popCD <- ddply(brCD, ~ Population, summarize, 
               Tot.pop = sum(Tot.pop))

popAB_CD <- popCD$Tot.pop/sum(popCD$Tot.pop)
names(popAB_CD) <- popCD$Population
relAB[2] <- list(popAB_CD)

#### For areas ####

areaCD <- ddply(brCD, ~ Area, summarize, 
                Tot.pop = sum(Tot.pop))
areaAB_CD <- areaCD$Tot.pop/sum(areaCD$Tot.pop)
names(areaAB_CD) <- areaCD$Area
relAB[3] <- list(areaAB_CD)
names(relAB) <- c("Colony", "Population", "Area")

CALDIO[2] <- list(relAB)

names(CALDIO)[2] <- "originRelAbund"

### Breeding sites ####

#' This is the geographic definition of sites in the "release" (breeding) season. Must be a spatial polygons object
#' 
#### For colonies ####
bre_coloniesCD <- as_Spatial(bre_colonies[bre_colonies$Sp == "CALDIO",])
bre_coloniesCD <- bre_coloniesCD[bre_coloniesCD$Colony != "Chafarinas",]
bre_coloniesCD <- bre_coloniesCD[bre_coloniesCD$Colony %in%
                                   representative_colonies$Colony,]
bre_coloniesCD@data <- droplevels(bre_coloniesCD@data)

orSites <- list()
orSites[1] <- list(bre_coloniesCD)

#### For populations ####

bre_popsCD <- as_Spatial(bre_pops[bre_pops$Species == "CALDIO",])
bre_popsCD <- bre_popsCD[bre_popsCD$Population != "NWAfrica",]
bre_popsCD@data <- droplevels(bre_popsCD@data)

orSites[2] <- list(bre_popsCD)

#### For areas ####

arPols_CD <- bre_pops_AREA

orSites[3] <- list(arPols_CD)
names(orSites) <- c("Colony", "Population", "Area")
CALDIO[3] <- list(orSites)

names(CALDIO)[3] <- "originSites"

### Assignment of individuals to breeding sites ####

#' Vector of length = n animals tracked assigning each one of them to a breeding site
#' 
#' We'll generate a dataframe with columns Colony, Population and Area and use the appropriate one at each level. 

Origins_CD <- all_complete %>% 
  filter(Species == "CALDIO") %>% 
  filter(Colony %in% representative_colonies$Colony) %>% 
  filter(Colony != "Chafarinas") %>% 
  dplyr::select(unique_ID_Bird, Colony , Population , Area) 

Origins_CD <- Origins_CD[!duplicated(Origins_CD),]

Origins_CD <- Origins_CD %>% 
  dplyr::arrange(unique_ID_Bird)

CD_BRAssignment <- list(as.character(Origins_CD$Colony),
                        as.character(Origins_CD$Population),
                        as.character(Origins_CD$Area))

names(CD_BRAssignment) <- c("Colony", "Population", "Area")

CALDIO[4] <- list(CD_BRAssignment)
names(CALDIO)[4] <- "originAssignment"

### Wintering (target) sites ####

#' Geographic definition of the sites in the non-release (wintering) season. These are the same whatever the level of organisation of breeding sites (Colony/Population/Area). However we must remove, from the object with all the possible areas, those not used by the species at hand. If there are "empty" wintering areas the function crashes. Must be a SpatialPolygons object
#' 

targetSites <- as_Spatial(win_areas[win_areas$ECOREGION %!in% 
                                      c("Agulhas_Benguela_Confluence", 
                                        "Namaqua",
                                        "North_Atlantic","Central_South_Atlantic",
                                        "PatagonianShelf_Falklands",
                                        "Uruguay_BuenosAires_Shelf",
                                        "SouthBrazilCurrent_RioGrande",
                                        "EasternBrazil"),])

targetSites@data$ECOREGION <- targetSites@data$ECOREGION
targetSites@data <- targetSites@data[,c(1,2)]

CALDIO[5] <- list(targetSites)
names(CALDIO)[5] <- "targetSites"

### Distance matrix between wintering locations ####

centroids <- data.frame(name = as.character(targetSites$ECOREGION),
                        lon = rgeos::gCentroid(targetSites, byid = T)@coords[,1],
                        lat = rgeos::gCentroid(targetSites, byid = T)@coords[,2])


distNBR <- GeoDistanceInMetresMatrix(centroids)
levelplot((distNBR)/1000, col.regions = viridis(100), contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

CALDIO[6] <- list(distNBR)
names(CALDIO)[6] <- "targetDist"

### Breeding points ####

#' SpatialPoints object with length = n animals tracked. Each point are the coordinates of the release (breeding) location (Colony/Population/Area) 
#' 

all_birds <- all_complete %>%
  filter(Species == "CALDIO") %>% 
  filter(Colony %in% representative_colonies$Colony) %>% 
  filter(Colony != "Chafarinas") %>% 
  dplyr::select(c("Colony", "unique_ID_Bird"))

all_birds <- all_birds[!duplicated(all_birds),]

pop_coordsCB <- breeding_info %>%
  filter(Species == "CALDIO") %>%
  dplyr::select(c("Colony", "Longitude", "Latitude"))

originPoints <- left_join(all_birds, pop_coordsCB)
originPoints$Colony <- as.factor(originPoints$Colony)

# important to order everytime so the order of these points corresponds to the order of the names in originAssignment
originPoints <- originPoints %>% 
  arrange(unique_ID_Bird)

originPoints <- SpatialPoints(coords = originPoints[,c(3,4)], 
                              proj4string = CRS(projections$WGS84))

CALDIO[7] <- list(originPoints)
names(CALDIO)[7] <- "originPoints"

 ### Wintering points ####

#' SpatialPoints object with length = n animals tracked. Each point are the coordinates of the estimated location in non-release (wintering) season. For our data we'll calculate the centroid of the 5% kernel of non-breeding locations
#' 
#' The function centroid_calculate is in the auxiliary script. As it takes a lot of time to run, we have run it once and stored the output in the Centroids.CALDIO object.  

## all_win <- all_complete %>%
##   filter(fen == "WIN" & Species == "CALDIO") %>%
##   filter(Colony %in% representative_colonies$Colony) %>%
##   filter(Colony != "Chafarinas") %>%
##   dplyr::select(c("unique_ID_Bird", "Population", "Colony", "Lon", "Lat", "fen"))
## 
## all_win <- droplevels(all_win)
## 
## centroids.CALDIO <- plyr::ddply(all_win, ~ unique_ID_Bird,
##                                 .fun = centroid_calculate,
##                                  Lon = "Lon", Lat = "Lat", group_ID = "Colony",
##                                 ind_ID = "unique_ID_Bird")
## 
## # save(centroids.CALDIO, file = "data/Centroids.CALDIO.2023.Rdata")

#' 
#' We use the generated object to generate the targetPoints object

load("data/Centroids.CALDIO.2023.RData")
centroids.CALDIO <- centroids.CALDIO[!duplicated(centroids.CALDIO),]

targetPoints <- centroids.CALDIO %>% arrange(unique_ID_Bird)
head(targetPoints)

targetPoints <- SpatialPoints(coords = targetPoints[,c(3,4)], 
                              proj4string = CRS(projections$WGS84))

CALDIO[8] <- list(targetPoints)
names(CALDIO)[8] <- "targetPoints"

### Geolocator bias  ####

#' For GLS data, a vector of length 2 indicating expected bias in geographic location of targetPoints in the units of the projection. This comes from a model run in the auxiliary script, and gives the coefficients and the variance/covariance matrix for Latitude and Longitude. We have run the model to calculate geo bias using the calibration datasets, where devices were groundtruthed at a known location. 

load("data/geo.error.model.Rdata")

geoBias <- coef(geo.error.model)
geoVCov <- vcov(geo.error.model)

CALDIO[9] <- list(geoBias)
CALDIO[10] <- list(geoVCov)
names(CALDIO)[9:10] <- c("geoBias", "geoVCov")

### Check all data is correct ####

plot(world_simple, col = "lightgray", border = "darkgray",
          xlim = c(CALDIO$targetSites@bbox[1,]), 
     ylim = c(CALDIO$targetSites@bbox[2,1], CALDIO$originPoints@bbox[2,2] + 10))
plot(CALDIO$targetSites, add = T, border = "firebrick3")
plot(CALDIO$originSites$Population, add = T, border = "dodgerblue3")
plot(CALDIO$targetPoints, add = T, col = "firebrick2", pch = 20)
plot(CALDIO$originPoints, add = T, col = "dodgerblue2", pch = 18, cex = 0.7)

# 3. Save all ####

saveRDS(CALDIO, file = "data/migcon_ready_data/CALDIO.RDS")

## the process repeats for calbor and caledw exactly the same
### For _Calonectris borealis_
#' 
#' #### Distance matrix between breeding locations

CALBOR <- list()
# For colonies
brCB <- breeding_info %>% 
  filter(Colony %in% representative_colonies$Colony) %>% 
  filter(Species == "CALBOR")

locsCB <- brCB %>% 
  dplyr::select(Colony, Longitude, Latitude) %>% 
  dplyr::rename_at(c("Colony", "Longitude", "Latitude"), ~c("name", "lon", "lat")) 

distBR_CB_col <- GeoDistanceInMetresMatrix(locsCB)

levelplot((distBR_CB_col)/1000, col.regions = viridis(100), 
          contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

# Populations
bre_pops_CB <- as_Spatial(bre_pops[bre_pops$Species == "CALBOR",])
bre_pops_CB <- bre_pops_CB[bre_pops_CB$Population %!in% 
                             c("WestAzores", "Terreros", "Galicia"), ]
bre_pops_CB@data <- droplevels(bre_pops_CB@data)

centroidsPopCB <- data.frame(name = bre_pops_CB$Population,
                             lon = rgeos::gCentroid(bre_pops_CB, 
                                                    byid = T)@coords[,1],
                             lat = rgeos::gCentroid(bre_pops_CB, 
                                                    byid = T)@coords[,2])

distBR_CB_pop <- GeoDistanceInMetresMatrix(centroidsPopCB)

levelplot((distBR_CB_pop)/1000, col.regions = viridis(100), 
          contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

# Areas
# convert to sf objects to do the "disolve" and convert back to sp 
bre_pops_sf <- st_as_sf(bre_pops_CB)

br_areas_sf_merged <- bre_pops_sf %>% 
  filter(Species == "CALBOR") %>% 
  group_by(Area) %>%
  summarise_all(first) %>% 
  dplyr::select("Area")

br_areas_sf_merged$Area <- as.factor(br_areas_sf_merged$Area)
br_areas_sf_merged$Area <- droplevels(br_areas_sf_merged$Area)

bre_pops_AREA <- as(br_areas_sf_merged, "Spatial")

centroidsArCB <- data.frame(name = droplevels(bre_pops_AREA$Area), 
                            lon = rgeos::gCentroid(bre_pops_AREA, 
                                                   byid = T)@coords[,1],
                            lat = rgeos::gCentroid(bre_pops_AREA, 
                                                   byid = T)@coords[,2])

distBR_CB_Ar <- GeoDistanceInMetresMatrix(centroidsArCB)
levelplot((distBR_CB_Ar)/1000, col.regions = viridis(100), 
          contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

distBR_CB <- list(distBR_CB_col, distBR_CB_pop, distBR_CB_Ar)
names(distBR_CB) <- c("Colony", "Population", "Area")
CALBOR[1] <- list(distBR_CB)
names(CALBOR) <- c("originDist")

#' 
#' 
#' #### Relative abundance of the breeding sites
# for colony
colAB_CB <- brCB$Tot.pop/sum(brCB$Tot.pop)
format(sum(colAB_CB), nsmall = 5)

names(colAB_CB) <- brCB$Colony

relAB <- list()
relAB[1] <- list(colAB_CB)


# for pop
popCB <- ddply(brCB, ~ Population, summarize, 
               Tot.pop = sum(Tot.pop))

popAB_CB <- popCB$Tot.pop/sum(popCB$Tot.pop)
names(popAB_CB) <- popCB$Population
sum(popAB_CB)

relAB[2] <- list(popAB_CB)

#for area
areaCB <- ddply(brCB, ~ Area, summarize, 
                Tot.pop = sum(Tot.pop))
areaAB_CB <- areaCB$Tot.pop/sum(areaCB$Tot.pop)

sum(areaAB_CB)

relAB[3] <- list(areaAB_CB)
names(relAB) <- c("Colony", "Population", "Area")

CALBOR[2] <- list(relAB)

names(CALBOR)[2] <- "originRelAbund"

#' 
#' #### Breeding sites
# colonies
colPols_CB <- bre_colonies %>% 
  filter(Sp == "CALBOR") %>% 
  filter(Colony %in% representative_colonies$Colony) 

ggplot(colPols_CB) + geom_sf(aes(fill = Colony))
st_is_valid(colPols_CB)
st_intersection(colPols_CB)$n.overlaps # if this is all 1s we're good

colPols_CB <- as_Spatial(colPols_CB)

orSites <- list()
orSites[1] <- list(colPols_CB)

# populations
popPols_CB <- bre_pops %>% 
  filter(Species == "CALBOR") %>% 
  filter(Population %!in% c("WestAzores", "Terreros", "Galicia")) 

ggplot(popPols_CB) + geom_sf(aes(fill = Population))
st_is_valid(popPols_CB)
st_intersection(popPols_CB)$n.overlaps

popPols_CB <- as_Spatial(popPols_CB)

orSites[2] <- list(popPols_CB)

# areas
arPols_CB <- bre_pops_AREA

orSites[3] <- list(arPols_CB)
names(orSites) <- c("Colony", "Population", "Area")
CALBOR[3] <- list(orSites)

names(CALBOR)[3] <- "originSites"

#' 
#' #### Assignment of individuals to breeding sites
Origins_CB <- all_complete %>% 
  filter(fen == "WIN") %>% 
  filter(Species == "CALBOR") %>% 
  filter(Colony %in% representative_colonies$Colony) %>% 
  dplyr::select(c("unique_ID_Bird", "Colony", "Population", "Area"))

Origins_CB <- Origins_CB[!duplicated(Origins_CB),]

Origins_CB <- Origins_CB %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) %>% 
  dplyr::arrange(unique_ID_Bird)

nrow(Origins_CB)

CB_BRAssignment <- list(as.character(Origins_CB$Colony),
                        as.character(Origins_CB$Population),
                        as.character(Origins_CB$Area))

names(CB_BRAssignment) <- c("Colony", "Population", "Area")

CALBOR[4] <- list(CB_BRAssignment)
names(CALBOR)[4] <- "originAssignment"

#' 
#' #### Wintering (target) sites
targetSitesCB <- win_areas %>% 
  filter(ECOREGION %!in% c("Equatorial_Atlantic", "Gulf_of_Guinea_West")) %>% 
  dplyr::select(ECOREGION)
  
st_is_valid(targetSitesCB)
st_intersection(targetSitesCB)$n.overlaps # apparently a few of these overlap which is not good

sf_use_s2(FALSE) # this turns off spherical geometry (shrug, but it makes it work)

overlapping <- st_intersection(targetSitesCB)

nonoverlapping <- overlapping %>% 
  filter(n.overlaps < 2)

ggplot(nonoverlapping) + 
  geom_sf(aes(fill = factor(ECOREGION)))

st_is_valid(nonoverlapping)

targetSitesCB <- as_Spatial(nonoverlapping)
CALBOR[5] <- list(targetSitesCB)
names(CALBOR)[5] <- "targetSites"

#' 
#' #### Distance matrix between wintering locations
centroids <- coordinates(targetSitesCB)

centroids <- data.frame(name = as.character(targetSitesCB$ECOREGION),
                        lon = centroids[,1],
                        lat = centroids[,2])


distNBR <- GeoDistanceInMetresMatrix(centroids)
levelplot((distNBR)/1000, col.regions = viridis(100), contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

CALBOR[6] <- list(distNBR)
names(CALBOR)[6] <- "targetDist"

#' 
#' #### Breeding points
all_win <- all_complete %>%
  mutate(Colony = if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) %>% 
  filter(Species == "CALBOR") %>% 
  filter(fen == "WIN") %>% 
  filter(Colony %in% representative_colonies$Colony) %>% 
  dplyr::select(c("unique_ID_Bird", "Colony", "fen")) 

all_birds <- all_win %>%
  dplyr::select(unique_ID_Bird, Colony)

all_birds <- all_birds[!duplicated(all_birds),]

pop_coordsCB <- breeding_info %>%
  filter(Species %in% c("CALBOR")) %>%
  filter(Colony %in% representative_colonies$Colony) %>% 
  dplyr::select(c("Colony", "Longitude", "Latitude")) 


originPoints <- left_join(all_birds, pop_coordsCB)
originPoints$Colony <- as.factor(originPoints$Colony)
unique(originPoints$Colony)

# important to order everytime so the order of these points corresponds to the order
# of the names in originAssignment

originPoints <- originPoints %>% 
  arrange(unique_ID_Bird)

nrow(originPoints)
names(originPoints)

originPoints <- SpatialPoints(coords = originPoints[,c(3,4)], 
                              proj4string = CRS(projections$WGS84))

CALBOR[7] <- list(originPoints)
names(CALBOR)[7] <- "originPoints"

#' 
#' #### Wintering points
#' 
#' We don't run this because it takes long time, we load the object instead 
## all_win <- all_complete %>%
##   filter(fen == "WIN" & Species == "CALBOR") %>%
##   dplyr::select(c("unique_ID_Bird", "Population",
##                   "Colony", "Lon", "Lat", "fen")) %>%
##   mutate(Colony = if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony))
## 
## centroids.CALBOR <- plyr::ddply(all_win, ~ unique_ID_Bird,
##                                 .fun = centroid_calculate,
##                                 Lon = "Lon", Lat = "Lat", group_ID = "Colony",
##                                 ind_ID = "unique_ID_Bird")
## 
## # save(centroids.CALBOR, file = "data/Centroids.CALBOR.2023.Rdata")

#' 
load("data/Centroids.CALBOR.2023.RData")
nrow(centroids.CALBOR)

targetPoints <- centroids.CALBOR %>% arrange(unique_ID_Bird)
targetPoints <- SpatialPoints(coords = targetPoints[,c(3, 4)], 
                              proj4string = CRS(projections$WGS84))

CALBOR[8] <- list(targetPoints)
names(CALBOR)[8] <- "targetPoints"

#' 
#' #### Geolocator bias 
load("data/geo.error.model.Rdata")

geoBias <- coef(geo.error.model)
geoVCov <- vcov(geo.error.model)

CALBOR[9] <- list(geoBias)
CALBOR[10] <- list(geoVCov)
names(CALBOR)[9:10] <- c("geoBias", "geoVCov")

#' 
#' #### Check all data is correct
plot(world_simple, col = "lightgray", border = "darkgray",
     xlim = c(CALBOR$targetSites@bbox[1,]), 
     ylim = c(CALBOR$targetSites@bbox[2,1], CALBOR$originPoints@bbox[2,2] + 10))
plot(CALBOR$targetSites, add = T, border = "firebrick3")
plot(CALBOR$originSites$Population, add = T, border = "dodgerblue3")
plot(CALBOR$targetPoints, add = T, col = "firebrick2", pch = 20)
plot(CALBOR$originPoints, add = T, col = "dodgerblue2", pch = 20)

#' 
#' #### Save all
#' 
saveRDS(CALBOR, file = "data/migcon_ready_data/CALBOR.RDS")

#' 
#' ### For _Calonectris edwardsii_ 
#' 
#' #### Distance matrix between breeding locations
CALEDW <- list()

brCE <- breeding_info %>% 
  filter(Species == "CALEDW") 

locsCE <- brCE %>% 
  dplyr::select(Colony, Longitude, Latitude) %>% 
  dplyr::rename_at(c("Colony", "Longitude", "Latitude"), ~c("name", "lon", "lat"))

distBR_CE_col <- GeoDistanceInMetresMatrix(locsCE)


CALEDW[1] <- list(distBR_CE_col)
names(CALEDW)[1] <- c("originDist")

#' 
#' #### Relative abundance of the breeding sites
colAB_CE <- brCE$Tot.pop/sum(brCE$Tot.pop)
names(colAB_CE) <- brCE$Colony
CALEDW[2] <- list(colAB_CE)
names(CALEDW)[2] <- "originRelAbund"

#' 
#' #### Breeding sites
colPols_CE <- as_Spatial(bre_colonies[bre_colonies$Sp == "CALEDW",])
plot(colPols_CE)

CALEDW[3] <- list(colPols_CE)
names(CALEDW)[3] <- "originSites"

#' 
#' #### Assignment of individuals to breeding sites
Origins_CE <- all_complete %>% 
  dplyr::filter(all_complete$Sp == "CALEDW") %>% 
  dplyr::select(c("unique_ID_Bird", "Colony", "Population", "Area"))

Origins_CE <- Origins_CE[!duplicated(Origins_CE),]

Origins_CE <- Origins_CE %>% 
  dplyr::arrange(unique_ID_Bird)

CE_BRAssignment <- as.character(Origins_CE$Colony)

CALEDW[4] <- list(CE_BRAssignment)
names(CALEDW)[4] <- "originAssignment"

#' 
#' #### Wintering sites
targetSitesCE <- as_Spatial(win_areas[win_areas$ECOREGION %in%
                                        c("NorthBrazilCurrent_SEBrazil",
                                          "SouthBrazilCurrent_RioGrande",
                                          "Uruguay_BuenosAires_Shelf"),])
                            
targetSitesCE@data <- targetSitesCE@data[,c(1,2)]

CALEDW[5] <- list(targetSitesCE)
names(CALEDW)[5] <- "targetSites"

#' 
#' #### Distance matrix between wintering locations
centroids <- coordinates(targetSitesCE)

centroids <- data.frame(name = as.character(targetSitesCE$ECOREGION),
                        lon = centroids[,1],
                        lat = centroids[,2])


distNBR <- GeoDistanceInMetresMatrix(centroids)
CALEDW[6] <- list(distNBR)
names(CALEDW)[6] <- "targetDist"

#' 
#' #### Breeding points
all_birds <- all_complete %>%
  filter(Species == "CALEDW") %>% 
  dplyr::select(c("Colony", "unique_ID_Bird"))


all_birds <- all_birds[!duplicated(all_birds),]

pop_coordsCE <- breeding_info %>%
  filter(Species == "CALEDW") %>%
  dplyr::select(c("Colony", "Longitude", "Latitude"))

originPoints <- left_join(all_birds, pop_coordsCE)
originPoints$Colony <- as.factor(originPoints$Colony)

originPoints <- originPoints %>% 
  arrange(unique_ID_Bird)

originPoints <- SpatialPoints(coords = originPoints[,c(3,4)], 
                              proj4string = CRS(projections$WGS84))


CALEDW[7] <- list(originPoints)
names(CALEDW)[7] <- "originPoints"

#' 
#' #### Wintering points
## all_win <- all_complete %>%
##   filter(fen == "WIN" & Species == "CALEDW") %>%
##   dplyr::select(c("unique_ID_Bird", "Population", "Colony", "Lon", "Lat",
##                   "fen"))
## 
## centroids.CALEDW <- plyr::ddply(all_win, ~ unique_ID_Bird,
##                                 .fun = centroid_calculate,
##                                 Lon = "Lon", Lat = "Lat", group_ID = "Colony",
##                                 ind_ID = "unique_ID_Bird")
## 
## # save(centroids.CALEDW, file = "data/Centroids.CALEDW.2023.Rdata")

#' 
load("data/centroids.CALEDW.RData")

targetPoints <- centroids.CALEDW[!duplicated(centroids.CALEDW),]
targetPoints <- SpatialPoints(coords = targetPoints[,c(3,4)], 
                              proj4string = CRS(projections$WGS84))


CALEDW[8] <- list(targetPoints)
names(CALEDW)[8] <- "targetPoints"

#' 
#' #### Geolocator bias
load("data/geo.error.model.Rdata")

geoBias <- coef(geo.error.model)
geoVCov <- vcov(geo.error.model)

CALEDW[9] <- list(geoBias)
CALEDW[10] <- list(geoVCov)
names(CALEDW)[9:10] <- c("geoBias", "geoVCov")

#' 
#' #### Check all data is correct
plot(world_simple, col = "lightgray", border = "darkgray",
     xlim = c(CALBOR$targetSites@bbox[1,]), 
     ylim = c(CALBOR$targetSites@bbox[2,1],CALBOR$originPoints@bbox[2,2] + 10))
plot(CALEDW$targetSites, add = T, border = "firebrick3")
plot(CALEDW$originSites, add = T, border = "dodgerblue3")
plot(CALEDW$targetPoints, add = T, col = "firebrick2", pch = 20)
plot(CALEDW$originPoints, add = T, col = "dodgerblue2", pch = 20)

#' 
#' #### Save all
#' 
saveRDS(CALEDW, file = "data/migcon_ready_data/CALEDW.RDS")

#' 
#' 
#' ### For all three species together
#' 
#' #### Distance matrix between breeding locations
CALONECTRIS <- list()

bre_pops_sf <- st_as_sf(bre_pops)

br_areas_sf_merged <- bre_pops_sf %>% 
  filter(Population %!in% c("WestAzores", "Terreros", "Galicia", "NWAfrica")) %>%
  group_by(Species) %>% 
  summarise_all(first) %>% 
  dplyr::select(Species)

plot(br_areas_sf_merged)

bre_pops_sp <- as(br_areas_sf_merged, "Spatial")

locsALL <- data.frame(names = sp_breeding$ECOREGION,
                      lon  = rgeos::gCentroid(as_Spatial(sp_breeding),
                                              byid = T)@coords[,1],
                      lat  = rgeos::gCentroid(as_Spatial(sp_breeding), 
                                              byid = T)@coords[,2])

distBR <- GeoDistanceInMetresMatrix(locsALL)
levelplot((distBR)/1000, col.regions = viridis(100), contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

CALONECTRIS[1] <- list(distBR)
names(CALONECTRIS) <- c("originDist")

#' 
#' #### Relative abundance of the breeding sites
Abund <- ddply(breeding_info, ~Species, summarise, 
               Tot.pop = sum(Tot.pop/sum(breeding_info$Tot.pop)))

Ab_ALL <- Abund$Tot.pop
format(sum(Ab_ALL), nsmall = 20)

names(Ab_ALL) <- Abund$Species
CALONECTRIS[2] <- list(Ab_ALL)
names(CALONECTRIS)[2] <- "originRelAbund"

#' 
#' #### Breeding sites
# CALONECTRIS[3] <- list(bre_pops_sp)
st_intersection(sp_breeding)$n.overlaps # apparently a few of these overlap which is not good

overlapping <- st_intersection(sp_breeding)

nonoverlapping <- overlapping %>% 
  filter(n.overlaps < 2)

ggplot(nonoverlapping) + 
  geom_sf(aes(fill = factor(ECOREGION)))

CALONECTRIS[3] <- list(as_Spatial(nonoverlapping))

names(CALONECTRIS)[3] <- "originSites"

#' 
#' #### Assignment of individuals to breeding sites
cdchafa <- all_complete %>% 
  filter(Species == "CALDIO" & Colony == "Chafarinas") %>% 
  dplyr::select(unique_ID_Bird) %>%
  distinct(.)

cdchafa <- droplevels(cdchafa)
all_complete <- droplevels(all_complete)
cdchafa$unique_ID_Bird %in% all_complete$unique_ID_Bird

Origins_ALL <- all_complete %>% 
  filter(fen == "WIN") %>% 
  filter(Colony %in% representative_colonies$Colony) %>% 
  filter(unique_ID_Bird %!in% cdchafa$unique_ID_Bird) %>% 
  dplyr::select(unique_ID_Bird, Colony, Species) %>% 
  distinct(.) %>% 
  arrange(unique_ID_Bird)


ALL_BRAssignment <- as.character(Origins_ALL$Sp)

CALONECTRIS[4] <- list(ALL_BRAssignment)
names(CALONECTRIS)[4] <- "originAssignment"

#' 
#' #### Wintering sites

st_is_valid(win_areas)

st_intersection(win_areas)$n.overlaps # apparently a few of these overlap which is not good

sf_use_s2(FALSE) # this turns off spherical geometry (shrug, but it makes it work)

overlapping <- st_intersection(win_areas)

nonoverlapping <- overlapping %>% 
  filter(n.overlaps < 2)

ggplot(nonoverlapping) + 
  geom_sf(aes(fill = factor(ECOREGION)))

st_is_valid(nonoverlapping)

targetSites <- as_Spatial(nonoverlapping)

CALONECTRIS[5] <- list(targetSites)
names(CALONECTRIS)[5] <- "targetSites"

#' 
#' #### Distance matrix between wintering sites
centroids <- coordinates(targetSites)
centroids <- data.frame(name = as.character(targetSites$ECOREGION),
                        lon = centroids[,1],
                        lat = centroids[,2])

distNBR <- GeoDistanceInMetresMatrix(centroids)
levelplot((distNBR)/1000, col.regions = viridis(100), contour = F, margin = F, 
          scales = list(x = list(rot = 45)))

CALONECTRIS[6] <- list(distNBR)
names(CALONECTRIS)[6] <- "targetDist"

#' 
#' #### Breeding points
all_birds <- all_complete %>%
  filter(fen == "WIN") %>% 
  filter(Colony %in% representative_colonies$Colony) %>% 
  filter(unique_ID_Bird %!in% cdchafa$unique_ID_Bird) %>% 
  dplyr::select(unique_ID_Bird, Colony) %>% 
  distinct(.) %>% 
  mutate(Colony = if_else(Colony == "Monta\xf1aClara", 
                          "MontañaClara", Colony)) %>% 
  arrange(unique_ID_Bird) 
  

pop_coords <- breeding_info %>%
  dplyr::select(Colony, Longitude, Latitude)

pop_coords$Colony <- as.character(pop_coords$Colony)

pop_coords <- pop_coords[!duplicated(pop_coords),] #chafarinas is twice and otherwise all chafarinas birds appear twice

originPoints <- left_join(all_birds, pop_coords, by = "Colony")

originPoints <- originPoints %>% 
  arrange(unique_ID_Bird)
head(originPoints)

originPoints <- SpatialPoints(coords = originPoints[,c(3,4)], 
                              proj4string = CRS(projections$WGS84))
plot(originPoints)

CALONECTRIS[7] <- list(originPoints)
names(CALONECTRIS)[7] <- "originPoints"

#' 
#' #### Witnering points 
## all_win <- all_complete %>%
##   filter(fen == "WIN") %>%
##   dplyr::select(c("unique_ID_Bird",
##                   "Colony", "Lon", "Lat", "fen"))
## 
## 
## centroids.ALL <- plyr::ddply(all_win, ~ unique_ID_Bird,
##                              .fun = centroid_calculate,
##                              Lon = "Lon", Lat = "Lat", group_ID = "Colony",
##                              ind_ID = "unique_ID_Bird")
## 
## # save(centroids.ALL, file = "data/Centroids.ALL.2023.Rdata")

#' 
load("data/Centroids.ALL.Rdata")

targetPoints <- centroids.ALL %>% 
  filter(unique_ID_Bird %!in% cdchafa$unique_ID_Bird) %>% 
  arrange(unique_ID_Bird)

head(targetPoints)
targetPoints <- SpatialPoints(coords = targetPoints[,c(3,4)], 
                              proj4string = CRS(projections$WGS84))


CALONECTRIS[8] <- list(targetPoints)
names(CALONECTRIS)[8] <- "targetPoints"

#' 
#' #### Geolocator bias
load("data/geo.error.model.Rdata")

geoBias <- coef(geo.error.model)
geoVCov <- vcov(geo.error.model)

CALONECTRIS[9] <- list(geoBias)
CALONECTRIS[10] <- list(geoVCov)
names(CALONECTRIS)[9:10] <- c("geoBias", "geoVCov")

#' 
#' #### Check all data is correct
plot(world_simple, col = "lightgray", border = "darkgray",
     xlim = c(CALONECTRIS$targetSites@bbox[1,]), 
     ylim = c(CALONECTRIS$targetSites@bbox[2,1], 
              CALONECTRIS$originPoints@bbox[2,2] + 10))
plot(CALONECTRIS$targetSites, add = T, border = "firebrick3")
plot(CALONECTRIS$originSites, add = T, border = "dodgerblue3")
plot(CALONECTRIS$targetPoints, add = T, col = "firebrick2", pch = 20)
plot(CALONECTRIS$originPoints, add = T, col = "dodgerblue2", pch = 20)

#' 
#' #### Save all
saveRDS(CALONECTRIS, file = "data/migcon_ready_data/CALONECTRIS.RDS")

