 
# 0. House keeping ####

setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")
# devtools::install_github("SMBC-NZP/MigConnectivity")

pacman::p_load("MigConnectivity", "sf")

# source functions
source("scripts/migcon_AUX.R")

# 1. Load data packages from previous script ####

CALBOR <- readRDS("data/migcon_ready_data/CALBOR.RDS")
CALDIO <- readRDS("data/migcon_ready_data/CALDIO.RDS")
CALEDW <- readRDS("data/migcon_ready_data/CALEDW.RDS")
CALONECTRIS <- readRDS("data/migcon_ready_data/CALONECTRIS.RDS")
BORDIO <- readRDS("data/migcon_ready_data/BORDIO.RDS")


# 2. Run all functions ####

## 2.1 For Calonectris diomedea ####

### Colony ####

M_CD_Col <- estMC(isGL = T, # Logical vector: light-level geolocator(T)/GPS(F)
                  geoBias = CALDIO$geoBias, # Light-level geolocator bias
                  geoVCov = CALDIO$geoVCov, #Light-level geolocator cova matrix
                  nSim = 1000, # tunning parameter, leave 1000 for GLS data
                  
                  # Target location dist matrix
                  targetDist = CALDIO$targetDist, 
                  # Non-breeding/target sites
                  targetSites = st_as_sf(CALDIO$targetSites), 
                  # Target locs from GLS
                  targetPoints = st_as_sf(CALDIO$targetPoints), 
                  
                  #Origin loc dist matrix
                  originDist = CALDIO$originDist$Colony, 
                  # Breeding sites 
                  originSites = st_as_sf(CALDIO$originSites$Colony),
                  # Capture Locations 
                  originPoints = st_as_sf(CALDIO$originPoints),
                  # Origin loc rel abun
                  originRelAbund = CALDIO$originRelAbund$Colony, 
                  
                  # originAssignment = CALDIO$originAssignment$Colony, 
                  # resampleProjection = CALDIO$originPoints@proj4string, 
                  
                  maintainLegacyOutput = T,
                  verbose = 1,   # output options - see help ??estMC
                  nSamples = 1000, # This is set low for example 
                  maxTries = 300) 

saveRDS(M_CD_Col, file = "outputs/M_CD_Col.RDS")

### Population ####

M_CD_Pop <- estMC(isGL = T, # Logical vector: light-level geolocator(T)/GPS(F)
                  geoBias = CALDIO$geoBias, # Light-level geolocator bias
                  geoVCov = CALDIO$geoVCov, #Light-level geolocator cova matrix
                  nSim = 1000, # tunning parameter, leave 1000 for GLS data
                  
                  # Target location dist matrix
                  targetDist = CALDIO$targetDist, 
                  # Non-breeding/target sites
                  targetSites = st_as_sf(CALDIO$targetSites), 
                  # Target locs from GLS
                  targetPoints = st_as_sf(CALDIO$targetPoints), 
                  
                  #Origin loc dist matrix
                  originDist = CALDIO$originDist$Population, 
                  # Breeding sites 
                  originSites = st_as_sf(CALDIO$originSites$Population),
                  # Capture Locations 
                  originPoints = st_as_sf(CALDIO$originPoints),
                  # Origin loc rel abun
                  originRelAbund = CALDIO$originRelAbund$Population, 
                  
                  # originAssignment = CALDIO$originAssignment$Colony, 
                  # resampleProjection = CALDIO$originPoints@proj4string, 
                  
                  maintainLegacyOutput = T,
                  verbose = 1,   # output options - see help ??estMC
                  nSamples = 1000, # This is set low for example 
                  maxTries = 300) 

saveRDS(M_CD_Pop, file = "outputs/M_CD_Pop.RDS")

## 2.2. For Calonectris borealis ####

### Colony ####

M_CB_Col <- estMC(isGL = T, # Logical vector: light-level geolocator(T)/GPS(F)
                  geoBias = CALBOR$geoBias, # Light-level geolocator bias
                  geoVCov = CALBOR$geoVCov, #Light-level geolocator cova matrix
                  nSim = 1000, # tunning parameter, leave 1000 for GLS data
                  
                  # Target location dist matrix
                  targetDist = CALBOR$targetDist, 
                  # Non-breeding/target sites
                  targetSites = st_as_sf(CALBOR$targetSites),
                  # Target locs from GLS
                  targetPoints = st_as_sf(CALBOR$targetPoints), 
                  
                  #Origin loc dist matrix
                  originDist = CALBOR$originDist$Colony, 
                  # Breeding sites 
                  originSites = st_as_sf(CALBOR$originSites$Colony),
                  # Capture Locations 
                  originPoints = st_as_sf(CALBOR$originPoints),
                  # Origin loc rel abun
                  originRelAbund = CALBOR$originRelAbund$Colony, 
                  
                  # originAssignment = CALBOR$originAssignment$Colony, 
                  # resampleProjection = CALBOR$originPoints@proj4string, 
                  
                  maintainLegacyOutput = T,
                  verbose = 1,   # output options - see help ??estMC
                  nSamples = 1000, # This is set low for example 
                  maxTries = 300) 

saveRDS(M_CB_Col, file = "outputs/M_CB_Col.RDS")


### Population ####

M_CB_Pop <- estMC(isGL = T, # Logical vector: light-level geolocator(T)/GPS(F)
                  geoBias = CALBOR$geoBias, # Light-level geolocator bias
                  geoVCov = CALBOR$geoVCov, #Light-level geolocator cova matrix
                  nSim = 1000, # tunning parameter, leave 1000 for GLS data
                  
                  # Target location dist matrix
                  targetDist = CALBOR$targetDist, 
                  # Non-breeding/target sites
                  targetSites = st_as_sf(CALBOR$targetSites), 
                  # Target locs from GLS
                  targetPoints = st_as_sf(CALBOR$targetPoints), 
                  
                  #Origin loc dist matrix
                  originDist = CALBOR$originDist$Population, 
                  # Breeding sites 
                  originSites = st_as_sf(CALBOR$originSites$Population),
                  # Capture Locations 
                  originPoints = st_as_sf(CALBOR$originPoints),
                  # Origin loc rel abun
                  originRelAbund = CALBOR$originRelAbund$Population, 
                  
                  # originAssignment = CALBOR$originAssignment$Colony, 
                  # resampleProjection = CALBOR$originPoints@proj4string, 
                  
                  maintainLegacyOutput = T,
                  verbose = 1,   # output options - see help ??estMC
                  nSamples = 1000, # This is set low for example 
                  maxTries = 300) 

saveRDS(M_CB_Pop, file = "outputs/M_CB_Pop.RDS")

## 2.3 For Calonectris edwardsii ####
 
### Colony ####

M_CE_Col <- estMC(isGL = T, # Logical vector: light-level geolocator(T)/GPS(F)
                  geoBias = CALEDW$geoBias, # Light-level geolocator bias
                  geoVCov = CALEDW$geoVCov, #Light-level geolocator cova matrix
                  nSim = 1000, # tunning parameter, leave 1000 for GLS data
                  
                  # Target location dist matrix
                  targetDist = CALEDW$targetDist, 
                  # Non-breeding/target sites
                  targetSites = st_as_sf(CALEDW$targetSites),
                  # Target locs from GLS
                  targetPoints = st_as_sf(CALEDW$targetPoints), 
                  
                  #Origin loc dist matrix
                  originDist = CALEDW$originDist, 
                  # Breeding sites 
                  originSites = st_as_sf(CALEDW$originSites),
                  # Capture Locations 
                  originPoints = st_as_sf(CALEDW$originPoints),
                  # Origin loc rel abun
                  originRelAbund = CALEDW$originRelAbund, 
                  
                  # originAssignment = CALEDW$originAssignment$Colony, 
                  # resampleProjection = CALEDW$originPoints@proj4string, 
                  
                  maintainLegacyOutput = T,
                  verbose = 1,   # output options - see help ??estMC
                  nSamples = 1000, # This is set low for example 
                  maxTries = 300) 

saveRDS(M_CE_Col, file = "outputs/M_CE_Col.RDS")


## 2.4 For all three species together ####

M_all <- estMC(isGL = T, # Logical vector: light-level geolocator(T)/GPS(F)
               geoBias = CALONECTRIS$geoBias, # Light-level geolocator location bias
               geoVCov = CALONECTRIS$geoVCov, #Light-level geolocator covariance matrix
               nSim = 1000, # tunning parameter, leave 1000 for GLS data
               
               # Target location distance matrix
               targetDist = CALONECTRIS$targetDist, 
               # Non-breeding / target sites
               targetSites = st_as_sf(CALONECTRIS$targetSites), 
               # Target locations from devices
               targetPoints = st_as_sf(CALONECTRIS$targetPoints),
               
               # Origin location distance matrix
               originDist = CALONECTRIS$originDist, 
               # Breeding / origin sites 
               originSites = st_as_sf(CALONECTRIS$originSites),
               # Capture Locations 
               originPoints = st_as_sf(CALONECTRIS$originPoints), 
               # Origin relative abundances
               originRelAbund = CALONECTRIS$originRelAbund, 
               # originAssignment = CALONECTRIS$originAssignment$Colony, 
               # resampleProjection = CALONECTRIS$originPoints@proj4string, 
               
               maintainLegacyOutput = T,
               verbose = 1,   # output options - see help ??estMC
               nSamples = 1000, # This is set low for example 
               maxTries = 300) 

save(M_all, file = "data/MigCon_CALONECTRIS.RData")

# 3. Store output ####

M_CD_Col <- readRDS("outputs/M_CD_Col.RDS")
M_CD_Pop <- readRDS("outputs/M_CD_Pop.RDS")
M_CB_Col <- readRDS("outputs/M_CB_Col.RDS")
M_CB_Pop <- readRDS("outputs/M_CB_Pop.RDS")

load("data/MigCon_CALEDW.Rdata")
load("data/MigCon_CALONECTRIS.RData")

## Migratory connectivity values ####

migconCALONECTRIS <- data.frame(Sp = "Calonectris spp.",
                                level = "Species",
                                lowerCI = M_all$simpleCI[1],
                                higherCI = M_all$simpleCI[2],
                                median = M_all$medianMC,
                                mean = M_all$meanMC)

migconCALBOR_Pop <- data.frame(Sp = "Cory's shearwater",
                               level = "Population",
                               lowerCI = M_CB_Pop$simpleCI[1],
                               higherCI = M_CB_Pop$simpleCI[2],
                               median = M_CB_Pop$medianMC,
                               mean = M_CB_Pop$meanMC)

migconCALBOR_Col <- data.frame(Sp = "Cory's shearwater",
                               level = "Colony",
                               lowerCI = M_CB_Col$simpleCI[1],
                               higherCI = M_CB_Col$simpleCI[2],
                               median = M_CB_Col$medianMC,
                               mean = M_CB_Col$meanMC)

migconCALDIO_Pop <- data.frame(Sp = "Scopoli's shearwater",
                               level = "Population",
                               lowerCI = M_CD_Pop$simpleCI[1],
                               higherCI = M_CD_Pop$simpleCI[2],
                               median = M_CD_Pop$medianMC,
                               mean = M_CD_Pop$meanMC)

migconCALDIO_Col <- data.frame(Sp = "Scopoli's shearwater",
                               level = "Colony",
                               lowerCI = M_CD_Col$simpleCI[1],
                               higherCI = M_CD_Col$simpleCI[2],
                               median = M_CD_Col$medianMC,
                               mean = M_CD_Col$meanMC)


migconCALEDW_Col <- data.frame(Sp = "Cabo Verde shearwater",
                               level = "Colony",
                               lowerCI = M_CE_Col$simpleCI[1],
                               higherCI = M_CE_Col$simpleCI[2],
                               median = M_CE_Col$medianMC,
                               mean = M_CE_Col$meanMC)


migcons <- rbind(migconCALONECTRIS, 
                 # migconBORDIO, 
                 migconCALBOR_Pop, migconCALBOR_Col, 
                 migconCALDIO_Pop, migconCALDIO_Col, migconCALEDW_Col)

write.csv(migcons, file = "outputs/migcon_values.csv", row.names = F)


## Migratory connectivity plot ####

ggplot(migcons) +
  geom_point(aes(x = level, y = median, col = Sp), 
             position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(x = level, ymin = lowerCI, ymax = higherCI, col = Sp), 
                position = position_dodge(width = 0.5), 
                width = 0.25) + 
  geom_hline(aes(yintercept = 0), col = "black", lty = 2) + 
  labs(x = "Aggregation level", y = "Migratory connectivity", col = "Species") +
  theme_bw()

ggsave("outputs/migconplot.pdf")

# 4. Psi matrices ####

#' Psi values represent the probability of transition between each breeding and wintering site. The results are given in an array of dimensions $X\times Y\times Z$ (X = no. of breeding sites, Y = no. of wintering sites, Z = no. of samples of the `estMC` function). From the Z dimension we calculate median and cred interval width and lower and higher limits, and save them outside R
## For Calonectris diomedea ####

median <- round(apply(M_CD_Col$samplePsi, c(2,3), median), 3)
upper <- apply(M_CD_Col$samplePsi, c(2,3), quantile, (probs = c(0.975)))
lower <- apply(M_CD_Col$samplePsi, c(2,3), quantile, (probs = c(0.025)))
conf <- upper - lower

medianCALDIO <- as.data.frame(t(round(median, 3)))
rownames(medianCALDIO) <- CALDIO$targetSites$ECOREGION
colnames(medianCALDIO) <- CALDIO$originSites$Colony$Colony
medianCALDIO %>%
  kable() %>%
  kable_styling()
# write.csv(medianCALDIO, file = "data/medianPsiCALDIO.csv")


upperCALDIO <- as.data.frame(t(round(upper,3)))
rownames(upperCALDIO) <- CALDIO$targetSites$ECOREGION
colnames(upperCALDIO) <- CALDIO$originSites$Colony$Colony
# write.csv(upperCALDIO, file = "data/upperPsiCALDIO.csv")

lowerCALDIO <- as.data.frame(t(round(lower, 3)))
rownames(lowerCALDIO) <- CALDIO$targetSites$ECOREGION
colnames(lowerCALDIO) <- CALDIO$originSites$Colony$Colony
# write.csv(lowerCALDIO, file = "data/lowerPsiCALDIO.csv")

confCALDIO <- as.data.frame(t(round(conf, 3)))
rownames(confCALDIO) <- CALDIO$targetSites$ECOREGION
colnames(confCALDIO) <- CALDIO$originSites$Colony$Colony
# write.csv(confCALDIO, file = "data/confPsiCALDIO.csv")

## For Calonectris borealis ####

median <- round(apply(M_CB_Col$samplePsi, c(2,3), median), 3)
upper <- apply(M_CB_Col$samplePsi, c(2,3), quantile, (probs = c(0.975)))
lower <- apply(M_CB_Col$samplePsi, c(2,3), quantile, (probs = c(0.025)))
conf <- upper - lower

medianCALBOR <- as.data.frame(t(round(median, 3)))
rownames(medianCALBOR) <- CALBOR$targetSites$ECOREGION
colnames(medianCALBOR) <- CALBOR$originSites$Colony$Colony
# write.csv(medianCALBOR, file = "data/medianPsiCALBOR.csv")
medianCALBOR %>%
  kable() %>%
  kable_styling()

upperCALBOR <- as.data.frame(t(round(upper,3)))
rownames(upperCALBOR) <- CALBOR$targetSites$ECOREGION
colnames(upperCALBOR) <- CALBOR$originSites$Colony$Colony
# write.csv(upperCALBOR, file = "data/upperPsiCALBOR.csv")

lowerCALBOR <- as.data.frame(t(round(lower, 3)))
rownames(lowerCALBOR) <- CALBOR$targetSites$ECOREGION
colnames(lowerCALBOR) <- CALBOR$originSites$Colony$Colony
# write.csv(lowerCALBOR, file = "data/lowerPsiCALBOR.csv")

confCALBOR <- as.data.frame(t(round(conf, 3)))
rownames(confCALBOR) <- CALBOR$targetSites$ECOREGION
colnames(confCALBOR) <- CALBOR$originSites$Colony$Colony
# write.csv(confCALBOR, file = "data/confPsiCALBOR.csv")


## For Calonectris edwardsii ####

median <- round(apply(M_CE_Col$samplePsi, c(2,3), median), 3)
upper <- apply(M_CE_Col$samplePsi, c(2,3), quantile, (probs = c(0.975)))
lower <- apply(M_CE_Col$samplePsi, c(2,3), quantile, (probs = c(0.025)))
conf <- upper - lower

medianCALEDW <- as.data.frame(t(round(median, 3)))
rownames(medianCALEDW) <- CALEDW$targetSites$ECOREGION
colnames(medianCALEDW) <- CALEDW$originSites$Colony
# write.csv(medianCALEDW, file = "data/medianPsiCALEDW.csv")
medianCALEDW %>%
  kable() %>%
  kable_styling()

upperCALEDW <- as.data.frame(t(round(upper,3)))
rownames(upperCALEDW) <- CALEDW$targetSites$ECOREGION
colnames(upperCALEDW) <- CALEDW$originSites$Colony
# write.csv(upperCALEDW, file = "data/upperPsiCALEDW.csv")

lowerCALEDW <- as.data.frame(t(round(lower, 3)))
rownames(lowerCALEDW) <- CALEDW$targetSites$ECOREGION
colnames(lowerCALEDW) <- CALEDW$originSites$Colony
# write.csv(lowerCALEDW, file = "data/lowerPsiCALEDW.csv")

confCALEDW <- as.data.frame(t(round(conf, 3)))
rownames(confCALEDW) <- CALEDW$targetSites$ECOREGION
colnames(confCALEDW) <- CALEDW$originSites$Colony
# write.csv(confCALEDW, file = "data/confPsiCALEDW.csv")


## For all three species ####

median <- round(apply(M_all$samplePsi, c(2,3), median), 3)

upper <- apply(M_all$samplePsi, c(2,3), quantile, (probs = c(0.975)))
lower <- apply(M_all$samplePsi, c(2,3), quantile, (probs = c(0.025)))
conf <- upper - lower

medianCALONECTRIS <- as.data.frame(t(round(median, 3)))
rownames(medianCALONECTRIS) <- CALONECTRIS$targetSites$ECOREGION
colnames(medianCALONECTRIS) <- CALONECTRIS$originSites$Species
# write.csv(medianCALONECTRIS, file = "data/medianPsiCalonectris.csv")
medianCALONECTRIS %>%
  kable() %>%
  kable_styling()

upperCALONECTRIS <- as.data.frame(t(round(upper,3)))
rownames(upperCALONECTRIS) <- CALONECTRIS$targetSites$ECOREGION
colnames(upperCALONECTRIS) <- CALONECTRIS$originSites$Species
# write.csv(upperCALONECTRIS, file = "data/upperPsiCalonectris.csv")

lowerCALONECTRIS <- as.data.frame(t(round(lower, 3)))
rownames(lowerCALONECTRIS) <- CALONECTRIS$targetSites$ECOREGION
colnames(lowerCALONECTRIS) <- CALONECTRIS$originSites$Species
# write.csv(lowerCALONECTRIS, file = "data/lowerPsiCalonectris.csv")

confCALONECTRIS <- as.data.frame(t(round(conf, 3)))
rownames(confCALONECTRIS) <- CALONECTRIS$targetSites$ECOREGION
colnames(confCALONECTRIS) <- CALONECTRIS$originSites$Species
# write.csv(confCALONECTRIS, file = "data/confPsiCalonectris.csv")

save(medianCALBOR, medianCALDIO, medianCALEDW, medianCALONECTRIS, file = "data/medianPsis.Rdata")

# 5. Distances between breeding colonies and wintering areas ####
 
#' We are going to calculate distances between breeding locations and centroids of wintering areas (avoiding land) to construct a matrix equal to the psi matrices but with distances
#' 
#' We first generate SpatialPoints for wintering and breeding locations, and project to a planar coordinate system

breeding_pops <- read.csv("data/colonies_info_100km_for_ArcGis.csv", 
                          row.names = NULL)

br_points <- breeding_pops %>% 
  filter(!is.na(Pop.Min)) %>% # remove colonies without pop.size
  ddply(~Radius100, summarize, Tot.pairs = sum(Pop.Min)) %>% # sum pop.size of all colonies inside the radius ("population")
  full_join(breeding_pops, by = "Radius100") %>%  # merge the sum with the original dataset
  filter(!is.na(Sample.Size)) %>% 
  rename_at(c("Radius100", "SampledColony"), ~c("Population", "Colony")) %>% 
  mutate(Tot.pop = Tot.pairs*2)  %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) 


coordinates(br_points) <- ~ Longitude + Latitude
br_points@proj4string <- CRS(projections$WGS84)
br_points <- spTransform(br_points, CRS(projections$EquidistConic))

win_areas <- as_Spatial(read_sf(dsn = "data/shapefiles", 
                                layer = "wintering_areas_1"))

win_points <- data.frame(Area = win_areas$ECOREGION, 
                         Longitude = coordinates(win_areas)[,1],
                         Latitude = coordinates(win_areas)[,2])
coordinates(win_points) <- ~Longitude + Latitude
win_points@proj4string <- CRS(projections$WGS84)
win_points <- spTransform(win_points, CRS(projections$EquidistConic))


#' We now calculate distances avoiding water. To do that, we use the functions of the package gdistance following instructions [here](https://stackoverflow.com/questions/28595414/calculate-distance-between-2-lon-lats-but-avoid-going-through-a-coastline-in-r):

load("data/Transition_layer.RData")

# initialise matrix with win_areas as rows and bre_colonies as columns
dist <- matrix(nrow = nrow(win_points), ncol = nrow(br_points))

#calculate shprtest path and distance between pairs of win and bre points, and fill in matrix
for (i in 1:nrow(win_points)) {
  # i = 2
  for (j in 1:nrow(br_points)) {
    # j = 1
    ds <- gdistance::shortestPath(Trans, coordinates(win_points)[i,], coordinates(br_points)[j,], output = "SpatialLines")
    # plot(water, main = paste("Distance between", as.character(win_points@data[i,1]), "and", as.character(br_points@data[j,5])), 
         # cex.main = 0.8)
    # plot(ds, lwd = 3, add = T)
    dis <- SpatialLinesLengths(ds, longlat = F)
    dist[i,j] <- dis
  }
}

rownames(dist) <- win_points$Area
colnames(dist) <- br_points$Colony
dist <- round(dist/1000, 0)
knitr::kable(dist, digits = 0)
# write.csv(dist, "outputs/distances_br_win.csv")

#' This matrix is saved as csv an, outside R, separated by species and reordered so it has the same structure as the psi matrices. We now load those matrices

# 6. Calculate distances between wintering centroids and breeding locations ####

#' We are going to use the same method explained above to calculate the shortest distance, avoiding land, between each individual's wintering centroid and breeding location, and between each individual's wintering centroid and the strait of Gibraltar (these we are going to use only for _Calonectris diomedea_, but we do it inside the loop so it gets calculated for everyone)

#' We first load and prepare the data

load("C:/Users/morer/Dropbox/Virginia/GLS_multi/data/Centroids.ALL.2023.Rdata")

#merge with colony locations
br_locs <- data.frame(Colony  = breeding_info$Colony,
                      Lon     = breeding_info$Longitude,
                      Lat     = breeding_info$Latitude, 
                      Species = breeding_info$Species)

for_distance <- centroids.ALL %>% left_join(br_locs)
for_distance <- for_distance[!duplicated(for_distance),]

for_distance <- for_distance[order(for_distance$Colony),]


#' Now we actually calculate the distances inside a loop. This is very time-consuming, so it was calculated once and then saved, sowe are able to avoid evaulating the next chunk

## # this is a shitty loop, should be a function but I can't seem to figure out what apply to use to apply it row-wise, and split to apply to each element of a list changes the order of the column unique_ID_Bird, so it messes up the output
## 
## # shortest_distance <- function(bre, win, Trans, proj = F) {
## distances <- rep(NA, nrow(for_distance))
## dist_from_strait <- rep(NA, nrow(for_distance))
## for (i in 1:nrow(for_distance)) {
##   # i = 1
##   # proj = F
##   x <- for_distance[i,]
##   bre <- x[,5:6]; win <- x[,3:4]
##   strait <- SpatialPoints(coords = data.frame(-5.5, 36))
##   bre = SpatialPoints(coords = bre)
##   win = SpatialPoints(coords = win)
##   if (proj == T) {
##     bre@proj4string <- CRS(projections$EquidistConic)
##     win@proj4string <- CRS(projections$EquidistConic)
##     strait@proj4string <- CRS(projections$EquidistConic)} else {
##       bre@proj4string <- CRS(projections$WGS84)
##       bre <- spTransform(bre, CRS(projections$EquidistConic))
##       win@proj4string <- CRS(projections$WGS84)
##       win <- spTransform(win, CRS(projections$EquidistConic))
##       strait@proj4string <- CRS(projections$WGS84)
##       strait <- spTransform(strait, CRS(projections$EquidistConic))}
##   ds <- gdistance::shortestPath(Trans, bre, win, output = "SpatialLines")
##   dG <- gdistance::shortestPath(Trans, strait, win, output = "SpatialLines")
##   # plot(water)
##   # plot(ds, lwd = 3, add = T)
##   # plot(dG, lwd = 3, lty = 6, col = "red", add = T)
##   dis <- SpatialLinesLengths(ds, longlat = F)
##   diG <- SpatialLinesLengths(dG, longlat = F)
##   print(paste("Distance between", x$unique_ID_Bird, "from", x$Colony, "and its wintering area is", round(dis/1000, 0), "km", sep = " "))
##   distances[i] <- dis
##   dist_from_strait[i] <- diG
##   # return(dis)
## }
## 
## for_distance$dist <- distances
## for_distance$distStrait <- dist_from_strait
## 
## save(for_distance, file = "outputs/distances_WinBre.Rdata")

