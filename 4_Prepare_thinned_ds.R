
# House keeping ####
setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")

pacman::p_load(spatstat, tidyverse, rgeos, sf, raster)
                    
'%!in%' <- function(x,y)!('%in%'(x,y))

source("scripts/migcon_AUX.R")
sourceDir("scripts/to_source")

# 0. Load all datasets ####

#representative colonies
boot_out_summary <- read.csv("embarcadero_SDMs/data/colony_bootstrap_summary.csv", row.names = NULL)
representative_colonies <- boot_out_summary %>% 
  filter(SampleSize > 5 & Representativity > 70) 

#data
all_complete <- read.csv("data/complete_fenology.csv", row.names = NULL)
all_complete <- all_complete %>% 
  mutate(Colsp = paste(Colony, Species, sep = "_")) %>% 
  filter(Colsp != "Chafarinas_CALDIO") %>% 
  filter(unique_ID_Bird != "6143080_2012") %>% #<6 win positions -> trouble with the centroid calculation
  filter(unique_ID_Bird %!in% c("L38927_first", "6195227_2014")) %>% 
  filter(Colony %in% representative_colonies$Colony) %>% 
  dplyr::select(Species, Colony, Population, unique_ID_Bird, Date_Time, Lat, Lon, fen) %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony))  

# population sizes info
breeding_pops <- read.csv("data/Pop_Col_info.csv", row.names = NULL)

colony_info <- breeding_pops %>% 
  filter(!is.na(n)) %>% 
  mutate(Colsp = paste(Colony, Species, sep = "_")) %>% 
  filter(Colsp != "Chafarinas_CALDIO") %>% 
  dplyr::select(Species, Population, Colony, Pop.Min = Min.pop, Longitude = Lon, Latitude = Lat) %>% 
  mutate(Tot.pop = Pop.Min*2) %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) 


# projections
load("data/projections.Rdata")

# worldmap
world_simple <- rworldmap::getMap(resolution = "high")
world_simple <- gBuffer(world_simple, byid = F, width = 0)
map_sf <- st_as_sf(world_simple)

predictors_original  <- stack("data/env_data/wintering/predictors_win_stack.grd")

# 1. Prepare presences dataset ####
 
#' Before actually running the SDM we'll
#' - Select only stationary non-breeding positions
#' - Remove points on land
 
all_nbr <- all_complete %>% 
  filter(fen %in% c("STAGE", "WIN")) %>% 
  dplyr::select(Species, Colony, unique_ID_Bird, Latitude = Lat, Longitude = Lon, fen) %>% 
  droplevels() %>% 
  st_as_sf(coords = c("Longitude", "Latitude")) %>%
  st_set_crs(st_crs(map_sf)) 

ovr <- sapply(st_intersects(all_nbr, map_sf), function(z) if (length(z) == 0) NA_integer_ else z[1])

all_nbr <- all_nbr %>%
  mutate(ovr = ovr) %>% 
  filter(is.na(ovr)) 

# 2. Prepare datasets for colony sdms ####
#' 
#' - Simulate 2000 points for each colony (could just subsample in the larger colonies but would need to simulate in the smaller ones so we subsample in everyone)
 
all_nbr_list <- split(all_nbr, all_nbr$Colony)
# for (i in seq_along(all_nbr_list)) {
for (i in seq_along(all_nbr_list)) {
  # i = 1
  col1 <- all_nbr_list[[i]] %>%
    rename(Population = Colony) %>% 
    dplyr::select(-ovr) %>% 
    sfc_as_cols(names = c("Longitude", "Latitude")) %>%
    st_set_geometry(NULL)
  
  # ColInfo must be a data frame with at least Population, Species and Pairs
  ColInfo <- colony_info %>%
    filter(Colony == as.character(unique(col1$Population))) %>% 
    dplyr::select(Species, Population = Colony, Pairs = Pop.Min)
  
  sim1 <- simulateDistribution(col1, ColInfo, scale = 186, multiFactor = 50)
  col1_sim <- as.data.frame(sim1)
  names(col1_sim) <- c("Colony", "Species", "Longitude", "Latitude")
  
  colony_unthinned_spatial[i] <- list(col1_sim)
  names(colony_unthinned_spatial)[i] <- as.character(unique(col1$Population))
  
  # thin
  ppp <- ppp(sim1@coords[,1], sim1@coords[,2], sim1@bbox[1,], sim1@bbox[2,])
  ppp <- ppp[!duplicated.ppp(ppp)]
  th <- rthin(ppp, 2000/ppp$n)
  col1_thinned <- as.data.frame(th)
  
  col1_thinned$Species <- unique(col1$Species)
  col1_thinned$Colony <- unique(col1$Population)
  names(col1_thinned) <- c("Longitude", "Latitude", "Species", "Colony")
  colony_thinned_df[i] <- list(col1_thinned)
  names(colony_thinned_df)[i] <- as.character(unique(col1$Population))
}

names(colony_thinned_df)
thinned_df <- do.call(rbind.data.frame, colony_thinned_df)
saveRDS(thinned_df, file = "embarcadero_SDMs/data/thinned_df_all_cols.RDS")

# 3. Prepare datasets for population sdms so the number of points for each colony is proportional to colony size ####

all_nbr <- all_nbr %>% 
  left_join(colony_info %>% dplyr::select(Colony, Population))

pop_list <- split(all_nbr, all_nbr$Population)

for(j in seq_along(pop_list)) {
  # j = 1
  temp_df <- pop_list[[j]]
  temp_list <- split(temp_df, temp_df$Colony)
  popj <- unique(temp_df$Population)
  col_list <- list()
  
  for (i in seq_along(temp_list)) {
    # i = 1
    col1 <- temp_list[[i]] %>%
      dplyr::select(-Population) %>% 
      rename(Population = Colony) %>% 
      sfc_as_cols(names = c("Longitude", "Latitude")) %>%
      st_set_geometry(NULL)
    
    # ColInfo must be a data frame with at least Population, Species and Pairs
    ColInfo <- colony_info %>%
      dplyr::select(-Population) %>% 
      filter(Colony == as.character(unique(col1$Population))) %>% 
      dplyr::select(Species, Population = Colony, Pairs = Pop.Min)
    
    sim1 <- simulateDistribution(col1, ColInfo, scale = 186, multiFactor = 50)
    col1_sim <- as.data.frame(sim1)
    names(col1_sim) <- c("Colony", "Species", "Longitude", "Latitude")

    col_list[i] <- list(col1_sim)
    names(pop_list)[i] <- as.character(unique(col1$Population))
  }
  
  # names(col_list)
  pop_df <- do.call(rbind.data.frame, col_list)
  pop_df$Population <- popj
  spj <- unique(pop_df$Species)
  pop_sp <- st_as_sf(pop_df, coords = c("Longitude", "Latitude"))
  
  # thin
  ppp <- ppp(st_coordinates(pop_sp)[,1], st_coordinates(pop_sp)[,2], 
             st_bbox(pop_sp)[c(1,3)], st_bbox(pop_sp)[c(2,4)])
  
  ppp <- ppp[!duplicated.ppp(ppp)]
  th <- rthin(ppp, 2000/ppp$n)
  pop_df_thinned <- as.data.frame(th)
  
  
  saveRDS(pop_df_thinned, file = paste0("embarcadero_SDMs/data/",spj, "_", popj, ".RDS"))
}

# 4. Prepare datasets for species sdms so the number of points for each colony is proportional to colony size ####

sp_list <- split(all_nbr, all_nbr$Species)

for(j in seq_along(sp_list)) {
  # j = 1
  temp_df <- sp_list[[j]]
  temp_list <- split(temp_df, temp_df$Colony)
  col_list <- list()
  
  for (i in seq_along(temp_list)) {
    # i = 1
    col1 <- temp_list[[i]] %>%
      rename(Population = Colony) %>% 
      sfc_as_cols(names = c("Longitude", "Latitude")) %>%
      st_set_geometry(NULL)
    
    # ColInfo must be a data frame with at least Population, Species and Pairs
    ColInfo <- colony_info %>%
      filter(Colony == as.character(unique(col1$Population))) %>% 
      dplyr::select(Species, Population = Colony, Pairs = Pop.Min)
    
    sim1 <- simulateDistribution(col1, ColInfo, scale = 186, multiFactor = 10)
    col1_sim <- as.data.frame(sim1)
    names(col1_sim) <- c("Colony", "Species", "Longitude", "Latitude")

    col_list[i] <- list(col1_sim)
    names(col_list)[i] <- as.character(unique(col1$Population))
  }
  
  # names(col_list)
  sp_df <- do.call(rbind.data.frame, col_list)
  sp_sp <- st_as_sf(sp_df, coords = c("Longitude", "Latitude"))
  spj <- unique(sp_df$Species)
  # thin
  ppp <- ppp(st_coordinates(sp_sp)[,1], st_coordinates(sp_sp)[,2], 
             st_bbox(sp_sp)[c(1,3)], st_bbox(sp_sp)[c(2,4)])
  
  ppp <- ppp[!duplicated.ppp(ppp)]
  th <- rthin(ppp, 2000/ppp$n)
  sp_df_thinned <- as.data.frame(th)
  sp_df_thinned$Species <- spj

  sp_df_thinned <- sp_df_thinned %>% 
    dplyr::select(Species, Longitude = x, Latitude = y)
  
  saveRDS(sp_df_thinned, file = paste0("embarcadero_SDMs/data/",spj, ".RDS"))
}

