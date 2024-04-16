 
# 0. House keeping ####
setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")

pacman::p_load(embarcadero, tidyverse, sf, rgeos, raster, dismo)
                    
'%!in%' <- function(x,y)!('%in%'(x,y))

# Load all datasets
predictors_original <- stack("data/env_data/wintering/all_vars.grd")

predictors_original <- stack(predictors_original, init(predictors_original, 'y'))

names(predictors_original)[2] <- "chla_var"
names(predictors_original)[8] <- "latitude"

colony_ds <- readRDS("embarcadero_SDMs/data/thinned_df_all_cols.RDS")

colony_sf <- colony_ds %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

# worldmap
world_simple <- rworldmap::getMap(resolution = "high")
world_simple <- gBuffer(world_simple, byid = F, width = 0)
map_sf <- st_as_sf(world_simple)

ovr <- sapply(st_intersects(colony_sf, map_sf), function(z) if (length(z) == 0) NA_integer_ else z[1])

colony_sf <- colony_sf %>%
  mutate(ovr = ovr) %>% 
  filter(is.na(ovr)) %>% 
  dplyr::select(-ovr)

colony_list <- colony_sf %>% 
  mutate(Longitude = st_coordinates(.)[,1], 
         Latitude = st_coordinates(.)[,2]) %>% 
  st_set_geometry(NULL) %>% 
  group_by(Colony) %>% 
  group_split()


# 1. Prepare all pres/abs ####

pres_abs_list <- list()

for(i in seq_along(colony_list)) {
  # i = 1
  temp <- colony_list[[i]]
  temp_sf <- temp %>% st_as_sf(coords = c("Longitude", "Latitude"))
  pres.cov <- raster::extract(predictors_original, temp_sf)
  
  absence <- randomPoints(predictors_original, nrow(temp_sf))
  abs.cov <- raster::extract(predictors_original, absence)
  
  pres.cov <- data.frame(pres.cov)
  pres.cov$pres <- 1
  
  abs.cov <- data.frame(abs.cov)
  abs.cov$pres <- 0
  
  # And one to bind them
  all.cov <-rbind(pres.cov, abs.cov)
  all.cov %>% drop_na -> all.cov

  pres_abs_list[i] <- list(all.cov) 
}

names(pres_abs_list) <- unique(colony_sf$Colony)

# 2. Run loop with all models, produce all predictions ####

model_list <- list()

for(i in 6:length(pres_abs_list)) {
  # i = 1
  print(i)
  all.cov <- pres_abs_list[[i]]
  
  xvars <- names(predictors_original)[-8]
  # remove lat lon for now and decide if we include it later
  
  sdm <-bart.step(x.data = all.cov[,xvars],
                  y.data = all.cov[,'pres'],
                  full = TRUE,
                  quiet = F, 
                  iter.step = 100, 
                  iter.plot = 100, 
                  tree.step = 10)
  saveRDS(sdm, file = paste0("embarcadero_SDMs/final_outputs/colony_models/", 
                             unique(colony_sf$Colony)[i] ,"_model.RDS"))
  
  model_list[i] <- list(sdm)
  
  prediction <- predict(object = sdm,
                        x.layers = predictors_original,
                        quantiles =c(0.025, 0.975),
                        splitby = 20,
                        quiet = TRUE)
  plot(prediction)
  names(prediction) <- c("Mean", "q2.5%", "q97.5%")
  writeRaster(prediction, 
              filename = paste0("embarcadero_SDMs/final_outputs/colony_models/", 
                                unique(colony_sf$Colony)[i] ,"_pred.grd"), 
              overwrite = T)
}

names(model_list) <- unique(colony_sf$Colony)
saveRDS(model_list, file = "embarcadero_SDMs/final_outputs/colony_models/col_model_list.RDS")
model_list <- readRDS("embarcadero_SDMs/final_outputs/colony_models/col_model_list.RDS")


# 3. Obtain outputs  ####

## Variable importance ####

#' We can now generate a plot of the importance of the included variables

varimp_list <- list()
for(i in seq_along(model_list)) {
  # i = 1
  
  varimp <- varimp(model_list[[i]], plots = F)
  
  varimp_ordered <- varimp %>% 
    arrange(-varimps) %>% 
    mutate(Sp_Population = names(model_list)[i])
  
  varimp_list[i] <- list(varimp_ordered)
}

varimp_pops <- map_df(varimp_list, bind_rows)
row.names(varimp_pops) <- NULL

varimp_pops <- varimp_pops %>% 
  dplyr::select(Colony = Sp_Population, Variable = names, Importance = varimps) %>% 
  mutate(Variable = recode(Variable, sst = "Sea Surface Temperature", chla_var = "Chlorophyll A Variability", 
                           logchla_lag3 = "Chlorophyll A", bati = "Bathimetry", 
                           sst_grad = "Sea Surface Temperature gradient", sal = "Salinity", 
                           slope = "Slope2"))

write.csv(varimp_pops, 
            file = "embarcadero_SDMs/final_outputs/colony_models/colony_varimp.csv")
  


## Partial response ####
#' These plots show the mean effect of each variable when measured at all other values of the other predictors, for all the values of the variable of interest. This will show us at which values of the covariate the dependence is higher.
#' With CI you can add credible intervals, and with trace = TRUE you can see the dependence plot for each of the trees (remember BARTs are regression trees), with the thicker line representing the mean of all the trees. 

for(i in seq_along(model_list)) {
  # i = 1
  sdm <- model_list[[i]]
  
  tmp_partial <- partial(sdm,
                         trace = FALSE,
                         ci = TRUE,
                         equal = TRUE,
                         smooth = 10, 
                         panels = F)
  
  saveRDS(tmp_partial, 
          file = paste0("embarcadero_SDMs/outputs/partial_plots_", names(model_list)[i], ".RDS"))
}


