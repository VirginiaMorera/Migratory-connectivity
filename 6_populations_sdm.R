#0. House keeping ####
setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")

pacman::p_load(embarcadero, tidyverse, sf, raster, dismo)
                    
'%!in%' <- function(x,y)!('%in%'(x,y))

# Load all datasets

predictors_original <- stack("data/env_data/wintering/all_vars.grd")

predictors_original <- stack(predictors_original, init(predictors_original, 'y'))

names(predictors_original)[2] <- "chla_var"
names(predictors_original)[8] <- "latitude"

pops_names <- str_replace(list.files("embarcadero_SDMs/data/populations_ds"), 
                          ".RDS", "")

pops <- list.files("embarcadero_SDMs/data/populations_ds", full.names = T)

pops_list <- purrr::map(pops, readRDS)

names(pops_list) <- pops_names


# 1. Prepare all pres/abs ####

pres_abs_list <- list()

for(i in seq_along(pops_list)) {
  # i = 1
  temp <- pops_list[[i]]
  temp_sf <- temp %>% st_as_sf(coords = c("x", "y"))
  
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
names(pres_abs_list) <- pops_names

# 2. Run loop with all models, produce all predictions ####

model_list <- list()

# for(i in seq_along(pres_abs_list)) {
for(i in 1:17) {
  # i = 1
  all.cov <- pres_abs_list[[i]]
  
  xvars <-names(predictors_original)[-8]
  # remove lat lon for now and decide if we include it later
  sdm_list <- list()
  for (j in 1:10) {
    print(paste("population", i, "out of 17 and iteration", j, "out of 10", sep = " "))
    sdm <-bart.step(x.data = all.cov[,xvars],
                    y.data = all.cov[,'pres'],
                    full = TRUE,
                    quiet = F, 
                    iter.step = 100, 
                    iter.plot = 100, 
                    tree.step = 10)
  sdm_list[j] <- list(sdm)}
  
  saveRDS(sdm_list, file = paste0("embarcadero_SDMs/iter_outputs/", pops_names[i] ,"_model_list.RDS"))
}
          
          
#   model_list[i] <- list(sdm)
#   
#   prediction <- predict(object = sdm,
#                         x.layers = predictors_original,
#                         quantiles =c(0.025, 0.975),
#                         splitby = 20,
#                         quiet = TRUE)
#   
#   names(prediction) <- c("Mean", "q2.5%", "q97.5%")
#   plot(prediction)
#   writeRaster(prediction, filename = paste0("embarcadero_SDMs/final_outputs/population_models/", pops_names[i] ,"_pred.grd"), 
#               overwrite = T)
# }

names(model_list) <- pops_names

saveRDS(model_list, file = "embarcadero_SDMs/final_outputs/population_models/model_list.RDS")
model_list <- readRDS("embarcadero_SDMs/final_outputs/population_models/model_list.RDS")


# 3. Obtain outputs ####

### Variable importance ####

pop_models <- list.files(path = "embarcadero_SDMs/iter_outputs/", pattern = "model_list.RDS", 
                         full.names = T)
pop_names <- str_sub(list.files(path = "embarcadero_SDMs/iter_outputs/", pattern = "model_list.RDS", 
                         full.names = F), end = -16)


pop_list <- list()
for (i in seq_along(pop_models)){
  # i = 1
  
  mod_list <- readRDS(pop_models[i]) 
  iter_list <- list()
  for (j in seq_along(mod_list))   {
    varimp <- varimp(mod_list[[j]]) %>% 
      arrange(-varimps) %>% 
      mutate(iter = paste0("iter_", j), 
             pop = pop_names[i])
    
    iter_list[j] <- list(varimp)
  }
  
  pop_list[i] <- list(bind_rows(iter_list, id = NULL))

}
  

all_varimps <- bind_rows(pop_list)

varimps_sum <- all_varimps %>% 
  separate(pop, into = c("Species", "Population"), sep = "_") %>% 
  group_by(Species, Population, names) %>% 
  summarise(mean = mean(varimps), 
            sd = sd(varimps), 
            iters = n())
  

ggplot(varimps_sum) + 
  geom_point(aes(x = names, y = mean, alpha = iters)) + 
  geom_errorbar(aes(x = names, ymin = mean-sd, ymax = mean+sd, alpha = iters)) +
  facet_grid(Population~., scales = "free") + 
  theme_bw()

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
  separate(Sp_Population, into = c("Species", "Population")) %>% 
  dplyr::select(Species, Population, Variable = names, Importance = varimps) %>% 
  mutate(Variable = recode(Variable, sst = "Sea Surface Temperature", chla_var = "Chlorophyll A Variability", 
                           logchla_lag3 = "Chlorophyll A", bati = "Bathimetry", 
                           sst_grad = "Sea Surface Temperature gradient", sal = "Salinity", 
                           slope = "Slope2"))

write.csv(varimp_pops, 
            file = "embarcadero_SDMs/final_outputs/population_models/populations_varimp.csv")
  

### Partial response ####

for(i in seq_along(model_list)) {
  print(i)
  # i = 1
  sdm <- model_list[[i]]
  
  tmp_partial <- partial(sdm,
                         trace = FALSE,
                         ci = TRUE,
                         equal = TRUE,
                         smooth = 10, 
                         panels = F)
  
  saveRDS(tmp_partial, 
          file = paste0("embarcadero_SDMs/final_outputs/population_models/partial_plots_", 
                        names(model_list)[i], ".RDS"))
}

