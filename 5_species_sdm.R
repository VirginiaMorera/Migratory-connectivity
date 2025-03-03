
# 0. House Keeping ####
rm(list = ls())
setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")

pacman::p_load(embarcadero, tidyverse, sf, raster, dismo, spatstat)
                    
'%!in%' <- function(x,y)!('%in%'(x,y))

source("scripts/0_migcon_AUX.R")

# 0. Load all datasets ####
predictors_original <- stack("data/env_data/wintering/all_vars.grd")

predictors_original <- stack(predictors_original, init(predictors_original, 'y'))

names(predictors_original)[2] <- "chla_var"
names(predictors_original)[8] <- "latitude"

predictors_scaled <- scale(predictors_original)

CALBOR <- readRDS("embarcadero_SDMs/data/iterated_ds/CALBOR_10_iters_list.RDS")
CALDIO <- readRDS("embarcadero_SDMs/data/iterated_ds/CALDIO_10_iters_list.RDS")
CALEDW <- readRDS("embarcadero_SDMs/data/iterated_ds/CALEDW_10_iters_list.RDS")

# 1. CALBOR ####

## 1.1 Prepare dataset ####

### Extract covariate values at presences ####

all.cov_list <- list()
for (i in seq_along(CALBOR)) {
  # i = 1
  CALBOR_sub <- CALBOR[[i]]
  
  CALBOR_sf <- CALBOR_sub %>% 
    st_as_sf(coords = c("Longitude", "Latitude"))
  
  pres.cov <- raster::extract(predictors_scaled, CALBOR_sf)
  
  ### Generate absences and covariate values at absences ####
  
  # From the vignette: 
  # BARTs are like boosted regression trees (BRTs) in that they are sensitive 
  # to assumed prevalence; anecdotally,I strongly suggest using an equal number 
  # of presences and absences in your training data. You can experiment with 
  # the demo data by changing “nrow(ticks)” to “5000” below if you want to see 
  # some odd model behavior
  
  absence <- randomPoints(predictors_scaled$bati, nrow(CALBOR_sf))
  abs.cov <- raster::extract(predictors_scaled, absence)
  
  absence <- randomPoints(predictors_scaled$bati, nrow(CALBOR_sf))
  abs.cov <- raster::extract(predictors_scaled, absence)
  
  ### Code the response ####
  
  #' Generate a vector with 1 for presences and 0 for absences (the response)
  pres.cov <- data.frame(pres.cov)
  pres.cov$pres <- 1
  
  abs.cov <- data.frame(abs.cov)
  abs.cov$pres <- 0
  
  # And one to bind them
  all.cov <-rbind(pres.cov, abs.cov)
  head(all.cov)
  table(all.cov$pres)
  
  # drop na's
  all.cov <- all.cov %>% 
    drop_na()
  
  # check that n is still more or less even between pres and abs
  table(all.cov$pres)
  all.cov_list[i] <- list(all.cov)
  
}



## 1.2. Run BART model ####

# Run a function that automatically does variable selection and runs the 
# fine tunned model

# covariate names without latitude
CB_stack <- stack()
CB_models <- list()
for(i in seq_along(all.cov_list)){
  # i = 1
  print(paste("Run model iter", i, sep = "_"))
  
  xvars <- names(predictors_scaled)[-8]
  
  all.cov_sub <- all.cov_list[[i]]
  calbor.sdm <-bart.step(
    x.data = all.cov_sub[,xvars],
    y.data = all.cov_sub[,'pres'],
    full = TRUE,
    quiet = F 
    # iter.step = 10,
    # iter.plot = 2,
    # tree.step = 5
  )
  # necessary step to save properly (according to Carlson's github but it doesn't work so do all this in one session)
  invisible(calbor.sdm$fit$state) 
  CB_models[i] <- list(calbor.sdm)
  # saveRDS(calbor.sdm, file = "embarcadero_SDMs/final_outputs/CALBOR.BART.sdm.RDS")
  # savedCALBOR.sdm <- readRDS("embarcadero_SDMs/outputs/CALBOR.BART.sdm.RDS")
  
  ## 1.3. Obtain outputs ####
  
  ### Make a definitive prediction ####
  
  # First we obtain the predicted mean and quantiles
  CB_prediction <- predict(object = calbor.sdm,
                           x.layers = predictors_scaled,
                           quantiles =c(0.5),
                           splitby = 20,
                           quiet = TRUE)
  plot(CB_prediction)
  
  # save for definitive plotting
  CB_stack <- stack(CB_stack, CB_prediction)
}


CB_stack2 <- dropLayer(CB_stack, c(1,3,5,7,9,11,13,15,17,19))
plot(CB_stack2)

writeRaster(CB_stack2, 
            filename = "embarcadero_SDMs/iter_outputs/CALBOR_pred_stack_2024.grd", 
            overwrite = T)

saveRDS(CB_models, file = "embarcadero_SDMs/iter_outputs/CALBOR_models.RDS")
CB_models <- readRDS("embarcadero_SDMs/iter_outputs/CALBOR_models.RDS")

summary(CB_models[[1]])

# We check how the mean prediction and the two quantiles look

# CB_prediction <- stack("embarcadero_SDMs/outputs/CALBOR_pred_stack.grd")
# par(mfrow = c(2,2))
# plot(CB_prediction[[1]],
#      box = FALSE,
#      axes = FALSE,
#      main ='CALBOR prediction mean',
#      zlim =c(0,1),
#      axis.args = list(at = pretty(0:1), 
#                       labels = pretty(0:1)),
#      legend.args = list(text = 'Probability', 
#                         side = 2, line = 1.3))
# 
# plot(CB_prediction[[2]],
#      box = FALSE,
#      axes = FALSE,
#      main ='CALBOR prediction 2.5%',
#      zlim =c(0,1),
#      axis.args = list(at = pretty(0:1), 
#                       labels = pretty(0:1)),
#      legend.args = list(text = 'Probability', 
#                         side = 2, line = 1.3))
# 
# plot(CB_prediction[[3]],
#      box = FALSE,
#      axes = FALSE,
#      main ='CALBOR prediction 97.5%',
#      zlim =c(0,1),
#      axis.args = list(at = pretty(0:1), 
#                       labels = pretty(0:1)),
#      legend.args = list(text = 'Probability', 
#                         side = 2, line = 1.3))
# 
# plot(CB_prediction[[3]] - CB_prediction[[2]],
#      box = FALSE,
#      axes = FALSE,
#      main ='CALBOR prediction CI width',
#      zlim =c(0,0.75),
#      axis.args = list(at = pretty(0:1), 
#                       labels = pretty(0:1)),
#      legend.args = list(text = 'CI width', 
#                         side = 2, line = 1.3))
# par(mfrow = c(1,1))
# 
# # We now produce a binary map using as cutoff the value suggested by the model summary
# 
# plot(CB_prediction[[1]] > 0.5104638  , # cutoff obtained from summary above
#      box = FALSE,
#      axes = FALSE, 
#      main = 'Binary prediction', 
#      axis.args = list(at = pretty(0:1), 
#                       labels = pretty(0:1)),
#      legend.args = list(text = 'Presence - absence', 
#                         side = 2, line = 1.3))
# 
# # We now produce a map of high uncertainty, plotting a binary map of the areas where the uncertainty (CI width) is in the >75% quantile
# 
# quant <- quantile(values(CB_prediction[[3]] - CB_prediction[[2]]), 0.75,
#                   na.rm = TRUE)
# 
# plot((CB_prediction[[3]] - CB_prediction[[2]]) > quant,
#      box = FALSE,
#      axes = FALSE,
#      main = "Highest uncertainty zones",
#      axis.args = list(at = pretty(0:1), 
#                       labels = pretty(0:1)),
#      legend.args = list(text = 'CI width', 
#                         side = 2, line = 1.3))

### Variable importance ####

# We can now generate a plot of the importance of the included variables
varimp_list <- list()

for(i in seq_along(CB_models)) {
  # i = 1
  mod_i <- CB_models[[i]]
  x <- varimp(mod_i, plots = F)
  varimp_i <- x %>% 
    mutate(iter = paste0("iter_", i))
  varimp_list[i] <- list(varimp_i) 
}


varimp_all <- bind_rows(varimp_list)
varimp_sum <- varimp_all %>% 
  group_by(names) %>% 
  summarise(mean = mean(varimps), 
            sd = sd(varimps)) %>% 
  ungroup()


ggplot(varimp_sum, aes(x = names)) + 
  geom_point(aes(y = mean)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean + sd)) + 
  theme_bw()

write.csv(varimp_sum, file = "embarcadero_SDMs/iter_outputs/CALBOR_varimp_iter.csv")

### Partial response ####

# These plots show the mean effect of each variable when measured at all other 
# values of the other predictors, for all the values of the variable of interest. 
# This will show us at which values of the covariate the dependence is higher.

# With CI you can add credible intervals, and with trace = TRUE you can see the 
# dependence plot for each of the trees (remember BARTs are regression trees), 
# with the thicker line representing the mean of all the trees. 
all_partials <- list()

for (i in seq_along(CB_models)) {
  # i = 1
  print(paste("Run model iter", i, sep = "_"))
  calbor.sdm <- CB_models[[i]]
  
  chla_partial_cb <- partial(calbor.sdm,
                             x.vars = "logchla_lag3",
                             trace = FALSE,
                             ci = TRUE,
                             equal = TRUE,
                             smooth = 10, 
                             panels = F)
  
  chla_i <- chla_partial_cb[[1]]$data %>% 
    mutate(variable = "chla", 
           iter = paste0("iter_", i))
  
  chla_var_partial_cb <- partial(calbor.sdm,
                                 x.vars = "chla_var",
                                 trace = FALSE,
                                 ci = TRUE,
                                 equal = TRUE,
                                 smooth = 10, 
                                 panels = F)
  
  chla_var_i <- chla_var_partial_cb[[1]]$data %>% 
    mutate(variable = "chla_var", 
           iter = paste0("iter_", i))
  
  sst_partial_cb <- partial(calbor.sdm,
                            x.vars = "sst",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)
  
  sst_i <- sst_partial_cb[[1]]$data %>% 
    mutate(variable = "sst", 
           iter = paste0("iter_", i))
  
  
  sal_partial_cb <- partial(calbor.sdm,
                            x.vars = "sal",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)
  
  
  sal_i <- sal_partial_cb[[1]]$data %>% 
    mutate(variable = "sal", 
           iter = paste0("iter_", i))
  
  all_i <- bind_rows(
    chla_i,
    chla_var_i,
    sst_i,
    sal_i
  )
  
  all_partials[i] <- list(all_i)
}


all_partials_df <- bind_rows(all_partials)


scaled_chla_vals <- scale(predictors_original$logchla_lag3[])
scaled_chlaVar_vals <- scale(predictors_original$chla_var[])
scaled_sst_vals <- scale(predictors_original$sst[])
scaled_sal_vals <- scale(predictors_original$sal[])


all_partials_df <- all_partials_df %>% 
  mutate(
    x2 = if_else(variable == "chla", x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center'), 
                 if_else(variable == "chla_var", x * attr(scaled_chlaVar_vals, 'scaled:scale') + attr(scaled_chlaVar_vals, 'scaled:center'), 
                         if_else(variable == "sst", x * attr(scaled_sst_vals, 'scaled:scale') + attr(scaled_sst_vals, 'scaled:center'), 
                                 x * attr(scaled_sal_vals, 'scaled:scale') + attr(scaled_sal_vals, 'scaled:center')))), 
    Species = "CALBOR")



ggplot(all_partials_df, aes(x = x2)) + 
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "lightblue") + 
  geom_line(aes(y = med)) + 
  theme_bw() + 
  facet_grid(iter ~ variable, scales = "free")

write.csv(all_partials_df, 
        file = "embarcadero_SDMs/iter_outputs/calbor_partial_df_iter.csv")


# 2. CALDIO ####

## 2.1 Prepare dataset ####

### Extract covariate values at presences ####

all.cov_list <- list()
for (i in seq_along(CALDIO)) {
  # i = 1
  CALDIO_sub <- CALDIO[[i]]
  
  CALDIO_sf <- CALDIO_sub %>% 
    st_as_sf(coords = c("Longitude", "Latitude"))
  
  pres.cov <- raster::extract(predictors_scaled, CALDIO_sf)
  
  ### Generate absences and covariate values at absences ####
  
  # From the vignette: 
  # BARTs are like boosted regression trees (BRTs) in that they are sensitive 
  # to assumed prevalence; anecdotally,I strongly suggest using an equal number 
  # of presences and absences in your training data. You can experiment with 
  # the demo data by changing “nrow(ticks)” to “5000” below if you want to see 
  # some odd model behavior
  
  absence <- randomPoints(predictors_scaled$bati, nrow(CALDIO_sf))
  abs.cov <- raster::extract(predictors_scaled, absence)
  
  absence <- randomPoints(predictors_scaled$bati, nrow(CALDIO_sf))
  abs.cov <- raster::extract(predictors_scaled, absence)
  
  ### Code the response ####
  
  #' Generate a vector with 1 for presences and 0 for absences (the response)
  pres.cov <- data.frame(pres.cov)
  pres.cov$pres <- 1
  
  abs.cov <- data.frame(abs.cov)
  abs.cov$pres <- 0
  
  # And one to bind them
  all.cov <-rbind(pres.cov, abs.cov)
  head(all.cov)
  table(all.cov$pres)
  
  # drop na's
  all.cov <- all.cov %>% 
    drop_na()
  
  # check that n is still more or less even between pres and abs
  table(all.cov$pres)
  all.cov_list[i] <- list(all.cov)
  
}



## 2.2. Run BART model ####

# Run a function that automatically does variable selection and runs the 
# fine tunned model

# covariate names without latitude
CD_stack <- stack()
CD_models <- list()
for(i in seq_along(all.cov_list)){
  # i = 1
  xvars <- names(predictors_scaled)[-8]
  
  all.cov_sub <- all.cov_list[[i]]
  caldio.sdm <-bart.step(
    x.data = all.cov_sub[,xvars],
    y.data = all.cov_sub[,'pres'],
    full = TRUE,
    quiet = F 
    # iter.step = 10,
    # iter.plot = 2,
    # tree.step = 5
  )
  # necessary step to save properly (according to Carlson's github but it doesn't work so do all this in one session)
  invisible(caldio.sdm$fit$state) 
  CD_models[i] <- list(caldio.sdm)
  # saveRDS(calbor.sdm, file = "embarcadero_SDMs/final_outputs/CALBOR.BART.sdm.RDS")
  # savedCALBOR.sdm <- readRDS("embarcadero_SDMs/outputs/CALBOR.BART.sdm.RDS")
  
  ## 2.3. Obtain outputs ####
  
  ### Make a definitive prediction ####
  
  # First we obtain the predicted mean and quantiles
  CD_prediction <- predict(object = caldio.sdm,
                           x.layers = predictors_scaled,
                           quantiles =c(0.5),
                           splitby = 20,
                           quiet = TRUE)
  plot(CD_prediction)
  
  # save for definitive plotting
  CD_stack <- stack(CD_stack, CD_prediction)
}


CD_stack2 <- dropLayer(CD_stack, c(1,3,5,7,9,11,13,15,17,19))
plot(CD_stack2)

writeRaster(CD_stack2, 
            filename = "embarcadero_SDMs/iter_outputs/CALDIO_pred_stack_2024.grd", 
            overwrite = T)

saveRDS(CD_models, file = "embarcadero_SDMs/iter_outputs/CALDIO_models.RDS")
CD_models <- readRDS("embarcadero_SDMs/iter_outputs/CALDIO_models.RDS")

### Variable importance ####

# We can now generate a plot of the importance of the included variables
varimp_list <- list()

for(i in seq_along(CD_models)) {
  # i = 1
  mod_i <- CD_models[[i]]
  x <- varimp(mod_i, plots = F)
  varimp_i <- x %>% 
    mutate(iter = paste0("iter_", i))
  varimp_list[i] <- list(varimp_i) 
}


varimp_all <- bind_rows(varimp_list)
varimp_sum <- varimp_all %>% 
  group_by(names) %>% 
  summarise(mean = mean(varimps), 
            sd = sd(varimps)) %>% 
  ungroup()


ggplot(varimp_sum, aes(x = names)) + 
  geom_point(aes(y = mean)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean + sd)) + 
  theme_bw()

write.csv(varimp_sum, file = "embarcadero_SDMs/iter_outputs/CALDIO_varimp_iter.csv")

### Partial response ####

# These plots show the mean effect of each variable when measured at all other 
# values of the other predictors, for all the values of the variable of interest. 
# This will show us at which values of the covariate the dependence is higher.

# With CI you can add credible intervals, and with trace = TRUE you can see the 
# dependence plot for each of the trees (remember BARTs are regression trees), 
# with the thicker line representing the mean of all the trees. 
all_partials <- list()

for (i in seq_along(CD_models)) {
  # i = 1
  print(paste("iter", i, sep = "_"))
  caldio.sdm <- CD_models[[i]]
  
  chla_partial_cd <- partial(caldio.sdm,
                             x.vars = "logchla_lag3",
                             trace = FALSE,
                             ci = TRUE,
                             equal = TRUE,
                             smooth = 10, 
                             panels = F)
  
  chla_i <- chla_partial_cd[[1]]$data %>% 
    mutate(variable = "chla", 
           iter = paste0("iter_", i))
  
  chla_var_partial_cd <- partial(caldio.sdm,
                                 x.vars = "chla_var",
                                 trace = FALSE,
                                 ci = TRUE,
                                 equal = TRUE,
                                 smooth = 10, 
                                 panels = F)
  
  chla_var_i <- chla_var_partial_cd[[1]]$data %>% 
    mutate(variable = "chla_var", 
           iter = paste0("iter_", i))
  
  sst_partial_cd <- partial(caldio.sdm,
                            x.vars = "sst",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)
  
  sst_i <- sst_partial_cd[[1]]$data %>% 
    mutate(variable = "sst", 
           iter = paste0("iter_", i))
  
  
  sal_partial_cd <- partial(caldio.sdm,
                            x.vars = "sal",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)
  
  
  sal_i <- sal_partial_cd[[1]]$data %>% 
    mutate(variable = "sal", 
           iter = paste0("iter_", i))
  
  all_i <- bind_rows(
    chla_i,
    chla_var_i,
    sst_i,
    sal_i
  )
  
  all_partials[i] <- list(all_i)
}


 all_partials_df <- bind_rows(all_partials)


scaled_chla_vals <- scale(predictors_original$logchla_lag3[])
scaled_chlaVar_vals <- scale(predictors_original$chla_var[])
scaled_sst_vals <- scale(predictors_original$sst[])
scaled_sal_vals <- scale(predictors_original$sal[])


all_partials_df <- all_partials_df %>% 
  mutate(
    x2 = if_else(variable == "chla", x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center'), 
                 if_else(variable == "chla_var", x * attr(scaled_chlaVar_vals, 'scaled:scale') + attr(scaled_chlaVar_vals, 'scaled:center'), 
                         if_else(variable == "sst", x * attr(scaled_sst_vals, 'scaled:scale') + attr(scaled_sst_vals, 'scaled:center'), 
                                 x * attr(scaled_sal_vals, 'scaled:scale') + attr(scaled_sal_vals, 'scaled:center')))), 
    Species = "CALBOR")



ggplot(all_partials_df, aes(x = x2)) + 
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "lightblue") + 
  geom_line(aes(y = med)) + 
  theme_bw() + 
  facet_grid(iter ~ variable, scales = "free")

write.csv(all_partials_df, 
          file = "embarcadero_SDMs/iter_outputs/caldio_partial_df_iter.csv")


# 3. CALEDW ####

## 3.1 Prepare dataset ####

### Extract covariate values at presences ####

all.cov_list <- list()
for (i in seq_along(CALEDW)) {
  # i = 1
  CALEDW_sub <- CALEDW[[i]]
  
  CALEDW_sf <- CALEDW_sub %>% 
    st_as_sf(coords = c("Longitude", "Latitude"))
  
  pres.cov <- raster::extract(predictors_scaled, CALEDW_sf)
  
  ### Generate absences and covariate values at absences ####
  
  # From the vignette: 
  # BARTs are like boosted regression trees (BRTs) in that they are sensitive 
  # to assumed prevalence; anecdotally,I strongly suggest using an equal number 
  # of presences and absences in your training data. You can experiment with 
  # the demo data by changing “nrow(ticks)” to “5000” below if you want to see 
  # some odd model behavior
  
  absence <- randomPoints(predictors_scaled$bati, nrow(CALEDW_sf))
  abs.cov <- raster::extract(predictors_scaled, absence)
  
  absence <- randomPoints(predictors_scaled$bati, nrow(CALEDW_sf))
  abs.cov <- raster::extract(predictors_scaled, absence)
  
  ### Code the response ####
  
  #' Generate a vector with 1 for presences and 0 for absences (the response)
  pres.cov <- data.frame(pres.cov)
  pres.cov$pres <- 1
  
  abs.cov <- data.frame(abs.cov)
  abs.cov$pres <- 0
  
  # And one to bind them
  all.cov <-rbind(pres.cov, abs.cov)
  head(all.cov)
  table(all.cov$pres)
  
  # drop na's
  all.cov <- all.cov %>% 
    drop_na()
  
  # check that n is still more or less even between pres and abs
  table(all.cov$pres)
  all.cov_list[i] <- list(all.cov)
  
}



## 3.2. Run BART model ####

# Run a function that automatically does variable selection and runs the 
# fine tunned model

# covariate names without latitude
CE_stack <- stack()
CE_models <- list()
for(i in seq_along(all.cov_list)){
  # i = 1
  xvars <- names(predictors_scaled)[-8]
  
  all.cov_sub <- all.cov_list[[i]]
  CALEDW.sdm <-bart.step(
    x.data = all.cov_sub[,xvars],
    y.data = all.cov_sub[,'pres'],
    full = TRUE,
    quiet = F 
    # iter.step = 10,
    # iter.plot = 2,
    # tree.step = 5
  )
  # necessary step to save properly (according to Carlson's github but it doesn't work so do all this in one session)
  invisible(CALEDW.sdm$fit$state) 
  CE_models[i] <- list(CALEDW.sdm)
  # saveRDS(calbor.sdm, file = "embarcadero_SDMs/final_outputs/CALBOR.BART.sdm.RDS")
  # savedCALBOR.sdm <- readRDS("embarcadero_SDMs/outputs/CALBOR.BART.sdm.RDS")
  
  ## 3.3. Obtain outputs ####
  
  ### Make a definitive prediction ####
  
  # First we obtain the predicted mean and quantiles
  CE_prediction <- predict(object = CALEDW.sdm,
                           x.layers = predictors_scaled,
                           quantiles =c(0.5),
                           splitby = 20,
                           quiet = TRUE)
  plot(CE_prediction)
  
  # save for definitive plotting
  CE_stack <- stack(CE_stack, CE_prediction)
}


CE_stack2 <- dropLayer(CE_stack, c(1,3,5,7,9,11,13,15,17,19))
plot(CE_stack2)

writeRaster(CE_stack2, 
            filename = "embarcadero_SDMs/iter_outputs/CALEDW_pred_stack_2024.grd", 
            overwrite = T)

saveRDS(CE_models, file = "embarcadero_SDMs/iter_outputs/CALEDW_models.RDS")
CE_models <- readRDS("embarcadero_SDMs/iter_outputs/CALEDW_models.RDS")

### Variable importance ####

# We can now generate a plot of the importance of the included variables
varimp_list <- list()

for(i in seq_along(CE_models)) {
  # i = 1
  mod_i <- CE_models[[i]]
  x <- varimp(mod_i, plots = F)
  varimp_i <- x %>% 
    mutate(iter = paste0("iter_", i))
  varimp_list[i] <- list(varimp_i) 
}


varimp_all <- bind_rows(varimp_list)
varimp_sum <- varimp_all %>% 
  group_by(names) %>% 
  summarise(mean = mean(varimps), 
            sd = sd(varimps)) %>% 
  ungroup()


ggplot(varimp_sum, aes(x = names)) + 
  geom_point(aes(y = mean)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean + sd)) + 
  theme_bw()

write.csv(varimp_sum, file = "embarcadero_SDMs/iter_outputs/CALEDW_varimp_iter.csv")

### Partial response ####

# These plots show the mean effect of each variable when measured at all other 
# values of the other predictors, for all the values of the variable of interest. 
# This will show us at which values of the covariate the dependence is higher.

# With CI you can add credible intervals, and with trace = TRUE you can see the 
# dependence plot for each of the trees (remember BARTs are regression trees), 
# with the thicker line representing the mean of all the trees. 
all_partials <- list()

for (i in seq_along(CE_models)) {
  # i = 1
  print(paste("iter", i, sep = "_"))
  CALEDW.sdm <- CE_models[[i]]
  
  chla_partial_CE <- partial(CALEDW.sdm,
                             x.vars = "logchla_lag3",
                             trace = FALSE,
                             ci = TRUE,
                             equal = TRUE,
                             smooth = 10, 
                             panels = F)
  
  chla_i <- chla_partial_CE[[1]]$data %>% 
    mutate(variable = "chla", 
           iter = paste0("iter_", i))
  
  chla_var_partial_CE <- partial(CALEDW.sdm,
                                 x.vars = "chla_var",
                                 trace = FALSE,
                                 ci = TRUE,
                                 equal = TRUE,
                                 smooth = 10, 
                                 panels = F)
  
  chla_var_i <- chla_var_partial_CE[[1]]$data %>% 
    mutate(variable = "chla_var", 
           iter = paste0("iter_", i))
  
  sst_partial_CE <- partial(CALEDW.sdm,
                            x.vars = "sst",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)
  
  sst_i <- sst_partial_CE[[1]]$data %>% 
    mutate(variable = "sst", 
           iter = paste0("iter_", i))
  
  
  sal_partial_CE <- partial(CALEDW.sdm,
                            x.vars = "sal",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)
  
  
  sal_i <- sal_partial_CE[[1]]$data %>% 
    mutate(variable = "sal", 
           iter = paste0("iter_", i))
  
  all_i <- bind_rows(
    chla_i,
    chla_var_i,
    sst_i,
    sal_i
  )
  
  all_partials[i] <- list(all_i)
}


all_partials_df <- bind_rows(all_partials)


scaled_chla_vals <- scale(predictors_original$logchla_lag3[])
scaled_chlaVar_vals <- scale(predictors_original$chla_var[])
scaled_sst_vals <- scale(predictors_original$sst[])
scaled_sal_vals <- scale(predictors_original$sal[])


all_partials_df <- all_partials_df %>% 
  mutate(
    x2 = if_else(variable == "chla", x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center'), 
                 if_else(variable == "chla_var", x * attr(scaled_chlaVar_vals, 'scaled:scale') + attr(scaled_chlaVar_vals, 'scaled:center'), 
                         if_else(variable == "sst", x * attr(scaled_sst_vals, 'scaled:scale') + attr(scaled_sst_vals, 'scaled:center'), 
                                 x * attr(scaled_sal_vals, 'scaled:scale') + attr(scaled_sal_vals, 'scaled:center')))), 
    Species = "CALBOR")



ggplot(all_partials_df, aes(x = x2)) + 
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "lightblue") + 
  geom_line(aes(y = med)) + 
  theme_bw() + 
  facet_grid(iter ~ variable, scales = "free")

write.csv(all_partials_df, 
          file = "embarcadero_SDMs/iter_outputs/CALEDW_partial_df_iter.csv")

# 4. CALONECTRIS ####

## 4.1 Prepare dataset ####

### Extract covariate values at presences ####

all.cov_list <- list()
for (i in seq_along(CALEDW)) {
  # i = 1
  CALBOR_sub <- CALBOR[[i]]
  CALDIO_sub <- CALDIO[[i]]
  CALEDW_sub <- CALEDW[[i]]
  
  CAL <- bind_rows(CALBOR_sub, CALDIO_sub, CALEDW_sub)
  CAL_sf <- CAL %>% 
    st_as_sf(coords = c("Longitude", "Latitude"))
  CAL_sf <- CAL_sf %>% 
    sample_frac(size = 1/3)
  
  pres.cov <- raster::extract(predictors_scaled, CAL_sf)
  
  ### Generate absences and covariate values at absences ####
  
  # From the vignette: 
  # BARTs are like boosted regression trees (BRTs) in that they are sensitive 
  # to assumed prevalence; anecdotally,I strongly suggest using an equal number 
  # of presences and absences in your training data. You can experiment with 
  # the demo data by changing “nrow(ticks)” to “5000” below if you want to see 
  # some odd model behavior
  
  absence <- randomPoints(predictors_scaled$bati, nrow(CAL_sf))
  abs.cov <- raster::extract(predictors_scaled, absence)
  
  absence <- randomPoints(predictors_scaled$bati, nrow(CAL_sf))
  abs.cov <- raster::extract(predictors_scaled, absence)
  
  ### Code the response ####
  
  #' Generate a vector with 1 for presences and 0 for absences (the response)
  pres.cov <- data.frame(pres.cov)
  pres.cov$pres <- 1
  
  abs.cov <- data.frame(abs.cov)
  abs.cov$pres <- 0
  
  # And one to bind them
  all.cov <-rbind(pres.cov, abs.cov)
  head(all.cov)
  table(all.cov$pres)
  
  # drop na's
  all.cov <- all.cov %>% 
    drop_na()
  
  # check that n is still more or less even between pres and abs
  table(all.cov$pres)
  all.cov_list[i] <- list(all.cov)
  
}



## 4.2. Run BART model ####

# Run a function that automatically does variable selection and runs the 
# fine tunned model

# covariate names without latitude
C_stack <- stack()
C_models <- list()
for(i in seq_along(all.cov_list)){
  # i = 1
  print(paste("iter", i, sep = " "))
  xvars <- names(predictors_scaled)[-8]
  
  all.cov_sub <- all.cov_list[[i]]
  CAL.sdm <-bart.step(
    x.data = all.cov_sub[,xvars],
    y.data = all.cov_sub[,'pres'],
    full = TRUE,
    quiet = F 
    # iter.step = 10,
    # iter.plot = 2,
    # tree.step = 5
  )
  # necessary step to save properly (according to Carlson's github but it doesn't work so do all this in one session)
  invisible(CAL.sdm$fit$state) 
  C_models[i] <- list(CAL.sdm)
  # saveRDS(calbor.sdm, file = "embarcadero_SDMs/final_outputs/CALBOR.BART.sdm.RDS")
  # savedCALBOR.sdm <- readRDS("embarcadero_SDMs/outputs/CALBOR.BART.sdm.RDS")
  
  ## 3.3. Obtain outputs ####
  
  ### Make a definitive prediction ####
  
  # First we obtain the predicted mean and quantiles
  C_prediction <- predict(object = CAL.sdm,
                          x.layers = predictors_scaled,
                          quantiles =c(0.5),
                          splitby = 20,
                          quiet = TRUE)
  
  # save for definitive plotting
  C_stack <- stack(C_stack, C_prediction)
}


C_stack2 <- dropLayer(C_stack, c(1,3,5,7,9,11,13,15,17,19))
plot(C_stack2)

writeRaster(C_stack2, 
            filename = "embarcadero_SDMs/iter_outputs/CAL_pred_stack_2024.grd", 
            overwrite = T)

saveRDS(C_models, file = "embarcadero_SDMs/iter_outputs/CAL_models.RDS")
C_models <- readRDS("embarcadero_SDMs/iter_outputs/CAL_models.RDS")

### Variable importance ####

# We can now generate a plot of the importance of the included variables
varimp_list <- list()

for(i in seq_along(C_models)) {
  # i = 1
  mod_i <- C_models[[i]]
  x <- varimp(mod_i, plots = F)
  varimp_i <- x %>% 
    mutate(iter = paste0("iter_", i))
  varimp_list[i] <- list(varimp_i) 
}


varimp_all <- bind_rows(varimp_list)
varimp_sum <- varimp_all %>% 
  group_by(names) %>% 
  summarise(mean = mean(varimps), 
            sd = sd(varimps)) %>% 
  ungroup()


ggplot(varimp_sum, aes(x = names)) + 
  geom_point(aes(y = mean)) + 
  geom_errorbar(aes(ymin = mean-sd, ymax = mean + sd)) + 
  theme_bw()

write.csv(varimp_sum, file = "embarcadero_SDMs/iter_outputs/CAL_varimp_iter.csv")

### Partial response ####

# These plots show the mean effect of each variable when measured at all other 
# values of the other predictors, for all the values of the variable of interest. 
# This will show us at which values of the covariate the dependence is higher.

# With CI you can add credible intervals, and with trace = TRUE you can see the 
# dependence plot for each of the trees (remember BARTs are regression trees), 
# with the thicker line representing the mean of all the trees. 
all_partials <- list()

for (i in seq_along(C_models)) {
  # i = 1
  print(paste("iter", i, sep = "_"))
  CAL.sdm <- C_models[[i]]
  
  chla_partial_C <- partial(CAL.sdm,
                             x.vars = "logchla_lag3",
                             trace = FALSE,
                             ci = TRUE,
                             equal = TRUE,
                             smooth = 10, 
                             panels = F)
  
  chla_i <- chla_partial_C[[1]]$data %>% 
    mutate(variable = "chla", 
           iter = paste0("iter_", i))
  
  chla_var_partial_C <- partial(CAL.sdm,
                                 x.vars = "chla_var",
                                 trace = FALSE,
                                 ci = TRUE,
                                 equal = TRUE,
                                 smooth = 10, 
                                 panels = F)
  
  chla_var_i <- chla_var_partial_C[[1]]$data %>% 
    mutate(variable = "chla_var", 
           iter = paste0("iter_", i))
  
  sst_partial_C <- partial(CAL.sdm,
                            x.vars = "sst",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)
  
  sst_i <- sst_partial_C[[1]]$data %>% 
    mutate(variable = "sst", 
           iter = paste0("iter_", i))
  
  
  sal_partial_C <- partial(CAL.sdm,
                            x.vars = "sal",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)
  
  
  sal_i <- sal_partial_C[[1]]$data %>% 
    mutate(variable = "sal", 
           iter = paste0("iter_", i))
  
  all_i <- bind_rows(
    chla_i,
    chla_var_i,
    sst_i,
    sal_i
  )
  
  all_partials[i] <- list(all_i)
}


all_partials_df <- bind_rows(all_partials)


scaled_chla_vals <- scale(predictors_original$logchla_lag3[])
scaled_chlaVar_vals <- scale(predictors_original$chla_var[])
scaled_sst_vals <- scale(predictors_original$sst[])
scaled_sal_vals <- scale(predictors_original$sal[])


all_partials_df <- all_partials_df %>% 
  mutate(
    x2 = if_else(variable == "chla", x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center'), 
                 if_else(variable == "chla_var", x * attr(scaled_chlaVar_vals, 'scaled:scale') + attr(scaled_chlaVar_vals, 'scaled:center'), 
                         if_else(variable == "sst", x * attr(scaled_sst_vals, 'scaled:scale') + attr(scaled_sst_vals, 'scaled:center'), 
                                 x * attr(scaled_sal_vals, 'scaled:scale') + attr(scaled_sal_vals, 'scaled:center')))))



ggplot(all_partials_df, aes(x = x2)) + 
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "lightblue") + 
  geom_line(aes(y = med)) + 
  theme_bw() + 
  facet_grid(iter ~ variable, scales = "free")

write.csv(all_partials_df, 
          file = "embarcadero_SDMs/iter_outputs/CAL_partial_df_iter.csv")
