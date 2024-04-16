
# 0. House Keeping ####
rm(list = ls())
setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")

pacman::p_load(embarcadero, tidyverse, sf, raster, dismo, spatstat)
                    
'%!in%' <- function(x,y)!('%in%'(x,y))

source("scripts/migcon_AUX.R")

# 0. Load all datasets ####
predictors_original <- stack("data/env_data/wintering/all_vars.grd")

predictors_original <- stack(predictors_original, init(predictors_original, 'y'))

names(predictors_original)[2] <- "chla_var"
names(predictors_original)[8] <- "latitude"

predictors_scaled <- scale(predictors_original)

CALBOR <- readRDS("embarcadero_SDMs/data/CALBOR.RDS")
CALDIO <- readRDS("embarcadero_SDMs/data/CALDIO.RDS")
CALEDW <- readRDS("embarcadero_SDMs/data/CALEDW.RDS")

# 1. CALBOR ####

## 1.1 Prepare dataset ####

### Extract covariate values at presences ####

CALBOR_sf <- CALBOR %>% 
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

## 1.2. Run BART model ####

# Run a function that automatically does variable selection and runs the 
# fine tunned model

# covariate names without latitude
xvars <- names(predictors_scaled)[-8]

calbor.sdm <-bart.step(
  x.data = all.cov[,xvars],
  y.data = all.cov[,'pres'],
  full = TRUE,
  quiet = F 
  # iter.step = 10,
  # iter.plot = 2,
  # tree.step = 5
)

# necessary step to save properly (according to Carlson's github but it doesn't work so do all this in one session)
invisible(calbor.sdm$fit$state) 
# saveRDS(calbor.sdm, file = "embarcadero_SDMs/final_outputs/CALBOR.BART.sdm.RDS")
# savedCALBOR.sdm <- readRDS("embarcadero_SDMs/outputs/CALBOR.BART.sdm.RDS")

## 1.3. Obtain outputs ####

### Make a definitive prediction ####

# First we obtain the predicted mean and quantiles
CB_prediction <- predict(object = calbor.sdm,
                         x.layers = predictors_scaled,
                         quantiles =c(0.025, 0.975),
                         splitby = 20,
                         quiet = TRUE)
plot(CB_prediction)
names(CB_prediction) <- c("Mean", "q2.5%", "q97.5%")

# save for definitive plotting
writeRaster(CB_prediction, 
            filename = "embarcadero_SDMs/final_outputs/CALBOR_pred_stack_2024.grd", 
            overwrite = T)

# We check how the mean prediction and the two quantiles look

# CB_prediction <- stack("embarcadero_SDMs/outputs/CALBOR_pred_stack.grd")
par(mfrow = c(2,2))
plot(CB_prediction[[1]],
     box = FALSE,
     axes = FALSE,
     main ='CALBOR prediction mean',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CB_prediction[[2]],
     box = FALSE,
     axes = FALSE,
     main ='CALBOR prediction 2.5%',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CB_prediction[[3]],
     box = FALSE,
     axes = FALSE,
     main ='CALBOR prediction 97.5%',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CB_prediction[[3]] - CB_prediction[[2]],
     box = FALSE,
     axes = FALSE,
     main ='CALBOR prediction CI width',
     zlim =c(0,0.75),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'CI width', 
                        side = 2, line = 1.3))
par(mfrow = c(1,1))

# We now produce a binary map using as cutoff the value suggested by the model summary

plot(CB_prediction[[1]] > 0.5104638  , # cutoff obtained from summary above
     box = FALSE,
     axes = FALSE, 
     main = 'Binary prediction', 
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Presence - absence', 
                        side = 2, line = 1.3))

# We now produce a map of high uncertainty, plotting a binary map of the areas where the uncertainty (CI width) is in the >75% quantile

quant <- quantile(values(CB_prediction[[3]] - CB_prediction[[2]]), 0.75,
                  na.rm = TRUE)

plot((CB_prediction[[3]] - CB_prediction[[2]]) > quant,
     box = FALSE,
     axes = FALSE,
     main = "Highest uncertainty zones",
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'CI width', 
                        side = 2, line = 1.3))

### Variable importance ####

# We can now generate a plot of the importance of the included variables
varimp(calbor.sdm, plots = TRUE)

varimp <- varimp(calbor.sdm, plots = F)

varimp_ordered <- varimp %>% 
  arrange(-varimps)

write.csv(varimp_ordered, file = "embarcadero_SDMs/final_outputs/CALBOR_varimp.csv")

### Partial response ####

# These plots show the mean effect of each variable when measured at all other 
# values of the other predictors, for all the values of the variable of interest. 
# This will show us at which values of the covariate the dependence is higher.

# With CI you can add credible intervals, and with trace = TRUE you can see the 
# dependence plot for each of the trees (remember BARTs are regression trees), 
# with the thicker line representing the mean of all the trees. 

chla_partial_cb <- partial(calbor.sdm,
                        x.vars = "logchla_lag3",
                        trace = FALSE,
                        ci = TRUE,
                        equal = TRUE,
                        smooth = 10, 
                        panels = F)

scaled_chla_vals <- scale(predictors_original$logchla_lag3[])

chla_part_df_cb <- chla_partial_cb[[1]]$data %>% 
  mutate(
    x2_log = x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center'),
    x2 = exp(x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center')), 
    Species = "CALBOR", 
    Variable = "Chlorophyll A")


chla_var_partial_cb <- partial(calbor.sdm,
                            x.vars = "chla_var",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)

scaled_chlaVar_vals <- scale(predictors_original$chla_var[])

chlaVar_part_df_cb <- chla_var_partial_cb[[1]]$data %>% 
  mutate(
    x2 = x * attr(scaled_chlaVar_vals, 'scaled:scale') + attr(scaled_chlaVar_vals, 'scaled:center'), 
    Species = "CALBOR", 
    Variable = "Chlorophyll A variability")

sst_partial_cb <- partial(calbor.sdm,
                            x.vars = "sst",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)

scaled_sst_vals <- scale(predictors_original$sst[])

sst_part_df_cb <- sst_partial_cb[[1]]$data %>% 
  mutate(x2 = x * attr(scaled_sst_vals, 'scaled:scale') + attr(scaled_sst_vals, 'scaled:center'), 
         Species = "CALBOR", 
         Variable = "Sea surface temperature")



calbor_partials <- bind_rows(sst_part_df_cb, chla_part_df_cb, chlaVar_part_df_cb)

write.csv(calbor_partials, 
        file = "embarcadero_SDMs/final_outputs/calbor_partial_df.csv")

# 2. CALDIO ####

## 2.1 Prepare dataset ####

### Extract covariate values at presences ####

CALDIO_sf <- CALDIO %>% 
  st_as_sf(coords = c("Longitude", "Latitude"))

pres.cov <- raster::extract(predictors_scaled, CALDIO_sf)

### Generate absences ####

absence <- randomPoints(predictors_scaled, nrow(CALDIO_sf))
abs.cov <- raster::extract(predictors_scaled, absence)

### Code the response ####
pres.cov <- data.frame(pres.cov)
pres.cov$pres <- 1

abs.cov <- data.frame(abs.cov)
abs.cov$pres <- 0

all.cov <-rbind(pres.cov, abs.cov)
head(all.cov)
table(all.cov$pres)

all.cov <- all.cov %>% 
  drop_na()

table(all.cov$pres)

## 2.2. Run BART model ####

xvars <-names(predictors_scaled)[-8]

caldio.sdm <-bart.step(x.data = all.cov[,xvars],
                y.data = all.cov[,'pres'],
                full = TRUE,
                quiet = T)

## 2.3. Obtain outputs ####

### Make a definitive prediction ####

CD_prediction <- predict(object = caldio.sdm,
                         x.layers = predictors_scaled,
                         quantiles =c(0.025, 0.975),
                         splitby = 20,
                         quiet = TRUE)
plot(CD_prediction)
names(CD_prediction) <- c("Mean", "q2.5%", "q97.5%")
# writeRaster(CD_prediction, filename = "embarcadero_SDMs/final_outputs/CALDIO_pred_stack_2024.grd", 
#             overwrite = T)

# We check how the mean prediction and the two quantiles look

# CD_prediction <- stack("embarcadero_SDMs/outputs/CALBOR_pred_stack_2024.grd")
par(mfrow = c(2,2))
plot(CD_prediction[[1]],
     box = FALSE,
     axes = FALSE,
     main ='CALDIO prediction mean',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CD_prediction[[2]],
     box = FALSE,
     axes = FALSE,
     main ='CALDIO prediction 2.5%',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CD_prediction[[3]],
     box = FALSE,
     axes = FALSE,
     main ='CALDIO prediction 97.5%',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CD_prediction[[3]] - CD_prediction[[2]],
     box = FALSE,
     axes = FALSE,
     main ='CALDIO prediction CI width',
     zlim =c(0,0.75),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'CI width', 
                        side = 2, line = 1.3))
par(mfrow = c(1,1))

# We now produce a binary map using as cutoff the value suggested by the model summary

plot(CD_prediction[[1]] > 0.4820289, # cutoff obtained from summary above
     box = FALSE,
     axes = FALSE, 
     main = 'Binary prediction', 
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Presence - absence', 
                        side = 2, line = 1.3))

# We now produce a map of high uncertainty, plotting a binary map of the areas where the uncertainty (CI width) is in the >75% quantile

quant <- quantile(values(CD_prediction[[3]] - CD_prediction[[2]]), 0.75,
                  na.rm = TRUE)

plot((CD_prediction[[3]] - CD_prediction[[2]]) > quant,
     box = FALSE,
     axes = FALSE,
     main = "Highest uncertainty zones",
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'CI width', 
                        side = 2, line = 1.3))

### Variable importance ####
#' 
# We can now generate a plot of the importance of the included variables

varimp(caldio.sdm, plots = TRUE)

varimp <- varimp(caldio.sdm, plots = F)

varimp_ordered <- varimp %>% 
  arrange(-varimps)

# write.csv(varimp_ordered, file = "embarcadero_SDMs/final_outputs/CALDIO_varimp.csv")

### Partial response ####

chla_var_partial_cd <- partial(caldio.sdm,
                            x.vars = "chla_var",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)

chlaVar_part_df_cd <- chla_var_partial_cd[[1]]$data %>% 
  mutate(
    x2 = x * attr(scaled_chlaVar_vals, 'scaled:scale') + attr(scaled_chlaVar_vals, 'scaled:center'), 
    Species = "CALDIO", 
    Variable = "Chlorophyll A variability")

chla_partial_cd <- partial(caldio.sdm,
                        x.vars = "logchla_lag3",
                        trace = FALSE,
                        ci = TRUE,
                        equal = TRUE,
                        smooth = 10, 
                        panels = F)

chla_part_df_cd <- chla_partial_cd[[1]]$data %>% 
  mutate(
    x2_log = x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center'),
    x2 = exp(x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center')), 
    Species = "CALDIO", 
    Variable = "Chlorophyll A")

sst_partial_cd <- partial(caldio.sdm,
                            x.vars = "sst",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)

sst_part_df_cd <- sst_partial_cd[[1]]$data %>% 
  mutate(
    x2 = x * attr(scaled_sst_vals, 'scaled:scale') + attr(scaled_sst_vals, 'scaled:center'), 
    Species = "CALDIO", 
    Variable = "Sea surface temperature")

caldio_partials <- bind_rows(sst_part_df_cd, chla_part_df_cd, chlaVar_part_df_cd)
write.csv(caldio_partials, 
          file = "embarcadero_SDMs/final_outputs/caldio_partial_df.csv")


# 3. CALEDW ####

## 3.1 Prepare dataset ####

### Extract covariate values at presences ####

CALEDW_sf <- CALEDW %>% 
  st_as_sf(coords = c("Longitude", "Latitude"))

pres.cov <- raster::extract(predictors_scaled, CALEDW_sf)

### Generate absences ####

absence <- randomPoints(predictors_scaled, nrow(CALEDW_sf))
abs.cov <- raster::extract(predictors_scaled, absence)

### Code the response ####

pres.cov <- data.frame(pres.cov)
pres.cov$pres <- 1

abs.cov <- data.frame(abs.cov)
abs.cov$pres <- 0

all.cov <-rbind(pres.cov, abs.cov)

table(all.cov$pres)

all.cov %>% 
  drop_na -> all.cov
table(all.cov$pres)

## 3.2. Run BART model ####
xvars <- names(predictors_scaled)[-8]
caledw.sdm <-bart.step(x.data = all.cov[,xvars],
                y.data = all.cov[,'pres'],
                full = TRUE,
                quiet = F, 
                iter.step = 100, 
                iter.plot = 100, 
                tree.step = 10)

# saveRDS(caledw.sdm, file = "embarcadero_SDMs/final_outputs/CALEDW.BART.sdm.RDS")

## 3.3. Obtain outputs ####
### Make a definitive prediction ####

CE_prediction <- predict(object = caledw.sdm,
                         x.layers = predictors_scaled,
                         quantiles =c(0.025, 0.975),
                         splitby = 20,
                         quiet = TRUE)
plot(CE_prediction)
names(CE_prediction) <- c("Mean", "q2.5%", "q97.5%")
# writeRaster(CE_prediction, filename = "embarcadero_SDMs/final_outputs/CALEDW_pred_stack_2024.grd", 
#             overwrite = T)

# CE_prediction <- stack("embarcadero_SDMs/final_outputs/CALBOR_pred_stack_2024.grd")
par(mfrow = c(2,2))
plot(CE_prediction[[1]],
     box = FALSE,
     axes = FALSE,
     main ='CALEDW prediction mean',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CE_prediction[[2]],
     box = FALSE,
     axes = FALSE,
     main ='CALEDW prediction 2.5%',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CE_prediction[[3]],
     box = FALSE,
     axes = FALSE,
     main ='CALEDW prediction 97.5%',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability', 
                        side = 2, line = 1.3))

plot(CE_prediction[[3]] - CE_prediction[[2]],
     box = FALSE,
     axes = FALSE,
     main ='CALEDW prediction CI width',
     zlim =c(0,0.75),
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'CI width', 
                        side = 2, line = 1.3))
par(mfrow = c(1,1))

plot(CE_prediction[[1]] > 0.444, # cutoff obtained from summary above
     box = FALSE,
     axes = FALSE, 
     main = 'Binary prediction', 
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'Presence - absence', 
                        side = 2, line = 1.3))

quant <- quantile(values(CE_prediction[[3]] - CE_prediction[[2]]), 0.75,
                  na.rm = TRUE)

plot((CE_prediction[[3]] - CE_prediction[[2]]) > quant,
     box = FALSE,
     axes = FALSE,
     main = "Highest uncertainty zones",
     axis.args = list(at = pretty(0:1), 
                      labels = pretty(0:1)),
     legend.args = list(text = 'CI width', 
                        side = 2, line = 1.3))

### Variable importance ####
varimp(caledw.sdm, plots = TRUE)

varimp <- varimp(caledw.sdm, plots = F)

varimp_ordered <- varimp %>% 
  arrange(-varimps)

write.csv(varimp_ordered, file = "embarcadero_SDMs/final_outputs/CALEDW_varimp.csv")

### Partial response ####

chla_var_partial_ce <- partial(caledw.sdm,
                            x.vars = "chla_var",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10, 
                            panels = F)

chlaVar_part_df_ce <- chla_var_partial_ce[[1]]$data %>% 
  mutate(
    x2 = x * attr(scaled_chlaVar_vals, 'scaled:scale') + attr(scaled_chlaVar_vals, 'scaled:center'), 
    Species = "CALEDW", 
    Variable = "Chlorophyll A variability")

chla_partial_ce <- partial(caledw.sdm,
                        x.vars = "logchla_lag3",
                        trace = FALSE,
                        ci = TRUE,
                        equal = TRUE,
                        smooth = 10, 
                        panels = F)

chla_part_df_ce <- chla_partial_ce[[1]]$data %>% 
  mutate(
    x2_log = x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center'),
    x2 = exp(x * attr(scaled_chla_vals, 'scaled:scale') + attr(scaled_chla_vals, 'scaled:center')), 
    Species = "CALEDW", 
    Variable = "Chlorophyll A")

sst_partial_ce <- partial(caledw.sdm,
                       x.vars = "sst",
                       trace = FALSE,
                       ci = TRUE,
                       equal = TRUE,
                       smooth = 10, 
                       panels = F)

sst_part_df_ce <- sst_partial_ce[[1]]$data %>% 
  mutate(
    x2 = x * attr(scaled_sst_vals, 'scaled:scale') + attr(scaled_sst_vals, 'scaled:center'),
    Species = "CALEDW", 
    Variable = "Sea surface temperature")

caledw_partials <- bind_rows(sst_part_df_ce, chla_part_df_ce, chlaVar_part_df_ce)

write.csv(caledw_partials, 
          file = "embarcadero_SDMs/final_outputs/caledw_partial_df.csv")




# 4. All Calonectris ####

## 4.1 Prepare dataset ####

### Extract covariate values at presences ####

Calonectris <- readRDS("embarcadero_SDMs/data/all_calonectris.RDS")
Calonectris_sf <- Calonectris %>%
  st_as_sf(coords = c("Longitude", "Latitude"))


pres.cov <- raster::extract(predictors_scaled, Calonectris_sf)

### Generate absences ####
absence <- randomPoints(predictors_scaled, nrow(Calonectris_sf))
abs.cov <- raster::extract(predictors_scaled, absence)

### Code the response ####

pres.cov <- data.frame(pres.cov)
pres.cov$pres <- 1

abs.cov <- data.frame(abs.cov)
abs.cov$pres <- 0

all.cov <-rbind(pres.cov, abs.cov)
head(all.cov)
table(all.cov$pres)

all.cov %>%
  drop_na -> all.cov
table(all.cov$pres)

## 4.2. Run BART model ####

xvars <-names(predictors_scaled)[-8]
calonectris.sdm <-bart.step(x.data = all.cov[,xvars],
                y.data = all.cov[,'pres'],
                full = TRUE,
                quiet = F,
                iter.step = 100,
                iter.plot = 100,
                tree.step = 10)

saveRDS(calonectris.sdm, file = "embarcadero_SDMs/final_outputs/Calonectris.BART.sdm.RDS")

## 4.3. Obtain outputs ####

### Make a definitive prediction ####

C_prediction <- predict(object = calonectris.sdm,
                        x.layers = predictors_scaled,
                        quantiles =c(0.025, 0.975),
                        splitby = 20,
                        quiet = TRUE)
plot(C_prediction)
names(C_prediction) <- c("Mean", "q2.5%", "q97.5%")
writeRaster(C_prediction, filename = "embarcadero_SDMs/final_outputs/Calonectris_pred_stack.grd",
            overwrite = T)

C_prediction <- stack("embarcadero_SDMs/outputs/Calonectris_pred_stack.grd")
par(mfrow = c(2,2))
plot(C_prediction[[1]],
     box = FALSE,
     axes = FALSE,
     main ='Calonectris prediction mean',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1),
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability',
                        side = 2, line = 1.3))

plot(C_prediction[[2]],
     box = FALSE,
     axes = FALSE,
     main ='Calonectris prediction 2.5%',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1),
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability',
                        side = 2, line = 1.3))

plot(C_prediction[[3]],
     box = FALSE,
     axes = FALSE,
     main ='Calonectris prediction 97.5%',
     zlim =c(0,1),
     axis.args = list(at = pretty(0:1),
                      labels = pretty(0:1)),
     legend.args = list(text = 'Probability',
                        side = 2, line = 1.3))

plot(C_prediction[[3]] - C_prediction[[2]],
     box = FALSE,
     axes = FALSE,
     main ='Calonectris prediction CI width',
     zlim =c(0,0.75),
     axis.args = list(at = pretty(0:1),
                      labels = pretty(0:1)),
     legend.args = list(text = 'CI width',
                        side = 2, line = 1.3))
par(mfrow = c(1,1))

plot(C_prediction[[1]] > 0.5366691 , # cutoff obtained from summary above
     box = FALSE,
     axes = FALSE,
     main = 'Binary prediction',
     axis.args = list(at = pretty(0:1),
                      labels = pretty(0:1)),
     legend.args = list(text = 'Presence - absence',
                        side = 2, line = 1.3))

quant <- quantile(values(C_prediction[[3]] - C_prediction[[2]]), 0.75,
                  na.rm = TRUE)

plot((C_prediction[[3]] - C_prediction[[2]]) > quant,
     box = FALSE,
     axes = FALSE,
     main = "Highest uncertainty zones",
     axis.args = list(at = pretty(0:1),
                      labels = pretty(0:1)),
     legend.args = list(text = 'CI width',
                        side = 2, line = 1.3))

### Variable importance ####
x <- varimp(calonectris.sdm, plots = TRUE)

varimp <- varimp(calonectris.sdm, plots = F)

varimp_ordered <- varimp %>%
  arrange(-varimps)

write.csv(varimp_ordered, file = "embarcadero_SDMs/final_outputs/Calonectris_varimp.csv")

### Partial response ####

bati_partial <- partial(calonectris.sdm,
                        x.vars = "bati",
                        trace = FALSE,
                        ci = TRUE,
                        equal = TRUE,
                        smooth = 10,
                        panels = F)
p5 <- bati_partial[[1]] + ggtitle("Bathimetry")

chla_var_partial <- partial(calonectris.sdm,
                            x.vars = "chla_var",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10,
                            panels = F)
p2 <- chla_var_partial[[1]] + ggtitle("Chlorophyl A variation")

chla_partial <- partial(calonectris.sdm,
                        x.vars = "logchla_lag3",
                        trace = FALSE,
                        ci = TRUE,
                        equal = TRUE,
                        smooth = 10,
                        panels = F)
p3 <- chla_partial[[1]] + ggtitle("Chlorophyl A")


sst_partial <- partial(calonectris.sdm,
                            x.vars = "sst",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10,
                            panels = F)
p1 <- sst_partial[[1]] + ggtitle("Sea surface temperature")

sst_grad_partial <- partial(calonectris.sdm,
                            x.vars = "sst_grad",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10,
                            panels = F)
p6 <- sst_grad_partial[[1]] + ggtitle("Sea surface temperature gradient")

sal_partial <- partial(calonectris.sdm,
                            x.vars = "sal",
                            trace = FALSE,
                            ci = TRUE,
                            equal = TRUE,
                            smooth = 10,
                            panels = F)
p4 <- sal_partial[[1]] + ggtitle("Salinity")

all_calonectris_partials <- list(p1, p2, p3, p4, p5, p6)


ggarrange(plotlist = all_calonectris_partials)

saveRDS(all_calonectris_partials, file = "embarcadero_SDMs/final_outputs/Calonectris_partial_plots.RDS")


