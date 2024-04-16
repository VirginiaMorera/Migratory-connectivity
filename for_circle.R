
# 0. House keeping####

setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")
# devtools::install_github("SMBC-NZP/MigConnectivity")

pacman::p_load(
  tidyverse, rgeos, sp, sf, purrr, circlize, Cairo, mapplots
)
  
# source functions
source("scripts/migcon_AUX.R")


# 1. Load data ####

## representative bird tracks ####

#representative colonies
boot_out_summary <- read.csv("embarcadero_SDMs/data/colony_bootstrap_summary.csv",
                             row.names = NULL)

representative_colonies <- boot_out_summary %>%
  dplyr::filter(SampleSize > 4 & Representativity > 70) %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) 


# tracks
all_complete <- read.csv("data/complete_fenology.csv", row.names = NULL)

# tracks from representative colonies
all_complete %<>%
  dplyr::filter(unique_ID_Bird %!in% c("L38927_first", "6195227_2014",
                                       "6143080_2012")) %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) %>% 
  dplyr::filter(Colony %in% representative_colonies$Colony) %>% 
  mutate(col_sp = paste(Colony, Species, sep = "_")) %>%
  filter(col_sp != "Chafarinas_CALDIO") %>% 
  dplyr::select(Species, Colony, unique_ID_Bird, Date_Time, Lat, Lon, fen) 


## world map ####
load("data/projections.Rdata")

# worldmap
world_simple <- rworldmap::getMap(resolution = "high")
world_simple <- gBuffer(world_simple, byid = F, width = 0)
world_simplekM <- spTransform(world_simple, CRS(projections$EquidistConic))
world_simple_sf <- st_as_sf(world_simple)

## breeding colonies data ####

breeding_pops <- read.csv("data/Pop_Col_info.csv", 
                          row.names = NULL)

population_sizes <- breeding_pops %>% 
  group_by(Species, Population) %>% 
  dplyr::summarise(Tot.pop = sum(Min.pop))
  
colony_sizes <- breeding_pops %>% 
  filter(!is.na(n)) %>% 
  mutate(Tot.pop = 2*Min.pop) %>% 
  dplyr::select(Species, Population, Colony, Tot.pop, n_col = n) %>% 
  mutate(Colony= if_else(Colony == "Monta\xf1aClara", "MontañaClara", Colony)) 

colony_locs <- breeding_pops %>% 
  dplyr::select(Species, Colony, Lon, Lat)
  
population_locs <-  breeding_pops %>% 
  group_by(Species, Population) %>% 
  dplyr::summarise(
    mean.lat = mean(Lat),
    mean.lon = mean(Lon)
  )
                      

## wintering areas data ####

win_areas <- as_Spatial(read_sf(dsn = "data/shapefiles", 
                                layer = "wintering_areas_1"))


plot(win_areas)
invisible(text(coordinates(win_areas), labels = as.character(win_areas$ECOREGION), cex = 0.5))



# 2. Prepare data for circular plots ####

load("data/centroids.ALL.2023.Rdata")

## Assign win area based on wintering centroids ####
for_circle <- centroids.ALL %>% 
  pmap(data.frame) %>% 
  map_dfr(assign_win_centroid, polys = win_areas, level = 2) 

# remove individuals assigned to no area
table(for_circle$Win_area, useNA = "ifany")
for_circle <- for_circle[!is.na(for_circle$Win_area), ]

## Assign species to individual id in the for_circle datasets ####
birds_species <- all_complete %>% 
  dplyr::select(unique_ID_Bird, Species) %>% 
  distinct()

for_circle2 <- for_circle %>% 
  inner_join(birds_species) %>%
  left_join(colony_sizes, by = c("Colony", "Species")) %>%
  group_by(Species, Population, Colony, Win_area, Tot.pop, n_col) %>%
  dplyr::summarise(n_col_win = n()) %>% 
  mutate(ind_weight = Tot.pop/n_col) %>% 
  ungroup()

## create for_circle datasets for each level ####

# for_circle for colonies
for_circle_col <- for_circle2 %>% 
  dplyr::select(Species, Colony, Win_area, n_col_win, ind_weight) %>% 
  mutate(circle_arch_value = n_col_win*ind_weight)

# for_circle for populations
for_circle_pop <- for_circle2 %>% 
  group_by(Species, Population, Win_area) %>% 
  dplyr::summarise(
    Tot.pop = sum(Tot.pop), 
    n_pop = sum(n_col),
    n_pop_win = sum(n_col_win),
    ind_weight = Tot.pop/n_pop) %>% 
  ungroup() %>% 
  mutate(circle_arch_value = n_pop_win*ind_weight) %>% 
  dplyr::select(Species, Population, Win_area, n_pop_win, ind_weight, circle_arch_value)
  
# for circle for sp
for_circle_sp <- for_circle2 %>% 
  group_by(Species, Win_area) %>% 
  dplyr::summarise(
    Tot.pop = sum(Tot.pop), 
    n_pop = sum(n_col),
    n_pop_win = sum(n_col_win),
    ind_weight = Tot.pop/n_pop) %>% 
  ungroup() %>% 
  mutate(circle_arch_value = n_pop_win*ind_weight) %>% 
  dplyr::select(Species, Win_area, n_pop_win, ind_weight, circle_arch_value)


saveRDS(for_circle_col, file = "outputs/for_circle_col.2023.RDS")
saveRDS(for_circle_pop, file = "outputs/for_circle_pop.2023.RDS")
saveRDS(for_circle_sp, file = "outputs/for_circle_sp.2023.RDS")

for_circle_col <- readRDS("outputs/for_circle_col.2023.RDS")
for_circle_sp <- readRDS("outputs/for_circle_sp.2023.RDS")

# 3. By species #####

## prepare dataset for circular plot ####

for_circle_sp <- for_circle_sp %>% 
  ungroup() %>% 
  # mutate(circle_arch_value = sqrt(circle_arch_value)) %>%
  dplyr::select(-n_pop_win, -ind_weight) 
  



# default order
union(for_circle_sp[[1]], for_circle_sp[[2]])
length(union(for_circle_sp[[1]], for_circle_sp[[2]]))

## set order in which sectors will be plotted ####
new_order <- c("CALBOR", "CALEDW", "CALDIO", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Equatorial_Atlantic", "Gulf_of_Guinea_West", "Gulf_of_Guinea_Central",  "Angolan_Current", "Namib", "Namaqua",
               "Agulhas_Benguela_Confluence", "Agulhas_Current", 
               "Central_South_Atlantic", "North_Atlantic", 
               "PatagonianShelf_Falklands", "Uruguay_BuenosAires_Shelf", "SouthBrazilCurrent_RioGrande", "NorthBrazilCurrent_SEBrazil",  "EasternBrazil")
length(new_order)

 ## set colours for each sector ####
grid.col <-  c(CALDIO = "#401a7f", CALBOR = "#F7D203", CALEDW = "#e500e5", 
               Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
               Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
               Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
               North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", Central_South_Atlantic = "#004C73", 
               EasternBrazil = "#E9FFBE",  NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
               Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")
length(grid.col)

## plot to pdf device ####
# png("circular_plots/Circular_by_species.png", height = 25, width = 25, units = "cm", res = 300)
Cairo(file = "outputs/circular_plots/Circular_by_species_no_names.png",
      type = "png",
      units = "in",
      width = 20,
      height = 20,
      pointsize = 32,
      dpi = 96,
      bg = "transparent")
# set plot parameters (must be restarted at the end with circos.clear())
circos.par(start.degree = 170, clock.wise = TRUE,  # define starting point and direction
           gap.after = c(rep(2, length(unique(for_circle_sp[[1]])) - 1), 7, # set gap between origin and destination
                         rep(2, length(unique(for_circle_sp[[2]])) - 1), 7))

chordDiagram(for_circle_sp, order = new_order, grid.col = grid.col, 
             transparency = 0.25,
             direction.type = c("diffHeight", "arrows"),
             annotationTrack = c("grid"), annotationTrackHeight = c(0.05),
             link.sort = T, link.decreasing = T, 
             link.zindex = rank(for_circle_sp[[3]]),
             directional = 1, link.arr.type = "big.arrow", 
             preAllocateTracks = list(track.height = 0.5))

# add labels perpendicular to arch
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
}, bg.border = NA) # here set bg.border to NA is important

circos.clear() # reset plotting parameters

dev.off()



# 4. For CALBOR at colony level ####

## prepare dataset for plotting ####
for_plot <- for_circle_col %>% 
  filter(Species == "CALBOR") %>% 
  dplyr::select(Colony, Win_area, circle_arch_value)

union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted ####
new_order <- c("Faial", "Pico", "Vila", "Chafarinas","Berlenga", "Selvagem", 
               "MontañaClara", "Timanfaya", "Veneguera", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Gulf_of_Guinea_Central", "Angolan_Current", "Namib", "Namaqua",
               "Agulhas_Benguela_Confluence", "Agulhas_Current", 
               "Central_South_Atlantic", "North_Atlantic", 
               "PatagonianShelf_Falklands", "Uruguay_BuenosAires_Shelf", "SouthBrazilCurrent_RioGrande", "NorthBrazilCurrent_SEBrazil",  "EasternBrazil")
length(new_order)


cbcol <- colorRampPalette(c("#401A7F", "#FFF400"), bias = 1)

## set colours for each sector ####
grid.col = c(Faial = cbcol(9)[1], Pico = cbcol(9)[2], Vila = cbcol(9)[3], 
             Chafarinas = cbcol(9)[4], Berlenga = cbcol(9)[5], 
             Selvagem = cbcol(9)[6], MontañaClara = cbcol(9)[7], Timanfaya = cbcol(9)[8], Veneguera = cbcol(9)[9],
             Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
             Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
             Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
             North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", Central_South_Atlantic = "#004C73", 
             EasternBrazil = "#E9FFBE",  NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
             Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

## plot to pdf device ####
# png("circular_plots/CALBOR_colonies.png", height = 25, width = 25, units = "cm", res = 300)
Cairo(file = "outputs/circular_plots/CALBOR_colonies_nonames.png",
      type = "png",
      units = "in",
      width = 20,
      height = 20,
      pointsize = 32,
      dpi = 96,
      bg = "transparent")# set plot parameters (must be restarted at the end with circos.clear())
circos.par(start.degree = 170, clock.wise = TRUE,  # define starting point and direction
           gap.after = c(rep(3, length(unique(for_plot[[1]])) - 1), 7, # set gap between origin and destination
                         rep(3, length(unique(for_plot[[2]])) - 1), 7))

chordDiagram(for_plot, order = new_order, grid.col = grid.col, 
             link.sort = T, link.decreasing = T, link.zindex = rank(for_plot[[3]]),
             directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", 
             annotationTrack = c("grid"),annotationTrackHeight = c(0.05), 
             preAllocateTracks = list(track.height = 0.5))

# add labels perpendicular to arch
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
}, bg.border = NA) # here set bg.border to NA is important


# reset plotting parameters
circos.clear()

dev.off()

# 5. For CALBOR at a population level ####

## prepare dataset for plotting ####
for_plot <- for_circle_pop %>% 
  filter(Species == "CALBOR") %>% 
  dplyr::select(Population, Win_area, circle_arch_value)


union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted ####
new_order <- c("CentralAzores", "EastAzores", "Chafarinas", "Berlengas", 
               "Selvagens", "WestCanaryIslands", "EastCanaryIslands", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Gulf_of_Guinea_Central", "Angolan_Current", "Namib", "Namaqua",
               "Agulhas_Benguela_Confluence", "Agulhas_Current", 
               "Central_South_Atlantic", "North_Atlantic", 
               "PatagonianShelf_Falklands", "Uruguay_BuenosAires_Shelf", "SouthBrazilCurrent_RioGrande", "NorthBrazilCurrent_SEBrazil",  "EasternBrazil")
length(new_order)

## set colours for each sector ####
grid.col = c(CentralAzores = cbcol(7)[1], EastAzores = cbcol(7)[2], Berlengas = cbcol(7)[3],
             Chafarinas = cbcol(7)[4], Selvagens = cbcol(7)[5], WestCanaryIslands = cbcol(7)[6], 
             EastCanaryIslands = cbcol(7)[7],
             Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
             Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", 
             Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
             Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
             North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", 
             Central_South_Atlantic = "#004C73", 
             EasternBrazil = "#E9FFBE",  NorthBrazilCurrent_SEBrazil = "#D1FF73", 
             SouthBrazilCurrent_RioGrande = "#98E600", 
             Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

## plot to pdf device ####
# png("circular_plots/CALBOR_populations.png", height = 25, width = 25, units = "cm", res = 300)
Cairo(file = "outputs/circular_plots/CALBOR_populations_nonames.png",
      type = "png",
      units = "in",
      width = 20,
      height = 20,
      pointsize = 32,
      dpi = 96,
      bg = "transparent")# set plot parameters (must be restarted at the end with circos.clear())# set plot parameters (must be restarted at the end with circos.clear())
circos.par(start.degree = 170, clock.wise = TRUE,  # define starting point and direction
           gap.after = c(rep(3, length(unique(for_plot[[1]])) - 1), 7, # set gap between origin and destination
                         rep(3, length(unique(for_plot[[2]])) - 1), 7))

chordDiagram(for_plot, order = new_order, grid.col = grid.col, 
             link.sort = T, link.decreasing = T, link.zindex = rank(for_plot[[3]]),
             directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", 
             annotationTrack = c("grid"),annotationTrackHeight = c(0.05), 
             preAllocateTracks = list(track.height = 0.5))
# # add labels perpendicular to arch
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
}, bg.border = NA) # here set bg.border to NA is important


# reset plotting parameters
circos.clear()

dev.off()


# 6. For CALDIO at a colony level ####

## prepare dataset for plotting ####
for_plot <- for_circle_col %>% 
  filter(Species == "CALDIO") %>% 
  dplyr::select(Colony, Win_area, circle_arch_value)

union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted ####
new_order <- c("Palomas", "Pantaleu", "CalaMorell", "Riou", "Frioul", "Giraglia", "Lavezzi", "Zembra", 
               "Linosa", "Filfla", "Malta","Tremiti", "Strofades", "Paximada", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Equatorial_Atlantic", "Gulf_of_Guinea_West", "Gulf_of_Guinea_Central",  "Angolan_Current", "Namib", 
               "Agulhas_Current", "NorthBrazilCurrent_SEBrazil")
               

length(new_order)

cdcol <- colorRampPalette(c("#0F4017", "#B8E502"), bias = 1)


## set colours for each sector ####
grid.col = c(Palomas = cdcol(14)[1], Pantaleu = cdcol(14)[2], CalaMorell = cdcol(14)[3], 
             Riou = cdcol(14)[4], Frioul = cdcol(14)[5], Giraglia = cdcol(14)[6], Lavezzi = cdcol(14)[7], Zembra = cdcol(14)[8], 
             Linosa = cdcol(14)[9], Filfla = cdcol(14)[10], Malta = cdcol(14)[11], Tremiti = cdcol(14)[12], Strofades = cdcol(14)[13], Paximada = cdcol(14)[14], 
             Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
             Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
             Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
             North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", Central_South_Atlantic = "#004C73", 
             EasternBrazil = "#E9FFBE",  NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
             Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

## plot to pdf device ####

# png("circular_plots/CALDIO_colonies.png", height = 25, width = 25, units = "cm", res = 300)
Cairo(file = "outputs/circular_plots/CALDIO_colonies_nonames.png",
      type = "png",
      units = "in",
      width = 20,
      height = 20,
      pointsize = 32,
      dpi = 96,
      bg = "transparent")# set plot parameters (must be restarted at the end with circos.clear())# set plot parameters (must be restarted at the end with circos.clear())

# set plot parameters (must be restarted at the end with circos.clear())
circos.par(start.degree = 170, clock.wise = TRUE,  # define starting point and direction
           gap.after = c(rep(3, length(unique(for_plot[[1]])) - 1), 7, # set gap between origin and destination
                         rep(3, length(unique(for_plot[[2]])) - 1), 7))

chordDiagram(for_plot, reduce = -1, order = new_order, grid.col = grid.col, 
             link.sort = T, link.decreasing = T, link.zindex = rank(for_plot[[3]]),
             directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", 
             annotationTrack = c("grid"),annotationTrackHeight = c(0.05), 
             preAllocateTracks = list(track.height = 0.5))

# add labels perpendicular to arch
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
}, bg.border = NA) # here set bg.border to NA is important


# reset plotting parameters
circos.clear()

dev.off()

# 7. For CALDIO at a population level ####

## prepare dataset for plotting ####
for_plot <- for_circle_pop %>% 
  filter(Species == "CALDIO") %>% 
  dplyr::select(Population, Win_area, circle_arch_value)

union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted ####
new_order <- c("Palomas", "BalearicIslands", "France", "CorsicaSardinia", 
               "SicilianChannel", "MiddleAdriatic", "IonianSea", "AegeanSea", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Equatorial_Atlantic",
               "Gulf_of_Guinea_West", "Gulf_of_Guinea_Central", 
               "Angolan_Current", "Namib", "Agulhas_Current", 
               "NorthBrazilCurrent_SEBrazil")

length(new_order)

## set colours for each sector ####
grid.col = c(Palomas = cdcol(14)[1], BalearicIslands = cdcol(14)[2], 
             France = cdcol(14)[3], CorsicaSardinia = cdcol(14)[4],
             SicilianChannel = cdcol(14)[8], MiddleAdriatic = cdcol(14)[12], 
             IonianSea = cdcol(14)[13], AegeanSea = cdcol(14)[14], 
             Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
             Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", 
             Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
             Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
             North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", 
             Central_South_Atlantic = "#004C73", EasternBrazil = "#E9FFBE",  
             NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
             Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

## plot to pdf device ####

# png("circular_plots/CALDIO_populations.png", height = 25, width = 25, units = "cm", res = 300)
Cairo(file = "outputs/circular_plots/CALDIO_populations.png",
      type = "png",
      units = "in",
      width = 20,
      height = 20,
      pointsize = 32,
      dpi = 96,
      bg = "transparent")# set plot parameters (must be restarted at the end with circos.clear())# set plot parameters (must be restarted at the end with circos.clear())

# set plot parameters (must be restarted at the end with circos.clear())
circos.par(start.degree = 170, clock.wise = TRUE,  # define starting point and direction
           gap.after = c(rep(3, length(unique(for_plot[[1]])) - 1), 7, # set gap between origin and destination
                         rep(3, length(unique(for_plot[[2]])) - 1), 7))

chordDiagram(for_plot, reduce = -1, order = new_order, grid.col = grid.col, 
             link.sort = T, link.decreasing = T, link.zindex = rank(for_plot[[3]]),
             directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", 
             annotationTrack = c("grid"),annotationTrackHeight = c(0.05), 
             preAllocateTracks = list(track.height = 0.5))
# add labels perpendicular to arch
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
}, bg.border = NA) # here set bg.border to NA is important

circos.clear()
dev.off()


# 8. For CALEDW at a colony level ####

## prepare dataset for plotting ####
for_plot <- for_circle_col %>% 
  filter(Species == "CALEDW") %>% 
  dplyr::select(Colony, Win_area, circle_arch_value)

union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted ####
new_order <- c("Raso", "CurralVelho", 
               "Uruguay_BuenosAires_Shelf", "SouthBrazilCurrent_RioGrande", "NorthBrazilCurrent_SEBrazil")

## set colours for each sector ####
grid.col = c(CurralVelho = "#401A7F", Raso = "#ffd14f", 
             NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
             Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

## plot to pdf device ####

# png("circular_plots/CALEDW_colonies.png", height = 25, width = 25, units = "cm", res = 300)
Cairo(file = "outputs/circular_plots/CALEDW_colonies_nonames.png",
      type = "png",
      units = "in",
      width = 20,
      height = 20,
      pointsize = 32,
      dpi = 96,
      bg = "transparent")# set plot parameters (must be restarted at the end with circos.clear())# set plot parameters (must be restarted at the end with circos.clear())


# set plot parameters (must be restarted at the end with circos.clear())
circos.par(start.degree = 170, clock.wise = TRUE,  # define starting point and direction
           gap.after = c(rep(1, length(unique(for_plot[[1]])) - 1), 3, # set gap between origin and destination
                         rep(1, length(unique(for_plot[[2]])) - 1), 3))

chordDiagram(for_plot, reduce = -1, order = new_order, grid.col = grid.col, 
             link.sort = T, link.decreasing = T, link.zindex = rank(for_plot[[3]]),
             directional = 1, direction.type = c("diffHeight", "arrows"), link.arr.type = "big.arrow", 
             annotationTrack = c("grid"), annotationTrackHeight = c(0.05), 
             preAllocateTracks = list(track.height = 0.5))


# add labels perpendicular to arch
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
}, bg.border = NA) # here set bg.border to NA is important


# reset plotting parameters
circos.clear()

dev.off()

