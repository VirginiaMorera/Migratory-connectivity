library(sf)
library(rnaturalearth)
library(ggplot2)
library(Cairo)
library(sp)
library(rgdal)
library(tidyverse)
library(circlize)
library(raster)

setwd("C:/Users/morer/Dropbox/Virginia/GLS_multi")


# Figure 1 ####

win <- read_sf("data/shapefiles/wintering_areas_1.shp")
bre <- read.csv("data/colonies_info_100km_for_ArcGis.csv")
pops <- read_sf("data/pops2024.kml")
r <- units::set_units(90, "km")

pops2 <- bre %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>% 
  st_buffer(bre, dist = r) %>% 
  st_union()

aoof <- read_sf("data/almeria_oran_Layer.shp") 
aoof_line <- aoof %>% 
  summarise(do_union=F) %>% 
  st_cast("LINESTRING") %>% 
  mutate(name = "AOOF")

  ### need to change so chafarinas appears as "both"


# worldmap
world_simple <- rworldmap::getMap(resolution = "high")
# world_simple <- gBuffer(world_simple, byid = F, width = 0)
wrld <- st_as_sf(world_simple)

grid.col <- c(Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
              Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
              Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
              North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", Central_South_Atlantic = "#004C73", 
              EasternBrazil = "#E9FFBE",  NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
              Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

sp.col <- c(CALDIO = "#9072ba", CALBOR = "#F7D203", CALEDW = "#e500e5", BOTH = "#59EB3F")


bre <- bre %>% 
  dplyr::select(Species, Population = Radius100, SampledColony, Sample.Size, Pop.Min, Pop.Max, Longitude, Latitude) %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)


plot_b <- ggplot() + 
  geom_sf(data = win, aes(fill = ECOREGION)) +
  geom_sf(data = wrld, col = "darkgray", fill = "lightgray") +
  scale_fill_manual(name = "",
                    values = grid.col,
                    labels = c(
                      "Agulhas current", 
                      "Agulhas-Benguela confluence", 
                      "Namaqua", 
                      "Namib", 
                      "Angolan current", 
                      "Central gulf of Guinea", 
                      "West gulf of Guinea",
                      "Sahelian upwelling", 
                      "Saharan upwelling",
                      "North Atlantic", 
                      "Equatorial Atlantic", 
                      "Central south Atlantic", 
                      "Eastern Brazil", 
                      "North Brazil current", 
                      "South Brazil current", 
                      "Uruguay/B. Aires shelf", 
                      "Patagonian shelf"
                    ), 
                    breaks = names(grid.col)) +
  guides(fill = guide_legend(
      ncol = 3, position = "bottom")) + 
  coord_sf(xlim = st_bbox(win)[c(1,3)], ylim = st_bbox(win)[c(2,4)]) +
  # geom_sf_text(data = win, aes(label = ECOREGION)) +
  theme_bw() + 
  # theme(legend.position="bottom") +
  NULL


plot_a <- ggplot() + 
  geom_sf(data = bre, aes(col = Species, size = Pop.Min)) +
  geom_sf(data = aoof_line, lwd = 1) +
  geom_sf(data = wrld, col = "darkgray", fill = "lightgray") +
  geom_sf(data = pops, col = "black", linetype = "dashed", fill = NA) +
  geom_sf(data = bre %>% dplyr::filter(!is.na(Sample.Size)), size = 0.7) +
  scale_colour_manual(name = "",
                    values = sp.col,
                    labels = c(
                      "C. diomedea",
                      "C. borealis",
                      "C. edwardsii",
                      "C. diomedea + C. borealis"
                      ),
                    breaks = names(sp.col)) +
  coord_sf(xlim = st_bbox(bre)[c(1,3)], ylim = st_bbox(bre)[c(2,4)]) +
  # geom_sf_text(data = win, aes(label = ECOREGION)) +
  theme_bw() + 
  theme(legend.position="bottom", legend.box="vertical", legend.margin=margin()) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  labs(size = "Population size") +
  scale_size(range = c(1, 10)) +
  NULL
  

# Cairo(file = "Drafts/Ecology/Figures/Fig1a.pdf",
#       type = "pdf",
#       units = "in",
#       width = 30*0.5,
#       height = 20*0.5,
#       # pointsize = 32,
#       # dpi = 96,
#       bg = "white")  
# plot(plot_a)
# dev.off()
# 
# Cairo(file = "Drafts/Ecology/Figures/Fig1b.pdf",
#       type = "pdf",
#       units = "in",
#       width = 30*0.5,
#       height = 20*0.5,
#       # pointsize = 32,
#       # dpi = 96,
#       bg = "white")  
# plot(plot_b)
# dev.off()

Cairo(file = "Drafts/Ecology/Figures/Fig1a.png",
      type = "png",
      units = "in",
      dpi = 700,
      width = 30*0.25,
      height = 20*0.25,
      bg = "transparent")  
plot(plot_a)
dev.off()

Cairo(file = "Drafts/Ecology/Figures/Fig1b.png",
      type = "png",
      units = "in",
      dpi = 700, 
      width = 30*0.2,
      height = 20*0.25,
      bg = "white")  
plot(plot_b)
dev.off()

rm(list=ls())   

# Figure 2 ####


for_circle_sp <- readRDS("outputs/for_circle_sp.2023.RDS")

## prepare dataset for circular plot 

for_circle_sp <- for_circle_sp %>% 
  ungroup() %>% 
  # mutate(circle_arch_value = sqrt(circle_arch_value)) %>%
  dplyr::select(-n_pop_win, -ind_weight) 




# default order
union(for_circle_sp[[1]], for_circle_sp[[2]])
length(union(for_circle_sp[[1]], for_circle_sp[[2]]))

## set order in which sectors will be plotted 
new_order <- c("CALBOR", "CALEDW", "CALDIO", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Equatorial_Atlantic", "Gulf_of_Guinea_West", "Gulf_of_Guinea_Central",  "Angolan_Current", "Namib", "Namaqua",
               "Agulhas_Benguela_Confluence", "Agulhas_Current", 
               "Central_South_Atlantic", "North_Atlantic", 
               "PatagonianShelf_Falklands", "Uruguay_BuenosAires_Shelf", "SouthBrazilCurrent_RioGrande", "NorthBrazilCurrent_SEBrazil",  "EasternBrazil")
length(new_order)

## set colours for each sector 
grid.col <-  c(CALDIO = "#401a7f", CALBOR = "#F7D203", CALEDW = "#e500e5", 
               Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
               Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
               Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
               North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", Central_South_Atlantic = "#004C73", 
               EasternBrazil = "#E9FFBE",  NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
               Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")
length(grid.col)

## plot to pdf device 

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
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
# }, bg.border = NA) # here set bg.border to NA is important

circos.clear() # reset plotting parameters

dev.off()

rm(list=ls())   

# Figure 3 ####

for_circle_col <- readRDS("outputs/for_circle_col.2023.RDS")
for_circle_pop <- readRDS("outputs/for_circle_pop.2023.RDS")

# For CALBOR colony

## prepare dataset for plotting 
for_plot <- for_circle_col %>% 
  filter(Species == "CALBOR") %>% 
  dplyr::select(Colony, Win_area, circle_arch_value)

union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted 
new_order <- c("Faial", "Pico", "Vila", "Chafarinas","Berlenga", "Selvagem", 
               "MontañaClara", "Timanfaya", "Veneguera", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Gulf_of_Guinea_Central", "Angolan_Current", "Namib", "Namaqua",
               "Agulhas_Benguela_Confluence", "Agulhas_Current", 
               "Central_South_Atlantic", "North_Atlantic", 
               "PatagonianShelf_Falklands", "Uruguay_BuenosAires_Shelf", "SouthBrazilCurrent_RioGrande", "NorthBrazilCurrent_SEBrazil",  "EasternBrazil")
length(new_order)


cbcol <- colorRampPalette(c("#401A7F", "#FFF400"), bias = 1)

## set colours for each sector¡
grid.col = c(Faial = cbcol(9)[1], Pico = cbcol(9)[2], Vila = cbcol(9)[3], 
             Chafarinas = cbcol(9)[4], Berlenga = cbcol(9)[5], 
             Selvagem = cbcol(9)[6], MontañaClara = cbcol(9)[7], Timanfaya = cbcol(9)[8], Veneguera = cbcol(9)[9],
             Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
             Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
             Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
             North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", Central_South_Atlantic = "#004C73", 
             EasternBrazil = "#E9FFBE",  NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
             Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

## plot to pdf device 
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

# # add labels perpendicular to arch
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
# }, bg.border = NA) # here set bg.border to NA is important


# reset plotting parameters
circos.clear()

dev.off()

# CALBOR population

## prepare dataset for plotting
for_plot <- for_circle_pop %>% 
  filter(Species == "CALBOR") %>% 
  dplyr::select(Population, Win_area, circle_arch_value)


union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted
new_order <- c("CentralAzores", "EastAzores", "Chafarinas", "Berlengas", 
               "Selvagens", "WestCanaryIslands", "EastCanaryIslands", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Gulf_of_Guinea_Central", "Angolan_Current", "Namib", "Namaqua",
               "Agulhas_Benguela_Confluence", "Agulhas_Current", 
               "Central_South_Atlantic", "North_Atlantic", 
               "PatagonianShelf_Falklands", "Uruguay_BuenosAires_Shelf", "SouthBrazilCurrent_RioGrande", "NorthBrazilCurrent_SEBrazil",  "EasternBrazil")
length(new_order)

## set colours for each sector
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

## plot to pdf device
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
# # # add labels perpendicular to arch
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
# }, bg.border = NA) # here set bg.border to NA is important
# 

# reset plotting parameters
circos.clear()

dev.off()


# CALDIO colony

## prepare dataset for plotting
for_plot <- for_circle_col %>% 
  filter(Species == "CALDIO") %>% 
  dplyr::select(Colony, Win_area, circle_arch_value)

union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted
new_order <- c("Palomas", "Pantaleu", "CalaMorell", "Riou", "Frioul", "Giraglia", "Lavezzi", "Zembra", 
               "Linosa", "Filfla", "Malta","Tremiti", "Strofades", "Paximada", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Equatorial_Atlantic", "Gulf_of_Guinea_West", "Gulf_of_Guinea_Central",  "Angolan_Current", "Namib", 
               "Agulhas_Current", "NorthBrazilCurrent_SEBrazil")


length(new_order)

cdcol <- colorRampPalette(c("#0F4017", "#B8E502"), bias = 1)


## set colours for each sector
grid.col = c(Palomas = cdcol(14)[1], Pantaleu = cdcol(14)[2], CalaMorell = cdcol(14)[3], 
             Riou = cdcol(14)[4], Frioul = cdcol(14)[5], Giraglia = cdcol(14)[6], Lavezzi = cdcol(14)[7], Zembra = cdcol(14)[8], 
             Linosa = cdcol(14)[9], Filfla = cdcol(14)[10], Malta = cdcol(14)[11], Tremiti = cdcol(14)[12], Strofades = cdcol(14)[13], Paximada = cdcol(14)[14], 
             Agulhas_Current = "#FFBDBD", Agulhas_Benguela_Confluence = "#FF7F7F", 
             Namaqua = "#E60000", Namib = "#A80000", Angolan_Current = "#730000", Gulf_of_Guinea_Central = "#832200", Gulf_of_Guinea_West = "#a83800", 
             Sahelian_Upwelling = "#FF5500", Saharan_Upwelling = "#FF8121", 
             North_Atlantic = "#BEE8FF", Equatorial_Atlantic = "#A35A39", Central_South_Atlantic = "#004C73", 
             EasternBrazil = "#E9FFBE",  NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
             Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

## plot to pdf device

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
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
# }, bg.border = NA) # here set bg.border to NA is important
# 
# 
# reset plotting parameters
circos.clear()

dev.off()

# CALDIO population

## prepare dataset for plotting
for_plot <- for_circle_pop %>% 
  filter(Species == "CALDIO") %>% 
  dplyr::select(Population, Win_area, circle_arch_value)

union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted
new_order <- c("Palomas", "BalearicIslands", "France", "CorsicaSardinia", 
               "SicilianChannel", "MiddleAdriatic", "IonianSea", "AegeanSea", 
               "Saharan_Upwelling", "Sahelian_Upwelling", "Equatorial_Atlantic",
               "Gulf_of_Guinea_West", "Gulf_of_Guinea_Central", 
               "Angolan_Current", "Namib", "Agulhas_Current", 
               "NorthBrazilCurrent_SEBrazil")

length(new_order)

## set colours for each sector
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

## plot to pdf device

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
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
# }, bg.border = NA) # here set bg.border to NA is important

circos.clear()
dev.off()


# 8 CALEDW colony 

## prepare dataset for plotting
for_plot <- for_circle_col %>% 
  filter(Species == "CALEDW") %>% 
  dplyr::select(Colony, Win_area, circle_arch_value)

union(for_plot[[1]], for_plot[[2]])

length(union(for_plot[[1]], for_plot[[2]]))

## set order in which sectors will be plotted
new_order <- c("Raso", "CurralVelho", 
               "Uruguay_BuenosAires_Shelf", "SouthBrazilCurrent_RioGrande", "NorthBrazilCurrent_SEBrazil")

## set colours for each sector
grid.col = c(CurralVelho = "#401A7F", Raso = "#ffd14f", 
             NorthBrazilCurrent_SEBrazil = "#D1FF73", SouthBrazilCurrent_RioGrande = "#98E600", 
             Uruguay_BuenosAires_Shelf = "#70A800", PatagonianShelf_Falklands = "#4C7300")

## plot to pdf device

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
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
#               facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.05))
# }, bg.border = NA) # here set bg.border to NA is important
# 

# reset plotting parameters
circos.clear()

dev.off()

rm(list=ls())   


#Figure 4 ####

wrld <- st_union(ne_countries(returnclass = "sf"))


CB_pred <- as(terra::rast("embarcadero_SDMs/final_outputs/CALBOR_pred_stack_2024.grd"), "Raster")
CD_pred <- as(terra::rast("embarcadero_SDMs/final_outputs/CALDIO_pred_stack_2024.grd"), "Raster")
CE_pred <- as(terra::rast("embarcadero_SDMs/final_outputs/CALEDW_pred_stack_2024.grd"), "Raster")

average_preds <- raster::stack(CB_pred$Mean, CD_pred$Mean, CE_pred$Mean)
names(average_preds) <- c("CALBOR", "CALDIO", "CALEDW")

plot(average_preds)

avg_preds_df <- as.data.frame(average_preds, xy = TRUE) %>% 
  pivot_longer(cols = c("CALBOR", "CALDIO", "CALEDW"), 
               names_to = "Species", values_to = "value") %>% 
  mutate(variable = "Mean suitability")

CI_width <- raster::stack(CB_pred$q97.5.-CB_pred$q2.5.,
                  CD_pred$q97.5.-CD_pred$q2.5.,
                  CE_pred$q97.5.-CE_pred$q2.5.)

names(CI_width) <- c("CALBOR", "CALDIO", "CALEDW")

CI_width_df <- as.data.frame(CI_width, xy = TRUE) %>% 
  pivot_longer(cols = c("CALBOR", "CALDIO", "CALEDW"), 
               names_to = "Species", values_to = "value") %>% 
  mutate(variable = "CI width")

all_preds <- bind_rows(avg_preds_df, CI_width_df) %>% 
  mutate(variable = factor(variable, levels = c("Mean suitability", "CI width")))

f4 <- ggplot() +
  geom_raster(data = all_preds , aes(x = x, y = y, fill = value)) +
  geom_sf(data = wrld, col = "black", fill = "gray50") +
  scale_fill_viridis_c(na.value = "gray50", limits = c(0, 1)) + 
  facet_grid(variable ~ Species) + 
  coord_sf(xlim = c(min(all_preds$x+10), max(all_preds$x-10)),
           ylim = c(min(all_preds$y+15), max(all_preds$y-10))) +
  theme_bw() + 
  labs(x = "Longitude (º)", y = "Latitude (º)", fill = "") + 
  theme(strip.text = element_text(size = 12))



Cairo::CairoPDF(file = "Drafts/Ecology/Figures/Fig4.pdf", width = 8*1.5, height = 6*1.5)
print(f4)
dev.off()

# Figure 5 ####

calbor_partials <- read.csv("embarcadero_SDMs/final_outputs/calbor_partial_df.csv")
caldio_partials <- read.csv("embarcadero_SDMs/final_outputs/caldio_partial_df.csv")
caledw_partials <- read.csv("embarcadero_SDMs/final_outputs/caledw_partial_df.csv")

partials <- bind_rows(calbor_partials, caldio_partials, caledw_partials) %>% 
  mutate(x2 = if_else(Variable == "Chlorophyll A", x2_log, x2),
         Variable = recode(Variable,  "Chlorophyll A" = "Chlorophyll A (mg·m-3) log scale", 
                           "Chlorophyll A variability" = "Chlorophyll A variability (%)",
                           "Sea surface temperature" = "Sea Surface Temperature (ºC)"))

f5 <- ggplot(partials, aes(x = x2)) + 
  geom_ribbon(aes(ymin = q05, ymax = q95), fill = "deepskyblue", alpha = 0.3) +
  geom_line(aes(y = med), linewidth = 0.5) + 
  theme_bw() + 
  labs(x = "Environmental variable", y = "Partial effect") + 
  ggh4x::facet_grid2(Species ~ Variable, scales = "free_x", independent = "x")

Cairo::CairoPDF(file = "Drafts/Ecology/Figures/Fig5.pdf", width = 8*1.4, height = 6*1.4)
print(f5)
dev.off()

# Figure S1 ####

migcons <- read.csv(file = "outputs/migcon_values.csv")

migcons2 <- migcons %>% 
  mutate(Sp = factor(Sp, levels = c("CALBOR", "CALDIO", "CALEDW", "CALONECTRIS"))) 

fs1 <- ggplot(migcons) +
  geom_point(aes(x = level, y = median, col = Sp), 
             position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(x = level, ymin = lowerCI, ymax = higherCI, col = Sp), 
                position = position_dodge(width = 0.5), 
                width = 0.25) + 
  geom_hline(aes(yintercept = 0), col = "black", lty = 2) + 
  labs(x = "Aggregation level", y = "Migratory connectivity", col = "Species") +
  theme_bw()

Cairo::CairoPDF(file = "Drafts/Ecology/Figures/FigS1.pdf", width = 8*1.4, height = 6*1)
print(fs1)
dev.off()
