# auxiliary functions for the mig connectivity script


# little function to make a "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# functions to calculate geographic distance matrix 
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="m")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}




# function to calculate centroid of win positions for each individual (passed to cluster to speed up)

centroid_calculate <- function(dt, Lon, Lat, ind_ID, group_ID) {
  #kernel.5 non-breeding:
  # dt <- all_win[all_win$unique_ID_Bird == unique(all_win$unique_ID_Bird)[[1]],]
  coor.xx.sp <- SpatialPoints(coords = dt[,c(Lon, Lat)])
  kxx <- adehabitatHR::kernelUD(coor.xx.sp, h = "href", grid = 1000, extent = 1)
  v5 <- adehabitatHR::getverticeshr(kxx, 5)
  centk5 <- rgeos::gCentroid(v5)
  d <- data.frame(unique_ID_Bird = unique(as.character(dt[, ind_ID])),
                  Colony         = unique(as.character(unique(dt[,group_ID]))),
                  Lon.win        = unique(centk5@coords[,1]),
                  Lat.win        = unique(centk5@coords[,2]))
  
  return(d)
}


#### geoBias ####
# for GLS data, vector of length 2 indicating expected bias in lonlat of targetPoints in the units of the projection

#### geoVCov ####
# # for GLS data 2x2 matrix with expected var/covar in lonlat of target points in the units of the projection
# 
# calibLocations <- read.csv("data/locs.csv", row.names = NULL) # this should be lon and lat of calibration locations
# calibLocations <- calibLocations %>% 
#   dplyr::filter(locality == "Veneguera") %>%  # because we've only used geos from Veneguera in this case. Ideally, a sample from more calib locations
#   dplyr::select(c("locality", "LON", "LAT"))
# calibLocations$locality <- as.factor(as.character(calibLocations$locality))
# 
# # Convert capture locations into SpatialPoints #
# CapLocs <- calibLocations
# coordinates(CapLocs) <- ~LON+LAT
# CapLocs@proj4string <- CRS(projections$WGS84)
# 
# # Project Capture locations # 
# 
# CapLocsM<-sp::spTransform(CapLocs, sp::CRS(projections$EquidistConic))
# 
# # # Retrieve raw breeding locations
# 
# all_br <- all_data %>% 
#   dplyr::filter(all_data$fen == "BRE" & all_data$Sp == "CALDIO")
# 
# # reduce dataset to months of late inc and early chick rearing, which are more restricted to around the colony. 
# # we are going to use this to calculate error of the geo by comparing to colony location, so this is important
# all_br$month <- month(as.POSIXlt(all_br$Date))
# table(all_br$month)
# 
# all_br <- all_br %>% 
#   dplyr::filter(month > 5 & month < 9)

# use calibration data to calculate error from known location (colony)

# calib_data <- read.csv("data/calib_veneguera_locations.csv") # this are lat, lon, id and calib site of the estimated positions during calibration
# names(calib_data)
# 
# Breeding_files_df <- data.frame(Bird_ID         = calib_data$geo_deployment,
#                                 CaptureLocation = "Veneguera", # in this instance. ideally, calib_data$Location
#                                 Date            = as.Date(calib_data$date, format = "%d/%m/%Y"),
#                                 Longitude       = calib_data$long,
#                                 Latitude        = calib_data$lat)
# 
# Breeding_files_df$Bird_ID <- as.factor(as.character(Breeding_files_df$Bird_ID))
# Breeding_files_df$CaptureLocation <- as.factor(as.character(Breeding_files_df$CaptureLocation))
# str(Breeding_files_df)
# 
# unique(Breeding_files_df$CaptureLocation) %in% unique(calibLocations$locality)
# 
# 
# Breeding_files <- split(Breeding_files_df, Breeding_files_df$Bird_ID)
# 
# 
# for(i in 1:length(Breeding_files)){
#   x <- Breeding_files[[i]]
#   x$Bird_ID <- as.factor(as.character(x$Bird_ID))
#   x$CaptureLocation <- as.factor(as.character(x$CaptureLocation))
#   x$Date <- as.Date(x$Date)
#   Breeding_files[i] <- list(x)}
# 
# # Turn the locations into spatialpoints #
# 
# B_GL <- lapply(Breeding_files, 
#                FUN = function(x){
#                  sp::SpatialPoints(cbind(x$Longitude,x$Latitude),
#                                    sp::CRS(projections$WGS84))})
# # Project into UTM projection #
# 
# B_GLmeters <- lapply(B_GL,
#                      FUN = function(x){sp::spTransform(x,
#                                                        sp::CRS(projections$EquidistConic))})
# 
# # Process to determine geolocator bias and variance-covariance in meters #
# 
# # generate empty vectors to store data #
# LongError<-rep(NA,length(B_GLmeters))
# LatError<-rep(NA,length(B_GLmeters))
# 
# # Calculate the error derived from geolocators from the true capture location 
# 
# Bird_Capture <- Breeding_files_df %>% 
#   dplyr::select(c("Bird_ID", "CaptureLocation"))
# 
# 
# Bird_Capture <- Bird_Capture[!duplicated(Bird_Capture),]
# 
# # longitude 
# 
# for(i in 1:nrow(Bird_Capture)){
#   # i = 1
#   x <- B_GLmeters[[as.character(Bird_Capture[i,1])]]
#   y <- CapLocsM[CapLocsM$locality == Bird_Capture[i,2],]
#   z <- mean(x@coords[,1]-y@coords[,1])
#   LongError[i] <- z
# }
# 
# summary(LongError)
# 
# # latitude
# 
# for(i in 1:nrow(Bird_Capture)){
#   # i = 2
#   x <- B_GLmeters[[as.character(Bird_Capture[i,1])]]
#   y <- CapLocsM[CapLocsM$locality == Bird_Capture[i,2],]
#   z <- mean(x@coords[,2]-y@coords[,2])
#   LatError[i] <- z
# }
# 
# summary(LatError)
# 
# plot(density(LatError))
# plot(density(LongError))
# 
# 
# # Get co-variance matrix for error of 
# # known non-breeding deployment sites 
# 
# # lm does multivariate normal models if you give it a matrix dependent variable!
# 
# geo.error.model <- lm(cbind(LongError,LatError) ~ 1) 



##### assign wintering area ####
# function to assign a wintering area to each individual 
# Level: 2 = ecoregion, 4 = Province, 6 = Realm, 9 = Lat_Zone
# test <- list()
assign_win_area <- function(track, polys, level = 2){
# for (i in seq_along(unique(all_complete$unique_ID_Bird))){
#   # i = 1
#   track <- all_complete[all_complete$unique_ID_Bird == unique(all_complete$unique_ID_Bird)[i], ]
#   polys <-  win_areas
#   level <- 2
  win_area <- track[track$fen == "WIN",]
  win_area <- data.frame(as.matrix(win_area))
  win_area$Lat <- as.numeric(as.character(win_area$Lat))
  win_area$Lon <- as.numeric(as.character(win_area$Lon))
  if (nrow(win_area) == 0) {main_area <- "no_wint"; prop <- "NA"} else{
    coordinates(win_area) <- ~Lon+Lat
    win_area@proj4string <- CRS(projections$WGS84)
    win_area <- spTransform(win_area, CRS(projections$EquidistConic))
    polys <- spTransform(polys, win_area@proj4string)
    ovr <- over(win_area, polys)
    ovr[,level] <- as.factor(as.character(ovr[,level]))
    areas <- as.data.frame(dplyr::count(ovr, ovr[,level]))
    if (is.na(areas$`ovr[, level]`[which.max(areas[,2]/sum(areas[,2]))])) {main_area <- "outside"} else{
      main_area <- as.character(areas$`ovr[, level]`[which.max(areas[,2]/sum(areas[,2]))])}
    prop <- areas[which.max(areas[,2]),][,2]/sum(areas[,2])
  }
  out <- data.frame(unique_ID_Bird = unique(track$unique_ID_Bird),
                    Bird = unique(track$Bird_ID),
                    Colony = unique(track$Colony),
                    Year = unique(track$year_rec),
                    Win_area = main_area,
                    Proportion = as.numeric(prop), 
                    Species = unique(track$Sp))
  return(out)
  # print(i); print(out)
  # test[i] <- list(out)
}

# for (i in seq_along(unique(centroids.ALL$unique_ID_Bird))) {
#   # i = 469
#   cent <- centroids.ALL[centroids.ALL$unique_ID_Bird == unique(centroids.ALL$unique_ID_Bird)[i], ]
#   cent <- centroids.ALL[centroids.ALL$unique_ID_Bird == "L38481_first", ]
#   print(cent$unique_ID_Bird)
#   polys <-  win_areas
#   level <- 2
assign_win_centroid <- function(cent, polys, level) {
  coordinates(cent) <-  ~ Lon.win + Lat.win
  cent@proj4string <- polys@proj4string
  cent$Win_area <- over(cent, polys)[,level]
  plot(polys, main = cent$unique_ID_Bird)
  plot(cent, add = T)
  cent <- as.data.frame(cent)
  return(cent)
}


calc_perc <- function(x){
  tot <- sum(x$circle_arch_value)
  x$perc <- round((x$circle_arch_value/tot)*100, 2)
  return(x)
}

sfc_as_cols <- function(x, geometry, names = c("x","y")) {
  if (missing(geometry)) {
    geometry <- sf::st_geometry(x)
  } else {
    geometry <- rlang::eval_tidy(enquo(geometry), x)
  }
  stopifnot(inherits(x,"sf") && inherits(geometry,"sfc_POINT"))
  ret <- sf::st_coordinates(geometry)
  ret <- tibble::as_tibble(ret)
  stopifnot(length(names) == ncol(ret))
  x <- x[ , !names(x) %in% names]
  ret <- setNames(ret,names)
  dplyr::bind_cols(x,ret)
}

sourceDir <- function(path, trace = TRUE, ...) {
  op <- options(); on.exit(options(op)) # to reset after each 
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}

# Neal Fultz's very handy function at 
# https://stackoverflow.com/questions/10287545/backtransform-scale-for-plotting

unscale <- function(z, 
                    center = attr(z, "scaled:center"), 
                    scale = attr(z, "scaled:scale")) 
  {
  if(!is.null(scale))  z <- sweep(z, 2, scale, `*`)
  if(!is.null(center)) z <- sweep(z, 2, center, `+`)
  structure(z,
            "scaled:center"   = NULL,
            "scaled:scale"    = NULL,
            "unscaled:center" = center,
            "unscaled:scale"  = scale
  )
}
