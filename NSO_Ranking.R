########################################################################################
########################################################################################
##            Ranking Unburned Islands for Northern Spotted Owl Management            ##
##                               Anthony J.  Martinez                                 ##
##                                    Spring 2019                                     ##
########################################################################################
########################################################################################

#################################################
##    Load packages and data                   ##
#################################################
library(raster)
library(rgdal)
library(rgeos)
library(reshape2)
library(dplyr)
library(broom)
library(ggplot2)
library(ggridges)
library(cowplot)
library(ggmap)

setwd("/Users/robandrus/Google Drive/WSU postdoc/Projects/SpottedOwl/Data")

# Compile list of available RdNBR rasters 
if(!all(file.exists("Inputs/perims.shp", "Inputs/unbs.shp"))){
  fires <- unique(toupper(substr(list.files(path = "Inputs/MTBS", pattern = "WA|OR"), 1, 21)))
}

# Load fire perimeters for all available RdNBR rasters
if(file.exists("Inputs/perims.shp")){
  perims <- readOGR("Inputs/perims.shp")
}else{
  perims <- readOGR("Inputs/MTBS/mtbs_perims_1984_2014.shp")
  perims <- perims[perims@data$Fire_ID %in% fires,]
  perims <- spTransform(perims, "+proj=aea +lat_1=43 +lat_2=48 +lat_0=34 +lon_0=-120 +x_0=600000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  writeOGR(perims, driver = "ESRI Shapefile", layer = "perims.shp", dsn = "Inputs/perims.shp")
}

# Load unburned islands for all available RdNBR rasters
if(file.exists("Inputs/unbs.shp")){
  unbs <- readOGR("Inputs/unbs.shp")
}else{
  unbs <- readOGR("Inputs/MTBS/unburned_areas.shp")
  unbs <- unbs[unbs@data$fire_id %in% fires,]
  unbs <- spTransform(unbs, "+proj=aea +lat_1=43 +lat_2=48 +lat_0=34 +lon_0=-120 +x_0=600000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  writeOGR(unbs, driver = "ESRI Shapefile", layer = "unbs.shp", dsn = "Inputs/unbs.shp")
}

# Add fire names to Unburned Islands
if(!"FireName" %in% colnames(unbs@data)){
  a <- perims@data[, c("Fire_Name", "Fire_ID")]
  b <- unbs@data$fire_id
  FireName <- a$Fire_Name[match(b, a$Fire_ID)]
  unbs@data$FireName <- FireName
  unbs@data <- unbs@data[, c("ID", "fire_id", "FireName", "Area", "Perim")]
  writeOGR(unbs, driver = "ESRI Shapefile", layer = "unbs.shp", dsn = "Inputs/unbs.shp")
  rm(a, b)
}


# Perform EEMS framework
if(file.exists("Inputs/unbs_NSO.shp")){
  unbs <- readOGR("Inputs/unbs_NSO.shp")
}else{
  # Import RDNBR mask
  if(file.exists("Intermediates/RdNBRMask.tif")){
    mask <- raster("Intermediates/RdNBRMask.tif")
  }else{
    # Import all RdNBR rasters
    rdnbr.files <- list.files(path = "Inputs/MTBS", pattern = "rdnbr.tif", ignore.case = T, full.names = T)
    rdnbr.rasters <- lapply(rdnbr.files, raster)
    
    # Combine rasters
    names(rdnbr.rasters) <- NULL
    rdnbr.rasters$fun <- max
    rdnbr.rasters$na.rm <- TRUE
    mask <- do.call(raster::mosaic, rdnbr.rasters)
    
    # Apply threshold to create mask
    mask[mask < 703] <- NA
    mask[mask >= 703]  <- 1
    
    # Project and save raster
    mask <- projectRaster(mask, suit)
    writeRaster(mask, filename = "Intermediates/RdNBRMask.tif", overwrite = T)
    rm(rdnbr.files, rdnbr.rasters)
  }
  
  # Mask rasters
  if(file.exists("Intermediates/suit_mask.tif")){
    suit <- raster("Intermediates/suit_mask.tif")
  }else{
    suit <- raster("Inputs/suit2.tif")
    suit <- mask(suit, mask, inverse = T, maskvalue = 0, filename = "Intermediates/suit_mask.tif")
  }
  
  if(file.exists("Intermediates/high_suit_mask.tif")){
    high.suit <- raster("Intermediates/high_suit_mask.tif")
  }else{
    high.suit <- raster("Inputs/high_suit2.tif")
    high.suit <- mask(high.suit, mask, inverse = T, maskvalue = 0, filename = "Intermediates/high_suit_mask.tif")
  }
  
  if(file.exists("Intermediates/core_mask.tif")){
    core <- raster("Intermediates/core_mask.tif")
  }else{
    #I actually made an updated core raster in ArcGIS.  It took too long in R
    #core <- raster("Inputs/core1.tif")
    #core <- mask(core, mask, inverse = T, maskvalue = 0, filename = "Intermediates/core_mask.tif")
  }
  
  if(file.exists("Intermediates/coree_mask.tif")){ # Coree is core + core-edge
    coree <- raster("Intermediates/coree_mask.tif")
  }else{
    coree <- raster("Inputs/coree1.tif")
    coree <- mask(coree, mask, inverse = T, maskvalue = 0, filename = "Intermediates/coree_mask.tif")
  }
  
  # Make neighborhood rasters with 2.5 km radius
  # I actually did this in ArcGIS (Focal Statistics) because it took too long in R
  if(!all(file.exists("Intermediates/suit_2500.tif", "Intermediates/high_suit_2500.tif", "Intermediates/coree_2500.tif"))){
    w1 <- focalWeight(suit, 2500, type = 'circle')
    w1[w1 > 0] <- 1
  }
  
  if(file.exists("Intermediates/suit_2500.tif")){
    suit.2500 <- raster("Intermediates/suit_2500.tif")
  }else{
    suit.2500 <- focal(suit, w1, fun = mean, na.rm = T, filename = "Intermediates/suit_2500.tif")
  }
  
  if(file.exists("Intermediates/high_suit_2500.tif")){
    high.suit.2500 <- raster("Intermediates/high_suit_2500.tif")
  }else{
    high.suit.2500 <- focal(high.suit, w1, fun = mean, na.rm = T, filename = "Intermediates/high_suit_2500.tif")
  }
  
  if(file.exists("Intermediates/coree_2500.tif")){
    core.2500 <- raster("Intermediates/coree_2500.tif")
  }else{
    core.2500 <- focal(coree, w1, fun = mean, na.rm = T, filename = "Intermediates/coree_2500.tif")
  }
  
  #############################
  # Start here if you want to #
  #  change model parameters  #
  #############################
  # Extract values from rasters to unburned islands
  # This takes about 21 hours on Lee's server
  if(file.exists("Intermediates/unbs_values.csv")){
    df <- read.csv("Intermediates/unbs_values.csv")
  }else{
    criteria.stack <- stack(suit, high.suit, core, suit.2500, high.suit.2500, core.2500)
    df <- extract(criteria.stack, unbs, fun = mean, na.rm = T, df = T) * 100
  }
  
  # Prepare data.frame
  df <- df[,-c(1:2)]
  colnames(df) <- c("Suit", "H_Suit", "Core", "Suit25", "H_Suit25", "Core25")
  df[is.na(df)] <- 0
  df$Core[df$Core > 0] <- 1 # Convert from proportion to true/false
  frac <- (2 * log(.25 * unbs@data$Perim)) / log(unbs@data$Area * 10000)
  df <- cbind(unbs@data, frac, df) # Combine unburned islands and extracted values
  
  # Perform EEMS evaluation
  # Set parameters
  thresh <- c(418, 100, 100, 1, 
              75, 75, 50,
              0.18, 0, 0, 0, 
              0, 0, 0)
  weights <- c(1,1,1,1,1,1,1,1,1)
  
  # Convert to fuzzy values
  Cvt2Fz <- function(True, False, variable){
    x <- df[,variable]
    t <- thresh[True]
    f <- thresh[False]
    m <- 2/(t-f)
    b <- 1-m*t
    ret <- (m*x)+b
    ret <- ifelse(ret > 1, 1, ret)
    ret <- ifelse(ret < -1, -1, ret)
    return(ret)
  }
  
  df <- cbind(ID.1 = seq(1, nrow(df)), df)
  df.fz <- data.frame(df,
                      F_Area = Cvt2Fz(1, 8, "Area"),
                      F_Suit = Cvt2Fz(2, 9, "Suit"),
                      F_H_Suit = Cvt2Fz(3, 10, "H_Suit"),
                      F_Core = Cvt2Fz(4, 11, "Core"),
                      F_Suit25 = Cvt2Fz(5, 12, "Suit25"),
                      F_H_Suit25 = Cvt2Fz(6, 13, "H_Suit25"),
                      F_Core25 = Cvt2Fz(7, 14, "Core25"))
  df.fz <- replace(df.fz, is.na(df.fz), -1)
  df.fz$F_Value_WI <- (weights[3]*df.fz$F_Area + 
                         weights[4]*df.fz$F_Suit + 
                         weights[5]*df.fz$F_H_Suit + 
                         weights[6]*df.fz$F_Core) / 4 # Union (mean)
  df.fz$F_Value_Out <- (weights[7]*df.fz$F_Suit25 + 
                          weights[8]*df.fz$F_H_Suit25 + 
                          weights[9]*df.fz$F_Core25) / 3 # Union (mean)
  df.fz$F_Value_Ref <- (weights[1]*df.fz$F_Value_WI + 
                          weights[2]*df.fz$F_Value_Out) / 2 # Union (mean)

  ranking <- function(fire){
    ret <- subset(df.fz, FireName == unique(df.fz$FireName)[fire])
    ret$F_Rank <- rank(-ret$F_Value_Ref)
    return(ret)
  }
  
  out <- df.fz[0,]
  for(i in 1:length(unique(df.fz$FireName))){
    out <- rbind(out, ranking(i))
  }
  out <- out[order(out$ID.1),]
  
  unbs@data <- out
  writeOGR(unbs, driver = "ESRI Shapefile", layer = "Inputs/unbs_NSO.shp", dsn = "Inputs/unbs_NSO.shp")
  rm(Cvt2Fz, ranking, i, thresh, weights, df, df.fz, out)
}

#################################################
##    Add lidar data                           ##
#################################################
load("Inputs/LidarStacks.Rdata")
colnames(unbs@data) <- c("ID.1", "ID", "FireID", "FireName", "Area", "Perim", "frac", "Suit", "H_Suit", "Core", "Suit25", "H_Suit25", 
                        "Core25", "F_Area", "F_Suit", "F_H_Suit", "F_Core", "F_Suit25", "F_H_Suit25", 
                        "F_Core25", "F_Value_WI", "F_Value_Out", "F_Value_Ref", "F_Rank")

ExtractLidar <- function(stack, fire){
  unb <- unbs[unbs@data$FireName == fire,]
  unb@data$ID <- as.numeric(as.character(unb@data$ID))
  unb <- spTransform(unb, proj4string(stack))
  unb.ras <- rasterize(unb, stack, field = "ID")
  unb.stack <- mask(stack, unb.ras)
  unb.values <- addLayer(unb.ras, unb.stack)
  names(unb.values) <- c("ID", names(unb.values)[-1])
  valuetable <- as.data.frame(getValues(unb.values))
  valuetable <- valuetable[!is.na(valuetable$ID),]
  valuetable <- aggregate(valuetable, by = list(valuetable$ID), FUN = mean)[,-1]
  ret <- merge.data.frame(unb@data, valuetable, by.x = "ID", by.y = "ID")
  return(ret)
}

df.BB <- ExtractLidar(stack.bb, "B&B COMPLEX (BOOTH)")
df.davis <- ExtractLidar(stack.davis, "DAVIS")
df.poison <- ExtractLidar(stack.poison, "POISON")
df.pole <- ExtractLidar(stack.pole, "POLE CREEK")
df.analysis <- rbind(df.BB, df.davis, df.poison, df.pole)
df.analysis$TRASP <-  (1 - cos((pi/180) * (df.analysis$topo_aspect_147p636F_98p424FEET-30))) / 2

#################################################
##    Analyze data                             ##
#################################################
# Prepare data
percentile <- 0.1
df.all <- df.analysis
df.all$var <- "Bottom 90%"

top <- function(data, percentile, fire){
  df <- subset(data, FireName == unique(data$FireName)[fire])
  cutoff <- quantile(df$F_Value_Ref, probs = 1 - percentile)
  df[df$F_Value_Ref >= cutoff, "var"] <- paste0("Top ", percentile*100, "%")
  return(df)
}

df.all <- rbind(top(df.all, percentile = percentile, 1),
                top(df.all, percentile = percentile, 2),
                top(df.all, percentile = percentile, 3),
                top(df.all, percentile = percentile, 4))

# Convert feet to meters
df.all$topo_elev_m <- df.all$topo_elevation_147p636F_98p424FEET #/ 3.281 #already in meters
df.all$elev_ave_m <- df.all$FIRST_RETURNS_elev_ave_6p5616plus_98p424FEET / 3.281
df.all$elev_P95_m <- df.all$FIRST_RETURNS_elev_P95_6p5616plus_98p424FEET / 3.281

# Shape data
df.all.melt <- melt(df.all, c("FireName", "var"))
df.all.means <- summarize_all(group_by(df.all, FireName, var), mean, na.rm = T)
df.all.mean <- summarize_all(group_by(df.all, FireName), mean, na.rm = T)
cols <- as.logical((sapply(df.all, class) == "numeric") + (colnames(df.all) == "FireName"))
df.all.top <- summarize_all(group_by(df.all[, cols], FireName), quantile, probs = (1 - percentile), na.rm = T)
x <- unique(as.character(df.all$FireName))

ks.area.frac  <- c(sapply(x, FUN = function(x){
  ks.test(as.numeric(subset(df.all, var == "Bottom 90%" | FireName == x, select = Area)[,1]), 
          as.numeric(subset(df.all, var == "Top 10%" | FireName == x, select = Area)[,1]))$p.value}),
  sapply(x, FUN = function(x){
    ks.test(as.numeric(subset(df.all, var == "Bottom 90%" | FireName == x, select = frac)[,1]), 
            as.numeric(subset(df.all, var == "Top 10%" | FireName == x, select = frac)[,1]))$p.value}))
ks.area.frac <- ifelse(ks.area.frac > 0.001, paste0("P = ", round(ks.area.frac, 3)), "P < 0.001")

ks.topo  <- c(sapply(x, FUN = function(x){
  ks.test(as.numeric(subset(df.all, var == "Bottom 90%" | FireName == x, select = topo_elev_m)[,1]), 
          as.numeric(subset(df.all, var == "Top 10%" | FireName == x, select = topo_elev_m)[,1]))$p.value}),
  sapply(x, FUN = function(x){
    ks.test(as.numeric(subset(df.all, var == "Bottom 90%" | FireName == x, select = TRASP)[,1]), 
            as.numeric(subset(df.all, var == "Top 10%" | FireName == x, select = TRASP)[,1]))$p.value}),
  sapply(x, FUN = function(x){
    ks.test(as.numeric(subset(df.all, var == "Bottom 90%" | FireName == x, select = topo_slope_147p636F_98p424FEET)[,1]), 
            as.numeric(subset(df.all, var == "Top 10%" | FireName == x, select = topo_slope_147p636F_98p424FEET)[,1]))$p.value}))
ks.topo <- ifelse(ks.topo > 0.001, paste0("P = ", round(ks.topo, 3)), "P < 0.001")

ks.canopy  <- c(sapply(x, FUN = function(x){
  ks.test(as.numeric(subset(df.all, var == "Bottom 90%" | FireName == x, select = AGB)[,1]), 
          as.numeric(subset(df.all, var == "Top 10%" | FireName == x, select = AGB)[,1]))$p.value}),
  sapply(x, FUN = function(x){
    ks.test(as.numeric(subset(df.all, var == "Bottom 90%" | FireName == x, select = elev_ave_m)[,1]), 
            as.numeric(subset(df.all, var == "Top 10%" | FireName == x, select = elev_ave_m)[,1]))$p.value}),
  sapply(x, FUN = function(x){
    ks.test(as.numeric(subset(df.all, var == "Bottom 90%" | FireName == x, select = elev_P95_m)[,1]), 
            as.numeric(subset(df.all, var == "Top 10%" | FireName == x, select = elev_P95_m)[,1]))$p.value}))
ks.canopy <- ifelse(ks.canopy > 0.001, paste0("P = ", round(ks.canopy, 3)), "P < 0.001")
#################################################
##    Plot data                                ##
#################################################

##########################
# Figure 4-7: Map/density#
##########################
perim.proj <- spTransform(perims, CRS("+proj=longlat +datum=WGS84"))
unb.proj <- spTransform(unbs, CRS("+proj=longlat +datum=WGS84"))

score.min <- round(min(unbs@data$F_Value_Ref, na.rm = T), 3)
score.max <- round(max(unbs@data$F_Value_Ref, na.rm = T), 3)
df.col <- data.frame(x = seq(from = score.min, to = score.max, by = (score.max - score.min)/(1000-1)),
                     y = rep(1, 1000),
                     col = colorRampPalette(c("#f7fcb9", "#004529"))(1000))

api_key <- "####"

FireMap <- function(firename, zoom = 11, xmin = 0.02, xmax = 0.02, ymin = 0.01, ymax = 0.01) {
  perim <- perim.proj[perim.proj@data$Fire_Name == firename,]
  e <- SpatialPoints(t(bbox(perim)), proj4string = CRS(proj4string(perim)))
  cen <- gCentroid(e)
  g <- get_googlemap(center = cen@coords, force = T, zoom = zoom, key = api_key, maptype = "satellite")
  p <- suppressWarnings(tidy(perim))
  u <- suppressWarnings(as.data.frame(tidy(unb.proj[unb.proj@data$FireNam == firename,], region = "F_Value_Ref")))
  u$id <- as.numeric(u$id)
  map <- ggmap(g)+
    geom_polygon(data = p, aes(x = long, y = lat, group = group),
                 color = "white", fill = NA, size = 1)+
    geom_polygon(data = u, aes(x = long, y = lat, group = group, fill = id, color = id))+
    scale_fill_gradient(low = "#BAE4BC", high = "#253494", na.value = "#BAE4BC", limits = c(score.min, score.max))+
    scale_color_gradient(low = "#BAE4BC", high = "#253494", na.value = "#BAE4BC", limits = c(score.min, score.max))+
    coord_map(projection = "mercator", xlim = c(e@coords[1]-0.02, e@coords[2]+0.02), ylim = c(e@coords[3]-0.01, e@coords[4]+0.01))+
    theme_void()+
    guides(color = F)+
    labs(fill = "Refugia\nscore")
  return(map)
}

ScoreDensity <- function(firename){
  ggplot(data = data.frame(subset(df.all, FireName == firename), dummy = 1), aes(x = F_Value_Ref, y = dummy, fill = ..x..))+
    geom_density_ridges_gradient()+
    geom_vline(data = subset(df.all.mean, FireName == firename), aes(xintercept = F_Value_Ref), linetype = "dashed", size = .8)+
    geom_vline(data = subset(df.all.top, FireName == firename), aes(xintercept = F_Value_Ref), linetype = "dotted", size = .8)+
    scale_fill_gradient(low = "#BAE4BC", high = "#253494", na.value = "#BAE4BC", limits = c(score.min, score.max))+
    theme(axis.title = element_text(size = 12),
          legend.position = "none")+
    labs(y = "Probability density", x = "Refugia importance score")+
    coord_cartesian(xlim = c(score.min, score.max))
}

# Figure 4: B&B
map.BB <- FireMap("B&B COMPLEX (BOOTH)")
density.BB <- ScoreDensity("B&B COMPLEX (BOOTH)")

Fig.BB <- ggdraw()+
  draw_plot(map.BB + guides(fill = F, color = F), 0, 0, 0.5, 1)+
  draw_plot(density.BB, 0.5, 0, 0.5, 1)+
  draw_plot(get_legend(map.BB), 0.9, 0.75, 0.1, 0.1)+
  draw_label("a", 0.01, 0.95, colour = "white")+
  draw_label("b", 0.53, 0.95)
ggsave("Outputs/Figure4.jpg", Fig.BB, width = 6.5, height = 5, units = "in", dpi = 600)

# Figure 5: Davis
map.Davis <- FireMap("DAVIS", zoom = 12)
density.Davis <- ScoreDensity("DAVIS")

Fig.Davis <- ggdraw()+
  draw_plot(map.Davis + guides(fill = F, color = F), 0, 0, 0.5, 1)+
  draw_plot(density.Davis, 0.5, 0, 0.5, 1)+
  draw_plot(get_legend(map.Davis), 0.9, 0.75, 0.1, 0.1)+
  draw_label("a", 0.01, 0.95)+
  draw_label("b", 0.53, 0.95)
ggsave("Outputs/Figure5.jpg", Fig.Davis, width = 6.5, height = 5, units = "in", dpi = 600)

# Figure 6: Poison
map.Poison <- FireMap("POISON", zoom = 12)
density.Poison <- ScoreDensity("POISON")

Fig.Poison <- ggdraw()+
  draw_plot(map.Poison + guides(fill = F, color = F), 0, 0, 0.5, 1)+
  draw_plot(density.Poison, 0.5, 0, 0.5, 1)+
  draw_plot(get_legend(map.Poison), 0.9, 0.75, 0.1, 0.1)+
  draw_label("a", 0.01, 0.95)+
  draw_label("b", 0.53, 0.95)
ggsave("Outputs/Figure6.jpg", Fig.Poison, width = 6.5, height = 5, units = "in", dpi = 600)

# Figure 7: Pole Creek
map.Pole <- FireMap("POLE CREEK")
density.Pole <- ScoreDensity("POLE CREEK")

Fig.Pole <- ggdraw()+
  draw_plot(map.Pole + guides(fill = F, color = F), 0, 0, 0.5, 1)+
  draw_plot(density.Pole, 0.5, 0, 0.5, 1)+
  draw_plot(get_legend(map.Pole), 0.9, 0.75, 0.1, 0.1)+
  draw_label("a", 0.01, 0.95)+
  draw_label("b", 0.53, 0.95)
ggsave("Outputs/Figure7.jpg", Fig.Pole, width = 6.5, height = 5, units = "in", dpi = 600)


##########################
# Figure 8: Area/Frac    #
##########################
plot.area <- ggplot(data = df.all, aes(x = Area, color = var, fill = var))+
  geom_density(size = 1, alpha = .2, adjust = .1)+
  geom_vline(data = df.all.means, aes(xintercept = Area, color = var), linetype = "dashed", size = 1)+
  facet_wrap(.~FireName, ncol = 1)+
  labs(x = "Area (ha)", y = " ")+
  theme(legend.position = "bottom", legend.justification = "center", legend.title = element_blank())+
  coord_cartesian(xlim = c(0,3.5))

plot.frac <- ggplot(data = df.all, aes(x = frac, color = var, fill = var))+
  geom_density(size = 1, alpha = .2)+
  geom_vline(data = df.all.means, aes(xintercept = frac, color = var), linetype = "dashed", size = 1)+
  facet_wrap(.~FireName, ncol = 1)+
  labs(x = "FRAC", y = " ", fill = "Refugia score", color = "Refugia score")+
  theme(legend.position = "bottom", legend.justification = "center",
        #legend.spacing.x = unit(0.5, "cm"),
        legend.text = element_text(margin = margin(r = 6, unit = "pt")))+
  coord_cartesian(xlim = c(1, 1.25))
plot.legend.a <- get_legend(plot.frac)
plot.area.frac <- plot_grid(plot.area + theme(legend.position="none"), 
                            plot.frac + theme(legend.position="none"), 
                            nrow = 1,
                            labels = "auto")
ks.a.x <- rep(c(0.4, 0.9), each = 4)
ks.a.y <- rep(c(0.94, 0.72, 0.5, 0.28), times = 2)

Fig.area.frac <- ggdraw()+
  draw_plot(plot.area.frac, 0, 0.05, 1, 0.95)+
  draw_text("Probability density", x = 0.02, y = 0.5, angle = 90)+
  draw_plot(plot.legend.a, x = 0.04, y = -0.47)+
  draw_text(ks.area.frac, x = ks.a.x, y = ks.a.y, size = 10)

ggsave("Outputs/Figure8.jpg", Fig.area.frac, width = 6.5, height = 5, units = "in", dpi = 600)

##########################
# Figure 9: Topography   #
##########################

# Topographic metrics
plot.t.elev <- ggplot(data = df.all, aes(x = topo_elev_m, color = var, fill = var))+
  geom_density(size = 1, alpha = .2)+
  geom_vline(data = df.all.means, aes(xintercept = topo_elev_m, color = var), linetype = "dashed", size = 0.7)+
  facet_wrap(.~FireName, ncol = 1)+
  scale_x_continuous(breaks = c(1000, 2000), minor_breaks = c(500, 1500)) +
  scale_y_continuous(breaks = c(0, 0.004, 0.008)) +
  labs(x = "Elevation (m)", y = "")+
  theme(legend.position = "none", legend.justification = "center", legend.title = element_blank())

plot.t.trasp <- ggplot(data = df.all, aes(x = TRASP, color = var, fill = var))+
  geom_density(size = 1, alpha = .2)+
  geom_vline(data = df.all.means, aes(xintercept = TRASP, color = var), linetype = "dashed", size = 0.7)+
  facet_wrap(.~FireName, ncol = 1)+
  scale_y_continuous(breaks = c(0, 1, 2)) +
  labs(x = "TRASP", y = "")+
  scale_x_continuous(breaks = round(seq(0, 1, by = 0.5),1)) +
  theme(legend.position = "none", legend.justification = "center", legend.title = element_blank())

plot.t.slope <- ggplot(data = df.all, aes(x = topo_slope_147p636F_98p424FEET, color = var, fill = var))+
  geom_density(size = 1, alpha = .2)+
  geom_vline(data = df.all.means, aes(xintercept = topo_slope_147p636F_98p424FEET, color = var), linetype = "dashed", size = 0.7)+
  facet_wrap(.~FireName, ncol = 1)+
  labs(x = "Slope (Degrees) ", y = "")+
  theme(legend.position = "none", legend.justification = "center", legend.title = element_blank())

plot.topo <- plot_grid(plot.t.elev, plot.t.trasp, plot.t.slope,
                       rel_widths = c(3.45, 3.18, 3.15),
                       nrow = 1,
                       labels = "auto")


ks.t.x <- c(0.28, 0.19, 0.28, 0.19, 0.61, 0.61, 0.61, 0.61, 0.95, 0.95, 0.95, 0.95)
ks.t.y <- rep(c(0.935, 0.72, 0.51, 0.29), times = 3)

Fig.topo <- ggdraw()+
  draw_plot(plot.topo, 0, 0.03, 1, 0.97)+
  draw_text("Probability density", x = 0.02, y = 0.5, angle = 90)+
  draw_plot(plot.legend.a, x = 0.08, y = -0.47)+
  draw_text(ks.topo, x = ks.t.x, y = ks.t.y, size = 10)

ggsave("Outputs/Figure9.jpg", Fig.topo, width = 6.5, height = 5, units = "in", dpi = 600)

##########################
# Figure 10: Canopy      #
##########################
plot.c.agb <- ggplot(data = df.all, aes(x = AGB, color = var, fill = var))+
  geom_density(size = 1, alpha = .2)+
  geom_vline(data = df.all.means, aes(xintercept = AGB, color = var), linetype = "dashed", size = 0.7)+
  facet_wrap(.~FireName, ncol = 1)+
  labs(x = "AGB (Mg/ha)", y = "")+
  scale_x_continuous(breaks = c(0, 150, 300)) +
  scale_y_continuous(breaks = c(0, 0.01, 0.02)) +
  theme(legend.position = "none", legend.justification = "center", legend.title = element_blank())

plot.c.elev <- ggplot(data = df.all, aes(x = elev_ave_m, color = var, fill = var))+
  geom_density(size = 1, alpha = .2)+
  geom_vline(data = df.all.means, aes(xintercept = elev_ave_m, color = var), linetype = "dashed", size = 0.7)+
  facet_wrap(.~FireName, ncol = 1)+
  labs(x = "Mean point elev. (m)", y = "")+
  scale_y_continuous(breaks = c(0, 0.2, 0.4)) +
  #scale_x_continuous(breaks = round(seq(0, 1, by = 0.5),1)) +
  theme(legend.position = "none", legend.justification = "center", legend.title = element_blank())

plot.c.p95 <- ggplot(data = df.all, aes(x = elev_P95_m, color = var, fill = var))+
  geom_density(size = 1, alpha = .2)+
  geom_vline(data = df.all.means, aes(xintercept = elev_P95_m, color = var), linetype = "dashed", size = 0.7)+
  facet_wrap(.~FireName, ncol = 1)+
  labs(x = "95% point elev. (m)   ", y = "")+
  theme(legend.position = "none", legend.justification = "center", legend.title = element_blank())

plot.canopy <- plot_grid(plot.c.agb, plot.c.elev, plot.c.p95,
                         rel_widths = c(3.4, 3.3, 3.3),
                         nrow = 1,
                         labels = "auto")

ks.c.x <- c(0.28, 0.28, 0.28, 0.28, 0.61, 0.61, 0.61, 0.61, 0.95, 0.95, 0.95, 0.95)
ks.c.y <- rep(c(0.935, 0.72, 0.51, 0.29), times = 3)

Fig.canopy <- ggdraw()+
  draw_plot(plot.canopy, 0, 0.03, 1, 0.97)+
  draw_text("Probability density", x = 0.02, y = 0.5, angle = 90)+
  draw_plot(plot.legend.a, x = 0.08, y = -0.47)+
  draw_text(ks.canopy, x = ks.c.x, y = ks.c.y, size = 10)

ggsave("Outputs/Figure10.jpg", Fig.canopy, width = 6.5, height = 5, units = "in", dpi = 600)

##########################
# Table S1: Correlations #
##########################
df.used <- subset(df.all, select = c("Area", "frac", "topo_elev_m", "TRASP", 
                                     "topo_slope_147p636F_98p424FEET", "AGB", "elev_ave_m", "elev_P95_m"))
colnames(df.used) <- c("Area", "FRAC", "Elevation", "TRASP", "Slope", "AGB", "Mean return height", "95% return height")

df.used <- na.omit(df.used)
table.S2 <- round(cor(df.used), 3)
table.S2[upper.tri(table.S2)] <- NA
diag(table.S2) <- 1
View(table.S2)
write.table(table.S2, "clipboard", sep="\t", row.names=T, na = "", quote = F)
