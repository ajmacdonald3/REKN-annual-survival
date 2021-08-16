####################### PREPARE SURVIVAL COVARIATES ############################

library(tidyverse) # data manipulation
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(rgdal) # package for geospatial analysis
library(rasterVis) # plotting rasters
library(viridis) # colour palette
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

########################### SNOW COVER COVARIATE ###############################

####---- 2009 ----####

# load data in netCDF format

nc_2009 <- nc_open("./data/snowc.2009.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2009.nc_metadata.txt")
  print(nc_2009)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2009 <- brick("./data/snowc.2009.nc", varname="snowc")

# subset to key nest initiation period
brick_2009 <- subset(brick_2009, 166:181) # june 15-30

# load shapefile of preferred rufa REKN breeding habitat
shp <- shapefile("./data/Preferred_Habitat_and_Range_Intersect.shp")
rufa_range <- spTransform(shp,  crs(brick_2009)) # reproject to same CRS as rasterbrick
rufa_range@bbox <- as.matrix(extent(4923499, 7624742, 6175588, 7534838))
plot(rufa_range)

# keep only cells inside rufa range
brick_2009_rufa <- mask(brick_2009, rufa_range)

# calculate mean value per cell for june 15-30
mean_2009_rufa <- calc(brick_2009_rufa, fun = mean, na.rm = TRUE)

mean_2009_rufa <- crop(x = mean_2009_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2009 <- getValues(mean_2009_rufa)
mean_snowc_2009 <- mean(vals_2009, na.rm=TRUE)

####---- 2010 ----####

# load data in netCDF format

nc_2010 <- nc_open("./data/snowc.2010.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2010.nc_metadata.txt")
  print(nc_2010)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2010 <- brick("./data/snowc.2010.nc", varname="snowc")

# subset to key nest initiation period
brick_2010 <- subset(brick_2010, 166:181) # june 15-30

# keep only cells inside rufa range
brick_2010_rufa <- mask(brick_2010, rufa_range)

# calculate mean value per cell for june 15-30
mean_2010_rufa <- calc(brick_2010_rufa, fun = mean, na.rm = TRUE)

mean_2010_rufa <- crop(x = mean_2010_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2010 <- getValues(mean_2010_rufa)
mean_snowc_2010 <- mean(vals_2010, na.rm=TRUE)

####---- 2011 ----####

# load data in netCDF format

nc_2011 <- nc_open("./data/snowc.2011.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2011.nc_metadata.txt")
  print(nc_2011)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2011 <- brick("./data/snowc.2011.nc", varname="snowc")

# subset to key nest initiation period
brick_2011 <- subset(brick_2011, 166:181) # june 15-30

# keep only cells inside rufa range
brick_2011_rufa <- mask(brick_2011, rufa_range)

# calculate mean value per cell for june 15-30
mean_2011_rufa <- calc(brick_2011_rufa, fun = mean, na.rm = TRUE)

mean_2011_rufa <- crop(x = mean_2011_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2011 <- getValues(mean_2011_rufa)
mean_snowc_2011 <- mean(vals_2011, na.rm=TRUE)

####--- 2012 ----####

# load data in netCDF format

nc_2012 <- nc_open("./data/snowc.2012.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2012.nc_metadata.txt")
  print(nc_2012)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2012 <- brick("./data/snowc.2012.nc", varname="snowc")

# subset to key nest initiation period
brick_2012 <- subset(brick_2012, 167:182) # june 15-30

# keep only cells inside rufa range
brick_2012_rufa <- mask(brick_2012, rufa_range)

# calculate mean value per cell for june 15-30
mean_2012_rufa <- calc(brick_2012_rufa, fun = mean, na.rm = TRUE)

mean_2012_rufa <- crop(x = mean_2012_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2012 <- getValues(mean_2012_rufa)
mean_snowc_2012 <- mean(vals_2012, na.rm=TRUE)

####---- 2013 ----####

# load data in netCDF format

nc_2013 <- nc_open("./data/snowc.2013.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2013.nc_metadata.txt")
  print(nc_2013)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2013 <- brick("./data/snowc.2013.nc", varname="snowc")

# subset to key nest initiation period
brick_2013 <- subset(brick_2013, 166:181) # june 15-30

# keep only cells inside rufa range
brick_2013_rufa <- mask(brick_2013, rufa_range)

# calculate mean value per cell for june 15-30
mean_2013_rufa <- calc(brick_2013_rufa, fun = mean, na.rm = TRUE)

mean_2013_rufa <- crop(x = mean_2013_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2013 <- getValues(mean_2013_rufa)
mean_snowc_2013 <- mean(vals_2013, na.rm=TRUE)

####---- 2014 ----####

# load data in netCDF format

nc_2014 <- nc_open("./data/snowc.2014.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2014.nc_metadata.txt")
  print(nc_2014)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2014 <- brick("./data/snowc.2014.nc", varname="snowc")

# subset to key nest initiation period
brick_2014 <- subset(brick_2014, 166:181) # june 15-30

# keep only cells inside rufa range
brick_2014_rufa <- mask(brick_2014, rufa_range)

# calculate mean value per cell for june 15-30
mean_2014_rufa <- calc(brick_2014_rufa, fun = mean, na.rm = TRUE)

mean_2014_rufa <- crop(x = mean_2014_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2014 <- getValues(mean_2014_rufa)
mean_snowc_2014 <- mean(vals_2014, na.rm=TRUE)

####---- 2015 ----####

# load data in netCDF format

nc_2015 <- nc_open("./data/snowc.2015.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2015.nc_metadata.txt")
  print(nc_2015)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2015 <- brick("./data/snowc.2015.nc", varname="snowc")

# subset to key nest initiation period
brick_2015 <- subset(brick_2015, 166:181) # june 15-30

# keep only cells inside rufa range
brick_2015_rufa <- mask(brick_2015, rufa_range)

# calculate mean value per cell for june 15-30
mean_2015_rufa <- calc(brick_2015_rufa, fun = mean, na.rm = TRUE)

mean_2015_rufa <- crop(x = mean_2015_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2015 <- getValues(mean_2015_rufa)
mean_snowc_2015 <- mean(vals_2015, na.rm=TRUE)

####---- 2016 ----####

# load data in netCDF format

nc_2016 <- nc_open("./data/snowc.2016.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2016.nc_metadata.txt")
  print(nc_2016)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2016 <- brick("./data/snowc.2016.nc", varname="snowc")

# subset to key nest initiation period
brick_2016 <- subset(brick_2016, 167:182) # june 15-30

# keep only cells inside rufa range
brick_2016_rufa <- mask(brick_2016, rufa_range)

# calculate mean value per cell for june 15-30
mean_2016_rufa <- calc(brick_2016_rufa, fun = mean, na.rm = TRUE)

mean_2016_rufa <- crop(x = mean_2016_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2016 <- getValues(mean_2016_rufa)
mean_snowc_2016 <- mean(vals_2016, na.rm=TRUE)

####---- 2017 ----####

# load data in netCDF format

nc_2017 <- nc_open("./data/snowc.2017.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2017.nc_metadata.txt")
  print(nc_2017)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2017 <- brick("./data/snowc.2017.nc", varname="snowc")

# subset to key nest initiation period
brick_2017 <- subset(brick_2017, 166:181) # june 15-30

# keep only cells inside rufa range
brick_2017_rufa <- mask(brick_2017, rufa_range)

# calculate mean value per cell for june 15-30
mean_2017_rufa <- calc(brick_2017_rufa, fun = mean, na.rm = TRUE)

mean_2017_rufa <- crop(x = mean_2017_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2017 <- getValues(mean_2017_rufa)
mean_snowc_2017 <- mean(vals_2017, na.rm=TRUE)

####---- 2018 ----####

# load data in netCDF format

nc_2018 <- nc_open("./data/snowc.2018.nc")
# Save the print(nc) dump to a text file
{
  sink("snowc.2018.nc_metadata.txt")
  print(nc_2018)
  sink()
}

# convert to rasterbrick format (layered rasters - one per day)
brick_2018 <- brick("./data/snowc.2018.nc", varname="snowc")

# subset to key nest initiation period
brick_2018 <- subset(brick_2018, 166:181) # june 15-30

# keep only cells inside rufa range
brick_2018_rufa <- mask(brick_2018, rufa_range)

# calculate mean value per cell for june 15-30
mean_2018_rufa <- calc(brick_2018_rufa, fun = mean, na.rm = TRUE)

mean_2018_rufa <- crop(x = mean_2018_rufa, y = extent(rufa_range)) # crop to zoom in on breeding range

# calculate overall mean
vals_2018 <- getValues(mean_2018_rufa)
mean_snowc_2018 <- mean(vals_2018, na.rm=TRUE)

####---- all years ----####

# create dataframe of annual means

Year <- c(2009:2018)
mean_snowc <- c(mean_snowc_2009, mean_snowc_2010, mean_snowc_2011, mean_snowc_2012,
                mean_snowc_2013, mean_snowc_2014, mean_snowc_2015, mean_snowc_2016,
                mean_snowc_2017, mean_snowc_2018)

arctic_snow_cov <- as.data.frame(cbind(Year, mean_snowc))

arctic_snow_cov <- arctic_snow_cov %>% 
  filter(!Year == 2009) %>% 
  mutate(Period = 1:9)

# plot rasters
colr <- colorRampPalette(RColorBrewer::brewer.pal(9, 'GnBu'))

# 2009

snow_raster_2009 <- levelplot(mean_2009_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2009 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
                              layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2009.png", width = 700, height = 700, res = 110)
print(snow_raster_2009)
dev.off()

# # set up basic things for maps 
# theme_set(theme_bw())
# world <- ne_countries(scale = "medium", returnclass = "sf")
# #lakes <- ne_download(scale = "medium", type = 'lakes', category = 'physical',
# #returnclass = "sf", destdir = "./map-data/lakes") # only need this first time downloading, saves shape file in folder called map-data in your working directory
# lakes <- ne_load(type = "lakes", scale = "medium", category = 'physical',
#                  returnclass = "sf",
#                  destdir = paste0(getwd(), "./map-data/lakes")) # use this if already downloaded shapefiles
# 
# # get better resolution basemap
# region <- ne_states(country = "Canada", returnclass = "sf") %>%
#   filter(name %in% c("Northwest Territories", "Nunavut"))
# 
# # plot map
# mean_2009_rufa_spdf <- as(mean_2009_rufa, "SpatialPixelsDataFrame")
# mean_2009_rufa_df <- as.data.frame(mean_2009_rufa_spdf)
# rufa_range_sf <- st_as_sf(rufa_range)
# 
# region <- st_transform(region, crs = st_crs(rufa_range_sf))
# lakes <- st_transform(lakes, crs = st_crs(rufa_range_sf))
# 
# xmin <- min(region$longitude) - 1.25
# xmax <- max(region$longitude) + 0.75
# ymin <- min(region$latitude) - 0.5
# ymax <- max(region$latitude) + 0.75
# 
# 
# test <- ggplot() +
#   geom_sf(data = region, colour = NA) +
#   geom_sf(data = lakes, fill = "white", colour = NA) +
#   geom_raster(data = mean_2009_rufa_df, aes(x = x, y = y, fill = layer)) +
#   geom_sf(data = rufa_range_sf, colour = "black", fill = NA) +
#   coord_sf(xlim = c(xmin, xmax),
#            ylim = c(ymin, ymax),
#            expand = FALSE)
# 
# png("test.png", width = 8, height = 6, units = "in", res = 600)
# 
# print(test)
# 
# dev.off()

# 2010

snow_raster_2010 <- levelplot(mean_2010_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2010 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
  layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2010.png", width = 700, height = 700, res = 110)
print(snow_raster_2010)
dev.off()

# 2011

snow_raster_2011 <- levelplot(mean_2011_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2011 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
  layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2011.png", width = 700, height = 700, res = 110)
print(snow_raster_2011)
dev.off()

# 2012

snow_raster_2012 <- levelplot(mean_2012_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2012 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
  layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2012.png", width = 700, height = 700, res = 110)
print(snow_raster_2012)
dev.off()

# 2013

snow_raster_2013 <- levelplot(mean_2013_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2013 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
  layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2013.png", width = 700, height = 700, res = 110)
print(snow_raster_2013)
dev.off()

# 2014

snow_raster_2014 <- levelplot(mean_2014_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2014 mean snow cover'),
                                par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                                col.regions=colr,                   
                                at=seq(0, 1, len=101)) +            
                              layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2014.png", width = 700, height = 700, res = 110)
print(snow_raster_2014)
dev.off()

# 2015

snow_raster_2015 <- levelplot(mean_2015_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2015 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
  layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2015.png", width = 700, height = 700, res = 110)
print(snow_raster_2015)
dev.off()

# 2016

snow_raster_2016 <- levelplot(mean_2016_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2016 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
  layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2016.png", width = 700, height = 700, res = 110)
print(snow_raster_2016)
dev.off()

# 2017

snow_raster_2017 <- levelplot(mean_2017_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2017 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
  layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2017.png", width = 700, height = 700, res = 110)
print(snow_raster_2017)
dev.off()

# 2018

snow_raster_2018 <- levelplot(mean_2018_rufa, 
                              margin=FALSE,                       
                              colorkey=list(
                                space='right',
                                labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1)),
                                axis.line=list(col='black')       
                              ),    
                              main=list(label='2018 mean snow cover'),
                              par.settings=list(
                                axis.line=list(col='transparent') 
                              ),
                              scales=list(draw=FALSE),            
                              col.regions=colr,                   
                              at=seq(0, 1, len=101)) +            
  layer(sp.polygons(rufa_range, lwd=0.01))

png(file="figures/snow_raster_2018.png", width = 700, height = 700, res = 110)
print(snow_raster_2018)
dev.off()

# all in one plot
stack_mean_rufa <- stack(mean_2009_rufa, mean_2010_rufa, mean_2011_rufa, mean_2012_rufa,
                         mean_2013_rufa, mean_2014_rufa, mean_2015_rufa, mean_2016_rufa,
                         mean_2017_rufa)
snow_rasters <- levelplot(stack_mean_rufa, 
          layout = c(3,3),
          margin=FALSE,                       
          colorkey=list(
            space='bottom',  
            title = "Mean snow cover\n(proportion)",
            labels=list(at=c(0, 0.2, 0.4, 0.6, 0.8, 1), font=2),
            axis.line=list(col='black'),
            width=0.75
          ),    
          par.settings=list(
            strip.border=list(col='transparent'),
            strip.background=list(col='transparent'),
            axis.line=list(col='transparent')
          ),
          scales=list(draw=FALSE),            
          col.regions=colr,                   
          at=seq(0, 1, len=50),
          names.attr=c("2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017")) +           
  layer(sp.polygons(rufa_range, lwd=0.1))

png(file="figures/snow_rasters.png")
print(snow_rasters)
dev.off()

########################### TEMPERATURE COVARIATE ##############################

# load data

ymo_weather <- read_csv("data/moosonee_weather_data.csv")

# filter and summarize annual means/totals

ymo_weather <- ymo_weather %>% 
  dplyr::select(Year, Month, `Mean Temp (C)`, `Total Precip (mm)`) %>% 
  filter(Month %in% c(5,6,7,8))

ymo_temp_cov <- ymo_weather %>% 
  filter(Month %in% c(5,6)) %>% 
  group_by(Year) %>% 
  summarize(spring_temp = mean(`Mean Temp (C)`, na.rm = TRUE)) %>% 
  filter(!Year == 2009) %>% 
  mutate(Period = 1:9)

######################### PRECIPITATION COVARIATE ##############################

ymo_precip_cov <- ymo_weather %>% 
  filter(Month %in% c(7,8)) %>% 
  group_by(Year) %>% 
  summarize(summer_precip = sum(`Total Precip (mm)`, na.rm = TRUE)) %>% 
  filter(!Year == 2018) %>% 
  mutate(Period = 1:9)

######################## DB DEPARTURE MASS COVARIATE ###########################

# load data

db_departmass <- read_csv("data/DB_departuremass_data.csv")

# summarize annual means

db_mass_cov <- db_departmass %>% 
  group_by(Year) %>% 
  summarize(dep_mass = mean(depart_mass)) %>% 
  filter(!Year == 2009) %>% 
  mutate(Period = 1:8) 
  
mean_dm <- mean(db_mass_cov$dep_mass)

############################ MERGE ALL COVARIATES ##############################

surv_covs <- left_join(ymo_temp_cov %>% dplyr::select(-Year),
                       ymo_precip_cov %>% dplyr::select(-Year))

surv_covs <- left_join(surv_covs, db_mass_cov %>% dplyr::select(-Year))

surv_covs <- left_join(surv_covs, arctic_snow_cov %>% dplyr::select((-Year))) %>% 
  dplyr::select(Period, spring_temp, summer_precip, dep_mass, mean_snowc) %>%  
  replace_na(list(dep_mass = mean_dm))

saveRDS(surv_covs, file = "data/survival_covariates.rds")

# standardize covariates

surv_covs_std <- surv_covs %>% 
  mutate(spring_temp_std = (spring_temp - mean(spring_temp))/sd(spring_temp),
         summer_precip_std = (summer_precip - mean(summer_precip))/sd(summer_precip),
         dep_mass_std = (dep_mass - mean(dep_mass))/sd(dep_mass),
         mean_snowc_std = (mean_snowc - mean(mean_snowc))/sd(mean_snowc)) %>% 
  dplyr::select(-spring_temp, -summer_precip, -dep_mass, -mean_snowc)

saveRDS(surv_covs_std, file = "data/survival_covariates_standardized.rds")
