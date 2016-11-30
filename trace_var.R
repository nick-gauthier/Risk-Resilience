# A script to import decadal means from TraCE outputs and plot patterns over
# the last deglaciation at three points in the western Mediterranean

library(raster) # netcdf import and raster processing
library(magrittr) # pipes for code readability
library(forecast) # fit arima model
library(reshape2) # melt data frames
library(ggplot2) # plotting
library(ncdf4)
library(rasterVis)
library(Bchron)

### Import TraCE data and extract values at three locations
setwd("~/Dropbox/ASU/Dissertation/Environment_data/Trace Variability") # Setup the directory

#####################################################################
## This little script exports average of different periods as tiff ##
#####################################################################

precip = brick('trace.01-36.22000BP.cam2.PRECT.22000BP_decavg_400BCE.nc')
temp = brick('trace.01-36.22000BP.cam2.TREFHT.22000BP_decavg_400BCE.nc')
extent(precip) = extent(0, 360, -90, 90)
rotated.precip <- rotate(precip)

grid <- extent(-20,20,30,60)  
europe.precip <- crop(rotated.precip, grid)
magda.precip = europe.precip[[201:801,]]
magda.precip = magda.precip * 2.592e+09

lower_magda_precip_a = mean(magda.precip[[1:201]]) ## Roughly 20-18 ka BP
lower_magda_precip_b = mean(magda.precip[[201:251]]) ## Roughly 18-17.5 ka BP
middle_magda_precip = mean(magda.precip[[251:401]]) ## Roughly 17.5-16 ka BP
upper_magda_precip_a = mean(magda.precip[[401:525]]) ## 16-14.75 ka BP
upper_magda_precip_b = mean(magda.precip[[425:601]]) ## 14.75-14 ka BP

writeRaster(lower_magda_precip_a, "lower_magda_a_precip.tif")
writeRaster(lower_magda_precip_b, "lower_magda_b_precip.tif")
writeRaster(middle_magda_precip, "middle_magda_precip.tif")
writeRaster(upper_magda_precip_a, "upper_magda_a_precip.tif")
writeRaster(upper_magda_precip_b, "upper_magda_b_precip.tif")

## Creating max and min maps
lower_max_precip_a = max(magda.precip[[1:201]]) ## Roughly 20-18 ka BP
lower_max_precip_b = max(magda.precip[[201:251]]) ## Roughly 18-17.5 ka BP
middle_max_precip = max(magda.precip[[251:401]]) ## Roughly 17.5-16 ka BP
upper_max_precip_a = max(magda.precip[[401:525]]) ## 16-14.75 ka BP
upper_max_precip_b = max(magda.precip[[425:601]]) ## 14.75-14 ka BP

writeRaster(lower_max_precip_a, "lower_max_a_precip.tif")
writeRaster(lower_max_precip_b, "lower_max_b_precip.tif")
writeRaster(middle_max_precip, "middle_max_precip.tif")
writeRaster(upper_max_precip_a, "upper_max_a_precip.tif")
writeRaster(upper_max_precip_b, "upper_max_b_precip.tif")

## Creating max and min maps
lower_min_precip_a = min(magda.precip[[1:201]]) ## Roughly 20-18 ka BP
lower_min_precip_b = min(magda.precip[[201:251]]) ## Roughly 18-17.5 ka BP
middle_min_precip = min(magda.precip[[251:401]]) ## Roughly 17.5-16 ka BP
upper_min_precip_a = min(magda.precip[[401:525]]) ## 16-14.75 ka BP
upper_min_precip_b = min(magda.precip[[425:601]]) ## 14.75-14 ka BP

writeRaster(lower_min_precip_a, "lower_min_a_precip.tif")
writeRaster(lower_min_precip_b, "lower_min_b_precip.tif")
writeRaster(middle_min_precip, "middle_min_precip.tif")
writeRaster(upper_min_precip_a, "upper_min_a_precip.tif")
writeRaster(upper_min_precip_b, "upper_min_b_precip.tif")

## For temperature

extent(temp) = extent(0, 360, -90, 90)
rotated.temp <- rotate(temp)

grid <- extent(-20,20,30,60)  
europe.temp <- crop(rotated.temp, grid)
magda.temp = europe.temp[[201:801,]]
magda.temp = magda.temp - 273.15

lower_magda_temp_a = mean(magda.temp[[1:201]]) ## Roughly 20-18 ka BP
lower_magda_temp_b = mean(magda.temp[[201:251]]) ## Roughly 18-17.5 ka BP
middle_magda_temp = mean(magda.temp[[251:401]]) ## Roughly 17.5-16 ka BP
upper_magda_temp_a = mean(magda.temp[[401:525]]) ## 16-14.75 ka BP
upper_magda_temp_b = mean(magda.temp[[425:601]]) ## 14.75-14 ka BP

writeRaster(lower_magda_temp_a, "lower_magda_a_temp.tif")
writeRaster(lower_magda_temp_b, "lower_magda_b_temp.tif")
writeRaster(middle_magda_temp, "middle_magda_temp.tif")
writeRaster(upper_magda_temp_a, "upper_magda_a_temp.tif")
writeRaster(upper_magda_temp_b, "upper_magda_b_temp.tif")

## Creating max and min maps
lower_max_temp_a = max(magda.temp[[1:201]]) ## Roughly 20-18 ka BP
lower_max_temp_b = max(magda.temp[[201:251]]) ## Roughly 18-17.5 ka BP
middle_max_temp = max(magda.temp[[251:401]]) ## Roughly 17.5-16 ka BP
upper_max_temp_a = max(magda.temp[[401:525]]) ## 16-14.75 ka BP
upper_max_temp_b = max(magda.temp[[425:601]]) ## 14.75-14 ka BP

writeRaster(lower_max_temp_a, "lower_max_a_temp.tif")
writeRaster(lower_max_temp_b, "lower_max_b_temp.tif")
writeRaster(middle_max_temp, "middle_max_temp.tif")
writeRaster(upper_max_temp_a, "upper_max_a_temp.tif")
writeRaster(upper_max_temp_b, "upper_max_b_temp.tif")

## Creating max and min maps
lower_min_temp_a = min(magda.temp[[1:201]]) ## Roughly 20-18 ka BP
lower_min_temp_b = min(magda.temp[[201:251]]) ## Roughly 18-17.5 ka BP
middle_min_temp = min(magda.temp[[251:401]]) ## Roughly 17.5-16 ka BP
upper_min_temp_a = min(magda.temp[[401:525]]) ## 16-14.75 ka BP
upper_min_temp_b = min(magda.temp[[425:601]]) ## 14.75-14 ka BP

writeRaster(lower_min_temp_a, "lower_min_a_temp.tif")
writeRaster(lower_min_temp_b, "lower_min_b_temp.tif")
writeRaster(middle_min_temp, "middle_min_temp.tif")
writeRaster(upper_min_temp_a, "upper_min_a_temp.tif")
writeRaster(upper_min_temp_b, "upper_min_b_temp.tif")

#########################
#########################

## This script transforms the Ocean temp and precip data to geoTIFF.

setwd("~/Dropbox/ASU/Dissertation/Environment_data/Modern_Sea_Temp") # Setup the directory

sst = brick('Ocean_Jan2013.nc')
grid <- extent(-20,20,30,60)  
sst.reduced <- crop(sst[[1]], grid)
writeRaster(sst.reduced, "Ocean_Jan2013.tif")

precip = brick('precip.mon.ltm.nc')
rotated.precip <- rotate(precip)
grid <- extent(-20,20,30,60)  
precip.reduced <- crop(rotated.precip, grid)
precip.mean = mean(precip.reduced[[1:12]])

writeRaster(precip.mean, "Ocean_Precip.tif")

#########################
#########################

# The 3.1104e+10 is to transform the precipitation into ml per year. 
# I used the downscaling scripts to calculate the difference I need to add based on weather stations.

loc = SpatialPoints(cbind(356, 43.4644))
precip.cantabria <- data.frame(value=extract(precip, loc)[1,] * 3.1104e+10 + 579.77, variable = "Precipitation (mm)", Region = "Cantabria", year = floor(temp@z[[1]]*-1000))
temp.cantabria <- data.frame(value=extract(temp, loc)[1,] - 273.15 + 4.10, variable = "Temperature (°C)", Region = "Cantabria", year = floor(temp@z[[1]]*-1000))

loc = SpatialPoints(cbind(0.675, 44.8317))
precip.aquitaine <- data.frame(value=extract(precip, loc)[1,] * 3.1104e+10 + 703.66, variable = "Precipitation (mm)", Region = "Aquitaine", year = floor(temp@z[[1]]*-1000))
temp.aquitaine <- data.frame(value=extract(temp, loc)[1,] - 273.15 + 2.34, variable = "Temperature (°C)", Region = "Aquitaine", year = floor(temp@z[[1]]*-1000))

precip.cantabria <- precip.cantabria[201:801,]
temp.cantabria <- temp.cantabria[201:801,]
precip.aquitaine <- precip.aquitaine[201:801,]
temp.aquitaine <- temp.aquitaine[201:801,]

var.dat = rbind(temp.cantabria,temp.aquitaine,precip.cantabria,precip.aquitaine)


# #This will be for Cantabria
# trace.prc <- brick('trace.01-36.22000BP.cam2.PRECT.22000BP_decavg_400BCE.nc')  %>%
#   raster::extract(matrix(c(356,0.675,43.46,44.8317), ncol = 2)) %>% # extract values at these coordinates (354 is for Cantabria, as this is on a 0-360 scale)
#   t %>% # transpose
#   extract(201:801,) %>% # get the values 20-14ka
#   multiply_by(2.592e+09) %>% # convert to mm/year
#   set_colnames(c("Cantabria", "Dordogne"))
# 
# trace.tmp <- brick('trace.01-36.22000BP.cam2.TREFHT.22000BP_decavg_400BCE.nc') %>%
#   raster::extract(matrix(c(356,0.675,43.46,44.8317), ncol = 2)) %>% t %>%
#   extract(201:801,) %>%
#   subtract(273.15) %>% # convert from kelvin to C
#   set_colnames(c("Cantabria", "Dordogne"))
# 
# ### Calculate detrended variances for each location and variable

# create a data frame of both variables from above
# var.dat <- cbind(trace.prc, trace.tmp)

# subset data frame by time period (time for Dordogne)
lm <- var.dat[1:250,] # Lower Magda
mm <- var.dat[251:400,] # Middle Magda
um <- var.dat[401:601,] # Upper Magda

# define a function to detrend each series using an automatic arima fit and then
# calculate the variances.

#*TO DO* make the outputs more readible, currently just prints out the results
# en masse in the order p.sw, p.nc, p.ne, t.sw, t.nc, t.ne
detrend.var <- function(x){
  for(i in 1:4){
    ar <- auto.arima(x[,i], seasonal = F)
    ar$residuals %>% var %>% print
  }
}

# apply the function to each time period
detrend.var(lm)
detrend.var(mm)
detrend.var(um)

###  Plot the the time series
# dat <- data.frame(Year = trace.tmp %>% melt %>% extract(,1) %>%
#                     as.character %>% strsplit(split="X.") %>% # remove 'x' from year names
#                     unlist %>%
#                     extract(seq(2,length(.), by=2)) %>%
#                     as.numeric,
#                   Region = trace.tmp %>% melt %>% extract(,2),
#                   'Temperature (°C)' = trace.tmp %>% melt %>% extract(,3),
#                   'Precipitation (mm)' = trace.prc %>% melt %>% extract(,3),
#                   check.names = F) %>%
#   melt(c("Year", 'Region'), c("Temperature (°C)", "Precipitation (mm)"))

ggplot(data = var.dat, aes(x = year, y = value)) +
  geom_vline(xintercept = c(20000,19000,17500,16500,14750,14000), color = 'grey', lwd = 1.5, alpha = .8) +
  geom_point(alpha = .8, aes(color = variable)) +
  facet_grid(variable~Region, scales = 'free_y', switch = 'y') +
  geom_smooth(color = "black", alpha = .8, span = 0.2, method = "loess") +
  labs(title = 'TraCE-21K Paleoclimate Simulation Decadal Means', x = '\n Years BP', y = '') +
  guides(color = "none") +
  scale_x_reverse() +
  theme_grey(base_size = 20) %+replace% theme(strip.background  = element_blank())

#save the output
ggsave('~/Dropbox/ASU/Dissertation/Environment_data/Trace Variability/TraCE Decadal Mean.png', width = 16.9, height = 8.61)

##CANTABRIAN DATES
lm <- var.dat[1:300,] # Lower Magda
mm <- var.dat[301:425,] # Middle Magda
um <- var.dat[426:601,] # Upper Magda

# define a function to detrend each series using an automatic arima fit and then
# calculate the variances.

#*TO DO* make the outputs more readible, currently just prints out the results
# en masse in the order p.sw, p.nc, p.ne, t.sw, t.nc, t.ne
detrend.var <- function(x){
  for(i in 1:4){
    ar <- auto.arima(x[,i], seasonal = F)
    ar$residuals %>% var %>% print
  }
}

# apply the function to each time period
detrend.var(lm)
detrend.var(mm)
detrend.var(um)
#detrend.var(ehol)


###  Plot the the time series
dat <- data.frame(Year = trace.tmp %>% melt %>% extract(,1) %>%
                    as.character %>% strsplit(split="X.") %>% # remove 'x' from year names
                    unlist %>%
                    extract(seq(2,length(.), by=2)) %>%
                    as.numeric,
                  Region = trace.tmp %>% melt %>% extract(,2),
                  'Temperature (°C)' = trace.tmp %>% melt %>% extract(,3),
                  'Precipitation (mm)' = trace.prc %>% melt %>% extract(,3),
                  check.names = F) %>%
  melt(c("Year", 'Region'), c("Temperature (°C)", "Precipitation (mm)"))


ggplot(data = dat, aes(x = Year, y = value)) +
  geom_vline(xintercept = c(17,15.75), color = 'grey', lwd = 1.5, alpha = .8) +
  geom_point(alpha = .8, aes(color = variable)) +
  facet_grid(variable~Region, scales = 'free_y', switch = 'y') +
  geom_smooth(color = "black", alpha = .8) +
  labs(title = 'TraCE-21K Paleoclimate Simulation Decadal Means', x = '\n 1,000 Years BP', y = '') +
  guides(color = "none") +
  scale_x_reverse(breaks = seq(12,22,2)) +
  theme_grey(base_size = 20) %+replace% theme(strip.background  = element_blank())

#save the output
ggsave('~/Dropbox/ASU/Dissertation/Environment_data/Trace Variability/TraCE Decadal Mean_Cantabria.png', width = 16.9, height = 8.61)


## CALIBRATING THE NGRIP DATA AND PLOTTING IT

ngrip <- read.csv("~/Dropbox/ASU/Dissertation/Environment_data/Trace Variability/NGRIP.csv")

ngrip = ngrip[501:1000,]

cal_ngrip = BchronCalibrate(ages = ngrip$age_bf_1950, ageSds = ngrip$error, calCurves = rep('intcal13',length(ngrip$age)))

for (i in 1:length(ngrip$age))
{ 
  date = paste("cal_ngrip$date",i,"$ageGrid",sep="")
  d = eval(parse(text = date))
  ngrip$cal_mean[i] = mean(d)
  ngrip$cal_std = sd(d)
}

dat <- subset(ngrip, cal_mean>=14000 & cal_mean<=20000)
age <- dat$cal_mean
d18o <- dat$d18o

ggplot(data = dat, mapping= aes(x = dat$cal_mean, y = dat$d18o)) + 
  geom_vline(xintercept = c(17500,16000), color = 'grey', lwd = 1.5, alpha = .8) +
  geom_line(color = "lightpink2") + 
  geom_smooth(method="loess", span=0.5, color="black", fill="grey", lwd=I(1)) + 
  labs(y="d18O\n", x="\ndate (cal years BP)") +
  scale_x_reverse()

## Uncalibrated graph

dat <- subset(ngrip, age_bf_1950>=14000 & age_bf_1950<=20000)
age <- dat$age_bf_1950
d18o <- dat$d18o

ggplot(data = dat, mapping= aes(x = dat$age_bf_1950, y = dat$d18o)) + 
  geom_vline(xintercept = c(17500,16000), color = 'grey', lwd = 1.5, alpha = .8) +
  geom_line(color = "lightpink2") + 
  geom_smooth(method="loess", span=0.5, color="black", fill="grey", lwd=I(1)) + 
  labs(y="d18O\n", x="\ndate (years BP)") +
  scale_x_reverse()

