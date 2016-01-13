###############################################################################
#  A script for fitting a GAMM to observed WorldClim data and the outputs of  #
#  CCSM4 climate models over Europe. Includes code to test model fit.         #
#  The GAMM models are saved at the end for use in other code.                #
#                                                                             #
#  Author: Nicolas Gauthier                                                   #
###############################################################################


###############################################################################
# Load necessary R packages

library(raster)     # functions for managing and processing raster data
library(mgcv)       # functions to fit and analyze GAMs
library(dismo)      # functions to sample points weighted by latitude
library(magrittr)   # piping functions for code readability


###############################################################################
# Import the observed present-day climatologies which we'll use to calibrate the
# model. We use WorldClim data at 5min resolution here, which can be easily 
# downloaded using the *getData()* function in the *raster* package. 

importWC <- function(var){
  raster::getData('worldclim', var = var, res = 5, download = T) %>% 
    crop(extent(-15, 40, 30, 58)) %>%     
    set_names(month.name)
}

obs.tmp <- importWC('tmean') %>% divide_by(10) # convert from degrees Celsius * 10
obs.prc <- importWC('prec')


###############################################################################
# Import climate-model outputs from the CCSM4 20th century historical run. Subset
# to the last 50 years of the simulation, average each month to create a
# climatology, rotate the map so that western hemisphere longitudes are negative,
# reproject to 5 arc minute resolution, and mask out ocean pixels

importGCM <- function(dir, level = 1){
  gcm.in <- brick(dir, level = level) %>%
    extract2(1201:1800) %>% # years 1950-2000
    stackApply(indices = 1:12, fun = mean) # monthly averages
  
  extent(gcm.in) <- extent(0, 360, -90, 90) # adjustment needed for rotate command
  gcm.in %>% rotate %>% projectRaster(obs.prc) %>% mask(obs.prc)
}

hist.prcl <- importGCM('GCM/b40.20th.track1.1deg.005.cam2.h0.PRECL.185001-200512.nc') %>% multiply_by(2629743830)
hist.prcc <- importGCM('GCM/b40.20th.track1.1deg.005.cam2.h0.PRECC.185001-200512.nc') %>% multiply_by(2629743830)
hist.tmp <- importGCM('GCM/b40.20th.track1.1deg.005.cam2.h0.TREFHT.185001-200512.nc') %>% subtract(273.15)
hist.psl <- importGCM('GCM/b40.20th.track1.1deg.005.cam2.h0.PSL.185001-200512.nc') %>% divide_by(1000)
hist.q <- importGCM('GCM/b40.20th.track1.1deg.005.cam2.h0.QREFHT.185001-200512.nc')
hist.u <- importGCM('GCM/b40.20th.track1.1deg.005.cam2.h0.U.185001-200512.nc', level = 26) # ~60m
hist.v <- importGCM('GCM/b40.20th.track1.1deg.005.cam2.h0.V.185001-200512.nc', level = 26) # ~60m
hist.z <- importGCM('GCM/b40.20th.track1.1deg.005.cam2.h0.Z3.185001-200512.nc', level = 23) # ~850 geopotential


###############################################################################
# Import a DEM and use it to generate topographic predictors: elevation,
# distance to the ocean, and the intensity of orographic lifting.

elev <- raster('alt_5m_bil/alt.bil') %>% crop(extent(-15, 40, 30, 58)) # import WorldClim 5m dem

# Calculate diffusive continentality (DCO) or distance to ocean in km.
dco <- raster('alt_5m_bil/alt.bil') %>% # reimport WorldClim 5m dem
  crop(extent(-20, 45, 30, 65)) %>%  # start with a wider region to get accurate distances
  reclassify(c(-Inf, Inf, NA, NA, NA, 1)) %>% # reverse NA and non-NA cells
  distance(doEdge = T) %>% # calculate the distances
  crop(extent(-15, 40, 30, 58)) %>% # crop to study region
  mask(elev) %>% # mask out ocean cells  
  divide_by(1000) # convert to km

# Calculate the velocity of orographic lifting as a function of wind direction,
# velocity, and terrain slope and aspect

dir <- (270 - overlay(hist.v, hist.u, fun = atan2) * (180 / pi)) %% 360
vel <- sqrt(hist.u ^ 2 + hist.v ^ 2)

aspect <- elev %>% terrain(opt='aspect', unit = 'degrees')
slope <- elev %>% terrain(opt = 'slope', unit = 'degrees')

delta <- abs(dir - aspect) # distance between wind angle and slope orientation
values(delta) <- ifelse(values(delta) > 180, 360 - values(delta), values(delta))
oro <- sin(slope * 2 * pi / 180) * cos(delta * pi / 180) * vel * .5


###############################################################################
# We can expect some random effects based on location and month, calculate these 
# values.

month <- obs.prc %>% setValues(c(as.factor(rep(1:12, each = ncell(oro)))))

lat <- elev %>%
  coordinates %>%
  extract( ,2) %>% 
  setValues(elev, .)

lon <- elev %>%
  coordinates %>%
  extract( ,1) %>% 
  setValues(elev, .)


###############################################################################
# Put all the predictor and response variables together, month by month, and
# remove the original files from the workspace.

cal.vars <- sapply(1:12, function(x){ 
  brick(obs.tmp[[x]], obs.prc[[x]], hist.tmp[[x]], hist.q[[x]], hist.prcc[[x]],
        hist.prcl[[x]], hist.psl[[x]], hist.z[[x]], elev, dco, oro[[x]], 
        month[[x]], lat, lon) %>%
    setNames(c('obs.tmp', 'obs.prc','TREFHT', 'Q', 'PRCC', 'PRCL', 'PSL', 'Z3', 
               'elev','dco', 'oro', 'month', 'lat','lon'))
})

rm(hist.tmp, hist.q, hist.prcc, hist.prcl, hist.psl, hist.z, hist.v, hist.u, 
   slope, aspect, oro, dir, vel, delta)

# Sample the variables at random points, weighting for latitude
cal.data <- lapply(cal.vars, function(x) (raster::extract(x, randomPoints(elev, 20000)) %>% data.frame)) %>% 
  do.call(rbind, .)



###############################################################################
# Fit the GAM for temperature

fit.tmp <- gam(obs.tmp ~ s(TREFHT, bs = 'cr') +
                 s(Z3, bs = 'cr') +
                 s(elev, bs = 'cr') +
                 s(dco, bs = 'cr'), 
               method = 'REML', data = cal.data)


fit.tmp
summary(fit.tmp)
gam.check(fit.tmp)
plot(fit.tmp, shade=T, seWithMean = T, pages = 1)


###############################################################################
# Fit the GAMs for precipitation occurrence and amount

fit.prc.occur <- bam(factor(obs.prc >= 1) ~ s(PRCC) + 
                       s(PRCL) + 
                       s(Z3) + 
                       s(TREFHT),
                     family = binomial, method = 'REML', data = cal.data)

summary(fit.prc.occur)
gam.check(fit.prc.occur)
plot(fit.prc.occur, seWithMean = T, shade = T, pages = 1)
levelplot(predict(cal.vars[[7]], fit.prc.occur, type = 'response'))
levelplot(obs.prc[[7]] < 1)


fit.prc <- bam(obs.prc ~ s(PRCC, bs = 'cr') +
                 s(PRCL, bs = 'cr') +
                 s(PSL, bs = 'cr') +
                 s(Q, bs = 'cr') + 
                 s(oro, bs = 'cr') +
                 s(Z3, bs = 'cr') +
                 s(elev, bs = 'cr') +
                 s(dco, bs = 'cr') +
                 s(month, bs = 're'),
               family = Gamma(link = 'log'), method = 'REML', data = cal.data[cal.data$obs.prc >= 1, ])


## set up an automatic system to validate outside of calibration domain, but same lat
# just calculate the RMSE and try to minimize