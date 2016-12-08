# A script to import decadal means from TraCE outputs and plot patterns over
# the last deglaciation at three points in the western Mediterranean

library(raster) # netcdf import and raster processing
library(magrittr) # pipes for code readability
library(forecast) # fit arima model
library(reshape2) # melt data frames
library(ggplot2) # plotting

### Import TraCE data and extract values at three locations

trace.prc <- brick('trace.01-36.22000BP.cam2.PRECT.22000BP_decavg_400BCE.nc') %>%
  raster::extract(matrix(c(1,40, 4,42, 14,46), ncol = 2, byrow = T)) %>% # extract values at these coordinates
  t %>% # transpose
  extract(1:1600, ) %>% # get the values up to 6ka
  multiply_by(3.154e+10) %>% # convert to mm/year
  set_colnames(c("Southwest", "North Central", "Northeast"))

trace.tmp <- brick('trace.01-36.22000BP.cam2.TREFHT.22000BP_decavg_400BCE.nc') %>%
  raster::extract(matrix(c(1,40, 4,42, 14,46), ncol = 2, byrow = T)) %>% t %>%
  extract(1:1600, ) %>%
  subtract(273.15) %>% # convert from kelvin to C
  set_colnames(c("Southwest", "North Central", "Northeast"))

### Calculate detrended variances for each location and variable

# create a data frame of both variables from above
var.dat <- cbind(trace.prc, trace.tmp)

# subset data frame by time period
lgm <- var.dat[1:300,] # LGM
ltgl <- var.dat[301:800,] # Late Glacial
trns <- var.dat[801:1200,] # Transition
ehol <- var.dat[1201:1600,] # Early Holocene

# define a function to detrend each series using an automatic arima fit and then
# calculate the variances.

#*TO DO* make the outputs more readible, currently just prints out the results
# en masse in the order p.sw, p.nc, p.ne, t.sw, t.nc, t.ne
detrend.var <- function(x){
  for(i in 1:6){
    ar <- auto.arima(x[,i], seasonal = F)
    ar$residuals %>% var %>% print
  }
}

# apply the function to each time period
detrend.var(lgm)
detrend.var(ltgl)
detrend.var(trns)
detrend.var(ehol)


#now try doing the arima lines for the whole time span
ar <- auto.arima(var.dat[,1], seasonal = F)

  for(i in 1:6){
    ar <- auto.arima(x[,i], seasonal = F)
    ar$residuals %>% var %>% print
  }


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
  geom_vline(xintercept = c(22, 19, 14, 10, 6), color = 'yellow', lwd = 1.5, alpha = .8) +
  geom_point(alpha = .8, aes(color = variable)) +
  facet_grid(variable~Region, scales = 'free_y', switch = 'y') +
  geom_smooth(color = "black", alpha = .8) +
  labs(title = 'TraCE-21K Paleoclimate Simulation Decadal Means', x = '\n 1,000 Years BP', y = '') +
  guides(color = "none") +
  scale_x_reverse(breaks = seq(6,22,2)) +
  theme_grey(base_size = 20) %+replace% theme(strip.background  = element_blank())

#save the output
ggsave('TraCE Decadal Mean.png', width = 16.9, height = 8.61)

# black and white
ggplot(data = dat, aes(x = Year, y = value)) +
  geom_vline(xintercept = c(22, 19, 14, 10, 6), lty = 2) +
  geom_point(color = 'grey', alpha = .3) +
  facet_grid(variable ~ Region, scales = 'free_y', switch = 'y') +
  geom_smooth(color = "black", alpha = .8) +
  labs(x = '1,000 Years BP', y = '') +
  guides(color = "none") +
  scale_x_reverse(breaks = seq(6,22,4)) +
  theme_bw(base_size = 20) %+replace% theme(strip.background  = element_blank())

#save the output
ggsave('TraCE Decadal Mean - BW.png', width = 16.9, height = 8.61)
