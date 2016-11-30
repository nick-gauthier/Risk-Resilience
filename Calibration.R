library(raster)     # functions for managing and processing raster data
library(mgcv)       # functions to fit and analyze GAMs
library(magrittr)   # piping functions for code readability
library(rasterVis)


cal.data <- read.csv('cal_data.csv')
load('val_vars.RData')

importWC <- function(var){
  raster::getData('worldclim', var = var, res = 5, download = T) %>% 
    crop(extent(-125, -70, 30, 58)) %>%     
    set_names(month.name)
}

obs.tmp <- importWC('tmean') %>% divide_by(10) # convert from degrees Celsius * 10
obs.prc <- importWC('prec')


###############################################################################
# Fit the GAM for temperature

fit.tmp <- gam(obs.tmp ~ s(TREFHT, bs = 'cr') +
                 s(Z3, bs = 'cr') +
                 s(elev, bs = 'cr'),
               method = 'REML', data = cal.data)


lapply(1:12, function(x){predict(val.vars[[x]], fit.tmp, type = 'response')}) %>% 
  brick %>%
  subtract(obs.tmp) %>%
  raise_to_power(2) %>%
  cellStats(mean) %>%
  mean %>%
  sqrt
  
test <- lapply(1:12, function(x){predict(val.vars[[x]], fit.tmp, type = 'response')}) %>% 
  brick 
levelplot(brick(c(test[[c(1,7)]],obs.tmp[[c(1,7)]])))
levelplot(test - obs.tmp, par.settings = BuRdTheme)




fit.tmp
summary(fit.tmp)
gam.check(fit.tmp)
plot(fit.tmp, shade=T, seWithMean = T, pages = 1)




###############################################################################
# Fit the GAMs for precipitation occurrence and amount

fit.prc.occur <- bam(factor(obs.prc >= 1) ~ s(PRCC) + 
                       s(PRCL) + 
                       s(Z3),
                     family = binomial, method = 'REML', data = cal.data)

summary(fit.prc.occur)
gam.check(fit.prc.occur)
plot(fit.prc.occur, seWithMean = T, shade = T, pages = 1)
levelplot(predict(cal.vars[[7]], fit.prc.occur, type = 'response'))
levelplot(obs.prc[[7]] < 1)


fit.prc <- bam(obs.prc ~ s(PRCC, bs = 'cr') +
                 s(PRCL, bs = 'cr') +
                 s(Z3, bs = 'cr'),
               family = Gamma(link = 'log'), method = 'REML', data = cal.data[cal.data$obs.prc >= 1, ])


lapply(1:12, function(x){predict(val.vars[[x]], fit.prc, type = 'response')}) %>% 
  brick %>%
  subtract(obs.prc) %>%
  raise_to_power(2) %>%
  mean %>%
  cellStats(mean) %>%
  sqrt

levelplot(obs.prc[[1]])
levelplot((predict(val.vars[[6]], fit.prc, type = 'response') - obs.prc[[6]])/obs.prc[[6]])
levelplot(brick(c(predict(val.vars[[6]], fit.prc, type = 'response'), obs.prc[[6]])))

