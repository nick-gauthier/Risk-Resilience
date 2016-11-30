library(raster)
library(magrittr)
library(rasterVis)


importGCM <- function(dir, var, sim = 'lgm', level = 1){
  gcm.in <- brick(dir, var = var, level = level) 
  
  if(sim == 'control'){gcm.in <- gcm.in[[1200:1799]]}
  if(sim == 'lgm'){gcm.in <- gcm.in[[4213:4812]]}
  
  gcm.in <- stackApply(gcm.in, indices = 1:12, fun = mean)
  extent(gcm.in) <- extent(0, 360, -90, 90)
  
  gcm.in %>%    
    rotate %>%
    projectRaster(obs.prc) %>%
    mask(obs.prc)
}


control.u <- importGCM('~/Desktop/GCM/b40.20th.track1.1deg.005.cam2.h0.U.185001-200512.nc', 
                       var = "U", 
                       sim = 'control',
                       level = 26)

control.v <- importGCM('~/Desktop/GCM/b40.20th.track1.1deg.005.cam2.h0.V.185001-200512.nc', 
                       var = "V", 
                       sim = 'control',
                       level = 26)

streamplot(stack(control.u[[1]], control.v[[1]]), isField = 'dXY', main = "January Wind")
streamplot(stack(control.u[[2]], control.v[[2]]), isField = 'dXY', main = "February Wind")
streamplot(stack(control.u[[3]], control.v[[3]]), isField = 'dXY', main = "March Wind")
streamplot(stack(control.u[[4]], control.v[[4]]), isField = 'dXY', main = "April Wind")
streamplot(stack(control.u[[5]], control.v[[5]]), isField = 'dXY', main = "May Wind")
streamplot(stack(control.u[[6]], control.v[[6]]), isField = 'dXY', main = "June Wind")
streamplot(stack(control.u[[7]], control.v[[7]]), isField = 'dXY', main = "July Wind")
streamplot(stack(control.u[[8]], control.v[[8]]), isField = 'dXY', main = "August Wind")
streamplot(stack(control.u[[9]], control.v[[9]]), isField = 'dXY', main = "September Wind")
streamplot(stack(control.u[[10]], control.v[[10]]), isField = 'dXY', main = "October Wind")
streamplot(stack(control.u[[11]], control.v[[11]]), isField = 'dXY', main = "November Wind")
streamplot(stack(control.u[[12]], control.v[[12]]), isField = 'dXY', main = "December Wind")



#############
# CCSM PMIP2 outputs
region <- extent(c(-10, 50, 30, 50))
obs.prc <- getData('worldclim', var = 'tmin', res=5) %>% crop(region)

uas <- brick('~/Dropbox/uas_Aclim_GISS-E2-R_midHolocene_r1i1p1_250001-259912-clim.nc') %>% rotate %>% projectRaster(obs.prc) %>% mask(obs.prc)
vas <- brick('~/Dropbox/vas_Aclim_GISS-E2-R_midHolocene_r1i1p1_250001-259912-clim.nc') %>% rotate %>% projectRaster(obs.prc) %>% mask(obs.prc)

streamplot(stack(uas[[1]], vas[[1]]), isField = 'dXY', main = "January Wind")
streamplot(stack(uas[[2]], vas[[2]]), isField = 'dXY', main = "February Wind")
streamplot(stack(uas[[3]], vas[[3]]), isField = 'dXY', main = "March Wind")
streamplot(stack(uas[[4]], vas[[4]]), isField = 'dXY', main = "April Wind")
streamplot(stack(uas[[5]], vas[[5]]), isField = 'dXY', main = "May Wind")
streamplot(stack(uas[[6]], vas[[6]]), isField = 'dXY', main = "June Wind")
streamplot(stack(uas[[7]], vas[[7]]), isField = 'dXY', main = "July Wind")
streamplot(stack(uas[[8]], vas[[8]]), isField = 'dXY', main = "August Wind")
streamplot(stack(uas[[9]], vas[[9]]), isField = 'dXY', main = "September Wind")
streamplot(stack(uas[[10]], vas[[10]]), isField = 'dXY', main = "October Wind")
streamplot(stack(uas[[11]], vas[[11]]), isField = 'dXY', main = "November Wind")
streamplot(stack(uas[[12]], vas[[12]]), isField = 'dXY', main = "December Wind")
