library(raster)
library(magrittr)

t31 <- brick('~/Downloads/trace.01-36.22000BP.cam2.PRECT.22000BP_decavgDJF_400BCE.nc') %>% 
  extract2(1) %>%
  rotate %>%
  rasterToPolygons %T>%
  plot

library(maptools)

writeSpatialShape(t31, 't31grid')
