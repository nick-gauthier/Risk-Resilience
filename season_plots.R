library(raster)
library(tidyverse)
library(magrittr)
samp.pts <- matrix(c(0, 40), ncol = 2, byrow = T)

getTrace <- function(x){
  brick(x) %>%
    raster::extract(samp.pts) %>%
    multiply_by(3.154e+10) %>%
    t
}

library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(4, "Spectral")))



plt.dat <- data_frame(year = rownames(getTrace('~/Downloads/trace.01-36.22000BP.cam2.PRECT.22000BP_decavgDJF_400BCE.nc')),
                   a = getTrace('~/Downloads/trace.01-36.22000BP.cam2.PRECT.22000BP_decavgDJF_400BCE.nc')[,1],
                   b = getTrace('~/Downloads/trace.01-36.22000BP.cam2.PRECT.22000BP_decavgMAM_400BCE.nc')[,1],
                   c = getTrace('~/Downloads/trace.01-36.22000BP.cam2.PRECT.22000BP_decavgJJA_400BCE.nc')[,1],
                   d = getTrace('~/Downloads/trace.01-36.22000BP.cam2.PRECT.22000BP_decavgSON_400BCE.nc')[,1]) %>%
    mutate(year = as.numeric(substring(year, 3)) * -1) %>%#filter(year < -6) %>%
  mutate(century = round(year, 1)) %>%
  group_by(century) %>%
  summarise_each(funs(median)) %>%
  select(-century) %>%
  gather(key, value, -year) %>%
  mutate(key = as.numeric(as.factor(key)))

ggplot(data = plt.dat, aes(x = key, y = value, group = year, color = year)) + 
  geom_smooth(se = F) + 
  scale_color_gradientn(colors = myPalette(100)) +
  theme_minimal()


###### GAM
library(mgcv)

trace.gam 