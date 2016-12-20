# West Mediterranean during the Last Deglaciation
Nick Gauthier  
`r Sys.Date()`  



# Introduction

We'll be comparing paleoclimate model estimates of temperature and precipitation over three points in the west Mediterranean to global paleoclimate proxies.

## Setup
Load all the packages we'll need for this analysis. 

<label for="tufte-mn-" class="margin-toggle">&#8853;</label><input type="checkbox" id="tufte-mn-" class="margin-toggle"><span class="marginnote">You'll need to have the netCDF libraries already installed on your system for ncdf4 to work.</span>

```r
library(ncdf4) # import GCM data
library(rgdal) # read GCM data
library(raster) # process GCM data
library(rasterVis) # plotting GCM data
library(tidyverse) # data management and plotting
library(magrittr) # pipes for code readability
library(EMD) # calculate trends in the data
```


# Climate Model

## Sample Locations
<p><span class="marginnote shownote"><!--
<div class="figure">-->
<img src="proxymodelcomparison_files/figure-html/map-1.png" alt="Locations of 3 sample points."  />
<!--
<p class="caption marginnote">-->Locations of 3 sample points.<!--</p>-->
<!--</div>--></span></p>

Create a matrix with the coordinates for the three locations of interest in the west Mediterranean. We'll be focusing on large grid cell averages, so the points do not have to be directly over land.

```r
samp.pts <- matrix(c(1, 40, 4, 42, 14, 46), 
                   ncol = 2, byrow = T)
```

## TraCE-21k

<p><span class="marginnote shownote"><!--
<div class="figure">-->
<img src="proxymodelcomparison_files/figure-html/trace_map-1.png" alt="TraCE21-k global precipitation and temperature"  /><img src="proxymodelcomparison_files/figure-html/trace_map-2.png" alt="TraCE21-k global precipitation and temperature"  />
<!--
<p class="caption marginnote">-->TraCE21-k global precipitation and temperature<!--</p>-->
<!--</div>--></span></p>

First, import data from the TraCE-21k paleoclimate simulation. Then extract temperature and precipitation values at three locations in the west Mediterranean. Use the *brick* function from **raster** to import decadal averages from the simulation. Put the coordinates for the three locations in a matrix, and use that matrix to and **raster's** *extract* function to get the values from the climate model brick. Convert the precipitation values to mm/year and temperature values to degrees Celsius. Finally, name the columns for each region appropriately. 

Now pull all the TraCE data into one data frame, with one row per year, and one column per variable/location combination. First *rbind* the two sets of TraCE data and *transpose* the results, turning the 6 rows into 6 columns. Add a column for the Year (in ka BP), and use to select only the entries earlier than 6,000 BP.


```r
trace.dat <- rbind(
  brick('trace.01-36.22000BP.cam2.PRECT.22000BP_decavg_400BCE.nc') %>%
    raster::extract(samp.pts) %>% # extract data at these coordinates
    multiply_by(3.154e+10), # convert to mm/year
  brick('trace.01-36.22000BP.cam2.TREFHT.22000BP_decavg_400BCE.nc') %>%
    raster::extract(samp.pts) %>% 
    subtract(273.15)) %>% # convert from kelvin to C
  t %>% # transpose
  as.data.frame %>%
  set_colnames(c('prc,Southwest', 'prc,North Central', 'prc,Northeast', 
                 'tmp,Southwest', 'tmp,North Central', 'tmp,Northeast')) %>%
  rownames_to_column('Year') %>%
  mutate(Year = as.numeric(substring(Year, 3))) %>%
  filter(Year > 6) # get all the decades up to 6ka BP
```



## Trend Analysis

Let's use the **EMD** package to calculate actual trend lines using the empirical mode decomposition approach.
<p><span class="marginnote shownote"><!--
<div class="figure">-->
<img src="proxymodelcomparison_files/figure-html/unnamed-chunk-2-1.png" alt="Empirical mode decomposition"  />
<!--
<p class="caption marginnote">-->Empirical mode decomposition<!--</p>-->
<!--</div>--></span></p>

Now organize the temperature and precipitation data to make plotting easier using functions from **tidyr**.
<label for="tufte-mn-" class="margin-toggle">&#8853;</label><input type="checkbox" id="tufte-mn-" class="margin-toggle"><span class="marginnote">Replace the variable names to make facet naming easier too.</span>

```r
trace.plot <- trace.dat %>% 
  gather(key, value, - Year) %>%
  separate(key, c('Variable', 'Region'), ',') %>%
  mutate(Variable = ifelse(
    Variable == 'tmp', 'Temperature (°C)', 'Precipitation (mm)'))

emd.res <- function(x) emd(x)$residue
trace.emd <- trace.dat %>%
  mutate_at(vars(-Year), emd.res) %>%
  gather(key, value, - Year) %>%
  separate(key, c('Variable', 'Region'), ',') %>%
  mutate(Variable = ifelse(
    Variable == 'tmp', 'Temperature (°C)', 'Precipitation (mm)'))
```

## Plotting
Plot everything with **ggplot2**.
<p><span class="marginnote shownote"><!--
<div class="figure">-->
<img src="proxymodelcomparison_files/figure-html/bwplot-1.png" alt="Black and white version"  />
<!--
<p class="caption marginnote">-->Black and white version<!--</p>-->
<!--</div>--></span></p>


```r
ggplot(data = trace.plot, aes(x = Year, y = value)) +
  facet_grid(Variable ~ Region, switch = 'y', scale = 'free_y') +
  geom_vline(xintercept = c(22, 19, 14, 10, 6), lty = 2) +
  geom_point(aes(color = Variable), alpha = .3) +
  geom_line(data = trace.emd, size = 1.2, color = "black", alpha = .8) +
  scale_x_reverse(breaks = seq(6,22,4)) +
  labs(x = '1,000 Years BP', y = '') +
  guides(color = "none") +
  theme_bw(base_size = 20) +
  theme(strip.background  = element_blank())
```

<div class="figure fullwidth">
<img src="proxymodelcomparison_files/figure-html/mainplot-1.png" alt=" "  />
<p class="caption marginnote shownote"> </p>
</div>

Calculate the detrended variances.

```r
emd.dat <- trace.dat %>%
  mutate_at(vars(-Year), emd.res)

(trace.dat - emd.dat) %>%
  select(-Year) %>%
  cbind(Year = trace.dat$Year, .) %>%
  mutate(Period = cut(Year, c(22, 19, 14, 10, 6))) %>%
  group_by(Period) %>%
  summarise_each(funs(var)) %>%
  select(-Year)
```

```
## # A tibble: 4 × 7
##    Period `prc,Southwest` `prc,North Central` `prc,Northeast`
##    <fctr>           <dbl>               <dbl>           <dbl>
## 1  (6,10]        463.0764            725.7196       1155.0184
## 2 (10,14]        809.8606           1028.3916       1051.4376
## 3 (14,19]       1453.6697           1422.2968       2251.1048
## 4 (19,22]        322.6892            398.9325        644.3691
## # ... with 3 more variables: `tmp,Southwest` <dbl>, `tmp,North
## #   Central` <dbl>, `tmp,Northeast` <dbl>
```

## Spatial patterns

Now let's import previously-downscaled ensemble equilibrium simulations of the Last Glacial Maximum and Mid Holocene, to estimate how the spatial patterns of climate variability have changed over time, and to test for consistency with the transient TraCE simulation.

First import the downscaled data.

```r
bbox <- extent(c(-10, 20, 35, 47))

lgm.prc <- brick('Downscaled/ensemble_prc_lgm.tif') %>%
  crop(bbox)
mh.prc <- brick('Downscaled/ensemble_prc_mh6k.tif') %>%
  crop(bbox)

lgm.tmp <- brick('Downscaled/ensemble_tmn_lgm.tif') %>%
  crop(bbox) %>%
  add(brick('Downscaled/ensemble_tmx_lgm.tif') %>% crop(bbox)) %>%
  divide_by(2)
mh.tmp <- brick('Downscaled/ensemble_tmn_mh6k.tif') %>%
  crop(bbox) %>%
  add(brick('Downscaled/ensemble_tmx_mh6k.tif') %>% crop(bbox)) %>%
  divide_by(2)
```

Calculate changes in seasonal precipitation and temperature.

```r
bySeason <- function(x, season, var){
  if(season == 'djf') {ids <- c(1,2,12)}
  if(season == 'jja') {ids <- c(6,7,8)}
  
  if(var == 'tmp') return(mean(x[[ids]]))
  if(var == 'prc') return(sum(x[[ids]]))
}

prc.change.map <- brick(c(
  (bySeason(mh.prc, 'djf', 'prc') - bySeason(lgm.prc, 'djf', 'prc')) * 100 / bySeason(lgm.prc, 'djf', 'prc'),
  (bySeason(mh.prc, 'jja', 'prc') - bySeason(lgm.prc, 'jja', 'prc')) * 100 / bySeason(lgm.prc, 'jja', 'prc')))

tmp.change.map <- brick(c(
  bySeason(mh.tmp, 'djf', 'tmp') - bySeason(lgm.tmp, 'djf', 'tmp'),
  bySeason(mh.tmp, 'jja', 'tmp') - bySeason(lgm.tmp, 'jja', 'tmp')))
```

Plot the results

```r
levelplot(prc.change.map, margin = F, names.attr = c('Winter', 'Summer'),
          main = 'Precipitation Change (%)\n LGM to Mid Holocene',
          par.settings = PuOrTheme(), 
          at = seq(-100,100,10))
```

<img src="proxymodelcomparison_files/figure-html/unnamed-chunk-7-1.png"  />

```r
levelplot(tmp.change.map, margin = F, names.attr = c('Winter', 'Summer'), 
          main = 'Temperature Change\n LGM to Mid Holocene',
          par.settings = BuRdTheme(),
          at = seq(-20,20,2))
```

<img src="proxymodelcomparison_files/figure-html/unnamed-chunk-7-2.png"  />

Now we can calculate changes in seasonality. For temperature, this is just the standard deviation of all 12 monthly averages. For precipitation, we will use the coefficient of variation.

```r
tmp.seasonality <- calc(mh.tmp, sd) - calc(lgm.tmp, sd)
prc.seasonality <- cv(mh.prc) - cv(lgm.prc)
```
Plot the results.

```r
levelplot(tmp.seasonality, margin = F, 
          main = 'Change in temperature seasonality (SD)\n LGM to Mid Holocene', 
          par.settings = BuRdTheme(), 
          at = seq(-4, 4, .4))
```

<img src="proxymodelcomparison_files/figure-html/unnamed-chunk-9-1.png"  />

```r
levelplot(prc.seasonality, margin = F, 
          main = 'Change in precipitation seasonality (CV)\n LGM to Mid Holocene', 
          par.settings = BuRdTheme(), 
          at = seq(-50, 50, 5))
```

<img src="proxymodelcomparison_files/figure-html/unnamed-chunk-9-2.png"  />

What about changes in spatial hetergeneity?

```r
wts <- matrix(c(0,0,1,0,0,0,1,1,1,0,1,1,1,1,1,0,1,1,1,0,0,0,1,0,0), nrow = 5)

tmp.hetero <- brick(c(
  bySeason(mh.tmp, 'djf', 'tmp') %>%
    focal(w = wts, sd, na.rm = T) %>%
    subtract(
      bySeason(lgm.tmp, 'djf', 'tmp') %>% 
        focal(w = wts, sd, na.rm = T)),
  bySeason(mh.tmp, 'jja', 'tmp') %>%
    focal(w = wts, sd, na.rm = T) %>%
    subtract(
      bySeason(lgm.tmp, 'jja', 'tmp') %>% 
        focal(w = wts, sd, na.rm = T)))) %>%
  mask(mh.tmp[[1]]) # clip buffer added by window
  
levelplot(tmp.hetero, margin = F, names.attr = c('Winter', 'Summer'), 
          main = 'Temperature heterogeneity (SD in 25km radius) change\n LGM to Mid Holocene',
          par.settings = BuRdTheme(), at = seq(-10, 10, 1))
```

<img src="proxymodelcomparison_files/figure-html/unnamed-chunk-10-1.png"  />
Same for precipitaiton.


```r
prc.hetero <- brick(c(
  bySeason(mh.prc, 'djf', 'prc') %>%
    focal(w = wts, sd, na.rm = T) %>%
    subtract(
      bySeason(lgm.prc, 'djf', 'prc') %>% 
        focal(w = wts, sd, na.rm = T)),
  bySeason(mh.prc, 'jja', 'prc') %>%
    focal(w = wts, sd, na.rm = T) %>%
    subtract(
      bySeason(lgm.prc, 'jja', 'prc') %>% 
        focal(w = wts, sd, na.rm = T)))) %>%
  mask(mh.prc[[1]]) # clip buffer added by window
  
  
levelplot(prc.hetero, margin = F, names.attr = c('Winter', 'Summer'), 
          main = 'Precipitation heterogeneity (SD in 25km radius) change\n LGM to Mid Holocene',
          par.settings = BuRdTheme(), at = seq(-500, 500, 50))
```

<img src="proxymodelcomparison_files/figure-html/unnamed-chunk-11-1.png"  />

<label for="tufte-mn-" class="margin-toggle">&#8853;</label><input type="checkbox" id="tufte-mn-" class="margin-toggle"><span class="marginnote">Compare these with raw gcm outputs to check the added value of SDM</span>










