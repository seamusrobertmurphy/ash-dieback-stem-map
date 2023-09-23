---
title: "Ash-Dieback Stem Detection"
author: "Vertinetik-SM"
date: "`r Sys.Date()`"
output: 
  urlcolor: blue
  html_document:
    toc: TRUE
    toc_depth: 5
    number_sections: FALSE
    df_print: tibble
---


```{r setup, echo=FALSE, message=FALSE,warning=FALSE, error=FALSE}
library(lidR)
library(ForestTools)
library(rgl)
library(pandocfilters)
library(rmarkdown)
library(formatR)
library(gitignore)
library(tinytex)
library(knitr)
library(raster)
library(webdriver)
library(webshot)
library(webshot2)
library(terra)
library(matlab)
#webshot::install_phantomjs(force = TRUE)
knit_hooks$set(webgl = hook_webgl)
knit_hooks$set(rgl.static = hook_rgl)
knitr::opts_chunk$set(echo = TRUE, warning=FALSE, error=FALSE, message = FALSE)
set.seed(23)
```

## Classify Normalise Rasterize

```{r, rgl.static=TRUE, cache=TRUE, eval=FALSE, echo=TRUE}
#Visualize
defaultCRS<-CRS("+init=EPSG:3857")
las_tile = readLAS("/media/seamus/USB1/Shavington/clouda340d379a59ddadf.las", select = 'xyzcr', filter = '-drop_class 19')
#Classify ground & noise
las_tile_csf = classify_ground(las_tile, csf(sloop_smooth=TRUE, 0.5, 1))
#las_tile_csf_so = classify_noise(las_tile_csf, sor(k=10, m=3))
#Clean noise
las_tile_csf_norm = normalize_height(las_tile_csf, knnidw())
#Clean overlaps
#las_tile_chm_so_norm_clean = filter_duplicates(las_tile_csf_so_norm)
#Rasterize
las_tile_chm = grid_canopy(las_tile_csf_norm, 1, dsmtin(8))
#plot(las_tile_chm, col = height.colors(50))

writeRaster(las_tile_chm, filename = "/media/seamus/USB1/Shavington/lead_htop_raster.tif", overwrite=TRUE)
plot(las_tile_chm)
```


## Variable Window Function 

```{r, echo=TRUE, eval=FALSE}

kernel <- matrix(1,3,3)
wf_plowright<-function(x){ 
  a=0.05
  b=0.6 
  y<-a*x+b 
  return(y)}
heights <- seq(0,40,0.5)
window_plowright <- wf_plowright(heights)
plot(heights, window_plowright, type = "l", ylim = c(0,12), xlab="point elevation (m)", ylab="window diameter (m)", main='Plowright, 2018')

las_tile_chm_smooth = focal(las_tile_chm, w = kernel, fun = median, na.rm = TRUE) 
ttops_1.5mfloor_plowright = ForestTools::vwf(las_tile_chm_smooth, wf_plowright, 1.5)
#ttops_1.5mfloor_plowright = ForestTools::vwf(CHM = las_tile_chm_smooth, winFun = wf_plowright, minHeight = 2)
ttops_1.5mfloor_plowright_sp = as_Spatial(ttops_1.5mfloor_plowright)
writeOGR(ttops_1.5mfloor_plowright_sp, "/media/seamus/USB1/Shavington", "ttops_1.5mfloor_plowright_sp", driver = "ESRI Shapefile") 
```

```{r, echo=FALSE}
stem_count_sf = st_as_sf("/media/seamus/USB1/Shavington/ttops_1.5mfloor_plowright.shp")

las_tile_chm = raster::raster("/media/seamus/USB1/Shavington/lead_htop_raster.tif")
#stem_count_ha_sf = st_as_sf(ttops_1.5mfloor_plowright)
plot(st_geometry(ttops_1.5mfloor_plowright["treeID"]), cex = 0.2, pch="+", col = 'red', lwd=1, alpha=1, add=TRUE) 
plot(las_tile_chm)
```
