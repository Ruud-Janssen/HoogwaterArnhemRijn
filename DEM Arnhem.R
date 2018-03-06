# install.packages("raster")
# install.packages("rgdal")
# install.packages("ggplot2")
# install.packages("ggmap")
# install.packages("gdistance")
# install.packages("sf")
# install.packages("mapview")
# install.packages("XML")
# install.packages("gdalUtils")
# install.packages("stringr")
#library(rgdal)
library(gdalUtils)
library(XML)
library(mapview)
library(sf)
library(gdistance)
library(raster)
library(rgdal)
library(ggplot2)
library(ggmap)
library(dplyr)
library(leaflet)
library(sp)
library(leaflet.extras)
library(stringr)






# Load 5m Dem files:
# - m5_40bz1.tif dtm arnhem east 
# - m5_40bz2.tif dtm arnhem west



#We should use the DTM (Digital Terrain Model)
#https://www.pdok.nl/nl/ahn3-downloads
DEM_dtm_west <- raster("input/m5_40az2.tif")
DEM_dtm_east <- raster("input/m5_40bz1.tif")
#DEM_dtm_west <- raster("input/m_40az2.tif") #high resolution, for sure an improvement
#DEM_dtm_east <- raster("input/m_40bz1.tif") #high resolution, for sure an improvement
#The North side of Arnhem (above the river Rijn) is not yet available in Ahn3 (the newest DEM information)

#I like to see values when printing the raster layer, including the min en max
#Therefor these weird calls 
values(DEM_dtm_west) <- values(DEM_dtm_west)
values(DEM_dtm_east) <- values(DEM_dtm_east)

#Combine the SouthEast and SouthWest of Arnhem (first the extent)
em = merge(extent(DEM_dtm_west),extent(DEM_dtm_east))

# create a raster with that extent, and the number of rows and colums
emRaster <- raster(em, crs = DEM_dtm_east@crs, resolution = c(5, 5))

DEM <- merge(emRaster, DEM_dtm_east)
DEM <- merge(DEM, DEM_dtm_west)

#plot(DEM_dtm_east)
#plot(DEM_dtm_west)
#plot(DEM)


#Lets add a rougher DEM to get values for all values
#CIAT-CSI  SRTM website, http://srtm.csi.cgiar.org
DEM_srtm90 <- raster("input/srtm_38_02.tif")
#plot(DEM_srtm90)
DEM_srtm90 <- crop(DEM_srtm90, extent(5.70, 6.09, 51.82, 52.02))
DEM_srtm90 <- projectRaster(DEM_srtm90, DEM, res = c(5, 5),  crs = DEM@crs)
#plot(DEM_srtm90)

#Create a complete DEM without gaps
DEM_full <- merge(DEM, DEM_srtm90)

# par(mfrow=c(1,3))
# plot(DEM)
# plot(DEM_srtm90)
# plot(DEM_full)


#function return the flooded area depending on the waterheight
calcFloodedArea <- function(waterLevel, DEM_Arnhem, src) {
  DEM_Arnhem <- DEM_Arnhem > waterLevel
  DEM_Arnhem[DEM_Arnhem == 1] <- NA
  t4 <- transition(DEM_Arnhem, transitionFunction = function(x){1}, directions = 4)
  t4.acc <- accCost(t4,src)
  t4.acc[t4.acc == Inf] <- NA
  t4.acc[t4.acc > 0] <- 1
  return (t4.acc)
} 

#src pinpoints points in the Rijn to start the cumulative cost function
#For testing the points position spot on in the Rijn, plot the points
#plot(DEM_full)
#points(c(194350, 185900), c(439800, 442900))
src <- SpatialPoints(cbind(c(194350, 185900), c(439800, 442900)))







#Create waterheight from 9.00m above NAP to 14.00m with steps of 0.2m
waterHeight <- seq(from = 9, to = 14, by = 0.2)
rm(DEM_dtm_east, DEM_dtm_west, DEM_srtm90)

#Calculate flooded areas
floodedAreas <- sapply(waterHeight, FUN = calcFloodedArea, DEM_Arnhem = DEM_full, src = src)



#A quick visualization of flooded areas
#remove first entry for a 5 by 5 image plot
floodedAreas[[1]] <- NULL
waterHeight <- waterHeight[2:26]

par(mfrow=c(5,5))
names(floodedAreas) <- waterHeight
n <- names(floodedAreas)
lapply(setNames(n, n), function(nameindex) {
    plot(
      floodedAreas[[nameindex]],
      axes = FALSE,
      legend = FALSE,
      xlab = "",
      ylab = "",
      main = nameindex
    )
  }
)
par(mfrow=c(1,1))




# Below a single raster waterheight example to convert raster to vector and visualize this
# waterHeight[20]
# plot(floodedAreas[[20]])
floodedArea12.8 <- floodedAreas[[20]]
floodedArea12.8[floodedArea12.8 == 0] <- NA
# still doesn't work unfortunately
# pol <- rasterToPolygons(floodedArea12.8, fun=function(x){x == 1})
# pol is 2Gb+ (or simple doesn't finish)...



#First export raster
writeRaster(floodedArea12.8, "output/RiverRijnArnhem_12_8m.tif", format = "GTiff", overwrite = T)


#convert the TIF to a shapefile using python!
source("gdal_polygonizeR.R")

gdal_polygonizeR(
  "\\output\\RiverRijnArnhem_12_8m.tif", 
  outshape = "\\output\\RiverRijnArnhem_12_8m.shp", 
  pypath = file.path("C:\\gdal-2.2.3\\swig\\python\\scripts\\gdal_polygonize.py"),
  overwrite = T,
  directory = normalizePath('.')
)

#Lets use the new sf package!
RiverRijnArnhem_12.8m <- st_read("output\\RiverRijnArnhem_12_8m.shp")


#Not sure why it lost its reference system when converting with gdal_polygonize. As it is within the .tif
RiverRijnArnhem_12.8m <- st_transform(RiverRijnArnhem_12.8m, 28992)
mapview(RiverRijnArnhem_12.8m)







#The results are not particularly satisfying above the Rijn. The 90m SRTM is not useful for this data problem
#Can we use AHN2 data, the previous DEM data. There is no direct download, but they have a wms and a wcs service
leaflet() %>% 
  addTiles() %>% 
  setView(lat = 51.97, lng = 5.91, zoom = 13) %>%
  addWMSTiles(
    "http://geodata.nationaalgeoregister.nl/ahn2/wms?SERVICE=WMS&",
    layers = "ahn2_05m_int", #or ahn2_5m
    options = WMSTileOptions(format = "image/tiff", transparent = TRUE),
    attribution = ""
  ) %>%
 #legend definitely needs customization, but that is for another time   
 addWMSLegend(uri = paste0(
   "http://geodata.nationaalgeoregister.nl/ahn2/wms?SERVICE=WMS&",
   "REQUEST=GetLegendGraphic&VERSION=1.1.1", #Or 1.3.0
   "&FORMAT=image/tiff&LAYER=ahn2_05m_int"
   )
 )

#Or grab the image straight from the URL
#http://geodata.nationaalgeoregister.nl/ahn2/wms?SERVICE=WMS&REQUEST=GetMap&VERSION=1.1.1&LAYERS=ahn2_05m_int&STYLES=&FORMAT=image%2Ftiff&TRANSPARENT=true&HEIGHT=2048&WIDTH=2048&SRS=EPSG%3A3857&BBOX=650631.9847634202,6785162.126818524,665307.8941941741,6799838.036249279
#http://geodata.nationaalgeoregister.nl/ahn2/wms?SERVICE=WMS&REQUEST=GetMap&VERSION=1.3.0&LAYERS=ahn2_05m_int&STYLES=&FORMAT=image%2Ftiff&TRANSPARENT=false&HEIGHT=2048&WIDTH=2048&CRS=EPSG%3A3857&BBOX=650631,6785162,665307,6799838



# However can we use the wcs (similar as wfs for features, but now for raster data)

# wcs = "https://geodata.nationaalgeoregister.nl/ahn2/wcs?"
# wcsGdal <- newXMLNode("WCS_GDAL")
# wcsGdal.serviceURL <- newXMLNode("ServiceURL", wcs, parent=wcsGdal)
# wcsGdal.layer <- newXMLNode("CoverageName", "ahn2__ahn2_5m", parent=wcsGdal)
# wcsGdal.xml = "output\\ahn_2_5m.xml"
# saveXML(wcsGdal, file = wcsGdal.xml)
# 
# #unfortunately I get in both cases of using gdalinfo
# gdalinfo(xml.out, verbose = T)
# system("OSGeo4W", input = paste0("gdalinfo ",normalizePath('.'),"\\", xml.out))

# Checking gdal_installation...
# Scanning for GDAL installations...
# Checking the gdalUtils_gdalPath option...
# GDAL version 2.2.3
# GDAL command being used: "C:\OSGeo4W64\bin\gdalinfo.exe" "input\ahn_2_5m.xml"
# ERROR 1: SSL certificate problem: self signed certificate in certificate chaingdalinfo failed - unable to open 'input\ahn_2_5m.xml'.
# [1] "ERROR 1: SSL certificate problem: self signed certificate in certificate chain"
# [2] "gdalinfo failed - unable to open 'input\\ahn_2_5m.xml'."                       
# attr(,"status")
# [1] 1
# Warning message:
#   running command '"C:\OSGeo4W64\bin\gdalinfo.exe" "input\ahn_2_5m.xml"' had status 1 

# Normally when this works you can grab the tif with gdal_translate
# Detailed example: http://gsif.isric.org/doku.php/wiki:tutorial_soilgrids



# This is a setback, I am not sure how to solve this
# The following Dutch site explains some things and why it works in QGIS (which is in my case similar)
# https://www.pdok.nl/sites/default/files/forum_attachments/user_484/ogc_wms_en_wfs_beveiligd_met_pki-certificaten_v0.1.pdf

# Now we downloaded the file via QGIS manually (make sure you set the map extent), which is fine for now as it is a single file
DEM_ahn2 <- raster("input/wcs_Ahn2_5m.tif") #AHN2, DEM not dtm, downloaded using wcs in qgis, with a certain mapextent 
values(DEM_ahn2) <- values(DEM_ahn2)
DEM_Ahn2_reproj <- projectRaster(DEM_ahn2, DEM, res = c(5, 5),  crs = DEM@crs)
#plot(DEM_Ahn2_reproj)
#lets update some weird extremely negative values
DEM_Ahn2_reproj[DEM_Ahn2_reproj < 0] <- 0
#plot(DEM_Ahn2_reproj)
#Lets put the NA values on 0 which is fine for out problem
DEM_Ahn2_reproj[is.na(DEM_Ahn2_reproj)] <- 0
#plot(DEM_Ahn2_reproj)

#merge the DEM's, first the AHN3 data and if not available the AHN2 data
DEM2 <- merge(DEM, DEM_Ahn2_reproj)
#plot(DEM2)

#lets fix the bridges (elevation = 0m), so the water doesn't go inwards at the spot of the bridges
writeRaster(DEM2, "output\\DEM_merged.tif")

#edit manually bridge height at the shores in QGIS (and a few spots where buildings were too close on the north side of the Rijn in the center)
#using QGIS plugin serval
DEM2 <- raster("output\\DEM_merged_bridgesAdapted.tif")
#plot(DEM2)

waterHeight <- seq(from = 9, to = 14, by = 0.2)
floodedAreas <- sapply(waterHeight, FUN = calcFloodedArea, DEM_Arnhem = DEM2, src = src)
names(floodedAreas) <- waterHeight
n <- names(floodedAreas)

lapply(setNames(n, n), function(nameindex) {
  writeRaster(
    floodedAreas[[nameindex]],
    paste0("output\\RiverRijnArnhemHoogwater", ifelse(as.numeric(nameindex) %% 1 != 0, paste0(str_replace(nameindex, "[.]", "-"),"0"), paste0(nameindex,"-00") ), ".tif")
  )
})



#visualize data
#remove first entry for a 5 by 5 image plot
floodedAreas[[1]] <- NULL
waterHeight <- waterHeight[2:26]

par(mfrow=c(5,5))
n <- names(floodedAreas)
lapply(setNames(n, n), function(nameindex) {
  plot(
    floodedAreas[[nameindex]],
    axes = FALSE,
    legend = FALSE,
    xlab = "",
    ylab = "",
    main = nameindex
  )
}
)
par(mfrow=c(1,1))
