#' @title gdal_polygonizeR
#' @description R Wrapper for the gdal_polygonize.py python script (http://www.gdal.org/gdal_polygonize.html)
#' This utility creates vector polygons for all connected regions of pixels in the raster sharing
#' a common pixel value. Each polygon is created with an attribute indicating the pixel value of
#' that polygon. Can be userful for example to create fishnet polygons (see "create_fishnet")
#'
#' @param x          `character`filename of a raster file, or "R" raster object
#' @param outshape   `character`filename of the desired output polygon shapefile
#' @param gdalformat  defaults to ESRI - don't change
#' @param pypath     `character` path of python  - if `NULL` (the default) the script tries to 
#' automatically retrieve it
#' @param readpoly   `logical` If TRUE sends back the shapefile as a SpataialPolygonsDataFrame to R
#'     (defaults to FALSE)
#' @param quiet      `logical` if TRUE, limit the messages (defaults to TRUE)
#' @param overwrite  `logical` If TRUE overwrite shapefile if already existing (defaults to FALSE)
#'
#' @return NULL, or a SpataialPolygonsDataFrame (if readpoly = TRUE)
#' @export
#' @importFrom raster writeRaster
#' @importFrom methods is
#' @author Original code by Jonh Baumgartner
#' [https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/](), with
#' slight modifications by L.Busetto
#'
#'
#' Modifications by R.Janssen
#' It now only supports rasterfiles (not raster objects) as input
#' Readpoly option is disappeared
#' Directory is added as param
#' Quiet is removed
#' Now it calls first osgeo4w and then the python raster to vector convert command
#' Becareful to check whether the shapefile has the right coordinate system afterwards!

gdal_polygonizeR <- function(inputraster, 
                             outshape,
                             gdalformat = 'ESRI Shapefile',
                             pypath     = NULL,
                             pyexe      = 'python',
                             overwrite  = FALSE,
                             directory  = "") {
  
  inputraster <- paste0(directory, inputraster)
  outshape <- paste0(directory, outshape)
  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
  
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep = '.'))
    if (any(f.exists)) {
      if (overwrite == FALSE) {
        stop(sprintf('File already exists: %s',
                     toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                    sep = '.')[f.exists])), call.=FALSE)
      } else (
        unlink(paste(outshape, c('shp', 'shx', 'dbf'), sep = '.'))
      )
    }
    system("OSGeo4W", input = paste0(sprintf('%1$s "%2$s" "%3$s" -f "%4$s" "%5$s.shp"', pyexe, pypath, inputraster, gdalformat, outshape), " -fieldname id"))
  }
}