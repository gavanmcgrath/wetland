
#' Polygonize a raster
#'
#' @param rast a raster
#'
#' @return a spatialPolygonsDataFrame
#' 
#' @details This function takes a raster object and returns a spatialPolygonsDataFrame
#' it creates unique, merged polygons for each spatially contigous set of unique values. 
#' Best not to run on a raster with continuous values.
#'
#' @examples
#' r1 <- raster::raster(nrows=30, ncols=30, xmn=0, xmx=30,ymn=0, ymx=30)
#' r1[] <- sample(c(1,2,0),length(r1), replace = TRUE)  
#' p1 <- polygonizer2(r1)
#' par(mfrow=c(1,2))
#' plot(p1,col=p1@data$layer)
#' plot(r1)
polygonizer2 <- function(rast){
  r.to.poly <- sf::as_Spatial(
    sf::st_as_sf(stars::st_as_stars(rast), 
                 as_points = FALSE, 
                 merge = TRUE)
  ) 
}

#' Make spatially contigous patches of like-values from a discrete raster
#'
#' @param x a binary Raster layer (0 or NA for background, and 1 for areas to be clumped)
#' @param distance the neighbourhood distance. Patches that occur within this distance of 
#'  one another will be clumped. This should be in the units of the CRS given.
#'  in p4s. If this is zero, the function will identify patches defined by 
#'   8-connectivity (i.e. queen's case).
#' @param p4s the CRS (see sp package) specifying an equal-area projection appropriate for the region. 
#' This will be used when calculating inter-patch distances. The returned objects will be in the
#'  original coordinate system.
#'  
#' @param givedist should the distance matrix be returned? (logical). Distances are in the
#' units of p4s, and are shortest cartesian distances between patches.
#' 
#' @return a list containing the new raster, a SpatialPolygonsDataFrame and optionally a matrix 
#' of distances.
#' 
#' @details #From  https://gist.github.com/johnbaums/
#'
patchify <- function(x, distance, p4s, givedist=TRUE) {
  
  require(raster)
  require(rgeos)
  require(SDMTools)
  if(!is(x, 'Raster')) x <- raster(x)
  if(!is(p4s, 'CRS')) p4s <- CRS(p4s)
  if(is.na(proj4string(x))) stop(substitute(x), ' lacks a CRS.')
  x[x == 0] <- NA
  cc <- ConnCompLabel(x)
  p <- polygonizer2(cc)
  p <- p[p@data$layer == 1,]
  suppressWarnings(proj4string(p) <- proj4string(x))
  pproj <- spTransform(p, p4s)
  bproj <- gBuffer(pproj, width=distance)
  pproj$patch <- over(pproj, disaggregate(bproj))
  b <- spTransform(pproj, CRS(proj4string(x)))
  pout <- aggregate(b, by='patch')
  pout$patch <- as.factor(pout$patch)
  rout <- rasterize(pout, x)
  out <- list(patchrast=rout, patchpoly=pout)
  if(isTRUE(givedist)) {
    poutproj <- spTransform(pout, p4s)
    d <- gDistance(poutproj, poutproj, byid=TRUE)
    out <- c(out, list(distance=d))
  } 
  return(out)
}


#' Download a map image from the ArcGIS REST API
#' From https://github.com/zappingseb/rayshaderanimate
#'
#' @param bbox bounding box coordinates (list of 2 points with long/lat values)
#' @param map_type map type to download - options are World_Street_Map, World_Imagery, World_Topo_Map
#' @param file file path to save to. Default is NULL, which will create a temp file.
#' @param width image width (pixels)
#' @param height image height (pixels)
#' @param sr_bbox Spatial Reference code for bounding box
#' 
#' @details This function uses the ArcGIS REST API, specifically the 
#' "Execute Web Map Task" task. You can find links below to a web UI for this
#' rest endpoint and API documentation.
#' 
#' Web UI: https://utility.arcgisonline.com/arcgis/rest/services/Utilities/PrintingTools/GPServer/Export%20Web%20Map%20Task/execute
#' API docs: https://developers.arcgis.com/rest/services-reference/export-web-map-task.htm
#'
#' @return file path for the downloaded .png map image
#'
#' @examples
#' bbox <- list(
#'   p1 = list(long = -122.522, lat = 37.707),
#'   p2 = list(long = -122.354, lat = 37.84)
#' )
#' image_size <- define_image_size(bbox, 600)
#' overlay_file <- get_arcgis_map_image(bbox, width = image_size$width,
#'                                      height = image_size$height)
#' 
get_arcgis_map_image <- function(bbox, map_type = "World_Street_Map", file = NULL, 
                                 width = 400, height = 400, sr_bbox = 4326) {
  require(httr)
  require(glue) 
  require(jsonlite)
  
  url <- parse_url("https://utility.arcgisonline.com/arcgis/rest/services/Utilities/PrintingTools/GPServer/Export%20Web%20Map%20Task/execute")
  
  # define JSON query parameter
  web_map_param <- list(
    baseMap = list(
      baseMapLayers = list(
        list(url = jsonlite::unbox(glue("https://services.arcgisonline.com/ArcGIS/rest/services/{map_type}/MapServer",
                                        map_type = map_type)))
      )
    ),
    exportOptions = list(
      outputSize = c(width, height)
    ),
    mapOptions = list(
      extent = list(
        spatialReference = list(wkid = jsonlite::unbox(sr_bbox)),
        xmax = jsonlite::unbox(max(bbox$p1$long, bbox$p2$long)),
        xmin = jsonlite::unbox(min(bbox$p1$long, bbox$p2$long)),
        ymax = jsonlite::unbox(max(bbox$p1$lat, bbox$p2$lat)),
        ymin = jsonlite::unbox(min(bbox$p1$lat, bbox$p2$lat))
      )
    )
  )
  
  res <- GET(
    url, 
    query = list(
      f = "json",
      Format = "PNG32",
      Layout_Template = "MAP_ONLY",
      Web_Map_as_JSON = jsonlite::toJSON(web_map_param))
  )
  
  if (status_code(res) == 200) {
    body <- content(res, type = "application/json")
    message(jsonlite::toJSON(body, auto_unbox = TRUE, pretty = TRUE))
    if (is.null(file)) 
      file <- tempfile("overlay_img", fileext = ".png")
    
    img_res <- GET(body$results[[1]]$value$url)
    img_bin <- content(img_res, "raw")
    writeBin(img_bin, file)
    message(paste("image saved to file:", file))
  } else {
    message(res)
  }
  invisible(file)
}

# From https://www.r-bloggers.com/a-3d-tour-over-lake-geneva-with-rayshader/
#' Make a horizontal colour bar
#'
#' @param col a character vector of colours
#' @param data.min numeric value, minimum value on colour bar
#' @param data.max numeric value, maximum value on colour bar
#' @param n.tick numeric value, number of tick marks on colour bar
#'
#' @return a plot object of the colour bar
#'
#' @examples
#'  zs <- seq(-0.1, 2.4, by = 0.1)
#'  i <- 10
#'  pal <- c(rep("lightblue",i),rep("lightgrey",length(zs) - i + 1))
#'  col_bar_h(pal,-0.1, 2.4, n.tick = 6)
#'
col_bar_h <- function(col, data.min, data.max, n.tick = 5){
  plot(c(data.min, data.max), c(0, 10), type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  tk <- pretty(c(data.min, data.max), n.tick + 2)
  tk <- tk[c(-which(tk < data.min), -which(tk > data.max))]
  axis(1, at = tk, col = NA, col.ticks = 1)
  inc <- (data.max - data.min)/ length(col)
  for(i in seq(1, length(col))){
    rect(data.min + inc * (i- 1), 0, data.min + inc * i, 10, col = col[i], border = NA)
  }
}


#Workflow
library(raster)
library(rgdal)
library(rayshader)
library(sp)
library(dplyr)

#Read in bounding polygon and buffer by 200 m
aoi <- readOGR(getwd(), layer = "aoi")
aoib <- buffer(aoi, 200)

#Read in raster elevation data and crop to the aoi
dem <- raster("GroundSurface.tif")
dem <- crop(dem, aoib)

#Make a list of logical rasters containing the water mask
list.of.dems <- list()
zs <- seq(-0.05,2.4,by=0.05)
for (i in 1:length(zs)){
  print(zs[i])
  dem2 <- dem < zs[i]
  clst.dem <- patchify(dem2, distance = 10, p4s = CRS(proj4string(dem)), givedist=TRUE)
  list.of.dems[[i]] <- clst.dem$patchrast
}

#Create the overlay map
overlay_file <- "./sf-map.png"
bbox1 <- extent(dem)
bbox <- list(
  p1 = list(long =  bbox1[1], lat =  bbox1[3]),
  p2 = list(long =  bbox1[2], lat =  bbox1[4])
)
overlaymap <- get_arcgis_map_image(
  bbox = bbox,
  map_type = "World_Imagery", 
  file = overlay_file,
  width = 1339 , #dimensions of the hillshade object
  height = 1056, #dimensions of the hillshade object
  sr_bbox = 32750) #Specify here the epsg code of the elevation rater

overlay_img <- png::readPNG(overlay_file)
overlay_img2 <- array(NA, dim=dim(overlay_img2)[c(2,1,3)])
overlay_img2[,,1] <- matlab::rot90(overlay_img[,,1],3)
overlay_img2[,,2] <- matlab::rot90(overlay_img[,,2],3)
overlay_img2[,,3] <- matlab::rot90(overlay_img[,,3],3)
overlay_img2[,,4] <- matlab::rot90(overlay_img[,,4],3)


#Loop through list.of.dems rendering a new image for each water map
export_folder <- getwd()
elmat <- as.matrix(dem)
elmat <- elmat[dim(elmat)[1]:1,]
hillshade <- elmat %>%
  sphere_shade(texture = "desert", zscale = 0.05, sunangle = 15)
for (i in 1:length(zs)){
  
  watermap <- as.matrix(list.of.dems[[i]])
  
  #Make colour bar
  img_scale <- magick::image_graph(width = 400, height = 60)
  par(mar = c(3, 0, 0, 0))
  pal <- c(rep("lightblue",i),rep("lightgrey",length(zs) - i + 1))
  col_bar_h(pal,-0.1, 2.4, n.tick = 6)
  dev.off()
  
  #3d shaded relief + water + overlay
  hillshade %>%
    add_overlay(overlay_img2,alphalayer=0.7) %>%
    add_water(watermap, color = "desert") %>%
    add_shadow(ray_shade(elmat, zscale = 0.05, multicore = TRUE, sunangle = 15), 0.5) %>%
    add_shadow(ambient_shade(elmat, multicore = TRUE), 0) %>%
    plot_3d(elmat, zscale = 0.05, fov = 0, theta = 25, zoom = 0.75, phi = 45, 
            windowsize = c(1000, 800),wateralpha = 0.5, water = TRUE, waterlinecolor = "white")
  Sys.sleep(5)
  file_img_i <- file.path(export_folder, paste0("render_", sprintf("%04d", i), ".png"))
  render_snapshot(filename = file_img_i, clear = TRUE)
  
  #Merge the 3d image with the scale bar
  map <- magick::image_read(file_img_i)
  img <- c(img_scale,map)
  img2 <- image_composite(map, img_scale, offset = "+30+50")
  img3 <- image_draw(img2, xlim = c(0, image_info(img2)[["width"]]), ylim = c(0, image_info(img2)[["height"]]))
  text(230, image_info(img2)[["height"]] - 20, paste("River Water Elevation (mAHD)"), cex = 1.2)
  dev.off()
  magick::image_write(img3, path = file_img_i, format = "png")
  
  #for smooth looping animations save the same image with a mirrored number
  file_img_i <- file.path(export_folder, paste0("render_", sprintf("%04d", 2*length(zs) - i + 1), ".png"))
  magick::image_write(img3, path = file_img_i, format = "png")
}


#I paste this manually into my wondows command window after changing the directory
# to where the images are saved.
# It requires ffmpeg (https://ffmpeg.org/)to be installed and the executables to be in your 
# PATH environment variable 
#'ffmpeg -y -r 5 -i render_%04d.png -c:v libx264 -vf fps=5,format=yuv420p temp_movie2.mp4')



