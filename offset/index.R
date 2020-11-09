## offsets

if (!requireNamespace("QPAD")) {
  if (!requireNamespace("remotes"))
    install.packages("remotes")
  remotes::install_github("psolymos/QPAD")
}
if (!requireNamespace("sp"))
  install.packages("sp")
if (!requireNamespace("maptools"))
  install.packages("maptools")
if (!requireNamespace("raster"))
  install.packages("raster")
if (!requireNamespace("intrval"))
  install.packages("intrval")

library(QPAD)
library(maptools)
library(intrval)
library(raster)

load_BAM_QPAD(version = 3)
if (getBAMversion() != "3")
  stop("This script requires BAM version 3")

source("functions.R")

rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
rd1 <- raster("./data/seedgrow.tif")
crs <- proj4string(rtree)

spp <- "OVEN"
## https://en.wikipedia.org/wiki/ISO_8601
dt <- "2019-06-07" # ISO 8601 in YYYY-MM-DD (0-padded)
tm <- "05:20" # ISO 8601 in hh:mm (24 hr clock, 0-padded)
lon <- -113.4938 # longitude WGS84 (EPSG: 4326)
lat <- 53.5461 # latitude WGS84 (EPSG: 4326)
dur <- 10 # mins
dis <- 100 # meters


x <- make_x(dt, tm, lon, lat, dur, dis)
o <- make_off(spp, x)
x
o


library(openxlsx)

SPP <- read.xlsx("~/Desktop/OFFSETS.xlsx", 1)
SPP <- as.character(SPP$Species)

X <- read.xlsx("~/Desktop/OFFSETS.xlsx", 2)
## need to deal with NAs
X$RECORD_TIME[is.na(X$RECORD_TIME)] <- mean(X$RECORD_TIME, na.rm=TRUE)

X$DATE <- convertToDate(X$sDATE)
tmp <- round(X$RECORD_TIME * 24 * 60) # minutes
HR <- as.character(tmp %/% 60)
HR <- ifelse(nchar(HR) < 2, paste0("0", HR), HR)
MIN <- as.character(tmp %% 60)
MIN <- ifelse(nchar(MIN) < 2, paste0("0", MIN), MIN)
X$TIME <- paste0(HR, ":", MIN)


x <- make_x(dt=X$DATE, tm=X$TIME,
  lon=X$Longitude, lat=X$Latitude,
  dur=3, dis=Inf, key=X$location)

OFF <- matrix(0, nrow(x), length(SPP))
rownames(OFF) <- x$key
colnames(OFF) <- SPP

for (spp in SPP) {
  o <- make_off(spp, x)
  OFF[,spp] <- o$offset
}

OUT <- data.frame(X, x, Offset=OFF)
write.csv(OUT, row.names=FALSE, file="~/Desktop/OFFSETS-out.csv")

