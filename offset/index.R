## offsets, works only in AB
## to extend beyon AP, need time zone tweaks + TREE/LCC layers

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

load_BAM_QPAD(version = 3)
if (getBAMversion() != "3")
  stop("This script requires BAM version 3")

rlcc <- raster("./data/lcc.tif")
rtree <- raster("./data/tree.tif")
rtz <- raster("./data/utcoffset.tif")
crs <- proj4string(rtree)

spp <- "OVEN"
## https://en.wikipedia.org/wiki/ISO_8601
dt <- "2019-06-07" # ISO 8601 in YYYY-MM-DD (0-padded)
tm <- "05:20" # ISO 8601 in hh:mm (24 hr clock, 0-padded)
lon <- -113.4938 # longitude WGS84 (EPSG: 4326)
lat <- 53.5461 # latitude WGS84 (EPSG: 4326)
dur <- 10 # mins
dis <- 100 # meters

## checking lengths
if (length(spp) > 1L)
  stop("spp argument must be length 1")
nn <- c(dt=length(dt), tm=length(tm), lon=length(lon), lat=length(lat), dur=length(dur), dis=length(dis))
n1 <- nn[nn == 1L]
n2 <- nn[nn > 1L]
if (!all(n2 == n2[1L]))
  stop("input lengths must be equal or 1")
n <- unname(if (length(n2)) n2[1L] else n1[1L])
if (length(dt) == 1L)
  dt <- rep(dt, n)
if (length(tm) == 1L)
  tm <- rep(tm, n)
if (length(lon) == 1L)
  lon <- rep(lon, n)
if (length(lat) == 1L)
  lat <- rep(lat, n)
if (length(dur) == 1L)
  dur <- rep(dur, n)
if (length(dis) == 1L)
  dis <- rep(dis, n)

## types
spp <- as.character(spp)
lat <- as.numeric(lat)
lon <- as.numeric(lon)
dur <- as.numeric(dur)
dis <- as.numeric(dis)
## parse date+time into POSIXlt
dt <- as.character(dt)
tm <- as.character(tm)
dtm <- strptime(paste0(dt, " ", tm, ":00"),
  format="%Y-%m-%d %H:%M:%S", tz="America/Edmonton")
day <- as.integer(dtm$yday)
hour <- as.numeric(round(dtm$hour + dtm$min/60, 2))

## checks
if (!(spp %in% getBAMspecieslist()))
  stop(sprintf("Species %s has no QPAD estimate", spp))
checkfun <- function(x, name="", range=c(-Inf, Inf)) {
  if (any(x[!is.na(x)] %)(% range))
    stop(sprintf("Parameter %s is out of range [%.0f, %.0f]", name, range[1], range[2]))
  invisible(NULL)
}
## BCR 4:14 included
## crs: WGS84 (EPSG: 4326)
## "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
#         min       max
#x -163.89547 -52.66936
#y   39.66214  68.98741
checkfun(lon, "lon", c(-163.89547, -52.66936))
checkfun(lat, "lat", c(39.66214, 68.98741))
checkfun(day, "day", c(0, 360))
checkfun(hour, "hour", c(0, 24))
checkfun(dur, "dur", c(0, Inf))
checkfun(dis, "dis", c(0, Inf))
if (any(is.infinite(lon)))
  stop("Parameter lon must be finite")
if (any(is.infinite(lat)))
  stop("Parameter lat must be finite")

## intersect here
xy <- data.frame(x=lon, y=lat)
coordinates(xy) <- ~ x + y
proj4string(xy) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
xy <- spTransform(xy, crs)

## LCC4 and LCC2
vlcc <- extract(rlcc, xy)
# 0: No data (NA/NA)
# 1: Temperate or sub-polar needleleaf forest (Conif/Forest)
# 2: Sub-polar taiga needleleaf forest (Conif/Forest)
# 5: Temperate or sub-polar broadleaf deciduous (DecidMixed/Forest)
# 6:  Mixed Forest (DecidMixed/Forest)
# 8: Temperate or sub-polar shrubland (Open/OpenWet)
# 10: Temperate or sub-polar grassland (Open/OpenWet)
# 11: Sub-polar or polar shrubland-lichen-moss (Open/OpenWet)
# 12: Sub-polar or polar grassland-lichen-moss (Open/OpenWet)
# 13: Sub-polar or polar barren-lichen-moss (Open/OpenWet)
# 14: Wetland (Wet/OpenWet)
# 15: Cropland (Open/OpenWet)
# 16: Barren Lands (Open/OpenWet)
# 17: Urban and Built-up (Open/OpenWet)
# 18: Water (NA/NA)
# 19: Snow and Ice (NA/NA)
lcclevs <- c("0"="", "1"="Conif", "2"="Conif", "3"="", "4"="",
  "5"="DecidMixed", "6"="DecidMixed", "7"="", "8"="Open", "9"="",
  "10"="Open", "11"="Open", "12"="Open", "13"="Open", "14"="Wet",
  "15"="Open", "16"="Open", "17"="Open", "18"="", "19"="")
lcc4 <- factor(lcclevs[vlcc+1], c("DecidMixed", "Conif", "Open", "Wet"))
lcc2 <- lcc4
levels(lcc2) <- c("Forest", "Forest", "OpenWet", "OpenWet")

## TREE
vtree <- extract(rtree, xy)
TREE <- vtree / 100
TREE[TREE %)(% c(0, 1)] <- 0

## UTC offset + 7 makes Alberta 0 (MDT offset)
tz <- extract(rtz, xy) + 7

## transform the rest
JDAY <- round(day / 365, 4) # 0-365
TREE <- round(vtree / 100, 4)
MAXDIS <- round(dis / 100, 4)
MAXDUR <- round(dur, 4)

## sunrise time adjusted by offset
sr <- sunriset(cbind("X"=lon, "Y"=lat),
  as.POSIXct(dtm, tz="America/Edmonton"),
  direction="sunrise", POSIXct.out=FALSE) * 24
TSSR <- round(unname((hour - sr + tz) / 24), 4)

## constant for NA cases
cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
## best model
mi <- bestmodelBAMspecies(spp, type="BIC",
    model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)

## make Xp and Xq
#' Design matrices for singing rates (`Xp`) and for EDR (`Xq`)
Xp <- cbind(
  "(Intercept)"=1,
  "TSSR"=TSSR,
  "JDAY"=JDAY,
  "TSSR2"=TSSR^2,
  "JDAY2"=JDAY^2)
Xq <- cbind("(Intercept)"=1,
  "TREE"=TREE,
  "LCC2OpenWet"=ifelse(lcc4 %in% c("Open", "Wet"), 1, 0),
  "LCC4Conif"=ifelse(lcc4=="Conif", 1, 0),
  "LCC4Open"=ifelse(lcc4=="Open", 1, 0),
  "LCC4Wet"=ifelse(lcc4=="Wet", 1, 0))

p <- rep(NA, n)
A <- q <- p
## constant for NA cases
cf0 <- exp(unlist(coefBAMspecies(spp, 0, 0)))
## best model
mi <- bestmodelBAMspecies(spp, type="BIC",
    model.sra=names(getBAMmodellist()$sra)[!grepl("DSLS", getBAMmodellist()$sra)])
cfi <- coefBAMspecies(spp, mi$sra, mi$edr)
## design matrices matching the coefs
Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
OKp <- rowSums(is.na(Xp2)) == 0
Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
OKq <- rowSums(is.na(Xq2)) == 0
## calculate p, q, and A based on constant phi and tau for the respective NAs
p[!OKp] <- sra_fun(MAXDUR[!OKp], cf0[1])
unlim <- ifelse(MAXDIS[!OKq] == Inf, TRUE, FALSE)
A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * MAXDIS[!OKq]^2)
q[!OKq] <- ifelse(unlim, 1, edr_fun(MAXDIS[!OKq], cf0[2]))
## calculate time/lcc varying phi and tau for non-NA cases
phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
p[OKp] <- sra_fun(MAXDUR[OKp], phi1)
unlim <- ifelse(MAXDIS[OKq] == Inf, TRUE, FALSE)
A[OKq] <- ifelse(unlim, pi * tau1^2, pi * MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(MAXDIS[OKq], tau1))
## log(0) is not a good thing, apply constant instead
ii <- which(p == 0)
p[ii] <- sra_fun(MAXDUR[ii], cf0[1])

DIG <- 6
out <- list(
  p=round(p, DIG),
  q=round(q, DIG),
  A=round(A, DIG),
  correction=round(p*A*q, DIG),
  offset=round(log(p) + log(A) + log(q), DIG))
out




