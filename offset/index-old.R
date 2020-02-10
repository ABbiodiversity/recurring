## offsets, works only in AB
## to extend beyon AP, need time zone tweaks + TREE/LCC layers

if (!requireNamespace("QPAD")) {
  if (!requireNamespace("remotes"))
    install.packages("remotes")
  remotes::install_github("psolymos/QPAD")
}
if (!requireNamespace("maptools"))
  install.packages("maptools")
if (!requireNamespace("intrval"))
  install.packages("intrval")

library(QPAD)
library(maptools)
library(intrval)

load_BAM_QPAD(version = 3)
if (getBAMversion() != "3")
  stop("This script requires BAM version 3")

spp <- "OVEN"
## https://en.wikipedia.org/wiki/ISO_8601
dt <- "2019-06-07" # ISO 8601 in YYYY-MM-DD (0-padded)
tm <- "05:20" # ISO 8601 in hh:mm (24 hr clock, 0-padded)
lon <- -100 # longitude WGS84 (EPSG: 4326)
lat <- 50 # latitude WGS84 (EPSG: 4326)
dur <- 10 # mins
dis <- 100 # meters
tree <- 100 # 0-100 percent tree cover
lcc <- "Conif" # values


## reproject xy
## extract tree and lcc
## reclass lcc
if (FALSE) {
if (!requireNamespace("raster"))
  install.packages("raster")
library(raster)
library(rgdal)
library(sp)
library(rgeos)
rt <- raster("data/spatial/tree.tif")
rl1 <- stack("d:/bam/BAM_data_v2019/gnm/data/subunits/bcr60all_1km.gri")
rl2 <- stack("d:/bam/BAM_data_v2019/gnm/data/subunits/bcr10all_1km.gri")
rl3 <- stack("d:/bam/BAM_data_v2019/gnm/data/subunits/bcr11all_1km.gri")
rl <- mosaic(rl1[["nalc"]], rl2[["nalc"]], rl2[["nalc"]], fun=mean)
rl <- spTransform(rl, proj4string(rt))
od <- setwd(file.path("~/Dropbox/courses/st-johns-2017", "data", "NatRegAB"))
AB <- readOGR(".", "Natural_Regions_Subregions_of_Alberta") # rgdal
ABpr <- gUnaryUnion(AB, rep(1, nrow(AB))) # province
setwd(od)
}


## types
spp <- as.character(spp)
lcc <- as.character(lcc)
lat <- as.numeric(lat)
lon <- as.numeric(lon)
tree <- as.integer(round(as.numeric(tree)))
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
  if (is.na(x))
    stop(sprintf("Parameter %s is NA, check intpu type", name))
  if (x %)(% range)
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
checkfun(tree, "tree", c(0, 100))
checkfun(dur, "dur", c(0, Inf))
checkfun(dis, "dis", c(0, Inf))
if (is.infinite(lon))
  stop("Parameter lon must be finite")
if (is.infinite(lat))
  stop("Parameter lat must be finite")
if (!(lcc %in% c("Conif", "DecidMixed", "Open", "Wet")))
  stop("Parameter lcc must be one of [Conif, DecidMixed, Open, Wet]")
#lcc4 <- factor(lcc, c("DecidMixed", "Conif", "Open", "Wet"))
#lcc2 <- lcc4
#levels(lcc2) <- c("Forest", "Forest", "OpenWet", "OpenWet")

## transform
JDAY <- round(day / 365, 4) # 0-365
TREE <- round(tree / 100, 4)
MAXDIS <- round(dis / 100, 4)
MAXDUR <- round(dur, 4)

sr <- sunriset(cbind("X"=lon, "Y"=lat),
  as.POSIXct(dtm, tz="America/Edmonton"),
  direction="sunrise", POSIXct.out=FALSE) * 24
TSSR <- round(unname((hour - sr) / 24), 4)


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
  "LCC2OpenWet"=ifelse(lcc %in% c("Open", "Wet"), 1, 0),
  "LCC4Conif"=ifelse(lcc=="Conif", 1, 0),
  "LCC4Open"=ifelse(lcc=="Open", 1, 0),
  "LCC4Wet"=ifelse(lcc=="Wet", 1, 0))

## design matrices matching the coefs
Xp2 <- Xp[,names(cfi$sra),drop=FALSE]
OKp <- rowSums(is.na(Xp2)) == 0
Xq2 <- Xq[,names(cfi$edr),drop=FALSE]
OKq <- rowSums(is.na(Xq2)) == 0
## calculate p, q, and A based on constant phi and tau for the respective NAs
if (OKp) {
  phi1 <- exp(drop(Xp2 %*% cfi$sra))
  p <- sra_fun(MAXDUR, phi1)
  if (p == 0)
    p <- sra_fun(MAXDUR, cf0[1])
} else {
  p <- sra_fun(MAXDUR, cf0[1])
}
unlim <- is.infinite(MAXDIS)
if (OKq) {
  tau1 <- exp(drop(Xq2 %*% cfi$edr))
  A <- ifelse(unlim, pi * tau1^2, pi * MAXDIS^2)
  q <- ifelse(unlim, 1, edr_fun(MAXDIS, tau1))
} else {
  A <- ifelse(unlim, pi * cf0[2]^2, pi * MAXDIS^2)
  q <- ifelse(unlim, 1, edr_fun(MAXDIS, cf0[2]))
}

out <- list(
  settings=list(jday=JDAY, tssr=TSSR, tree=TREE, lcc=lcc, maxdur=MAXDUR, maxdis=MAXDIS, species=spp),
  results=list(p=p, q=q, A=A, correction=p*A*q, offset=log(p) + log(A) + log(q))
)
out

## todo
## - get NALCMS & TREE layers (full N Am)
## - add intersection capabilities
## - solve TZ issue (TZ offset?)
## - vectorized API: csv or json input


# getting NALC
library(raster)

fl <- list.files("d:/bam/BAM_data_v2019/gnm/data/subunits",
  pattern = ".gri")
rl <- list()
for (i in seq_along(fl)) {
  st <- stack(paste0("d:/bam/BAM_data_v2019/gnm/data/subunits/", fl[i]))
  rl[[i]] <- st[["nalc"]]
}
nalc <- mosaic(rl[[1]], rl[[2]], rl[[3]], rl[[4]], rl[[5]], rl[[6]], rl[[7]], rl[[8]],
  rl[[9]], rl[[10]], rl[[11]], rl[[12]], rl[[13]], rl[[14]], rl[[15]], rl[[16]],
  rl[[17]], rl[[18]], rl[[19]], rl[[20]], rl[[21]], rl[[22]], rl[[23]], rl[[24]],
  fun=mean)
writeRaster(nalc, file="d:/spatial/offset-layers/nalc.tif")

nalc <- raster("d:/spatial/offset-layers/nalc.tif")

# tree
rt <- raster("d:/spatial/offset-layers/treecoverlcc1.tif")
rt <- crop(rt, extent(nalc))
rt <- mask(rt, nalc)

# tz
library(rgdal)
library(rgeos)
tz <- readOGR("d:/spatial/offset-layers/wtz")
tzz <- spTransform(tz, proj4string(nalc))
r <- rasterize(tz, nalc, as.integer(tz@data$UTC_OFFSET), fun='last')

tz1 <- readOGR("d:/spatial/offset-layers/us")
tz2 <- readOGR("d:/spatial/offset-layers/canada")
#tz1 <- gSimplify(tz1, tol = 0.05)
#tz2 <- gSimplify(tz2, tol = 0.05)

r1 <- rasterize(tz1, nalc, as.integer(tz1$TZID))
r2 <- rasterize(tz2, nalc, as.integer(tz2$TZID))


library(raster)
nalc <- raster("d:/spatial/offset-layers/nalc.tif")

tree <- raster("d:/spatial/offset-layers/treecoverlcc1.tif")
tree <- crop(tree, extent(nalc))

tz <- raster("d:/spatial/offset-layers/tz.tif")
tz <- crop(tz, extent(nalc))

tree2 <- resample(tree, nalc, "ngb")
tree2 <- mask(tree2, nalc)

tz2 <- resample(tz, nalc, "ngb")
tz2 <- mask(tz2, nalc)
values(tz2)[!is.na(values(tz2)) & values(tz2) > 5] <- -10

plot(nalc)
plot(tree2)
plot(tz2)

writeRaster(nalc, file="d:/spatial/offset-layers/offset-lcc.tif")
writeRaster(tree2, file="d:/spatial/offset-layers/offset-tree.tif")
writeRaster(tz2, file="d:/spatial/offset-layers/offset-tz.tif")

library(rgdal)
library(raster)

## adjust these paths as needed
ROOT <- "d:/bam/BAM_data_v2019/gnm"
load(file.path(ROOT, "data", "BAMdb-GNMsubset-2020-01-08.RData"))
SUB <- readOGR(dsn=file.path(ROOT, "data", "sub-bound"), "BCRSubunits")
BCR <- readOGR(dsn=file.path(ROOT, "data", "bcr"), "bcrfinallcc")
PROV <- readOGR(dsn=file.path(ROOT, "data", "prov"), "province_state_line")
LAKES <- readOGR(dsn=file.path(ROOT, "data", "lakes"), "lakes_lcc")
LAKES <- spTransform(LAKES, proj4string(BCR))

plot(tz)
plot(PROV, col="#aaaaaa", add=TRUE)


    rast <- raster(file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-SD.tif")))
    q <- quantile(values(rast), 0.999, na.rm=TRUE)
    values(rast)[!is.na(values(rast)) & values(rast)>q] <- q
    jpeg(file.path(ROOT, "artifacts", spp, paste0("pred-", spp, "-CAN-SD.jpg")),
        height=4000, width=6000, res=300)
    op <- par(cex.main=3, mfcol=c(1,1), oma=c(0,0,0,0), mar=c(0,0,5,0))
    plot(rast, col=bluegreen.colors(15), axes=FALSE, legend=TRUE, main=paste(spp, "SD"), box=FALSE)
    plot(PROV, col="#aaaaaa", add=TRUE)
    plot(LAKES,col="#aaaaff", border=NA, add=TRUE)
    plot(BCR, add=TRUE)
    par(op)
    dev.off()
