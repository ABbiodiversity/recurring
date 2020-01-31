## offsets

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

spp <- "OVEN"
day <- 180 # 0-365
hour <- 5.2 # 0-24
lon <- 100 # longitude
lat <- 50 # latitude
tree <- 100 # 0-100 percent tree cover
lcc <- "Conif" # values
dur <- 10 # mins
dis <- 100 # meters

## types
spp <- as.character(spp)
lcc <- as.character(lcc)
day <- as.integer(round(as.numeric(day)))
hour <- as.numeric(hour)
lat <- as.numeric(lat)
lon <- as.numeric(lon)
tree <- as.integer(round(as.numeric(tree)))
dur <- as.numeric(dur)
dis <- as.numeric(dis)

## checks
if (!(spp %in% getBAMspecieslist()))
  stop(sprintf("Species %s has no QPAD estimate", spp))
checkfun <- function(x, name="", range=c(-Inf, Inf)) {
  if (is.na(x))
    stop(sprintf("Parameter %s is NA, check intpu type", name))
  if (x %)(% range)
    stop(sprintf("Parameter %s is out of range (%.0f - %.0f)", name, range[1], range[2]))
  invisible(NULL)
}
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
day <- day / 365
tree <- tree / 100
dis <- dis / 100

## figure out a date/time value
JL <- as.POSIXct(dd$DATI, tz="America/Edmonton")
JL <- Sys.time()
sr <- sunriset(cbind("X"=lon, "Y"=lat), JL, direction="sunrise", POSIXct.out=FALSE) * 24
tssr <- unname((hour - sr) / 24)


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
  "TSSR"=tssr,
  "JDAY"=day,
  "TSSR2"=tssr^2,
  "JDAY2"=day^2)
Xq <- cbind("(Intercept)"=1,
  "TREE"=tree,
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
p[!OKp] <- sra_fun(dd$MAXDUR[!OKp], cf0[1])
unlim <- ifelse(dd$MAXDIS[!OKq] == Inf, TRUE, FALSE)
A[!OKq] <- ifelse(unlim, pi * cf0[2]^2, pi * dd$MAXDIS[!OKq]^2)
q[!OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[!OKq], cf0[2]))
## calculate time/lcc varying phi and tau for non-NA cases
phi1 <- exp(drop(Xp2[OKp,,drop=FALSE] %*% cfi$sra))
tau1 <- exp(drop(Xq2[OKq,,drop=FALSE] %*% cfi$edr))
p[OKp] <- sra_fun(dd$MAXDUR[OKp], phi1)
unlim <- ifelse(dd$MAXDIS[OKq] == Inf, TRUE, FALSE)
A[OKq] <- ifelse(unlim, pi * tau1^2, pi * dd$MAXDIS[OKq]^2)
q[OKq] <- ifelse(unlim, 1, edr_fun(dd$MAXDIS[OKq], tau1))
## log(0) is not a good thing, apply constant instead
ii <- which(p == 0)
p[ii] <- sra_fun(dd$MAXDUR[ii], cf0[1])
## store, next
off[,spp] <- log(p) + log(A) + log(q)

