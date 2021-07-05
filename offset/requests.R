## 2021-07-05 Lisa Venier

## run setup scripts

## 5-min Unlimited (0-5 min, 0-100-Inf m)
pk1 <- read.csv("~/Desktop/x/LV_BAM_PKEY_AOU_Atlasonly.txt", stringsAsFactors = FALSE)
ss1 <- read.csv("~/Desktop/x/LV_BAM_XY_AOU_Atlasonly.txt", stringsAsFactors = FALSE)
pk1 <- data.frame(pk1, ss1[match(pk1$SS, ss1$SS),])
pk1$MAXDUR <- 5
pk1$MAXDIS <- Inf
## 3-min Unlimited (0-3 min, 0-Inf m)
pk2 <- read.csv("~/Desktop/x/LV_BBS_V3_PKEY_AOU.txt", stringsAsFactors = FALSE)
ss2 <- read.csv("~/Desktop/x/LV_BBS_V3_XY_AOU.txt", stringsAsFactors = FALSE)
pk2 <- data.frame(pk2, ss2[match(pk2$SS, ss2$SS),])
pk2$MAXDUR <- 3
pk2$MAXDIS <- Inf

cn <- intersect(colnames(pk1),colnames(pk2))
pk <- rbind(pk1[,cn], pk2[,cn])
data.frame(NAS=colSums(is.na(pk)))

pk$dt <- paste0(pk$YEAR, "-0", pk$MONTH, "-", ifelse(pk$DAY < 10, "0", ""), pk$DAY)
pk$tm <- paste0(ifelse(pk$HOUR < 10, "0", ""), pk$HOUR, ifelse(pk$MIN < 10, ":0", ":"), pk$MIN)
summary(as.Date(pk$dt))

x <- make_x(
  dt=pk$dt,
  tm=pk$tm,
  lon=pk$X,
  lat=pk$Y,
  dur=pk$MAXDUR,
  dis=pk$MAXDIS,
  key=pk$PKEY)
head(x)
summary(x)

SPP <- getBAMspecieslist()

OFF <- matrix(0, nrow(x), length(SPP))
rownames(OFF) <- x$key
colnames(OFF) <- SPP

for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  o <- make_off(spp, x)
  OFF[,spp] <- o$offset
}

OUT <- data.frame(x, Offset=OFF)
write.csv(OUT, row.names=FALSE,
  file="~/Desktop/x/LV_BAM_BBS_ON_2021-07-05-offsets.csv")



## 2021-01-21

library(readxl)
library(mefa4)

XX <- read_excel("~/Desktop/BAM MB-ON-BBATLAS data 20Jan20201.xlsx")
XX <- as.data.frame(XX)
X <- nonDuplicated(as.data.frame(XX), XX$PKEY, TRUE)
X$dt <- as.character(X$survey_date)
nc <- nchar(as.character(X$survey_time))
X$tm <- substr(as.character(X$survey_time), nc-7, nc-3)
X <- X[!is.na(X$latitude),]

unique(X$protocol_distance)
unique(X$protocol_duration)

x <- make_x(dt=X$dt, tm=X$tm,
  lon=X$longitude, lat=X$latitude,
  dur=5, dis=Inf, key=X$PKEY)
head(x)

SPPx <- unique(X$species_code)
SPPo <- getBAMspecieslist()
SPP <- sort(intersect(SPPx, SPPo))

OFF <- matrix(0, nrow(x), length(SPP))
rownames(OFF) <- x$key
colnames(OFF) <- SPP

for (spp in SPP) {
  cat(spp, "\n")
  flush.console()
  o <- make_off(spp, x)
  OFF[,spp] <- o$offset
}

OUT <- data.frame(x, Offset=OFF)
write.csv(OUT, row.names=FALSE,
  file="~/Desktop//BAM MB-ON-BBATLAS data 20Jan20201-offsets.csv")



## 2020

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

