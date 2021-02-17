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

