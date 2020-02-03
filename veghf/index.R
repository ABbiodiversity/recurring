library(mefa4)
library(DBI)
source("function.R")

UID_COL    = "ABMI_ID_Wi"
BASE_YR    = 2019 # or column name as character
FILE       = "sites-example.sqlite"
TABLE      = "Summary_Buffers" # only if sqlite
SUB_COL    = "Section" # or NULL
SUB_VAL    = c("NE","NW","SE","SW")

if (endsWith(FILE, ".csv")) {
  d <- read.csv(FILE)
} else {
  db <- dbConnect(RSQLite::SQLite(), FILE)
  dbListTables(db)
  d <- dbReadTable(db, TABLE)
  dbDisconnect(db)
  d <- make_char2fact(d)
}

if (!is.null(SUB_COL)) {
  d <- d[d[[SUB_COL]] %in% SUB_VAL,]
}

d_long <- make_vegHF_wide_v6(d,
    col.label=UID_COL,
    col.year=BASE_YR,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE)
d_wide0 <- make_vegHF_wide_v6(d_long,
    col.label=UID_COL,
    col.year=BASE_YR,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, widen_only=TRUE)

dx <- nonDuplicated(d, d[[UID_COL]], TRUE)[rownames(d_wide0[[1]]),]
d_wide <- fill_in_0ages_v6(d_wide0, dx$NSRNAME, ages_list)

stub <- strsplit(FILE, "\\.")[[1]]
stub <- paste0(stub[-length(stub)], collapse="_")
save(d_long, d_wide, file=paste0(stub, "_", Sys.Date(), ".RData"))


