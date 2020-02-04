## provide settings here when using in interactive mode
if (interactive()) {

  ## --- EDIT BELOW THIS LINE ---

    ## unique ID (can be all the same value)
  UID_COL    = "ABMI_ID_Wi"

  ## base year for surveys or HF inventory
  ## this is used to calculate years since last disturbence
  ## (i.e. base year - origin year)
  ## use a numeric value when it is the same for each record
  ## use a character value to indicate a field when it
  ## varies by record
  BASE_YR    = 2019 # or column name as character

  ## input file name, can contain the path (i.e. /dir/file.csv)
  ## file type can be .csv or .sqlite
  FILE       = "sites-example.sqlite"

  ## table name for SQLite database
  ## ignored for csv files
  TABLE      = "Summary_Buffers"

  ## optional, field name to be used for subsetting
  ## can be NULL (ignored)
  SUB_COL    = "Section" # or NULL
  ## values in the <SUB_COL> field to keep
  ## can be a single value or a character vector
  SUB_VAL    = c("NE","NW","SE","SW")

  ## optional, the name of the output file
  ## it can contain path as well (e.g. /dir/file.RData)
  ## if NULL, <FILE>_YYYY-MM-DD.RData is used
  OUTPUT     = NULL

  ## keep as TRUE when a Shape_Area field is present
  ## (long and wide format summaries can be calculated)
  ## set it to FALSE when e.g. doing point intersections
  ## (only long summary can be calculated)
  AREA       = TRUE

  ## add comments here, e.g. describing the characteristics
  ## of the input (backfilled v6.1 + 2017 HFI) when it is
  ## not trivial from file name, or describe purpose of
  ## the summaries as a reminder
  COMMENTS   = ""

  ## --- DON'T EDIT BELOW THIS LINE ---

} else {
  source(commandArgs(trailingOnly = TRUE)[1L])
}

.veghf_settings <- list(
  input=list(
    UID_COL=UID_COL,
    BASE_YR=BASE_YR,
    FILE=FILE,
    TABLE=TABLE,
    SUB_COL=SUB_COL,
    SUB_VAL=SUB_VAL,
    OUTPUT=OUTPUT,
    AREA=AREA,
    COMMENTS=COMMENTS),
  session=list(date=Sys.time(), info=sessionInfo())
)

## install required species if needed
if (!requireNamespace("mefa4"))
  install.packages("mefa4")
if (!requireNamespace("DBI"))
  install.packages("DBI")
if (!requireNamespace("RSQLite"))
  install.packages("RSQLite")

## load required packages
library(mefa4)
library(DBI)
library(RSQLite)

## source required functions
## this also loads lookup tables
source("function.R")

## read in data
if (endsWith(tolower(FILE), ".csv")) {
  d <- read.csv(FILE)
} else {
  db <- dbConnect(RSQLite::SQLite(), FILE)
  dbListTables(db)
  d <- dbReadTable(db, TABLE)
  dbDisconnect(db)
  d <- make_char2fact(d)
}

## take a subset if needed
if (!is.null(SUB_COL)) {
  d <- d[d[[SUB_COL]] %in% SUB_VAL,]
}

## long format summary
d_long <- make_vegHF_wide_v6(d,
    col.label=UID_COL,
    col.year=BASE_YR,
    col.HFyear="YEAR",
    col.HABIT="Combined_ChgByCWCS",
    col.SOIL="Soil_Type_1",
    sparse=TRUE, HF_fine=TRUE, wide=FALSE)
## do not do this when Shape_Area is missing
if (AREA) {
  ## wide format summary with >=0 unknown age area
  d_wide0 <- make_vegHF_wide_v6(d_long,
      col.label=UID_COL,
      col.year=BASE_YR,
      col.HFyear="YEAR",
      col.HABIT="Combined_ChgByCWCS",
      col.SOIL="Soil_Type_1",
      sparse=TRUE, HF_fine=TRUE, widen_only=TRUE)
  ## wide format summary with unknown age area redistributed
  dx <- nonDuplicated(d, d[[UID_COL]], TRUE)[rownames(d_wide0[[1]]),]
  d_wide <- fill_in_0ages_v6(d_wide0, dx$NSRNAME, ages_list)
} else {
  d_wide <- NULL
}

## output file name
if (is.null(OUTPUT)) {
  stub <- strsplit(FILE, "\\.")[[1]]
  stub <- paste0(stub[-length(stub)], collapse="_")
  OUTPUT <- paste0(stub, "_", Sys.Date(), ".RData")
}

## save output
save(d_long, d_wide, .veghf_settings, file=OUTPUT)
