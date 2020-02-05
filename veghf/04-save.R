## output file name
if (is.null(OUTPUT)) {
  stub <- strsplit(FILE, "\\.")[[1]]
  stub <- paste0(stub[-length(stub)], collapse="_")
  OUTPUT <- paste0(stub, "_", Sys.Date(), ".RData")
}

## session info etc
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

## save output
cat("Saving results:\n", OUTPUT)
save(d_long, d_wide, .veghf_settings, file=OUTPUT)
cat("\n\nDONE\n\n")
