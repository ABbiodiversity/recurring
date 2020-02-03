# veghf

## reclass for natural veg
recl <- structure(list(V6_COMB = structure(1:30, .Label = c("Alkali",
  "AlpineLarch", "Bare", "Conif", "Decid", "Fir", "GraminoidFen",
  "GrassHerb", "Marsh", "Mixedwood", "Pine", "Shrub", "ShrubbyBog",
  "ShrubbyFen", "ShrubbySwamp", "SnowIce", "Spruce", "TreedBog-BSpr",
  "TreedFen-BSpr", "TreedFen-Decid", "TreedFen-Larch", "TreedFen-Mixedwood",
  "TreedSwamp-Conif", "TreedSwamp-Decid", "TreedSwamp-Fir", "TreedSwamp-Forest",
  "TreedSwamp-Mixedwood", "TreedSwamp-Spruce", "TreedWetland-Mixedwood",
  "Water"), class = "factor"), MERGED = structure(c(4L, 13L, 1L,
  13L, 2L, 13L, 3L, 4L, 5L, 6L, 7L, 8L, 9L, 10L, 11L, 12L, 13L,
  14L, 15L, 15L, 15L, 15L, 16L, 16L, 16L, 16L, 16L, 16L, 15L, 17L
  ), .Label = c("Bare", "Decid", "GraminoidFen", "GrassHerb", "Marsh",
  "Mixedwood", "Pine", "Shrub", "ShrubbyBog", "ShrubbyFen", "ShrubbySwamp",
  "SnowIce", "Spruce", "TreedBog", "TreedFen", "TreedSwamp", "Water"
  ), class = "factor")), class = "data.frame", row.names = c(NA,
  -30L))

## HF types
hftypes <- read.csv("lookup-hf-type-v2014.csv")
hftypes <- droplevels(hftypes[!is.na(hftypes$HF_GROUP_COMB) &
  !duplicated(hftypes$FEATURE_TY),])
hfgroups <- read.csv("lookup-hf-class-v2014.csv")
hflt <- hfgroups[match(hftypes$HF_GROUP_COMB, hfgroups$HF_GROUP_COMB),]
rownames(hflt) <- hftypes$FEATURE_TY

## age distribution by NR, NSR
load("ages-by-nsr.Rdata")

## functions

#' SQL tables represent character not factor
make_char2fact <- function(x) {
    if (is.null(dim(x)))
        if (is.character(x))
            return(as.factor(x))
    for (i in seq_len(ncol(x)))
        if (is.character(x[,i]))
            x[,i] <- as.factor(x[,i])
        x
}

## HF_fine T=use 2014 classification, F=use previous
make_vegHF_wide_v6 <-
function(d, col.label, col.year=NULL, col.HFyear=NULL,
col.HABIT=NULL, col.SOIL=NULL, wide=TRUE, sparse=FALSE,
HF_fine=TRUE, widen_only=FALSE) {

    TreedClassesCC <- c("Decid", "Mixedwood", "Pine", "Spruce")
    TreedClasses   <- c(TreedClassesCC, "TreedBog", "TreedFen", "TreedSwamp")
    NontreedClasses <- c("GrassHerb", "Shrub",
        "GraminoidFen", "Marsh",
        "ShrubbyBog", "ShrubbyFen", "ShrubbySwamp",
        "Bare", "SnowIce", "Water")
    if (!HF_fine) {
        HFLab <- c("BorrowpitsDugoutsSumps", "Canals", "CultivationCropPastureBareground",
            "CutBlocks", "HighDensityLivestockOperation", "IndustrialSiteRural",
            "MineSite", "MunicipalWaterSewage", "OtherDisturbedVegetation",
            "PeatMine", "Pipeline", "RailHardSurface", "RailVegetatedVerge",
            "Reservoirs", "RoadHardSurface", "RoadTrailVegetated", "RoadVegetatedVerge",
            "RuralResidentialIndustrial", "SeismicLine", "TransmissionLine",
            "Urban", "WellSite", "WindGenerationFacility")
    } else {
        HFLab <- c("UrbanIndustrial", "UrbanResidence", "RuralResidentialIndustrial",
            "IndustrialSiteRural", "WindGenerationFacility", "OtherDisturbedVegetation",
            "MineSite", "PeatMine", "WellSite", "Pipeline", "TransmissionLine",
            "SeismicLineNarrow", "SeismicLineWide", "RoadHardSurface", "RailHardSurface",
            "RoadTrailVegetated", "RoadVegetatedVerge", "RailVegetatedVerge",
            "CultivationCrop", "CultivationAbandoned", "CultivationRoughPasture",
            "CultivationTamePasture", "HighDensityLivestockOperation",
            "CutBlocks", "BorrowpitsDugoutsSumps", "MunicipalWaterSewage",
            "Reservoirs", "Canals")
    }
    RfLab <- c(paste0(rep(TreedClasses, each=11),
        c("0","R","1","2","3","4","5","6","7","8","9")),
        NontreedClasses)
    ## use necessary CC-forest classes once evaluated
    CrOnlyLab <- c(HFLab,
        paste0("CC", paste0(rep(TreedClassesCC, each=5),
        c("R","1","2","3","4"))))
    HLEVS <- c(TreedClasses, NontreedClasses)
    AllLabels <- c(RfLab, CrOnlyLab)

    SoilLab <- c("UNK", "Water", #"Wetland",
        "BdL", "BlO", "CS", "Cy", "Gr", "LenA", "LenSP",
        "LenT", "LenS", "Li", "Lo", "LtcC", "LtcD", "LtcH", "LtcS", "Ov",
        "Sa", "Sb", "SL", "SwG", "Sy", "TB")

    ## designate a label column (there are different column names in use)
    d$LABEL <- d[,col.label]

    if (!widen_only) {

        if (is.null(d[,col.HABIT]))
            stop("Shoot -- check the damn HABIT column...")
        d <- droplevels(d[!is.na(d[,col.HABIT]) & d[,col.HABIT] != "",])
        d$HABIT <- reclass(d[,col.HABIT], as.matrix(recl), all=TRUE)

        d$HF_Year <- d[,col.HFyear]
        if (any(is.na(d$LABEL)))
            stop("missing LABEL")
        #    d <- d[!is.na(d$LABEL),]
        ## designate a year column
        if (is.null(col.year)) {
            THIS_YEAR <- as.POSIXlt(Sys.Date())$year + 1900
            d$SampleYear <- THIS_YEAR
        } else {
            if (is.numeric(col.year)) {
                if (length(col.year) > 1)
                    stop("length of col.year > 1")
                THIS_YEAR <- col.year
                d$SampleYear <- THIS_YEAR
            } else {
                THIS_YEAR <- NA
                d$SampleYear <- d[,col.year]
            }
        }
        ## use upper-case labels for FEATURE_TY
        levels(d$FEATURE_TY) <- toupper(levels(d$FEATURE_TY))

        d$ORIGIN_YEAR <- d$Origin_Year

        #### Footprint classes:
        ## check if we have all the feature types in the lookup table
        ## "" blank is for non-HF classes in current veg
        levels(d$FEATURE_TY)[levels(d$FEATURE_TY) == "''"] <- ""
        levels(d$FEATURE_TY)[levels(d$FEATURE_TY) == " "] <- ""
        if (!all(setdiff(levels(d$FEATURE_TY), rownames(hflt)) == "")) {
            print(setdiff(levels(d$FEATURE_TY), c("", rownames(hflt))))
            stop("HF diff found, see above")
        }
        ## classify feature types according to the chosen level of HF designation
        ## which comes from hf.level column of hflt (HF lookup table)
        if (HF_fine) {
            d$HFclass <- hflt$HF_GROUP_COMB[match(d$FEATURE_TY, rownames(hflt))]
        } else {
            d$HFclass <- hflt$HF_GROUP[match(d$FEATURE_TY, rownames(hflt))]
        }
        d$HFclass <- as.factor(d$HFclass)
        ## HFclass inherits all levels from hflt[,hf.level]
        ## need to add in the blank for further parsing
        levels(d$HFclass) <- c(levels(d$HFclass), "")
        d$HFclass[is.na(d$HFclass)] <- ""

        ## slivers (tiny polys with no veg info):
        #stopifnot(max(d$Shape_Area[d$VEGclass == ""]) < 1)
        if (any(d$HABIT == ""))
            warning(paste("blank HABIT:", sum(d$Shape_Area[d$HABIT == ""]), "m^2"))
        #d <- d[d$HABIT != "",]
        d$HABIT <- droplevels(d$HABIT)

        #### HABIT/EC classes:
        ## follow HABIT/EC classes, but there are few oddities when outside of AVI
        #d$VEGclass <- d$EC_Type
        d$VEGclass <- d$HABIT
    #    levels(d$VEGclass)[levels(d$VEGclass) == "Non-Veg"] <- "NonVeg"
    #    levels(d$VEGclass) <- gsub("/", "", levels(d$VEGclass))

        if (length(setdiff(d$VEGclass, HLEVS)) > 0)
            stop(paste("check HABIT classes", setdiff(d$VEGclass, HLEVS)))
        #tb <- cbind(c("", HLEVS), c("", HLEVS))
    #return(setdiff(levels(d$VEGclass), tb[,1]))
        #levels(d$VEGclass) <- tb[match(levels(d$VEGclass), tb[,1]),2]

        #tmp <- aggregate(d$Shape_Area, list(lcc=d$EC_Type), sum)
        #tmp$p <- round(100*tmp$x/sum(tmp$x),2)
        #tmp2 <- aggregate(d$Shape_Area, list(lcc=d$VEGclass), sum)
        #tmp2$p <- round(100*tmp2$x/sum(tmp2$x),2)


    #    NontreedClasses <- setdiff(levels(d$VEGclass), TreedClasses)
    #    NontreedClasses <- NontreedClasses[NontreedClasses != ""]

        #### Age info for backfilled (Rf) and current (Cr)
        ## reference age class 0=no age (either not forest or no info)
        ## 1=0-19, 2=20-39, etc.
        d$ORIGIN_YEAR[!is.na(d$ORIGIN_YEAR) & d$ORIGIN_YEAR == 9999] <- NA
        d$ORIGIN_YEAR[!is.na(d$ORIGIN_YEAR) & d$ORIGIN_YEAR > d$SampleYear] <- NA
        d$AgeRf <- as.integer(sign(d$ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$ORIGIN_YEAR) / 20)))
        ## truncate reference age classes at 9 = 160+
        d$AgeRf[d$AgeRf > 9L] <- 9L
        ## catching origin_year > sample_year instances: this defaults to old
        d$AgeRf[d$AgeRf < 1] <- 9L
        ## placeholder for recent burn (0-9 years)
        tmp <- as.integer(sign(d$ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$ORIGIN_YEAR) / 10)))
        d$AgeRf[tmp == 1L] <- 999L
        ## set 0 year in treed habitats as max (assumed old forest)
    #    d$AgeRf[d$AgeRf == 0L & d$VEGclass %in% TreedClasses] <- 9L
        ## unknown age is set to 0
        #table(d$AgeRf, d$VEGclass, useNA="a") # check NA ORIGIN_YEAR values
        #d$AgeRf[is.na(d$AgeRf)] <- 0L

        ## incorporate HF year for cutblocks
        d$CC_ORIGIN_YEAR <- d$ORIGIN_YEAR
        ii <- d$HFclass == "CutBlocks"
        ii[ii & !is.na(d$ORIGIN_YEAR) & d$HF_Year >= d$ORIGIN_YEAR] <- TRUE
        ii[ii & is.na(d$ORIGIN_YEAR)] <- TRUE
        d$CC_ORIGIN_YEAR[ii] <- d$HF_Year[ii]
        d$CC_ORIGIN_YEAR[!is.na(d$CC_ORIGIN_YEAR) & d$CC_ORIGIN_YEAR == 9999] <- NA
        ## age for current with cutblock ages
        d$AgeCr <- as.integer(sign(d$CC_ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$CC_ORIGIN_YEAR) / 20)))
        ## truncate current age classes at 9
        d$AgeCr[d$AgeCr > 9L] <- 9L
        ## catching origin_year > sample_year instances: this defaults to old
        d$AgeCr[d$AgeCr < 1] <- 9L
        ## placeholder for recent CC (0-9 years)
        tmp <- as.integer(sign(d$CC_ORIGIN_YEAR) * (1 + floor((d$SampleYear - d$CC_ORIGIN_YEAR) / 10)))
        d$AgeCr[tmp == 1L] <- 999L
        ## unknown age is set to 0
        #table(d$AgeCr, d$VEGclass, useNA="a") # check NA ORIGIN_YEAR values
        #d$AgeCr[is.na(d$AgeCr)] <- 0L
        #table(d$AgeCr,useNA="a")

        ## correcting reference age class based on cutblock info:
        ## these happened as a result of backfilling, so we accept HF age instead
        ## but this should be rare (ref age must be >= current)
        #ii <- !is.na(d$AgeCr) & d$AgeCr > d$AgeRf & d$AgeCr < 999L
        #if (sum(ii)>0)
        #    warning(paste("AgeCr > AgeRf for this many cases:", sum(ii)))
        #d$AgeRf[ii] <- d$AgeCr[ii]
        d$AgeRf[is.na(d$AgeRf)] <- 0L
        #table(rf=d$AgeRf,cr=d$AgeCr,useNA="a")
        ## turning age values into factor:
        ## 0=no age info,
        ## 1:9=valid age classes for treed veg types,
        ## ""=non-treed
        ## 999=placeholder for _R_ecent burn "R"
        d$AgeRf <- factor(d$AgeRf, levels=c(as.character(c(0:9, 999)), ""))
        ## NA --> "0" as unknown age class
        d$AgeRf[is.na(d$AgeRf)] <- "0"
        ## age is not relevant in non-treed veg types
        d$AgeRf[!(d$VEGclass %in% TreedClasses)] <- ""
        ## burn
        levels(d$AgeRf)[levels(d$AgeRf)=="999"] <- "R"

        ## making current age as factor
        d$AgeCr <- factor(d$AgeCr, levels=c(as.character(c(0:9, 999)), ""))
        ## NA --> "0" as unknown age class
        d$AgeCr[is.na(d$AgeCr)] <- "0"
        ## age is not relevant in non-treed veg types (no HF)
        d$AgeCr[d$VEGclass %in% NontreedClasses & d$HFclass == ""] <- ""
        ## age is not relevant outside of cutblocks
        d$AgeCr[!(d$HFclass %in% c("", "CutBlocks"))] <- ""
        ## recent CC
        levels(d$AgeCr)[levels(d$AgeCr)=="999"] <- "R"
        #table(current=d$AgeCr, reference=d$AgeRf)

        #### Combining VEG, HF and Age:
        ## reference VEG + Age labels:
        d$VEGAGEclass <- interaction(d$VEGclass, d$AgeRf, drop=TRUE, sep="", lex.order=TRUE)
        levels(d$VEGAGEclass) <- c(levels(d$VEGAGEclass),
            setdiff(RfLab, levels(d$VEGAGEclass)))

        ## manage CC labels
        ## current veg+hf
        d$VEGHFclass <- d$VEGclass
        #CClabels <- paste0("CC", levels(d$VEGclass)[levels(d$VEGclass) != ""])
        CClabels <- paste0("CC", levels(d$VEGclass))
        tmp <- setdiff(levels(d$HFclass), levels(d$VEGclass))
        tmp <- tmp[!(tmp %in% c("", "CutBlocks"))]
        levels(d$VEGHFclass) <- c(levels(d$VEGHFclass), tmp, CClabels)
        ## add non-CC HF types
        d$VEGHFclass[!(d$HFclass %in% c("", "CutBlocks"))] <- d$HFclass[!(d$HFclass %in% c("", "CutBlocks"))]
        ## should later the non-merchendisable forests with CC should be redistributed?
        ## e.g. after producing the wide format
        ## update CC labels obly for <= 80 yr CC (usually this does not happen
        ## just to make sure labels are OK)
        ## anything above age class >4 is turned into 4 to avoid labeling issues (shrubland)
        d$AgeCr[d$HFclass == "CutBlocks" & d$AgeCr %in% c("5","6","7","8","9")] <- "4"
        ii <- d$HFclass == "CutBlocks" & d$AgeCr %in% c("0","R","1","2","3","4")
        if (sum(ii) > 0)
            d$VEGHFclass[ii] <- paste0("CC", as.character(d$VEGclass[ii]))

        ## labels where backfilled cutblock label is not forested habitat
        ## right now I just collapse them to see % of the areas
        ## it is usually < 10% at this scale so it might be safe to ignore them
        ## usually young ages, but ranges R-1-2-3
        #ii <- unlist(lapply(paste0("CC", NontreedClasses), grep, x=levels(d$VEGHFclass)))
        #levels(d$VEGHFclass)[ii] <- "CCOpenTypes"
        ## unknown types under 'CC' considered as 'CCOpenTypes'
        #levels(d$VEGHFclass)[levels(d$VEGHFclass) == "CC"] <- "CCOpenTypes"
        ## treed wetlands
        #ii <- unlist(lapply(paste0("CC", TreedWetClasses), grep, x=levels(d$VEGHFclass)))
        #levels(d$VEGHFclass)[ii] <- "CCWetTypes"

        ## current VEG + HF + Age labels:
        d$VEGHFAGEclass <- interaction(d$VEGHFclass, d$AgeCr, drop=TRUE, sep="", lex.order=TRUE)
        ## labels with 0 age category are also to be fixed later ------> hard stuff
        #ii <- unlist(lapply(paste0(TreedClasses, 0), grep, x=levels(d$VEGHFAGEclass)))
        #levels(d$VEGHFAGEclass)[ii] <- "CCproblem"
        ## Labels for output columns
        levels(d$VEGHFAGEclass) <- c(levels(d$VEGHFAGEclass), setdiff(AllLabels, levels(d$VEGHFAGEclass)))

        #### soils:
        if (is.null(d[,col.SOIL]))
            stop("Shoot -- check the damn SOIL column...")
        d$SOILclass <- d[,col.SOIL]
        ## need to have the UNKnown class to be able to deal with NAs
        if (!is.factor(d$SOILclass))
            d$SOILclass <- as.factor(d$SOILclass)
        if (!any(levels(d$SOILclass) == ""))
            levels(d$SOILclass) <- c(levels(d$SOILclass), "")
        ## dealing with NAs
        d$SOILclass[is.na(d$SOILclass)] <- ""
        ## unknown soil type outside of GVI and Dry Mixedwood
        levels(d$SOILclass)[levels(d$SOILclass) == ""] <- "UNK"
        levels(d$SOILclass)[levels(d$SOILclass) == " "] <- "UNK"
        ## get rid of modifiers
        levels(d$SOILclass) <- sapply(strsplit(levels(d$SOILclass), "-"), function(z) z[1L])
        ## add in Water label
        levels(d$SOILclass) <- c(levels(d$SOILclass), "Water")

        ## treat these as Water or Wetland?
        levels(d$SOILclass)[levels(d$SOILclass) %in% c("Len","LenW","Ltc","LtcR")] <- "Water"
    #    levels(d$SOILclass)[levels(d$SOILclass) %in% c("Len","LenW","Ltc","LtcR")] <- "Wetland"
        ## DEM/EC based Water class overrides soil
        d$SOILclass[d$VEGclass == "Water"] <- "Water"

        levels(d$SOILclass) <- c(levels(d$SOILclass), setdiff(SoilLab, levels(d$SOILclass)))
        d$SOILHFclass <- d$SOILclass
        levels(d$SOILHFclass) <- c(levels(d$SOILHFclass), levels(d$HFclass)[levels(d$HFclass) != ""])
        d$SOILHFclass[d$HFclass != ""] <- d$HFclass[d$HFclass != ""]
        ## NOTE: current UNK can be smaller than reference UNK, it can be turned into HF
        ## currently this is not tracked

        ## for point intersection or transition matrix processing, etc.
        if (!wide)
            return(d)
    }

    SoilHFLab <- levels(d$SOILHFclass)
    #### crosstabs
    ## veg reference
    VegRf <- Xtab(Shape_Area ~ LABEL + VEGAGEclass, d)
    ## veg + HF current
    VegCr <- Xtab(Shape_Area ~ LABEL + VEGHFAGEclass, d)
    ## soils (`reference`)
    SoilRf <- Xtab(Shape_Area ~ LABEL + SOILclass, d)
    ## soils (`current`, soil + HF)
    SoilCr <- Xtab(Shape_Area ~ LABEL + SOILHFclass, d)

    if (!all(colnames(VegRf) %in% RfLab)) {
        cat(colnames(VegRf)[!(colnames(VegRf) %in% RfLab)], sep="\n")
        stop("Unexpected VegRf label")
    }
    if (!all(colnames(VegCr) %in% AllLabels)) {
        cat(colnames(VegCr)[!(colnames(VegCr) %in% AllLabels)], sep="\n")
        stop("Unexpected VegCr label")
    }
    if (!all(colnames(SoilRf) %in% SoilLab)) {
        cat(colnames(SoilRf)[!(colnames(SoilRf) %in% SoilLab)], sep="\n")
        stop("Unexpected SoilRf label")
    }
    if (!all(colnames(SoilCr) %in% SoilHFLab)) {
        cat(colnames(SoilCr)[!(colnames(SoilCr) %in% SoilHFLab)], sep="\n")
        stop("Unexpected SoilCr label")
    }

    rn <- rownames(VegRf) # make sure row labels are identical across tables
    VegRf <- VegRf[rn, RfLab, drop=FALSE]
    VegCr <- VegCr[rn, AllLabels, drop=FALSE]
    SoilRf <- SoilRf[rn, SoilLab, drop=FALSE]
    SoilCr <- SoilCr[rn, SoilHFLab, drop=FALSE]

    out <- list(veg_current=VegCr,
        veg_reference=VegRf,
        soil_current=SoilCr,
        soil_reference=SoilRf)
#        rs_veg_current=structure(as.integer(rsVegCr), names=names(rsVegCr)),
#        rs_veg_reference=structure(as.integer(rsVegRf), names=names(rsVegRf)),
#        rs_soil_current=structure(as.integer(rsSoilCr), names=names(rsSoilCr)),
#        rs_soil_reference=structure(as.integer(rsSoilRf), names=names(rsSoilRf)))
    if (!sparse)
        out <- lapply(out, as.matrix)

    ## year for each row
    tmp <- nonDuplicated(d, LABEL, TRUE)
    tmp <- tmp[rownames(VegCr),]
    #out$sample_year <- THIS_YEAR
    out$sample_year <- tmp$SampleYear

    out
}

fill_in_0ages_v6 <- function(x, NSR, ages_list, rm0=TRUE) {
    Target <- names(ages_list)

    Ages <- c("0", "R", as.character(1:9))
    NSR <- droplevels(as.factor(NSR))
    NSRs <- levels(NSR)
    for (current in c(TRUE, FALSE)) {
        xx <- if (current)
            x$veg_current else x$veg_reference
        xx <- as.matrix(xx)
        ag <- if (current)
            AvgAgesNSROld$current else AvgAgesNSROld$reference
        ag2 <- if (current)
            AvgAgesAllOld$current else AvgAgesAllOld$reference
        for (nsr in NSRs) {
            cat(ifelse(current, "current:", "reference:"), nsr)
            flush.console()
            for (i in Target) {
                Cols <- paste0(i, Ages)
                j <- NSR == nsr
                if (any(j)) {
                    p0 <- ag[[i]][nsr,]
                    if (sum(p0) == 0)
                        p0 <- ag2[[i]]
                    Mat <- xx[j, Cols, drop=FALSE]
                    Mat0 <- Mat
                    ## multiply Mat[,1] (unknown age) with this matrix
                    Unk <- Mat[,1] * t(matrix(p0, length(Ages), sum(j)))
                    Mat[,1] <- 0 # will be 0 and redistributed from Unk
                    Mat <- Mat + Unk
                    xx[j, Cols] <- Mat # ridiculously slow as sparse matrix
                    if (sum(Mat0)-sum(Mat) > 10^-6)
                        cat("\n\ttype:", i, "| diff =", round((sum(Mat0)-sum(Mat))/10^6))
                }
            }
            cat(" --- OK\n")
        }
        if (current) {
            x$veg_current <- as(xx, "dgCMatrix")
        } else {
            x$veg_reference <- as(xx, "dgCMatrix")
        }
    }
    if (rm0) {
        excl <- c("Decid0", "Mixedwood0", "Pine0", "Spruce0", "TreedBog0", "TreedFen0",
            "TreedSwamp0", "CutBlocks")
        x$veg_current <- x$veg_current[,!(colnames(x$veg_current) %in% excl)]
        x$veg_reference <- x$veg_reference[,!(colnames(x$veg_reference) %in% excl)]
    }
    x
}
