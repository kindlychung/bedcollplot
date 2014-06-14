Plotcoll = setRefClass("Plotcoll",
    fields = list(
    bedstem="character",
    shiftstemCommon = "character",
    bedpath="character",
    fampath="character",
    bimpath="character",
    shiftFilesStem="character",
    corruptFileList="character",
    nsnp="integer",
    nindiv="integer",
    bytesSnp="integer",
    bytesTotal="integer",
    largestNshift="integer"
    ))
Plotcoll$methods(
    sayhi = function(personname) {
        message("Hi ", personname, ", how are you? I miss you a lot!")
    }
)
Plotcoll$methods(
    getNshiftStr = function(shiftpath) {
        nshiftStr = gsub(".*_shift_(\\d{4}).*", "\\1", shiftpath)
        nshiftStr
    }
)
Plotcoll$methods(
    getNshift = function(shiftpath) {
        nshift = as.integer(getNshiftStr(shiftpath))
        nshift
    }
)
Plotcoll$methods(
    getShiftStem = function(shiftpath) {
        shiftstem = paste(shiftstemCommon, getNshiftStr(shiftpath), sep="")
        shiftstem
    }
)
Plotcoll$methods(
    getstem = function(pathname) {
        gsub("(.*?)(\\.[^./]+)+", "\\1", pathname)
    }
)
Plotcoll$methods(
    readout = function(tagname) {
        outfiles = sapply(shiftFilesStem, function(eachstem) {
            paste(eachstem, tagname, sep = "_")
        })
        print(outfiles)
    }
)
Plotcoll$methods(
    initialize = function(bed) {
        bedstem <<- getstem(bed)
        bedpath <<- paste(bedstem, "bed", sep=".")
        fampath <<- paste(bedstem, "fam", sep=".")
        bimpath <<- paste(bedstem, "bim", sep=".")
        shiftstemCommon <<- paste(bedstem, "_shift_", sep = "")
        shiftFilesStem <<- Sys.glob(paste(bedstem, "_shift_*.bed", sep=""))
        shiftFilesStem <<- sapply(shiftFilesStem, getstem)
        names(shiftFilesStem) <<- NULL
        nsnp <<- R.utils::countLines(bimpath)
        nindiv <<- R.utils::countLines(fampath)
        bytesSnp <<- as.integer(ceiling(nindiv / 4))
        bytesTotal <<- bytesSnp * nsnp
    }
)

plotobj = Plotcoll("../ssknEx.bed")
plotobj$readout("20140614-233528.assoc.linear")

plotobj$bedpath
plotobj$shiftstemCommon
plotobj$shiftFilesStem
plotobj$bimpath
plotobj$fampath
plotobj$nsnp
plotobj$nindiv
plotobj$bytesSnp
plotobj$bytesTotal
plotobj$sayhi("Kaiyin")
