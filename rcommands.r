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
        nshiftStr = as.integer(gsub(".*_shift_(\\d{4}).*", "\\1", shiftpath))
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
        gsub("(.*?)(\\..*)+", "\\1", pathname)
    }
)
Plotcoll$methods(
    initialize = function(bed) {
        bedstem <<- getstem(bed)
        shiftstemCommon <<- paste(bedstem, "_shift_", sep = "")
        bedpath <<- paste(bedstem, "bed", sep=".")
        fampath <<- paste(bedstem, "fam", sep=".")
        bimpath <<- paste(bedstem, "bim", sep=".")
        shiftFilesStem <<- Sys.glob(paste(bedstem, "_shift_*.bed", sep=""))
        shiftFilesStem <<- sapply(shiftFilesStem, getstem)
        names(shiftFilesStem) <<- NULL
        nsnp <<- R.utils::countLines(bimpath)
        nindiv <<- R.utils::countLines(fampath)
        bytesSnp <<- as.integer(ceiling(nindiv / 4))
        bytesTotal <<- bytesSnp * nsnp
    }
)

plotobj = Plotcoll("ssknEx.bed")
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

mEdit = setRefClass("mEdit", fields = list(data="matrix", edits="list"))
mEdit$methods(
    edit = function(i, j, value) {
        backup = list(i, j, data[i, j])
        data[i, j] <<- value
        edits <<- c(edits, list(backup))
        invisible(value)
    }
)
mEdit$methods(
undo = function() {
    prev = edits
    if(length(prev)) {
        prev = prev[[length(prev)]]
    }
    else {
        stop("No more edits to undo!")
    }
    edit(prev[[1]], prev[[2]], prev[[3]])
    length(edits) <<- length(edits) - 2
    invisible(prev)
}
)
mEdit$methods(
    show = function() {
        message("ClassName: ", classLabel(class(.self)))
        message("Data:")
        methods::show(data)
        message("Undo list length: ", length(edits))
    }
)
mEdit$methods(
    .DollarNames.mEdit = function(x, pattern) {
        grep(pattern, getRefClass(class(x))$methods(), value=TRUE)
    }
)

x = matrix(1:24, 3, 8)
xx = mEdit(data=x)
xx$edit(2,2,0)
xx$show()
xx$edit(3, 5, 1)
xx$show()
xx$undo()
xx$show()

mv = setRefClass(
"matrixViewer",
fields=c("viewerDevice", "viewerFile"),
contains="mEdit"
                 )
mv$methods(
    .DollarNames.mEdit = function(x, pattern) {
        grep(pattern, getRefClass(class(x))$methods(), value=TRUE)
    }
)

mv$methods(
    view = function() {
        ## dd = dev.cur();
        ## dev.set(viewerDevice)
        ## devAskNewPage(FALSE)
        image(
            data,
            main=paste("After", length(edits), "edits")
        )
        ## dev.set(dd)
    }
)
mv$methods(
    edit = function(i,j, value) {
        callSuper(i,j, value)
        view()
    }
)
mv$methods(
    initialize = function(file="./mv.pdf", ...) {
        viewerFile <<- file
        ## pdf(viewerFile)
        ## viewerDevice <<- dev.cur()
        ## dev.set(dev.prev())
        callSuper(...)
    }
)
mv$methods(
    finalize = function() {
        dev.off(viewerDevice)
    }
)


x = matrix(rnorm(64, 0, 34), 8, 8)
xx = mv(file="/tmp/x.pdf", data=x)
xx$edit(2,2,0)
xx$edit(3, 5, 1)
xx$edit(4, 4, 2.3)
xx$undo()
xx$view()

