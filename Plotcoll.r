require(Matrix)
pminNoNA = function(...) {
    pmin(..., na.rm = TRUE)
}

Plotcoll = setRefClass("Plotcoll",
    fields = list(
        bedstem="character",
        shiftstemCommon = "character",
        bedpath="character",
        fampath="character",
        bimpath="character",
        shiftFilesStem="character",
        nsnp="numeric",
        nindiv="numeric",
        nshiftMax = "numeric",
        nshiftStrs = "character",
        chr = "matrix",
        snp = "matrix",
        bp = "matrix",
        pvals = "matrix"
    ))
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
    updateShiftFilesStem = function() {
        shiftFilesStem <<- Sys.glob(paste(bedstem, "_shift_*.bed", sep=""))
        shiftFilesStem <<- getstem(shiftFilesStem)
        names(shiftFilesStem) <<- NULL
        nshiftStrs <<- as.character(getNshift(shiftFilesStem))
    }
)

Plotcoll$methods(
    readout = function(tagname) {
        outRdata = paste(bedstem, "RData", sep = ".")
        if(file.exists(outRdata)) {
            message("Reading from previously generated data...")
            load(outRdata)

            # make sure dimensions are right
            dim1 = dim(chr)
            dim2 = dim(snp)
            dim3 = dim(bp)
            dim4 = dim(pvals)
            if(
               length(unique(c(dim1[1], dim2[1], dim3[1], dim4[1]))) > 1 |
               length(unique(c(dim1[2], dim2[2], dim3[2], dim4[2]))) > 1 
               ) {
                stop(paste(outRdata, "is corrupt, you had better remove it!", sep = " "))
            }

            # figure out what files to read
            hasColnames = colnames(chr)
            targetColnames = c("0", nshiftStrs)
            addIdx = which(!(targetColnames %in% hasColnames))
            message("Index for adding: ")
            print(addIdx)

            selectIdx = which(hasColnames %in% targetColnames)
            removeIdx = ! selectIdx
            removeItems = hasColnames(removeIdx)


            message("Has cols: ")
            print(hasColnames)
            message("Should have cols: ")
            print(targetColnames)



            ## if(ncol(pvals) == nshiftMax + 1) {
            ##     message("Good, I have all the p values. I will just use them.")
            ##     pvals <<- pvals
            ## } else {
            ##     message("I only have a part of all p values, need to read them all in again.")
            ##     needReadFromPlinkOut = TRUE
            ## }
            needReadFromPlinkOut = FALSE
        } else {
            needReadFromPlinkOut = TRUE
        }

        if(needReadFromPlinkOut == TRUE) {
            # set up matrices
            pvalsNcols = length(nshiftStrs) + 1
            pvals <<- matrix(NA, nsnp, pvalsNcols)
            chr <<- pvals
            snp <<- pvals
            bp <<- pvals
            # read the principle (unshifted) results
            outfilePrincip = paste(bedstem, tagname, sep='.')
            message("Reading ", outfilePrincip, "...")
            datPrincip = readplinkout(outfilePrincip)
            chr[, 1] <<- datPrincip$CHR
            snp[, 1] <<- datPrincip$SNP
            bp[, 1] <<- datPrincip$BP
            pvals[, 1] <<- datPrincip$P

            # read shifted results
            outfiles = sapply(shiftFilesStem, function(eachstem) {
                paste(eachstem, tagname, sep = ".")
            })
            outfiles = setNames(outfiles, NULL)

            pcounter = 2
            for(outfile in outfiles) {
                message("Reading ", outfile, "...")
                datshift = readplinkout(outfile, "P")
                chr[, pcounter] <<- datPrincip$CHR
                snp[, pcounter] <<- datPrincip$SNP
                bp[, pcounter] <<- datPrincip$BP
                pvals[, pcounter] <<- datshift$P
                pcounter = pcounter + 1
            }

            matColnames = c("0", nshiftStrs)
            colnames(chr) <<- colnames(snp) <<- colnames(bp) <<- colnames(pvals) <<- matColnames
            message("Saving to ", outRdata, "...")
            save(chr, snp, bp, pvals, file=outRdata)
        }
    }
)

Plotcoll$methods(
    basemh = function() {
        baseplotObj = Mhplot(chr, bp, pvals[, 1])
        baseplot = baseplotObj$getmhplot()
        baseplot
    }
)

Plotcoll$methods(
    minpmh = function() {
        minpplotObj = Mhplot(chr, bp, do.call(pminNoNA, as.data.frame(pvals)))
        minpplot = minpplotObj$getmhplot()
        minpplot = minpplot + geom_hline(yintercept = -log10(5e-8 / nshiftMax), color="yellow", alpha=0.2)
        minpplot
    }
)

Plotcoll$methods(
    initialize = function(bed) {
        bedstem <<- getstem(bed)
        bedpath <<- paste(bedstem, "bed", sep=".")
        fampath <<- paste(bedstem, "fam", sep=".")
        bimpath <<- paste(bedstem, "bim", sep=".")
        shiftstemCommon <<- paste(bedstem, "_shift_", sep = "")
        updateShiftFilesStem()
        nsnp <<- R.utils::countLines(bimpath)
        nindiv <<- R.utils::countLines(fampath)
        nshiftMax <<- length(shiftFilesStem)
        pvals <<- matrix(1, 1, 1)
        chr <<- matrix(0)
        snp <<- matrix(0)
        bp <<- matrix(0)
    }
)

setwd("~/data/rs1exoseq/testReadout/")
plotobj = Plotcoll("ssknEx.bed")
plotobj$nshiftStrs
plotobj$shiftFilesStem
plotobj$readout("assoc.linear")
head(plotobj$snp)
head(plotobj$chr)
head(plotobj$bp)
head(plotobj$pvals)
require(ggplot2)
rs1exoBaseplot = plotobj$basemh()
rs1exoMinpplot = plotobj$minpmh()
ggsave(rs1exoBaseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/rs1exoBaseplotSskn.png", width = 9, height = 3)
ggsave(rs1exoMinpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/rs1exoMinpplotSskn.png", width = 9, height = 3)

##############

## setwd("~/data/rs1exoseq/AgeSexSskn")
## plotobj = Plotcoll("ssknEx.bed")
## plotobj$readout("assoc.linear")
## head(plotobj$snp)
## head(plotobj$chr)
## head(plotobj$bp)
## head(plotobj$pvals)
## require(ggplot2)
## rs1exoBaseplot = plotobj$basemh()
## rs1exoMinpplot = plotobj$minpmh()
## ggsave(rs1exoBaseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/rs1exoBaseplotSskn.png", width = 9, height = 3)
## ggsave(rs1exoMinpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/rs1exoMinpplotSskn.png", width = 9, height = 3)

## setwd("~/data/rs1exoseq/AgeSexHskn")
## plotobj = Plotcoll("ssknEx.bed")
## plotobj$readout("assoc.linear")
## head(plotobj$snp)
## head(plotobj$chr)
## head(plotobj$bp)
## head(plotobj$pvals)
## require(ggplot2)
## rs1exoBaseplot = plotobj$basemh()
## rs1exoMinpplot = plotobj$minpmh()
## ggsave(rs1exoBaseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/rs1exoBaseplotHskn.png", width = 9, height = 3)
## ggsave(rs1exoMinpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/rs1exoMinpplotHskn.png", width = 9, height = 3)

## setwd("~/data/sskn_regions_from_fan/AgeSexSskn/")
## plotobj = Plotcoll("sskn_reg.bed")
## plotobj$readout("assoc.linear")
## head(plotobj$snp)
## head(plotobj$chr)
## head(plotobj$bp)
## head(plotobj$pvals)
## require(ggplot2)
## rs1exoBaseplot = plotobj$basemh()
## rs1exoMinpplot = plotobj$minpmh()
## ggsave(rs1exoBaseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegBaseplotHskn.png", width = 9, height = 5)
## ggsave(rs1exoMinpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegMinpplotHskn.png", width = 9, height = 5)
## plotobj$bedpath
## plotobj$bedstem
## plotobj$shiftstemCommon
## plotobj$shiftFilesStem

## plotobj$bimpath
## plotobj$fampath
## plotobj$nsnp
## plotobj$nindiv
## plotobj$bytesSnp
## plotobj$bytesTotal
## plotobj$sayhi("Kaiyin")

