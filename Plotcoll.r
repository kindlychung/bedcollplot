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
        chr = "numeric",
        snp = "character",
        bp = "numeric",
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
    gettag = function(tagname) {
        tag = gsub("(\\d{8}-\\d{6}).*", "\\1", tagname)
    }
)
Plotcoll$methods(
    readout = function(tagname) {
        outRdata = paste(bedstem, "RData", sep = ".")
        if(file.exists(outRdata)) {
            message("Reading from previously generated data...")
            load(outRdata)
            chr <<- chr
            snp <<- snp
            bp <<- bp
            pvals <<- pvals
        } else {
            # clean up pvals
            pvals <<- matrix(NA, nsnp, nshiftMax+1)
            # read the principle (unshifted) results
            outfilePrincip = paste(bedstem, tagname, sep='.')
            message("Reading ", outfilePrincip, "...")
            if(chr[1] == 0 | snp[1] == "0" | bp[1] == 0) {
                datPrincip = readplinkout(outfilePrincip)
                chr <<- datPrincip$CHR
                snp <<- datPrincip$SNP
                bp <<- datPrincip$BP
            } else {
                datPrincip = readplinkout(outfilePrincip, "P")
            }

            # todo: use bigmatrix!
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
                pvals[, pcounter] <<- datshift$P
                pcounter = pcounter + 1
            }

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
        shiftFilesStem <<- Sys.glob(paste(bedstem, "_shift_*.bed", sep=""))
        shiftFilesStem <<- sapply(shiftFilesStem, getstem)
        names(shiftFilesStem) <<- NULL
        nsnp <<- R.utils::countLines(bimpath)
        nindiv <<- R.utils::countLines(fampath)
        nshiftMax <<- length(shiftFilesStem)
        pvals <<- matrix(1, 1, 1)
        chr <<- 0
        snp <<- "0"
        bp <<- 0
    }
)

setwd("~/data/rs1exoseq/AgeSexSskn")
plotobj = Plotcoll("ssknEx.bed")
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

setwd("~/data/rs1exoseq/AgeSexHskn")
plotobj = Plotcoll("ssknEx.bed")
plotobj$readout("assoc.linear")
head(plotobj$snp)
head(plotobj$chr)
head(plotobj$bp)
head(plotobj$pvals)
require(ggplot2)
rs1exoBaseplot = plotobj$basemh()
rs1exoMinpplot = plotobj$minpmh()
ggsave(rs1exoBaseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/rs1exoBaseplotHskn.png", width = 9, height = 3)
ggsave(rs1exoMinpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/rs1exoMinpplotHskn.png", width = 9, height = 3)

setwd("~/data/sskn_regions_from_fan/AgeSexSskn/")
plotobj = Plotcoll("sskn_reg.bed")
plotobj$readout("assoc.linear")
head(plotobj$snp)
head(plotobj$chr)
head(plotobj$bp)
head(plotobj$pvals)
require(ggplot2)
rs1exoBaseplot = plotobj$basemh()
rs1exoMinpplot = plotobj$minpmh()
ggsave(rs1exoBaseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegBaseplotHskn.png", width = 9, height = 5)
ggsave(rs1exoMinpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegMinpplotHskn.png", width = 9, height = 5)
plotobj$bedpath
plotobj$bedstem
plotobj$shiftstemCommon
plotobj$shiftFilesStem
## plotobj$bimpath
## plotobj$fampath
## plotobj$nsnp
## plotobj$nindiv
## plotobj$bytesSnp
## plotobj$bytesTotal
## plotobj$sayhi("Kaiyin")
