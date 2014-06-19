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
        chr2 = "matrix",
        snp2 = "matrix",
        bp2 = "matrix",
        pvals = "matrix",
        minPvals = "numeric"
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
        nshiftStrs <<- c("0", nshiftStrs)
        shiftFilesStem <<- c(bedstem, shiftFilesStem)
    }
)

Plotcoll$methods(
    shiftDataByCol = function(dat) {
        nr = nrow(dat)
        nc = ncol(dat)
        newdat = matrix(NA, nr, nc)
        for(i in 1:nc) {
            colShiftN = as.integer(colnames(dat)[i])
            newdat[1:(nr-colShiftN), i] = dat[(colShiftN+1):nr, 1]
        }
        colnames(newdat) = colnames(dat)
        newdat
    }
)

Plotcoll$methods(
    snp2info = function() {
        chr2 <<- shiftDataByCol(chr)
        snp2 <<- shiftDataByCol(snp)
        bp2 <<-shiftDataByCol(bp)
    }
)

Plotcoll$methods(
    readout = function(tagname) {
        outRdata = paste(bedstem, "RData", sep = ".")
        if(file.exists(outRdata)) {
            message("Reading from previously generated data...")
            load(outRdata)
            currentEnv = environment()

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
            # add new files
            hasColnames = colnames(chr)
            addIdx = which(!(nshiftStrs %in% hasColnames))
            addColnames = nshiftStrs[addIdx]
            addFiles = shiftFilesStem[addIdx]

            if(length(addFiles) > 0) {
                fileAdded = TRUE
                message("It looks like you have added some new files, let me read them: ")
                print(addFiles)

                pcounter = ncol(pvals) + 1
                for(i in 1:length(addFiles)) {
                    addf = addFiles[i]
                    addf = paste(addf, tagname, sep = ".") 
                    addcoln = addColnames[i]
                    message("Reading ", addf, "...")
                    datshift = readplinkout(addf)
                    currentEnv$chr = cbind(currentEnv$chr, datshift$CHR)
                    currentEnv$snp = cbind(currentEnv$snp, datshift$SNP)
                    currentEnv$bp = cbind(currentEnv$bp, datshift$BP)
                    currentEnv$pvals = cbind(currentEnv$pvals, datshift$P)
                    colnames(currentEnv$chr)[pcounter] = addcoln
                    colnames(currentEnv$snp)[pcounter] = addcoln
                    colnames(currentEnv$bp)[pcounter] = addcoln
                    colnames(currentEnv$pvals)[pcounter] = addcoln
                    pcounter = pcounter + 1
                }
            } else {
                fileAdded = FALSE
            }


            # remove colums if corresponding files are removed
            # e.g. if x_shift_0002.bed is removed, then remove column "2"
            # first you need to update hasColnames!
            hasColnames = colnames(chr)
            selectIdx = which(hasColnames %in% nshiftStrs)
            rmIdx = which(!(hasColnames %in% nshiftStrs))
            if(length(selectIdx) != length(hasColnames)) {
                fileRemoved = TRUE
                rmColname = hasColnames[rmIdx]
                rmFiles = sapply(rmColname, function(col) {
                    sprintf("%s%04d", shiftstemCommon, as.integer(col))
                })
                names(rmFiles) = NULL
                message("It looks like you have removed some file(s), let me delete corresponding data:")
                print(rmFiles)
                currentEnv$chr = currentEnv$chr[, selectIdx]
                currentEnv$snp = currentEnv$snp[, selectIdx]
                currentEnv$bp = currentEnv$bp[, selectIdx]
                currentEnv$pvals = currentEnv$pvals[, selectIdx]
            } else {
                fileRemoved = FALSE
            }

            # assign into the class!
            chr <<- chr
            snp <<- snp
            bp <<- bp
            pvals <<- pvals

            if(fileAdded | fileRemoved) {
                message("Saving changes back to RData file...")
                save(chr, snp, bp, pvals, file=outRdata)
            } else {
                message("I did not see any change, keep the RData file as it is.")
            }
        } else {
            # set up matrices
            pvalsNcols = length(nshiftStrs)
            pvals <<- matrix(NA, nsnp, pvalsNcols)
            chr <<- pvals
            snp <<- pvals
            bp <<- pvals

            # read shifted results
            outfiles = sapply(shiftFilesStem, function(eachstem) {
                paste(eachstem, tagname, sep = ".")
            })
            outfiles = setNames(outfiles, NULL)

            for(i in 1:length(outfiles)) {
                outfile = outfiles[i]
                coln = nshiftStrs[i]
                message("Reading ", outfile, "...")
                datshift = readplinkout(outfile)
                chr[, i] <<- datshift$CHR
                snp[, i] <<- datshift$SNP
                bp[, i] <<- datshift$BP
                pvals[, i] <<- datshift$P
            }

            colnames(chr) <<- colnames(snp) <<- colnames(bp) <<- colnames(pvals) <<- nshiftStrs
            message("Saving to ", outRdata, "...")
            save(chr, snp, bp, pvals, file=outRdata)
        }
        minPvals <<- do.call(pminNoNA, as.data.frame(pvals))
        message("Update info for SNP2...")
        snp2info()
    }
)

Plotcoll$methods(
    basemh = function() {
        baseplotObj = Mhplot(chr[, 1], bp[, 1], pvals[, 1])
        baseplot = baseplotObj$getmhplot()
        baseplot
    }
)

Plotcoll$methods(
    minpmh = function() {
        minpplotObj = Mhplot(chr[, 1], bp[, 1], minPvals)
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
        chr <<- matrix(, 0, 0)
        snp <<- matrix(, 0, 0)
        bp <<- matrix(, 0, 0)
        chr2 <<- matrix(, 0, 0)
        snp2 <<- matrix(, 0, 0)
        bp2 <<- matrix(, 0, 0)
    }
)



##1######################################
setwd("~/data/sskn_regions_from_fan/AgeSexSskn/")
plotobj = Plotcoll("sskn_reg.bed")
plotobj$readout("assoc.linear")
datForLocusZoom = data.frame(CHR=plotobj$chr[, 1], SNP=plotobj$snp[, 1], BP=plotobj$bp[, 1], P=plotobj$minPvals)

for(nchr in unique(datForLocusZoom$CHR)) {
    tmpdat = datForLocusZoom[which(datForLocusZoom$CHR == nchr), ]
    print(nchr)
    print(range(tmpdat$BP))
    print(range(tmpdat$P, na.rm = TRUE))
}

write.table(datForLocusZoom, file = "datForLocusZoom.txt", col.names=TRUE, quote=FALSE, row.names=FALSE)

require(ggplot2)
## Baseplot = plotobj$basemh() + geom_point(size=19)
## Minpplot = plotobj$minpmh() + geom_point(size=19)
Baseplot = plotobj$basemh()
Minpplot = plotobj$minpmh()
print(Baseplot)
print(Minpplot)
ggsave(rs1exoBaseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegSsknBaseplot.png", width = 9, height = 5)
ggsave(rs1exoMinpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegSsknMinpplot.png", width = 9, height = 5)
##2######################################

## setwd("~/data/sskn_regions_from_fan/AgeSexRed/")
## plotobj = Plotcoll("sskn_reg.bed")
## plotobj$readout("assoc.logistic")
## require(ggplot2)
## Baseplot = plotobj$basemh() + geom_point(size=15)
## Minpplot = plotobj$minpmh() + geom_point(size=15)
## ## png("/tmp/x.png", width = 5000, height = 2500)
## print(Baseplot)
## ## dev.off()
## print(Minpplot)
## ggsave(Baseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegRedBaseplot.png", width = 9, height = 5, scale=5)
## ggsave(Minpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegRedMinpplot.png", width = 9, height = 5, scale=5)

## plotobj$bimpath
## plotobj$fampath
## plotobj$nsnp
## plotobj$nindiv
## plotobj$bytesSnp
## plotobj$bytesTotal
## plotobj$sayhi("Kaiyin")


## ?png
## ?ggsave


## setwd("~/data/rs1exoseq/testReadout/")
## plotobj = Plotcoll("ssknEx.bed")
## plotobj$nshiftStrs
## plotobj$shiftFilesStem
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
