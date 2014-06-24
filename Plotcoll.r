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
        nshift = "numeric",
        nshiftMax = "numeric",
        nshiftStrs = "character",
        chr = "matrix",
        snp = "matrix",
        bp = "matrix",
        chr2 = "matrix",
        snp2 = "matrix",
        bp2 = "matrix",
        pvals = "matrix",
        minPvals = "numeric",
        minPvalsBonfer = "numeric",
        ntests = "numeric", 
        chrunique = "numeric"
    ))

Plotcoll$methods(
    getNshiftStr = function(shiftpath) {
        nshiftStr = gsub(".*_shift_(\\d{4}).*", "\\1", shiftpath)
        nshiftStr
    }
)
Plotcoll$methods(
    getNshift = function(shiftpath) {
        nshiftRes = as.integer(getNshiftStr(shiftpath))
        nshiftRes
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
        nshift <<- c(0, getNshift(shiftFilesStem))
        nshiftStrs <<- as.character(nshift)
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


        if(chrunique[1] == 0) {
            message("Nothing in the set of unique CHRs, calculating it...")
            chrunique <<- unique(chr[, 1])
        } else {
            message("Unique set of CHRs:")
            print(chrunique)
        }

        minPvals <<- do.call(pminNoNA, as.data.frame(pvals))
        
        ntests <<- rep(0, nsnp)
        # start with a reversed vector with inversed signs, see http://tinyurl.com/lnm48hj
        for(i in nshift) {
            ntests[1:i] <<- ntests[1:i] - 1
        }
        ntests <<- ntests - (ntests[1] - 1)
        # reverse it back
        ntests <<- rev(ntests)

        # Bonferroni correction for minimal p values
        minPvalsBonfer <<- minPvals * ntests
        minPvalsBonfer <<- ifelse(minPvalsBonfer > 1, 1, minPvalsBonfer)

        message("Update info for SNP2...")
        snp2info()
    }
)

Plotcoll$methods(
    basemh = function(chrfilter=NULL, bplower=NULL, bpupper=NULL, pvallower=NULL, pvalupper=NULL) {
        filter = rep(TRUE, nsnp)
        if(! is.null(chrfilter)) {
            filter = chr[, 1] %in% chrfilter
        }
        if(! is.null(bplower)) {
            filter = filter & bp[, 1] >= bplower
        }
        if(! is.null(bpupper)) {
            filter = filter & bp[, 1] <= bpupper
        }
        if(! is.null(pvallower)) {
            filter = filter & pvals[, 1] > pvallower
        }
        if(! is.null(pvalupper)) {
            filter = filter & pvals[, 1] < pvalupper
        }
        baseplotObj = Mhplot(chr[filter, 1], bp[filter, 1], pvals[filter, 1])
        baseplot = baseplotObj$getmhplot()
        baseplot
    }
)

Plotcoll$methods(
    minpmh = function(
        chrfilter=NULL, bplower=NULL, bpupper=NULL,
        pvallower=NULL, pvalupper=NULL, minpcorrect=TRUE
    ) {
        filter = rep(TRUE, nsnp)
        if(! is.null(chrfilter)) {
            filter = chr[, 1] %in% chrfilter
        }
        if(! is.null(bplower)) {
            filter = filter & bp[, 1] >= bplower
        }
        if(! is.null(bpupper)) {
            filter = filter & bp[, 1] <= bpupper
        }
        if(! is.null(pvallower)) {
            filter = filter & pvals[, 1] > pvallower
        }
        if(! is.null(pvalupper)) {
            filter = filter & pvals[, 1] < pvalupper
        }

        if(minpcorrect == TRUE) {
            minpplotObj = Mhplot(chr[filter, 1], bp[filter, 1], minPvalsBonfer[filter])
        } else {
            minpplotObj = Mhplot(chr[filter, 1], bp[filter, 1], minPvals[filter])
        }
        minpplot = minpplotObj$getmhplot()
        minpplot
    }
)

Plotcoll$methods(
    contrastplot = function(chrfilter=NULL, bplower=NULL,
        bpupper=NULL, pvallower=NULL, pvalupper=NULL, minpcorrect=TRUE
    ) {
        filter = rep(TRUE, nsnp)
        if(! is.null(chrfilter)) {
            filter = chr[, 1] %in% chrfilter
        }
        if(! is.null(bplower)) {
            filter = filter & bp[, 1] >= bplower
        }
        if(! is.null(bpupper)) {
            filter = filter & bp[, 1] <= bpupper
        }
        if(! is.null(pvallower)) {
            filter = filter & pvals[, 1] > pvallower
        }
        if(! is.null(pvalupper)) {
            filter = filter & pvals[, 1] < pvalupper
        }
        filter = which(filter)

        # combine the variables
        bothChr = c(chr[filter, 1], chr[filter, 1])
        bothBp = c(bp[filter, 1], bp[filter, 1])
        if(minpcorrect == TRUE) {
            bothP = c(pvals[filter, 1], minPvalsBonfer[filter])
        } else {
            bothP = c(pvals[filter, 1], minPvals[filter])
        }
        colorvec = rep(c("Single SNP", "QCDH"), each=length(filter))

        # get the order right
        posorder = order(bothChr, bothBp)
        bothChr = bothChr[posorder]
        bothBp  = bothBp[posorder]
        bothP   = bothP[posorder]
        colorvec = colorvec[posorder]

        Mhplot(bothChr, bothBp, bothP, colorvec=colorvec)$getmhplot()
    }
)

Plotcoll$methods(
    basepMinp = function(logpvals=TRUE, minpcorrect=TRUE) {
        df = data.frame(basep = pvals[, 1])
        if(minpcorrect == TRUE) {
            df$minp = minPvalsBonfer
        } else {
            df$minp = minPvals
        }

        if(logpvals == TRUE) {
            df$minp = -log10(df$minp)
            df$basep = -log10(df$basep)
        }

        bmplot = ggplot(df, aes(basep, minp)) + geom_point(alpha=0.4)
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
        chrunique <<- 0
        minPvalsBonfer <<- 1
        minPvals <<- 1
    }
)



##1######################################
datadirs = Sys.glob("~/data/sskn_regions_from_fan/rand*")
for(datadir in datadirs) {
    ## datadir = "~/data/sskn_regions_from_fan/rand.1"
    setwd(datadir)
    plotobj = Plotcoll("sskn_reg.bed")
    plotobj$readout("assoc.linear")
    outf = paste0("/tmp/", basename(datadir), ".png")
    message(sprintf("Plot base and min p values to %s", outf))
    bmplot = plotobj$basepMinp(minpcorrect = FALSE)
    ggsave(bmplot, filename = outf, width = 10, height = 5, dpi = 600)
}

setwd("~/data/sskn_regions_from_fan/AgeSexSskn/")
plotobj = Plotcoll("sskn_reg.bed")
plotobj$readout("assoc.linear")

plotobj$nsnp
plotobj$nshift
plotobj$nshiftMax
plotobj$nshiftStrs
plotobj$chrunique

x = plotobj$basepMinp(minpcorrect = FALSE)
print(x)

## datForLocusZoom = data.frame(CHR=plotobj$chr[, 1], SNP=plotobj$snp[, 1], BP=plotobj$bp[, 1], P=plotobj$minPvals)
## for(nchr in unique(datForLocusZoom$CHR)) {
##     tmpdat = datForLocusZoom[which(datForLocusZoom$CHR == nchr), ]
##     print(nchr)
##     print(range(tmpdat$BP))
##     print(range(tmpdat$P, na.rm = TRUE))
## }
## write.table(datForLocusZoom, file = "datForLocusZoom.txt", col.names=TRUE, quote=FALSE, row.names=FALSE)

require(ggplot2)

## for(i in plotobj$chrunique) {
##     myplot = plotobj$contrastplot(chrfilter = i)
##     outf = paste0("/tmp/", i, ".png")
##     ggsave(myplot, filename = outf, width = 10, height = 5, dpi = 600)
## }

ggsave(c1, filename = "/tmp/c1.png", width = 10, height = 5, dpi = 600)
ggsave(c2, filename = "/tmp/c2.png", width = 10, height = 5, dpi = 600)
Baseplot = plotobj$basemh()
print(Baseplot)
Minpplot = plotobj$minpmh()
print(Minpplot)
ggsave(Baseplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegSsknBaseplot.pdf", width = 10, height = 5)
ggsave(Minpplot, file="/home/kaiyin/projManuscripts/qcdh/figs/ssknRegSsknMinpplot.pdf", width = 10, height = 5)
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
