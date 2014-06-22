Mhplot = setRefClass(
    "Mhplot",
    fields = list(
    chr="numeric",
    snp="character",
    bp="numeric",
    pvals="numeric",
    nsnp="numeric",
    mlogp="numeric",
    nchr="numeric",
    chrunique="numeric",
    # apperant chr, see function below
    achr="numeric",
    # base-pair position scaled within chr
    sbp="numeric"
    # todo: base-pair position scaled to real length of chromosome!
    )
)

Mhplot$methods(
    # apparent chr number
    # e.g. chr(1, 1, 4, 4, 5, 5, 7, 7) ==> chr(1, 1, 2, 2, 3, 3, 4, 4)
    # usefu for plotting
    appChr = function() {
        achr <<- rep(0, nchr)
        for(i in chrunique) {
            achr <<- achr + (chr >= i)
        }
    }
)

Mhplot$methods(
    # sort out the chr stuf, scale by each chr
    # notice I used the apparent chr here
    scaleByChr = function() {
        allchrScaledPos = rep(NA, nsnp)
        for (chrI in chrunique) {
            chrcheck = (chrI == chr)
            chrNsnp = sum(chrcheck)
            firstPos = bp[which(chrcheck)[1]]
            firstPos = rep(firstPos, chrNsnp)
            posDiff = bp[chrcheck] - firstPos
            # make sure it's scaled to the range (0, 1)
            scaledPos = posDiff / (posDiff[chrNsnp] + 0.5)
            allchrScaledPos[chrcheck] = scaledPos
        }
        sbp <<- achr + allchrScaledPos
    }
)

Mhplot$methods(
    getdf = function() {
        df = data.frame(chr, achr, bp, sbp, mlogp)
        df
    }
)

Mhplot$methods(
    getmhplot = function(annotation=NULL, chrColor=TRUE) {
        mydat = getdf()
        maxlogp = ceiling(max(mlogp, na.rm = TRUE))
        minlogp = min(mlogp, na.rm = TRUE)
        gwthresh = -log10(5e-8)
        
        # if there are more than one chr, then x axis should be labeled by scaled BP (sbp)
        # if there is only one chr, then x axis should be labeled by BP
        if(length(unique(chr)) > 1) {
            myplot = ggplot(mydat, aes(sbp, mlogp)) + xlab("CHR") + 
                     scale_x_continuous(breaks=unique(achr),
                                        minor_breaks=NULL,
                                        labels=unique(chr))
        } else {
            myplot = ggplot(mydat, aes(bp, mlogp)) + xlab(paste("Position on CHR", chr[1]))
            chrColor = FALSE
        }


        if(chrColor == TRUE) {
            myplot = myplot + geom_point(aes(color=factor(achr%%2))) + 
                scale_color_manual(values = c("gray20", "gray60"), guide=FALSE) 
        } else if (chrColor == FALSE) {
            myplot = myplot + geom_point()
        } else {
            myplot = myplot + geom_point(color=chrColor)
        }

        myplot = myplot + scale_y_continuous(limits=c(minlogp, maxlogp), minor_breaks=NULL) +
                geom_hline(yintercept=gwthresh, alpha=.3, color="blue") + 
                ylab("-log P")

        # todo: annotation
        if(!is.null(annotation)) {
            if(snp == "0") {
                stop("You must initialize me with an vector of SNP names if you want to annotate!")
            } else {
            }
        }

        return(myplot)
    }
)


Mhplot$methods(
    initialize = function(chrinit, bpinit, pvalsinit, snpinit="0") {
        if(any(sort(chrinit) != chrinit)) {
            message("I will not sort CHR for you, since this might mess up with the other columns of your data!")
            stop("CHR must be in order (1,2,3...22, for example)!")
        }
        nsnp <<- length(bpinit)
        if(length(chrinit) != nsnp | length(pvalsinit) != nsnp) {
            stop("CHR, BP and P do not match in length!")
        }
        snp <<- snpinit
        chr <<- chrinit
        bp <<-bpinit
        pvals <<- pvalsinit
        mlogp <<- -log10(pvals)
        nchr <<- length(chr)
        chrunique <<- unique(chr)
        appChr()
        scaleByChr()
    }
)



## require(ggplot2)
## plinkout = readplinkout("~/data/sskn_regions_from_fan/AgeSexRed/sskn_reg.assoc.logistic")
## ## plinkout = plinkout[which(plinkout$CHR == 20), ]
## plinkout = plinkout[which(plinkout$CHR == 16), ]
## head(plinkout)
## plinkplotObj = Mhplot(plinkout$CHR, plinkout$BP, plinkout$P)
## plinkplot = plinkplotObj$getmhplot()
## print(plinkplot)
## plinkplot = plinkplot + geom_point(size=8)
## print(plinkplot)
