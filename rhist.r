testset = readplinkout("../ssknEx.assoc.linear")
smalltest = testset[which(testset$P < 1e-1 & testset$CHR != 1), ]
require(ggplot2)
x = Mhplot(smalltest$CHR, smalltest$BP, smalltest$P/500)$getmhplot()
print(x)
## head(x)


## x$chr
## x$bp
## x$achr
## x$sbp

## plot(x$chr, x$sbp %% 1)
## plot(x$chr, x$sbp)
## table(sort(smalltest$CHR) == smalltest$CHR)
## length(sort(smalltest$CHR))
## length(smalltest$CHR)
## table(smalltest$CHR, useNA = "always")

    

## testset = scaleByChr(testset)
## head(testset)
## pvalset = testset[testset$P < 1e-4, ]
## myplot = ggplot(pvalset, aes(x=SBP, y=BPDIST)) + geom_point(aes(size=P, alpha=1/P^2, color=BPDIST)) + scale_size_continuous(range=c(1, 10))
## ggsave(file = "/tmp/myplot.pdf", plot=myplot, width = 9, height = 3)


## head(testset)
## infodat = testset[, c("CHR", "SNP", "BP", "SBP")]
## class(infodat)
## head(infodat)
## nsnp = nrow(infodat)
## nsnp
## shift1fn = "ssknEx_shift_0001.assoc.linear"
## nshift = as.integer(gsub(".*_shift_(\\d{4}).*", "\\1", shift1fn))
## nshift
## shift1 = readplinkout(shift1fn, shift=TRUE)  
## head(shift1)
## dim(shift1)
## infoShift1.1 = infodat[1:(nsnp - nshift), ]
## head(infoShift1.1)
## dim(infoShift1.1)
## infoShift1.2 = infodat[(1+nshift):nsnp, c("SNP", "BP")]
## setnames(infoShift1.2, c("SNP2", "BP2"))
## head(infoShift1.2)
## dim(infoShift1.2)

## shift1 = cbind(infoShift1.1, infoShift1.2, shift1)
## shift1$BPDIST = shift1$BP2 - shift1$BP
## testset = rbind(testset, shift1)
## maxbpdist = testset$BPDIST[nrow(testset)]
## testset$BPDIST[testset$BPDIST == 0] = -maxbpdist
## testset$BPDIST[nrow(testset)]
## head(testset)

## pvalset = testset[testset$P < 1e-4, ]
## myplot = ggplot(pvalset, aes(x=SBP, y=BPDIST)) + geom_point(aes(size=P, alpha=1/P^2, color=BPDIST)) + scale_size_continuous(range=c(1, 10))
## ggsave(file = "/tmp/myplot.pdf", plot=myplot, width = 9, height = 3)

