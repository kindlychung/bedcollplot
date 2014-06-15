readplinkout = function(filename, colnameSelect=c("CHR", "SNP", "BP", "P")) {
    ## colnameSelect=c("CHR", "SNP", "BP", "P", "TEST")
    ## filename = "../ssknEx.assoc.linear"
    ## colClsAll = c("numeric", "character", "numeric", "character", "character", rep("numeric", 4))
    ## names(colClsAll) = c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P")

    tmp1 = read.table(filename, header=TRUE, nrows=2)
    cnames = colnames(tmp1)
    if(! all(colnameSelect %in% cnames)) {
        stop("Some requested columns unavailable!")
    }

    testIdx = which(cnames == "TEST")
    filtercmdTest = paste("sed -r 's/^ +//;s/ +/\t/g' ",
        filename, " | cut -f", 
        testIdx, " | tail -n +2",
        sep=""
    )
    filtercmdTest
    alltests = read.table(pipe(filtercmdTest), header = FALSE)
    ncovar = nrow(unique(alltests))

    colSelectIdx = which(cnames %in% colnameSelect)
    colSelString = paste(colSelectIdx, collapse = ",")

    filtercmd = paste("sed -r 's/^ +//;s/ +/\t/g' ",
        filename, " | cut -f", 
        colSelString, " | sed -n '2~", ncovar, "p'",
        sep=""
    )
    filtercmd

    ## filtercmd
    plinkRes0 = read.table(pipe(filtercmd), header = FALSE, stringsAsFactors = FALSE)
    plinkRes0 = setNames(plinkRes0, cnames[colSelectIdx])
    head(plinkRes0)
    plinkRes0
}
