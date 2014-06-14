Table = setRefClass("Table",
    fields = list(
        filename="character",
        ncol="numeric",
        nrow="numeric"
    )
)
Table$methods(
    read = function(selRow, header=TRUE, ...) {
        if(! ncol) {
            tmp = read.table(filename, header=header, nrow=2, ...)
            print(tmp)
            ncol <<- ncol(tmp)
        }
        if(! nrow) {
            nrow <<- R.utils::countLines(filename)
            # trailing newlines would make nrow incorrect, make sure they are not there
            lastLine = scan(filename, what="", skip=nrow-1, quiet=TRUE)
            if(! length(lastLine)) {
                stop("Trailing white space at the end of file!")
            }
            if(header) {
                nrow <<- nrow - 1
            }
        }

        selRowLen = length(selRow)
        if(selRowLen == nrow) {
            selRowFinal = selRow
        } else {
            if(nrow %% selRowLen == 0) {
                selRowFinal = rep(selRow, nrow/selRowLen)
            } else {
                stop("Number of rows is not a multiple of length of selRow!")
            }
        }

        con = file(filename)
        open(con)
        dfColnames = colnames(tmp)
        firstLine = scan(con, what="", nlines=1, quiet=TRUE)
        print(firstLine)
        df = lapply(selRowFinal, function(sel) {
            line = read.table(con, header=FALSE, nrow=1, ...)
            if(sel == TRUE) line else NULL
        })
        close(con)
        df = do.call(rbind, df)
        df = setNames(df, dfColnames)
        df
    }
)
Table$methods(
    initialize = function(filenameinit) {
        filename <<- filenameinit
        ncol <<- 0
        nrow <<- 0
    }
)

x = Table("/tmp/test.assoc")
x$read(selRow = c(1, 0, 0), colClasses=c("numeric", rep("NULL", 8)))




