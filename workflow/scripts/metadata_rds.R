#!/usr/bin/Rscript

intxt <- unlist(snakemake@input)
outrds <- unlist(snakemake@output)

### Get names for input and output file vectors
names(intxt) <-  gsub('metadata\\.([^\\.]+)\\.txt.gz', '\\1', basename(intxt))
names(outrds) <- gsub('metadata\\.([^\\.]+)\\.rds', '\\1', basename(outrds))
stopifnot(all(names(intxt) == names(outrds)))


### Read each table and add to list
###     Feature tables (gene_features, tx_features, etc.) contain headers;
###     otherwise it is a 2-column mapping table and the variable names are
###     indicated in the filename: i.e. tid_gid has columns for tid and gid.
annot.tables <- lapply(intxt, function(f) {
    tname <- basename(f)
    splitname <- unlist(strsplit(tname, '_'))
    is_ftable <- splitname[2] == 'features'
    tbl <- read.table(f, sep='\t', header=is_ftable)
    if(!is_ftable) names(tbl) <- splitname
    if(!any(duplicated(tbl[,1]))) {
        row.names(tbl) <- tbl[,1]
        tbl[,1] <- NULL
    }
    tbl
})

### Save RDS files
ret <- lapply(names(annot.tables), function(tname) {
    saveRDS(annot.tables[[tname]], file=outrds[tname])
})

