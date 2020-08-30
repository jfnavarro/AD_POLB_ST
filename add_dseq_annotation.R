#!/usr/bin/env Rscript
# give paths to DESeq output tables

library(org.Mm.eg.db)

addExtraAnnotation <- function(table) {
  table$ensembl = mapIds(org.Mm.eg.db, keys=rownames(table), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
  table$name = mapIds(org.Mm.eg.db, keys=rownames(table), column="GENENAME", keytype="SYMBOL", multiVals="first")
  table$entrez = mapIds(org.Mm.eg.db, keys=rownames(table), column="ENTREZID", keytype="SYMBOL", multiVals="first")
  return(table)
}

args = commandArgs(trailingOnly=T)

for (path in args) {
    x = read.table(path, sep="\t", header=T, row.names=1)
    x = addExtraAnnotation(x)
    write.table(x, file=paste("annotated_", path, sep=""), sep="\t", row.names=T, col.names=T)
}
