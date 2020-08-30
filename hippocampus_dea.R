#!/usr/bin/env Rscript
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(DESeq2))
suppressMessages(library(zinbwave))
suppressMessages(library(BiocParallel))
register(MulticoreParam(8))

# Load in-house functions
source("libraries.R")     

# Load full dataset and meta data
counts = read.delim("data/merged_counts.tsv", sep="\t", header=T, row.names=1)
meta = read.delim("data/meta.tsv", sep="\t", header=T, row.names=1)

# Genes as rows
counts = t(counts)

# keep only CA1, CA2-3 and DG
counts = counts[,meta$cluster %in% c("CA1", "CA2-3", "DG")]
meta = meta[colnames(counts),]
meta$cluster = droplevels(meta$cluster)

MIN_COUNT = 1
# Remove all genes with no/low expression (10% spots detected per cluster with a count >= 1)
#  CA1 CA2-3    DG  
#  238   200   179  
sub_ca1 = counts[,meta$cluster == "CA1"]
sub_ca23 = counts[,meta$cluster == "CA2-3"]
sub_dg = counts[,meta$cluster == "DG"]
to_keep_genes = unique(c(rownames(sub_ca1[rowSums(sub_ca1 >= MIN_COUNT) > 23,]),
                         rownames(sub_ca23[rowSums(sub_ca23 >= MIN_COUNT) > 20,]),
                         rownames(sub_dg[rowSums(sub_dg >= MIN_COUNT) > 17,])))
counts = counts[to_keep_genes,]
# Remove all spots with low number of genes (~1% genes detected with count >= 1)
counts = counts[,colSums(counts >= MIN_COUNT) > 200]

# Update the meta 
meta = meta[colnames(counts),]
meta = as.data.frame(meta)

REDUCE = FALSE
if (REDUCE) {
  # Keep only 600 spots per cluster (50 per animal) to reduce computational time
  to_keep_spots = list()
  num_spots = 50
  for(c in unique(meta$cluster)) {
    for(a in unique(meta$animal)) {
      rows = rownames(meta[meta$cluster == c & meta$animal == a,])
      to_keep_spots = c(to_keep_spots, sample(rows, min(num_spots, length(rows))))
    }
  }
  # Update counts and meta
  counts = counts[,unlist(to_keep_spots)]
  meta = meta[unlist(to_keep_spots),]
}

# Build models
dds = do_dea_deseq(counts, meta, formulaN="~ chip + condition", test="Wald")

# Test vs WT
write_deseq_results(dds, meta, "deseq", FDR_T=0.1, do_spatial=FALSE)

# Test vs AD
clusters = unique(meta$cluster)
for (c in clusters) {
  # ADP vs AD in cluster c
  res = results(dds, contrast=c("condition", paste(c, "ADPolB", sep="_"), paste(c, "AD", sep="_")),
                parallel=TRUE, alpha=0.1)
  res = lfcShrink(dds, contrast=c("condition", paste(c, "ADPolB", sep="_"), paste(c, "AD", sep="_")),
                  res=res, type="normal", parallel=TRUE)
  res = as.data.frame(res)
  res = na.omit(res[order(res$padj),])
  write.table(res, file=paste("deseq", c, "ADPolB_vs_AD.tsv", sep="_"), sep="\t")
}

# Save the workspace
save.image(file="hippocampus_dea.RData")
             
