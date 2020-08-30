#!/usr/bin/env Rscript
suppressMessages(library(scran))
suppressMessages(library(scater))
suppressMessages(library(DESeq2))
suppressMessages(library(BiocParallel))
register(MulticoreParam(10))

# Load in-house functions
source("libraries.R")

args = commandArgs(trailingOnly=TRUE)
name_prefix = args[1]

counts = read.delim("data/merged_counts.tsv", sep="\t", header=T, row.names=1)
meta = read.delim("data/meta.tsv", sep="\t", header=T, row.names=1)

# Genes as rows
counts = t(counts)

# Remove ONL spots since we are not interested in them 
counts = counts[,!meta$cluster == "ONL"]
meta = meta[colnames(counts),]
meta$cluster = droplevels(meta$cluster)

MIN_COUNT = 1
# Remove all genes with no/low expression (10% spots detected per cluster with a count >= 1)
# EPL  GCL   GL  MCL 
# 1631  992 1285  698 
sub_gcl = counts[,meta$cluster == "GCL"]
sub_mcl = counts[,meta$cluster == "MCL"]
sub_epl = counts[,meta$cluster == "EPL"]
sub_gl = counts[,meta$cluster == "GL"]
to_keep_genes = unique(c(rownames(sub_gcl[rowSums(sub_gcl >= MIN_COUNT) > 90,]),
                         rownames(sub_mcl[rowSums(sub_mcl >= MIN_COUNT) > 70,]),
                         rownames(sub_epl[rowSums(sub_epl >= MIN_COUNT) > 160,]),
                         rownames(sub_epl[rowSums(sub_gl >= MIN_COUNT) > 120,])))
counts = counts[to_keep_genes,]
# Remove all spots with low number of genes (~1% genes detected with count >= 1)
counts = counts[,colSums(counts >= MIN_COUNT) > 200]

# Update the meta 
meta = meta[colnames(counts),]
meta = as.data.frame(meta)

REDUCE = FALSE
if (REDUCE) {
  # Keep only 240 spots per cluster (10 per animal/genotype/cluster) to reduce computational time
  to_keep_spots = list()
  num_spots = 10
  for(c in unique(meta$cluster)) {
    for(g in unique(meta$genotype)) {
      for(a in unique(meta$animal)) {
        rows = rownames(meta[meta$cluster == c & meta$genotype == g & meta$animal == a,])
        if (length(rows) > 0) {
            to_keep_spots = c(to_keep_spots, sample(rows, min(num_spots, length(rows))))
        }
      }
    }
  }
  # Update counts and meta
  counts = counts[,unlist(to_keep_spots)]
  meta = meta[unlist(to_keep_spots),]
}


# Perform the DEA
dds = do_dea_deseq(counts, meta, test="Wald")
write_deseq_results(dds, meta, paste0(name_prefix, "_deseq"), 0.1)

# Save the workspace
save.image(file=paste0(name_prefix, ".RData"))
             
