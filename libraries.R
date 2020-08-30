# A set of functions to generate figures and analysis
library(enrichR)
library(pheatmap)
library(RColorBrewer)
library(VennDiagram)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DESeq2)
library(plyr)
library(dplyr)
library(AUCell)
library(GSEABase)
library(GSVA)
library(fgsea)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(stringr)
library(SingleR)
library(qdapTools)
library(data.table)
library(repr)
library(ggpubr)

# The EnrichR databases
dbs = c("GO_Biological_Process_2018", "KEGG_2019_Mouse", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")

doEnrichR <- function(genes, substring, min_fdr=0.1, min_pvalue=0.1, min_genes=3) {
  enriched = enrichr(genes, dbs)
  
  go_biological = enriched[[1]]
  go_biological = go_biological[which(go_biological$Adjusted.P.value <= min_fdr & 
                                        go_biological$P.value <= min_pvalue & 
                                        str_count(go_biological$Genes, ";") >= min_genes - 1), 
                                c(1,2,3,4,7,8,9)]
  if (dim(go_biological)[1] > 1) {
    write.table(go_biological, file=paste(substring, "GO_biological", ".tsv", sep="_"), 
                sep="\t", row.names=F, col.names=T)
  }
  
  kegg = enriched[[2]]
  kegg = kegg[which(kegg$Adjusted.P.value < min_fdr & 
                      kegg$P.value <= min_pvalue & 
                      str_count(kegg$Genes, ";") >= min_genes - 1), 
              c(1,2,3,4,8,9)]
  if (dim(kegg)[1] > 1) {
    write.table(kegg, file=paste(substring, "KEGG", ".tsv", sep="_"), 
                sep="\t", row.names=F, col.names=T)
  }
  
  go_molecular = enriched[[3]]
  go_molecular = go_molecular[which(go_molecular$Adjusted.P.value< min_fdr & 
                                      go_molecular$P.value <= min_pvalue & 
                                      str_count(go_molecular$Genes, ";") >= min_genes - 1),
                              c(1,2,3,4,7,8,9)]
  if (dim(go_molecular)[1] > 1) {
    write.table(go_molecular, file=paste(substring, "GO_molecular", ".tsv", sep="_"), 
                sep="\t", row.names=F, col.names=T)
  }
  
  go_cellular = enriched[[4]]
  go_cellular = go_cellular[which(go_cellular$Adjusted.P.value < min_fdr & 
                                    go_cellular$P.value <= min_pvalue & 
                                    str_count(go_cellular$Genes, ";") >= min_genes - 1),
                            c(1,2,3,4,7,8,9)]
  if (dim(go_cellular)[1] > 1) {
    write.table(go_cellular, file=paste(substring, "GO_cellular", ".tsv", sep="_"), 
                sep="\t", row.names=F, col.names=T)
  }
}

mixedToFloat <- function(x){
  is.integer  <- grepl("^\\d+$", x)
  is.fraction <- grepl("^\\d+\\/\\d+$", x)
  is.mixed    <- grepl("^\\d+ \\d+\\/\\d+$", x)
  stopifnot(all(is.integer | is.fraction | is.mixed))
  
  numbers <- strsplit(x, "[ /]")
  
  ifelse(is.integer,  as.numeric(sapply(numbers, `[`, 1)),
         ifelse(is.fraction, as.numeric(sapply(numbers, `[`, 1)) /
                  as.numeric(sapply(numbers, `[`, 2)),
                as.numeric(sapply(numbers, `[`, 1)) +
                  as.numeric(sapply(numbers, `[`, 2)) /
                  as.numeric(sapply(numbers, `[`, 3))))
}

heatmap_pathways <- function(pathways, name, max=25) {
  new_ids = unlist(lapply(pathways$Term, extract_go_id))
  new_terms = unlist(lapply(pathways$Term, extract_go_term))
  pathways$Term = new_terms
  pathways$ID = new_ids
  rownames(pathways) = pathways$Term
  
  all_pathways = pathways$Term[1:min(nrow(pathways),max)]
  all_genes = getPathwaysGenes(pathways, all_pathways)
  data = data.frame(matrix(ncol=length(all_pathways), nrow=length(all_genes)))
  rownames(data) = all_genes
  colnames(data) = all_pathways
  
  for(pathway in all_pathways) {
    genes = getPathwaysGenes(pathways, pathway)
    data[genes,pathway] = 1
  }
  
  ind = apply(data, 1, function(x) all(is.na(x)))
  data = data[!ind,] 
  data[is.na(data)] = 0
  
  pdf(file=paste0(name, ".pdf"), width=14, height=14)
  print(pheatmap(data, cluster_rows=FALSE, cellwidth=5, cellheight=5, show_rownames=TRUE, legend=FALSE, border_color=NA,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize_row=4, fontsize_col=4, angle_col=45,
                 color = colorRampPalette(brewer.pal(n=7, name ="Greys"))(100)))
  print(pheatmap(t(data), cluster_rows=FALSE, cellwidth=5, cellheight=5, show_rownames=TRUE, legend=FALSE, border_color=NA,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize_row=4, fontsize_col=4, angle_col=45,
                 color = colorRampPalette(brewer.pal(n=7, name ="Greys"))(100)))
  dev.off()
}

barplot_pathways <- function(pathways, name, max=25) {
  new_ids = unlist(lapply(pathways$Term, extract_go_id))
  new_terms = unlist(lapply(pathways$Term, extract_go_term))
  pathways$Term = new_terms
  pathways$ID = new_ids
  rownames(pathways) = pathways$Term
  
  all_genes = getPathwaysGenes(pathways, pathways$Term)
  
  pathways = pathways[1:min(nrow(pathways),max),]
  pathways$Count = str_count(pathways$Genes, ";") + 1
  pathways$Ratio = pathways$Count / length(all_genes)
  #pathways$Ratio = mixedToFloat(terms$Ratio)
  pathways$Term = factor(pathways$Term, levels=pathways$Term[order(pathways$Combined.Score)])
  
  p = ggplot(pathways, aes(x=Term, y=Combined.Score, fill=Ratio), position=position_stack(reverse=TRUE)) + 
    geom_bar(stat="identity") + 
    theme_classic() +
    scale_fill_continuous(low="red", high="blue", name="Ratio", guide=guide_colorbar(reverse=TRUE)) + 
    coord_flip() +
    ggtitle("Enriched terms") + xlab(NULL) + ylab("Score") +
    theme(axis.line=element_line(colour="black", size=1, linetype="solid"),
          axis.text.x=element_text(size=10, colour="black"),
          axis.text.y=element_text(size=10, colour="black"),
          axis.title.x=element_text(face="bold", size=12),
          axis.title.y=element_blank(),
          legend.text=element_text(size=10),
          legend.title=element_text(face="bold", size=10),
          plot.title=element_text(face="bold", size=14, colour="black"))
  pdf(file=paste0(name, ".pdf"), width=12, height=12)
  print(p)
  dev.off()
}

getPathwaysGenes <- function(enrichment_results, pathways) {
  pathway_genes = c()
  for(pathway in pathways) {
    genes = c()
    if (length(strsplit(enrichment_results[pathway,]$Genes, ";")) > 0) {
      genes = c(genes, strsplit(enrichment_results[pathway,]$Genes, ";")[[1]])
    }
    pathway_genes = c(pathway_genes, unique(genes))
  }
  pathway_genes = na.omit(pathway_genes)
  pathway_genes = unique(sapply(pathway_genes, tolower))
  firstup <- function(x) {
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    return(x)
  }
  pathway_genes = sapply(pathway_genes, firstup)
  return(pathway_genes)
}

extract_go_id <- function(go_term) {
  s = unlist(strsplit(go_term, " "))
  go = gsub("\\(", "", s[length(s)])
  go = gsub("\\)", "", go)
  return(go)
}

extract_go_term <- function(go_term) {
  s = unlist(strsplit(go_term, " "))
  term = do.call(paste, c(as.list(s[1:length(s)-1]), sep = " "))
  return(term)
}

convert_genes <- function(genes) {
  return(gsub(";", ", ", genes))
}

spot_distribution <- function(meta) {
  n_row = length(unique(meta$cluster)) * length(unique(meta$genotype))
  df = data.frame(matrix(vector(), n_row, 3,
                         dimnames=list(paste(rep(unique(meta$cluster), 
                                                 each=length(unique(meta$genotype))), 
                                             unique(meta$genotype), sep="_"), 
                                       c("Cluster", "Count", "Genotype"))),
                  stringsAsFactors=F)
  for (c in unique(meta$cluster)) {
    for (g in unique(meta$genotype)) {
      index = paste(c, g, sep="_")
      df[index, "Cluster"] = c
      df[index, "Genotype"] = g
      df[index, "Count"] = nrow(meta[meta$cluster == c & meta$genotype == g,])
    }
  }
  p = ggplot(df, aes(x=Cluster, y=Count, fill=Genotype), position=position_stack(reverse=TRUE)) + 
    geom_bar(stat="identity") + 
    theme_classic() +
    coord_flip() +
    ggtitle("Spots distribution") + xlab("Cluster") + ylab("#Spots") +
    scale_fill_manual("Genotype", values = c("WT" = "#ACDFE7", "PB" = "#F5BED7", "3xAD" = "#2484BC", "3xPB" = "#D53C44"))
  theme(axis.line=element_line(colour="black", size=1, linetype="solid"),
        axis.text.x=element_text(size=10, colour="black"),
        axis.text.y=element_text(size=10, colour="black"),
        axis.title.x=element_text(face="bold", size=12),
        axis.title.y=element_text(face="bold", size=12),
        legend.text=element_text(size=10),
        legend.title=element_text(face="bold", size=10),
        plot.title=element_text(face="bold", size=14, colour="black"))
  pdf(file="distribution_spots.pdf", width=6, height=6)
  print(p)
  dev.off()
}

do_dea_deseq <- function(counts, meta, test="LRT", formulaN="~condition", formulaR="~1", suffix="") {
  # Normalize and run DESeq2 (optimized for single cell)
  dds = SummarizedExperiment(assays=list(counts=counts))
  colData(dds)$animal = as.factor(meta$animal)
  colData(dds)$genotype = as.factor(meta$genotype)
  colData(dds)$cluster = as.factor(meta$cluster)
  colData(dds)$condition = factor(paste(meta$cluster, meta$genotype, sep="_"))
  if ("chip" %in% colnames(meta)) {
    colData(dds)$chip = as.factor(meta$chip)
  }
  dds = DESeqDataSet(dds, design=formula(formulaN))
  if (test == "LRT") {
    dds = DESeq(dds, parallel=TRUE, test=test, reduced=formula(formulaR), 
                useT=TRUE, minmu=1e-6, sfType="poscounts", minReplicatesForReplace=Inf)
  } else {
    dds = DESeq(dds, parallel=TRUE, test=test, 
                useT=TRUE, minmu=1e-6, sfType="poscounts", minReplicatesForReplace=Inf)   
  }
  pdf(paste0(suffix, "disp_est_deseq.pdf"))
  DESeq2::plotDispEsts(dds)
  DESeq2::plotMA(dds)
  plot(colSums(counts), sizeFactors(dds), xlab="sums of counts", ylab="size factors")
  plot(colMeans(counts), sizeFactors(dds), xlab="means of counts", ylab="size factors")
  dev.off()
  return(dds)
}

write_deseq_results <- function(dds, meta, prefix, FDR_T=0.1, do_shrink=TRUE, do_genotype=TRUE, do_spatial=TRUE) {
  i = 1
  clusters = unique(meta$cluster)
  for (c in clusters) {
    for (t in unique(meta$genotype)) {
      # AD vs WT in cluster c
      if (t != "WT" && do_genotype) {
        res = results(dds, contrast=c("condition", 
                                      paste(c, t, sep="_"), 
                                      paste(c, "WT", sep="_")), 
                      parallel=TRUE, alpha=FDR_T)
        if (do_shrink) {
          res = lfcShrink(dds, contrast=c("condition", 
                                          paste(c, t, sep="_"), 
                                          paste(c, "WT", sep="_")),
                          res=res, type="normal", parallel=TRUE)
        }
        res = as.data.frame(res)
        res = na.omit(res[order(res$padj),])
        write.table(res, file=paste(prefix, c, t, "vs_WT.tsv", sep="_"), sep="\t")
      }
      # Cluster c vs other clusters 
      for (c2 in clusters[i:length(clusters)]) {
        if (c != c2 && do_spatial) {
          res = results(dds, contrast=c("condition",
                                        paste(c, t, sep="_"),
                                        paste(c2, t, sep="_")), 
                        parallel=TRUE, alpha=FDR_T)
          if (do_shrink) {
            res = lfcShrink(dds, contrast=c("condition", 
                                            paste(c, t, sep="_"), 
                                            paste(c2, t, sep="_")),
                            res=res, type="normal", parallel=TRUE)
          }
          res = as.data.frame(res)
          res = na.omit(res[order(res$padj),])
          write.table(res, file=paste(prefix, t, c, "vs", c2, ".tsv", sep="_"), sep="\t")
        }
      }
    }
    i = i + 1
  }
}

enrichment_scores_gsva <- function(counts, geneSets, n_cores=1) {
  gbm = gsva(counts, geneSets, mx.diff=FALSE, verbose=FALSE, parallel.sz=n_cores, min.sz=2)
  return(gbm)
}

enrichment_scores_auc <- function(counts, geneSets, n_cores=1, ngenes=200) {
  cells_rankings = AUCell_buildRankings(counts, nCores=n_cores, plotStats=FALSE)
  cells_AUC = AUCell_calcAUC(geneSets, cells_rankings, nCores=n_cores, aucMaxRank=ngenes)
  cells_assignment = AUCell_exploreThresholds(cells_AUC, plotHist=FALSE, assign=TRUE)
  auc_table = getAUC(cells_AUC)
  return(auc_table)
}

enrichment_scores_naive <- function(counts, geneSets, USE_MEAN=FALSE, MIN_LENGTH=5) {
  scores = data.frame(matrix(vector(), 
                             length(colnames(counts)), length(names(geneSets)),
                             dimnames=list(colnames(counts), names(geneSets))),
                      stringsAsFactors=F)
  for (gene_set in names(geneSets)) {
    genes_in_set = unlist(geneIds(geneSets[gene_set]), use.names=FALSE)
    genes_in_set = intersect(genes_in_set, rownames(counts))
    if (length(genes_in_set) >= MIN_LENGTH) {
      slice = counts[genes_in_set,]
      num_genes = length(rownames(slice))
      for (spot in colnames(counts)) {
        if (USE_MEAN) {
          scores[spot,gene_set] = mean(slice[,spot])
        } else {
          scores[spot,gene_set] = sum(slice[,spot]) / sum(counts[,spot])
        }
      }
    }
  }
  scores[is.na(scores)] = 0
  return(t(scores))
}

enrichment_scores_fgsa <- function(counts, geneSets, n_cores=1, MAX_SIZE=500, MIN_SIZE=3) {
  scores = data.frame(matrix(vector(), 
                             length(colnames(counts)), length(names(geneSets)),
                             dimnames=list(colnames(counts), names(geneSets))),
                      stringsAsFactors=F)
  pathways = list()
  for (gene_set in names(geneSets)) {
    pathways[[gene_set]] = unlist(geneIds(geneSets[gene_set]), use.names=FALSE)
  }
  colnames(scores) = names(pathways)
  
  for (spot in colnames(counts)) {
    ranks = as.vector(counts[,spot])
    names(ranks) = rownames(counts)
    ranks = sort(ranks, decreasing=FALSE)
    fgseaRes = fgseaMultilevel(pathways, ranks, nproc=n_cores, scoreType="pos", minSize=MIN_SIZE, maxSize=MAX_SIZE)
    for (gene_set in fgseaRes$pathway) {
      scores[spot,gene_set] = fgseaRes[fgseaRes$pathway == gene_set]$ES
    }
  }
  scores[is.na(scores)] = 0
  return(t(scores))
}

enrichment_scores_singleR <- function(counts, ann) {
  imm = ImmGenData()
  mouse = MouseRNAseqData()
  
  pred = SingleR(test=counts, ref=imm, labels=imm$label.main)
  scores = as.data.frame(pred$scores)
  rownames(scores) = rownames(pred)
  scores = range01(scores)
  heatmap_counts(t(scores), ann, "singleCell_imm_scores.pdf", COL, ann_colors, 5)
  
  pred = SingleR(test=counts, ref=mouse, labels=mouse$label.main)
  scores = as.data.frame(pred$scores)
  rownames(scores) = rownames(pred)
  scores = range01(scores)
  heatmap_counts(t(scores), ann, "singleCell_mouse_scores.pdf", COL, ann_colors, 5)
}

volcano_plot <- function(markers_table, labels, cols, name, p_val_threshold=0.05, fc_threshold=0.5) {
  # Markers table must has the following columns:
  # p_val avg_logFC, p_val_adj, cluster, gene
  n_clusters = length(unique(markers_table$cluster))
  clusters = unique(markers_table$cluster)
  # Check if any cluster is missing and if so re-arrange
  if (n_clusters != length(labels)) {
    i = 0
    for (c in clusters) {
      markers_table[markers_table$cluster == c,]$cluster = i
      i = i + 1
    }
    labels = labels[clusters + 1]
  }
  # Create summary table
  markers_table_summary = markers_table %>% 
    group_by(cluster) %>% 
    arrange(-avg_logFC) %>%
    summarize(threshold=fc_threshold, 
              max_val=max(fc_threshold,avg_logFC[1]), 
              low.threshold=-fc_threshold, 
              min_val=min(-fc_threshold,avg_logFC[n()]))
  markers_table = markers_table %>% 
    group_by(cluster) %>%
    arrange(-avg_logFC) %>%
    mutate(DE=ifelse(p_val_adj < p_val_threshold, "yes", "no"),
           sign.flat=ifelse(avg_logFC >= fc_threshold | avg_logFC <= -fc_threshold, gene, "")) 
  
  # Add labels to clusters
  markers_table_summary$labels = labels
  
  # Flat DE plot
  p = ggplot() +
    theme_classic() +
    theme(axis.text.x=element_blank(), axis.line.x=element_blank(), 
          axis.ticks.x=element_blank())
  
  pos = position_jitter(width = 0.4, seed = 1)
  for (i in 1:n_clusters) {
    c = as.character(unique(markers_table$cluster))[i]
    markers_table_subset = subset(markers_table, cluster == c)
    
    p = p + 
      geom_point(data=markers_table_subset, 
                 aes(x=cluster, y=avg_logFC, color=DE), alpha=0.7, size=0.4, position=pos) +
      annotate("rect", xmin=i-1.4, xmax=i-.4, ymin=-0.1, ymax=0.1, fill=cols[i], color="black") +
      annotate("rect", xmin=i-1.4, xmax=i-.4, ymin=as.numeric(markers_table_summary[i, 2]), 
               ymax=as.numeric(markers_table_summary[i, 3]), alpha=.05, fill="black") +
      annotate("rect", xmin=i-1.4, xmax=i-.4, ymin=as.numeric(markers_table_summary[i, 5]), 
               ymax=as.numeric(markers_table_summary[i, 4]), alpha=.05, fill="black") +
      geom_text_repel(data=markers_table_subset, aes(x=cluster, y=avg_logFC, label=sign.flat), 
                      position=pos, size=3, segment.size=0.1)
  }
  p = p + scale_color_manual(values = c("no" = "black", "yes" = "red")) +
    geom_text(data=markers_table_summary, aes(x=0:(n_clusters-1), y=rep(0, n_clusters), label=labels), color="black") +
    ggtitle("DE genes per cluster") +
    labs(color=paste0("adj-pval < ", p_val_threshold))
  p
  ggsave(name, limitsize=FALSE, width=14, height=14)
}

heatmap_genes <- function(candidate_genes, 
                          dg_ad_vs_wt, dg_adpolb_vs_wt, dg_polb_vs_wt,
                          ca1_ad_vs_wt, ca1_adpolb_vs_wt, ca1_polb_vs_wt,
                          ca23_ad_vs_wt, ca23_adpolb_vs_wt, ca23_polb_vs_wt,
                          name, fontsize, width, height) {
  
  data = data.frame(ADP_CA1=numeric(), AD_CA1=numeric(), POLB_CA1=numeric(),
                    ADP_CA2_3=numeric(), AD_CA2_3=numeric(), POLB_CA2_3=numeric(),
                    ADP_DG=numeric(), AD_DG=numeric(), POLB_DG=numeric())
  
  for(gene in candidate_genes) {
    data[gene,] = c(ca1_adpolb_vs_wt[gene,]$log2FoldChange, ca1_ad_vs_wt[gene,]$log2FoldChange, ca1_polb_vs_wt[gene,]$log2FoldChange,
                    ca23_adpolb_vs_wt[gene,]$log2FoldChange, ca23_ad_vs_wt[gene,]$log2FoldChange, ca23_polb_vs_wt[gene,]$log2FoldChange,
                    dg_adpolb_vs_wt[gene,]$log2FoldChange, dg_ad_vs_wt[gene,]$log2FoldChange, dg_polb_vs_wt[gene,]$log2FoldChange)
    
  }
  
  ind = apply(data, 1, function(x) all(is.na(x)))
  data = data[!ind,]
  data[is.na(data)] = 0.0
  data$AD_DG = round(as.numeric(data$AD_DG), digits=2)
  data$ADP_DG  = round(as.numeric(data$ADP_DG), digits=2)
  data$POLB_DG  = round(as.numeric(data$POLB_DG), digits=2)
  data$AD_CA1 = round(as.numeric(data$AD_CA1), digits=2)
  data$ADP_CA1 = round(as.numeric(data$ADP_CA1), digits=2)
  data$POLB_CA1 = round(as.numeric(data$POLB_CA1), digits=2)
  data$AD_CA2_3 = round(as.numeric(data$AD_CA2_3), digits=2)
  data$ADP_CA2_3 = round(as.numeric(data$ADP_CA2_3), digits=2)
  data$POLB_CA2_3 = round(as.numeric(data$POLB_CA2_3), digits=2)
  
  data_nopb = data[,c("ADP_CA1", "AD_CA1", "ADP_CA2_3", "AD_CA2_3", "ADP_DG", "AD_DG")]
  
  COL = colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)

  pdf(paste0(name, ".pdf"), width=width, height=height)
  print(pheatmap(data, cluster_rows=FALSE, show_rownames=TRUE, cellwidth=4, cellheight=4,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize=fontsize,
                 angle_col=315, color=COL, border_color=NA, na_col="white"))
  print(pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cellwidth=4, cellheight=4,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize=fontsize,
                 angle_col=315, color=COL, border_color=NA, na_col="white"))
  print(pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cellwidth=4, cellheight=4,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize=fontsize,
                 angle_col=315, color=COL, border_color=NA, na_col="white"))
  print(pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cellwidth=4, cellheight=4,
                 cluster_cols=TRUE, show_colnames=TRUE, fontsize=fontsize,
                 angle_col=315, color=COL, border_color=NA, na_col="white"))
  dev.off()
}

heatmap_genes_MOB <- function(candidate_genes, 
                              epl_ad_vs_wt, epl_adpolb_vs_wt, epl_polb_vs_wt,
                              gcl_ad_vs_wt, gcl_adpolb_vs_wt, gcl_polb_vs_wt,
                              mcl_ad_vs_wt, mcl_adpolb_vs_wt, mcl_polb_vs_wt,
                              gl_ad_vs_wt, gl_adpolb_vs_wt, gl_polb_vs_wt,
                              name, fontsize, width, height) {
  data = data.frame(AD_EPL=numeric(), ADP_EPL=numeric(), POLB_EPL=numeric(),
                    AD_GCL=numeric(), ADP_GCL=numeric(), POLB_GCL=numeric(),
                    AD_MCL=numeric(), ADP_MCL=numeric(), POLB_MCL=numeric(),
                    AD_GL=numeric(), ADP_GL=numeric(), POLB_GL=numeric())
  
  for(gene in candidate_genes) {
    data[gene,] = c(epl_ad_vs_wt[gene,]$log2FoldChange, epl_adpolb_vs_wt[gene,]$log2FoldChange, epl_polb_vs_wt[gene,]$log2FoldChange,
                    gcl_ad_vs_wt[gene,]$log2FoldChange, gcl_adpolb_vs_wt[gene,]$log2FoldChange, gcl_polb_vs_wt[gene,]$log2FoldChange,
                    mcl_ad_vs_wt[gene,]$log2FoldChange, mcl_adpolb_vs_wt[gene,]$log2FoldChange, mcl_polb_vs_wt[gene,]$log2FoldChange,
                    gl_ad_vs_wt[gene,]$log2FoldChange, gl_adpolb_vs_wt[gene,]$log2FoldChange, gl_polb_vs_wt[gene,]$log2FoldChange)
    
  }
  
  ind = apply(data, 1, function(x) all(is.na(x)))
  data = data[!ind,]
  data[is.na(data)] = 0.0
  data$AD_EPL = round(as.numeric(data$AD_EPL), digits=2)
  data$ADP_EPL = round(as.numeric(data$ADP_EPL), digits=2)
  data$POLB_EPL = round(as.numeric(data$POLB_EPL), digits=2)
  data$AD_GCL = round(as.numeric(data$AD_GCL), digits=2)
  data$ADP_GCL = round(as.numeric(data$ADP_GCL), digits=2)
  data$POLB_GCL = round(as.numeric(data$POLB_GCL), digits=2)
  data$AD_MCL = round(as.numeric(data$AD_MCL), digits=2)
  data$ADP_MCL = round(as.numeric(data$ADP_MCL), digits=2)
  data$POLB_MCL = round(as.numeric(data$POLB_MCL), digits=2)
  data$AD_GL = round(as.numeric(data$AD_GL), digits=2)
  data$ADP_GL = round(as.numeric(data$ADP_GL), digits=2)
  data$POLB_GL = round(as.numeric(data$POLB_GL), digits=2)
  
  data_nopb = data[,c("AD_EPL", "ADP_EPL", "AD_GCL", "ADP_GCL", "AD_MCL", "ADP_MCL", "AD_GL", "ADP_GL")]
  
  COL = colorRampPalette(rev(brewer.pal(n=7, name="RdBu")))(100)
  
  pdf(paste0(name, ".pdf"), width=width, height=height)
  print(pheatmap(data, cluster_rows=FALSE, show_rownames=TRUE, cellwidth=4, cellheight=4,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize=fontsize,
                 angle_col=315, color=COL, border_color=NA, na_col="white"))
  print(pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cellwidth=4, cellheight=4,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize=fontsize,
                 angle_col=315, color=COL, border_color=NA, na_col="white"))
  print(pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cellwidth=4, cellheight=4,
                 cluster_cols=FALSE, show_colnames=TRUE, fontsize=fontsize,
                 angle_col=315, color=COL, border_color=NA, na_col="white"))
  print(pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, cellwidth=4, cellheight=4,
                 cluster_cols=TRUE, show_colnames=TRUE, fontsize=fontsize,
                 angle_col=315, color=COL, border_color=NA, na_col="white"))
  dev.off()
}

do_coepr <- function(data, name, h_value, linkage="average", fontsize=10) {
  corr = cor(t(data))
  corr[is.na(corr)] = 0
  dist = as.dist(1 - corr)
  clust = fastcluster::hclust(dist, method=linkage)
  cut = cutree(clust, h=h_value)
  names(cut) = rownames(data)
  cut_ann = as.data.frame(cut)
  colnames(cut_ann) = c("Module")
  rownames(cut_ann) = rownames(data)
  cut_ann$Module = as.factor(cut_ann$Module)
  
  pdf(paste0("dendogram_coepr_", name, ".pdf"))
  plot(clust, labels=FALSE, cex=0.1, hang=-1)
  rect.hclust(clust, h=h_value)
  dev.off()
  
  pdf(paste0("heatmap_coepr_", name, ".pdf"))
  print(pheatmap(corr, clustering_distance_cols=dist, clustering_distance_rows=dist,
                 cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=TRUE,
                 annotation_row=cut_ann, annotation_col=cut_ann,
                 show_colnames=FALSE, fontsize=fontsize, clustering_method=linkage))
  dev.off()
  
  data = data.frame(module=character(), genes=character(), stringsAsFactors=FALSE)
  entrez_genes = c()
  modules = c()
  for (m in unique(cut)) {
    data[m,"module"] = m
    genes = paste(names(cut[cut == m]), collapse=";")
    data[m,"genes"] = genes
    split_genes = unlist(strsplit(genes, ";"))
    modules = c(modules, rep(m, length(split_genes)))
    entrez_genes = c(entrez_genes, mapIds(org.Mm.eg.db, keys=split_genes, column="ENTREZID", keytype="SYMBOL", multiVals="first"))
    doEnrichR(split_genes, paste(paste0(name, "_coepr_module_"), m, "_pathways", sep=""), min_genes=3)
  }         
  write.table(data, paste0("coexpression_modules_", name, ".tsv"), sep="\t", row.names=FALSE)
  
  data_cluster = data.frame(Entrez=entrez_genes,  group=modules)
  res1 = compareCluster(Entrez~group, data=data_cluster, fun='enrichGO', organism="mmu", pvalueCutoff=0.05)
  pdf(paste0(name, "_coepr_enrichGO.pdf"))
  print(dotplot(res1, title="MODULES PATHWAYS GO", font.size=10))
  dev.off()
}

heatmap_counts <- function(data, ann, filename, col_pal, ann_colors, font=8) {
  pdf(filename)
  print(pheatmap(data, cluster_rows=FALSE, show_rownames=TRUE, color=col_pal,
                 cluster_cols=FALSE, annotation_col=ann, annotation_colors=ann_colors, show_colnames=FALSE, fontsize=font))
  print(pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, color=col_pal,
                 cluster_cols=FALSE, annotation_col=ann, annotation_colors=ann_colors, show_colnames=FALSE, fontsize=font))
  print(pheatmap(data, cluster_rows=FALSE, show_rownames=TRUE, color=col_pal,
                 cluster_cols=TRUE, annotation_col=ann, annotation_colors=ann_colors, show_colnames=FALSE, fontsize=font))
  print(pheatmap(data, cluster_rows=TRUE, show_rownames=TRUE, color=col_pal,
                 cluster_cols=TRUE, annotation_col=ann, annotation_colors=ann_colors, show_colnames=FALSE, fontsize=font))
  dev.off()
}

firstup <- function(x) {
  x = tolower(x)
  substr(x, 1, 1) = toupper(substr(x, 1, 1))
  return(x)
}

range01 <- function(x) {
  return((x-min(x))/(max(x)-min(x)))
}

parse_dea_results <- function(filename, PVALUE=0.1, FDR=0.1, FC=0.5) {
  dea = read.delim(filename, sep="\t", header=T, row.names=1)
  dea$log2FoldChange = round(dea$log2FoldChange, digits=3)
  dea$padj = round(dea$padj, digits=3)
  dea$pvalue = round(dea$pvalue, digits=3)
  dea = dea[abs(dea$log2FoldChange) >= FC & dea$padj <= FDR & dea$pvalue <= PVALUE,]
  dea = dea[rev(order(abs(dea$log2FoldChange), dea$padj, dea$pvalue)),]
  rownames(dea) = gsub("[.]", "-", rownames(dea))
  # Remove mt (mito) genes
  if (length(grep("mt-", rownames(dea))) > 0) {
    dea = dea[-grep("^mt-", rownames(dea)),]
  }
  return(dea)
}

fc_corr <- function(dea1, dea2, xlabel, ylabel, title) {
  shared = intersect(rownames(dea1), rownames(dea2))
  df = data.frame(X=dea1[shared, "log2FoldChange"], 
                  Y=dea2[shared, "log2FoldChange"])
  df$label = shared
  df[sign(df$X) == sign(df$Y),]$label = ""

  pdf(paste0(title, "_correlation_shared.pdf"))
  options(repr.plot.width=14, repr.plot.height=14)
  p = ggplot(df, aes(X, Y, label=label)) +
    geom_point(size=1, col=ifelse(sign(df$X) != sign(df$Y), "red", "black")) +
    geom_smooth(method=lm, linetype="dashed", color="darkred", fill="blue") +
    guides(fill=FALSE) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    geom_text_repel(size=2, segment.size=0.25, data=subset(df, label != "")) +
    labs(title=title, x=xlabel, y=ylabel) +
    stat_cor(method="pearson", label.x=-1, label.y=1)
  print(p)
  dev.off()
  
  both = unique(rownames(dea1), rownames(dea2))
  df = data.frame(X=dea1[both, "log2FoldChange"], 
                  Y=dea2[both, "log2FoldChange"])
  df$label = both
  df[is.na(df)] = 0
  df[sign(df$X) == sign(df$Y) & df$X != 0 & df$Y != 0,]$label = ""
  
  pdf(paste0(title, "_correlation_both.pdf"))
  options(repr.plot.width=14, repr.plot.height=14)
  p = ggplot(df, aes(X, Y, label=label)) +
    geom_point(size=1, col=ifelse(df$X == 0 | df$Y == 0 | sign(df$X) != sign(df$Y), "red", "black")) +
    geom_smooth(method=lm, linetype="dashed", color="darkred", fill="blue") +
    guides(fill=FALSE) +
    scale_x_continuous(expand=c(0, 0)) +
    scale_y_continuous(expand=c(0, 0)) +
    geom_text_repel(size=2, segment.size=0.25, data=subset(df, label != "")) +
    labs(title=title, x=xlabel, y=ylabel) +
    stat_cor(method="pearson", label.x=-1, label.y=1)
  print(p)
  dev.off()
}
