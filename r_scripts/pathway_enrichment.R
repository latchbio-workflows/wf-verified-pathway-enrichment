suppressMessages(suppressWarnings({
  library(clusterProfiler)
  library(KEGG.db)
  library(DOSE)
  library(ReactomePA)
  library(msigdbr)
  library(pathview)
  library(ggridges)
  library(ggplot2)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(org.Sc.sgd.db)
  library(org.Mmu.eg.db)
}))

args <- commandArgs(trailingOnly = TRUE)

dge_file <- args[1]
root_path <- args[2]
genome <- args[3]

# dge_file <- "/root/pe/deseq2_results.csv"
# root_path <- "/root/pe/output9"
# genome <- "homo_sapiens"

orgDb <- switch(genome,
                "homo_sapiens" = org.Hs.eg.db,
                "mus_musculus" = org.Mm.eg.db,
                "rattus_norvegicus" = org.Rn.eg.db,
                "rhesus_macaque" = org.Mmu.eg.db,
                "saccharomyces_cerevisiae" = org.Sc.sgd.db,
                stop("Invalid genome selected"))

KEGG_organism <- switch(genome,
               "homo_sapiens" = 'hsa',
               "mus_musculus" = 'mmu',
               "rattus_norvegicus" = 'rnu',
               "rhesus_macaque" = 'mcc',
               "saccharomyces_cerevisiae" = 'sce',
               stop("Invalid genome selected"))

WP_organism <- switch(genome,
                        "homo_sapiens" = 'Homo sapiens',
                        "mus_musculus" = 'Mus musculus',
                        "rattus_norvegicus" = 'Rattus norvegicus',
                        "rhesus_macaque" = 'None',
                        "saccharomyces_cerevisiae" = 'Saccharomyces cerevisiae',
                        stop("Invalid genome selected"))

msig_species <- switch(genome,
                       "homo_sapiens" = "Homo sapiens",
                       "mus_musculus" = "Mus musculus",
                       "rattus_norvegicus" = "Rattus norvegicus",
                       "rhesus_macaque" = "Macaca mulatta",
                       "saccharomyces_cerevisiae" = "Saccharomyces cerevisiae",
                       stop("Invalid genome selected for MSigDB"))

# Load Data
dge_data <- read.csv(dge_file, row.names = 1)

#------------------------- Data Preprocessing ------------------------#

# Create a ranked list of DGEs
ranked_dges <- dge_data$log2FoldChange
names(ranked_dges) <- rownames(dge_data)
ranked_dges <- sort(ranked_dges, decreasing = TRUE)  # Sorting for GSEA

# Identify up- and down-regulated DGEs
upregulated_dges <- rownames(subset(dge_data, log2FoldChange > 1 & padj < 0.05))
downregulated_dges <- rownames(subset(dge_data, log2FoldChange < -1 & padj < 0.05))

# Map gene symbols to Entrez IDs
upregulated_dges_entrez <- mapIds(orgDb, keys = upregulated_dges, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first") %>% na.omit()
downregulated_dges_entrez <- mapIds(orgDb, keys = downregulated_dges, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first") %>% na.omit()
ranked_dges_entrez <- mapIds(orgDb, keys = names(ranked_dges), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first") %>% na.omit()

# Create a new data frame combining ranked_list and entrez_ranked_list
ranked_dges_entrez_dataframe <- data.frame(SYMBOL = names(ranked_dges),
                                           log2FoldChange = ranked_dges,
                                           ENTREZID = as.numeric(ranked_dges_entrez[names(ranked_dges)]))
ranked_dges_entrez_dataframe <- na.omit(ranked_dges_entrez_dataframe)
ranked_dges_entrez_dataframe <- ranked_dges_entrez_dataframe[order(-ranked_dges_entrez_dataframe$log2FoldChange),]

create_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
}

create_dir(root_path)

#------------------------- Gene Ontology (GO) ------------------------#

go_path <- file.path(root_path, "GO")
create_dir(go_path)

# Run ORA for GO
upregulated_dges_ora_go <- enrichGO(gene = upregulated_dges,
                                    OrgDb = orgDb,
                                    keyType = 'SYMBOL', ont = "ALL", pool = FALSE,
                                    pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
upregulated_dges_ora_go@result <- upregulated_dges_ora_go@result[order(upregulated_dges_ora_go@result$p.adjust), ]
write.csv(upregulated_dges_ora_go@result, file.path(go_path, "upregulated_dges_go_ora.csv"), row.names = FALSE)

downregulated_dges_ora_go <- enrichGO(gene = downregulated_dges,
                                      OrgDb = orgDb,
                                      keyType = 'SYMBOL', ont = "ALL", pool = FALSE,
                                      pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
downregulated_dges_ora_go@result <- downregulated_dges_ora_go@result[order(downregulated_dges_ora_go@result$p.adjust), ]
write.csv(downregulated_dges_ora_go@result, file.path(go_path, "downregulated_dges_go_ora.csv"), row.names = FALSE)

# Run GSEA for GO
dges_gsea_go <- gseGO(geneList = ranked_dges,
                      OrgDb = orgDb,
                      keyType = 'SYMBOL', ont = "ALL",
                      pAdjustMethod = "BH", pvalueCutoff = 0.05,
                      by = "fgsea", verbose = TRUE, eps = 0)
dges_gsea_go@result <- dges_gsea_go@result[order(dges_gsea_go@result$p.adjust), ]
write.csv(dges_gsea_go@result, file.path(go_path, "dges_go_gsea.csv"), row.names = FALSE)


#------------------------- KEGG Pathways -------------------------#

kegg_path <- file.path(root_path, "KEGG")
create_dir(kegg_path)

# Run ORA for KEGG (tries to connect online first)
tryCatch({
  upregulated_dges_ora_kegg <- enrichKEGG(gene = upregulated_dges_entrez, organism = KEGG_organism,
                                          pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                          use_internal_data = FALSE)
}, error = function(e) {
  upregulated_dges_ora_kegg <- enrichKEGG(gene = upregulated_dges_entrez, organism = KEGG_organism,
                                          pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                          use_internal_data = TRUE)
})
upregulated_dges_ora_kegg@result <- upregulated_dges_ora_kegg@result[order(upregulated_dges_ora_kegg@result$p.adjust), ]
write.csv(upregulated_dges_ora_kegg@result, file.path(kegg_path, "upregulated_dges_kegg_ora.csv"), row.names = FALSE)

tryCatch({
  downregulated_dges_ora_kegg <- enrichKEGG(gene = downregulated_dges_entrez, organism = KEGG_organism,
                                            pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                            use_internal_data = FALSE)
}, error = function(e) {
  downregulated_dges_ora_kegg <- enrichKEGG(gene = downregulated_dges_entrez, organism = KEGG_organism,
                                            pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                            use_internal_data = TRUE)
})
downregulated_dges_ora_kegg@result <- downregulated_dges_ora_kegg@result[order(downregulated_dges_ora_kegg@result$p.adjust), ]
write.csv(downregulated_dges_ora_kegg@result, file.path(kegg_path, "downregulated_dges_kegg_ora.csv"), row.names = FALSE)

# Run GSEA for KEGG
tryCatch({
  dges_gsea_kegg <- gseKEGG(geneList = setNames(ranked_dges_entrez_dataframe$log2FoldChange, ranked_dges_entrez_dataframe$ENTREZID),
                            organism = KEGG_organism,
                            pAdjustMethod = "BH", pvalueCutoff = 0.05,
                            by = "fgsea", verbose = TRUE, eps = 0,
                            use_internal_data = FALSE)
}, error = function(e) {
  dges_gsea_kegg <- gseKEGG(geneList = setNames(ranked_dges_entrez_dataframe$log2FoldChange, ranked_dges_entrez_dataframe$ENTREZID),
                            organism = KEGG_organism,
                            pAdjustMethod = "BH", pvalueCutoff = 0.05,
                            by = "fgsea", verbose = TRUE, eps = 0,
                            use_internal_data = TRUE)
})

entrez_to_symbol <- setNames(ranked_dges_entrez_dataframe$SYMBOL, ranked_dges_entrez_dataframe$ENTREZID)

# Add gene names to the results
if (!is.null(dges_gsea_kegg) && nrow(dges_gsea_kegg@result) > 0) {
  dges_gsea_kegg@result <- dges_gsea_kegg@result %>%
    mutate(
      core_enrichment_symbols = sapply(strsplit(core_enrichment, "/"), function(ids) {
        symbols <- entrez_to_symbol[ids]
        paste(symbols[!is.na(symbols)], collapse = "/")
      })
    )
}

dges_gsea_kegg@result <- dges_gsea_kegg@result[order(dges_gsea_kegg@result$p.adjust), ]
write.csv(dges_gsea_kegg@result, file.path(kegg_path, "dges_kegg_gsea.csv"), row.names = FALSE)

#  KEGG Pathview results
pathview_dir <- file.path(kegg_path, "Pathview")
create_dir(pathview_dir)

original_dir <- getwd()
setwd(pathview_dir)

if (!is.null(dges_gsea_kegg) && nrow(dges_gsea_kegg@result) > 0) {

  # Select top pathways (e.g., top 5)
  top_pathways <- head(dges_gsea_kegg@result$ID, 20)

  for (pathway in top_pathways) {
    tryCatch({
      cat(paste("  Running pathview on", pathway, "\n"))

      gene_data <- setNames(ranked_dges_entrez_dataframe$log2FoldChange,
                            ranked_dges_entrez_dataframe$ENTREZID)

      pathview(gene.data = gene_data,
               pathway.id = pathway,
               species = KEGG_organism,
               out.suffix = "pathview",
               kegg.native = TRUE,
               same.layer = FALSE,
               gene.idtype = "ENTREZ",
               out.dir = pathview_dir)

      cat(paste("  Pathview files for", pathway, "saved in", pathview_dir, "\n"))

    }, error = function(e) {
      cat(paste("    Failed to generate KEGG pathview for", pathway, ":", e$message, "\n"))
    })
  }
} else {
  cat("No significant KEGG pathways found for pathview visualization.\n")
}
setwd(original_dir)

#------------------------- Disease Ontology (DO) -------------------------#

do_path <- file.path(root_path, "DO")
create_dir(do_path)

# Run ORA for Disease Ontology
upregulated_dges_ora_do <- enrichDO(gene = upregulated_dges_entrez, ont = "DO",
                                    pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
upregulated_dges_ora_do@result <- upregulated_dges_ora_do@result[order(upregulated_dges_ora_do@result$p.adjust), ]
write.csv(upregulated_dges_ora_do@result, file.path(do_path, "upregulated_dges_do_ora.csv"), row.names = FALSE)

downregulated_dges_ora_do <- enrichDO(gene = downregulated_dges_entrez, ont = "DO",
                                      pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
downregulated_dges_ora_do@result <- downregulated_dges_ora_do@result[order(downregulated_dges_ora_do@result$p.adjust), ]
write.csv(downregulated_dges_ora_do@result, file.path(do_path, "downregulated_dges_do_ora.csv"), row.names = FALSE)

# Run GSEA for Disease Ontology
dges_gsea_do <- gseDO(geneList = setNames(ranked_dges_entrez_dataframe$log2FoldChange, ranked_dges_entrez_dataframe$ENTREZID),
                      pAdjustMethod = "BH", pvalueCutoff = 0.05, eps = 0,
                      by = "fgsea", verbose = TRUE)
dges_gsea_do@result <- dges_gsea_do@result[order(dges_gsea_do@result$p.adjust), ]
write.csv(dges_gsea_do@result, file.path(do_path, "dges_do_gsea.csv"), row.names = FALSE)

#------------------------- Wikipathways -------------------------#

wp_path <- file.path(root_path, "WP")
create_dir(wp_path)

if (WP_organism != 'None') {
  # Run ORA for WP
  upregulated_dges_ora_wp <- enrichWP(gene = upregulated_dges_entrez, organism = WP_organism,
                                      pAdjustMethod = "BH", pvalueCutoff = 0.05)

  if (!is.null(upregulated_dges_ora_wp) && !is.null(upregulated_dges_ora_wp@result)) {
    upregulated_dges_ora_wp@result <- upregulated_dges_ora_wp@result[order(upregulated_dges_ora_wp@result$p.adjust), ]
    write.csv(upregulated_dges_ora_wp@result, file.path(wp_path, "upregulated_dges_wp_ora.csv"), row.names = FALSE)
  }

  downregulated_dges_ora_wp <- enrichWP(gene = downregulated_dges_entrez, organism = WP_organism,
                                        pAdjustMethod = "BH", pvalueCutoff = 0.05)

  if (!is.null(downregulated_dges_ora_wp) && !is.null(downregulated_dges_ora_wp@result)) {
    downregulated_dges_ora_wp@result <- downregulated_dges_ora_wp@result[order(downregulated_dges_ora_wp@result$p.adjust), ]
    write.csv(downregulated_dges_ora_wp@result, file.path(wp_path, "downregulated_dges_wp_ora.csv"), row.names = FALSE)
  }

  # Run GSEA for WP
  dges_gsea_wp <- gseWP(geneList = setNames(ranked_dges_entrez_dataframe$log2FoldChange, ranked_dges_entrez_dataframe$ENTREZID),
                        organism = WP_organism,
                        pAdjustMethod = "BH", pvalueCutoff = 0.05, eps = 0,
                        by = "fgsea", verbose = TRUE)

  if (!is.null(dges_gsea_wp) && !is.null(dges_gsea_wp@result)) {
    dges_gsea_wp@result <- dges_gsea_wp@result[order(dges_gsea_wp@result$p.adjust), ]
    write.csv(dges_gsea_wp@result, file.path(wp_path, "dges_wp_gsea.csv"), row.names = FALSE)
  }
}

#------------------------- MSigDB Analysis -------------------------#

msig_path <- file.path(root_path, "MSIG")
create_dir(msig_path)

msig_db <- msigdbr(species = msig_species) %>%
  dplyr::select(gs_name, entrez_gene)

# Run ORA for MSigDB
upregulated_msig_ora <- enricher(gene = upregulated_dges_entrez,
                                 TERM2GENE = msig_db,
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 0.05)

if (!is.null(upregulated_msig_ora) && !is.null(upregulated_msig_ora@result)) {
  upregulated_msig_ora@result <- upregulated_msig_ora@result[order(upregulated_msig_ora@result$p.adjust), ]
  write.csv(upregulated_msig_ora@result, file.path(msig_path, "upregulated_dges_msig_ora.csv"), row.names = FALSE)
}

downregulated_msig_ora <- enricher(gene = downregulated_dges_entrez,
                                   TERM2GENE = msig_db,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff = 0.05)

if (!is.null(downregulated_msig_ora) && !is.null(downregulated_msig_ora@result)) {
  downregulated_msig_ora@result <- downregulated_msig_ora@result[order(downregulated_msig_ora@result$p.adjust), ]
  write.csv(downregulated_msig_ora@result, file.path(msig_path, "downregulated_dges_msig_ora.csv"), row.names = FALSE)
}

# Run GSEA for MSigDB
msig_gsea <- GSEA(geneList = setNames(ranked_dges_entrez_dataframe$log2FoldChange, ranked_dges_entrez_dataframe$ENTREZID),
                  TERM2GENE = msig_db,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  eps = 0,
                  verbose = TRUE)

if (!is.null(msig_gsea) && !is.null(msig_gsea@result)) {
  msig_gsea@result <- msig_gsea@result[order(msig_gsea@result$p.adjust), ]
  write.csv(msig_gsea@result, file.path(msig_path, "dges_msig_gsea.csv"), row.names = FALSE)
}

#------------------------- Graph Functions -------------------------#

width_all = 2000
height_all = 1000

create_dir(file.path(root_path, "GO/graphs"))
create_dir(file.path(root_path, "KEGG/graphs"))
create_dir(file.path(root_path, "DO/graphs"))
create_dir(file.path(root_path, "WP/graphs"))
create_dir(file.path(root_path, "MSIG/graphs"))

# Dotplot - Done
generateDotPlot <- function(enrichment_object, title, file_path) {
  if (nrow(enrichment_object) == 0) {
    return(FALSE)
  } else {
    png(paste0(file_path, "_dot.png"), width = width_all, height = height_all, res = 100)
    p <- dotplot(enrichment_object, showCategory = 20, title = title)
    print(p)
    dev.off()

    # Save data to CSV
    write.csv(p$data, paste0(file_path, "_dot.csv"), row.names = FALSE)
  }
}

generateCnetPlot <- function(enrichment_object, title, file_path, fold_change) {
  if (enrichment_object@readable == FALSE) {
    enrichment_readable <- setReadable(enrichment_object, orgDb, 'ENTREZID')
  } else {
    enrichment_readable <- enrichment_object
  }

  png(paste0(file_path, "_cnet.png"), width = width_all, height = height_all, res = 100)
  p <- cnetplot(enrichment_readable, showCategory = 10, color.params = list(foldChange = fold_change))
  print(p)
  dev.off()

  # Save data to CSV
  edge_data <- p$data
  write.csv(edge_data, paste0(file_path, "_cnet.csv"), row.names = FALSE)
}

generateHeatPlot <- function(enrichment_object, title, file_path, fold_change) {
  if (enrichment_object@readable == FALSE) {
    enrichment_readable <- setReadable(enrichment_object, orgDb, 'ENTREZID')
  } else {
    enrichment_readable <- enrichment_object
  }

  png(paste0(file_path, "_heat.png"), width = width_all, height = height_all, res = 100)
  p <- heatplot(enrichment_readable, showCategory = 10, foldChange = fold_change)
  print(p)
  dev.off()

  # Save data to CSV
  heat_data <- p$data
  write.csv(heat_data, paste0(file_path, "_heat.csv"), row.names = FALSE)
}

generateTreePlot <- function(enrichment_object, title, file_path) {
  enrichment_object_pairwise <- pairwise_termsim(enrichment_object)
  png(paste0(file_path, "_tree.png"), width = width_all, height = height_all, res = 100)
  p <- treeplot(enrichment_object_pairwise, cluster.params = list(n = 5, method = "ward.D", label_words_n = 5))
  print(p)
  dev.off()

  # Save data to CSV
  tree_data <- p$data
  write.csv(tree_data, paste0(file_path, "_tree.csv"), row.names = FALSE)
}

generateEmapPlot <- function(enrichment_object, title, file_path) {
  enrichment_termsim <- pairwise_termsim(enrichment_object)
  png(paste0(file_path, "_emap.png"), width = width_all, height = height_all, res = 100)
  p <- emapplot(enrichment_termsim, showCategory = 30)
  print(p)
  dev.off()

  # Save data to CSV
  emap_data <- p$data
  write.csv(emap_data, paste0(file_path, "_emap.csv"), row.names = FALSE)
}


generateBarPlot <- function(enrichment_object, title, file_path) {
  if (nrow(enrichment_object) == 0 || inherits(enrichment_object, "gseaResult")) {
    cat("Skipping barplot for", title, "- empty result or GSEA result\n")
    return(FALSE)
  } else {
    png(paste0(file_path, "_bar.png"), width = width_all, height = height_all, res = 100)
    p <- barplot(enrichment_object, showCategory=20) +
      labs(title = title)
    print(p)
    dev.off()

    bar_data <- p$data
    write.csv(bar_data, paste0(file_path, "_bar.csv"), row.names = FALSE)

    # Generate and save an additional plot with -log10(p.adjust) as bar height
    png(paste0(file_path, "_bar_qscore.png"), width = width_all, height = height_all, res = 100)
    p_qscore <- enrichment_object %>%
      mutate(qscore = -log(p.adjust, base=10)) %>%
      barplot(x="qscore") +
      labs(title = paste(title, "- QScore"))
    print(p_qscore)
    dev.off()

    # Save qscore data to CSV
    qscore_data <- p_qscore$data
    write.csv(qscore_data, paste0(file_path, "_bar_qscore.csv"), row.names = FALSE)
  }
}

generateRidgePlot <- function(enrichment_object, title, file_path) {
  if (nrow(enrichment_object) == 0 || !inherits(enrichment_object, "gseaResult")) {
    cat("Skipping ridgeplot for", title, "- not a GSEA result or empty result\n")
    return(FALSE)
  } else {
    png(paste0(file_path, "_ridge.png"), width = width_all, height = height_all, res = 100)
    p <- ridgeplot(enrichment_object) + labs(x = "enrichment distribution", title = title)
    print(p)
    dev.off()

    # Save data to CSV
    ridge_data <- p$data
    write.csv(ridge_data, paste0(file_path, "_ridge.csv"), row.names = FALSE)
  }
}

generateAllPlots <- function(dge_data, title, output_path, ranked_dge) {
  if (nrow(dge_data) == 0) {
    cat("Skipping plots for", title, "due to empty results\n")
    return()
  }
  try(generateBarPlot(dge_data, title, output_path), silent = TRUE)
  try(generateDotPlot(dge_data, title, output_path), silent = TRUE)
  try(generateRidgePlot(dge_data, title, output_path), silent = TRUE)
  try(generateCnetPlot(dge_data, title, output_path, ranked_dge), silent = TRUE)
  try(generateHeatPlot(dge_data, title, output_path, ranked_dge), silent = TRUE)
  try(generateTreePlot(dge_data, title, output_path), silent = TRUE)
  try(generateEmapPlot(dge_data, title, output_path), silent = TRUE)

}

#------------------------- Plot Graphs -------------------------#

# GO
generateAllPlots(upregulated_dges_ora_go, "GO Upregulated Genes (ORA)", file.path(root_path, "GO/graphs/upregulated_dges_go_ora"), ranked_dges)
generateAllPlots(downregulated_dges_ora_go, "GO Downregulated Genes (ORA)", file.path(root_path, "GO/graphs/downregulated_dges_go_ora"), ranked_dges)
generateAllPlots(dges_gsea_go, "GO All Genes (GSEA)", file.path(root_path, "GO/graphs/dges_go_gsea"), ranked_dges)

# KEGG
generateAllPlots(upregulated_dges_ora_kegg, "KEGG Upregulated Genes (ORA)", file.path(root_path, "KEGG/graphs/upregulated_dges_kegg_ora"), ranked_dges)
generateAllPlots(downregulated_dges_ora_kegg, "KEGG Downregulated Genes (ORA)", file.path(root_path, "KEGG/graphs/downregulated_dges_kegg_ora"), ranked_dges)
generateAllPlots(dges_gsea_kegg, "KEGG All Genes (GSEA)", file.path(root_path, "KEGG/graphs/dges_kegg_gsea"), ranked_dges)

# DO
generateAllPlots(upregulated_dges_ora_do, "DO Upregulated Genes (ORA)", file.path(root_path, "DO/graphs/upregulated_dges_do_ora"), ranked_dges)
generateAllPlots(downregulated_dges_ora_do, "DO Downregulated Genes (ORA)", file.path(root_path, "DO/graphs/downregulated_dges_do_ora"), ranked_dges)
generateAllPlots(dges_gsea_do, "DO All Genes (GSEA)", file.path(root_path, "DO/graphs/dges_do_gsea"), ranked_dges)

# WP
generateAllPlots(upregulated_dges_ora_wp, "WP Upregulated Genes (ORA)", file.path(root_path, "WP/graphs/upregulated_dges_wp_ora"), ranked_dges)
generateAllPlots(downregulated_dges_ora_wp, "WP Downregulated Genes (ORA)", file.path(root_path, "WP/graphs/downregulated_dges_wp_ora"), ranked_dges)
generateAllPlots(dges_gsea_wp, "WP All Genes (GSEA)", file.path(root_path, "WP/graphs/dges_wp_gsea"), ranked_dges)

# MSigDB
generateAllPlots(upregulated_msig_ora, "MSigDB Upregulated Genes (ORA)", file.path(root_path, "MSIG/graphs/upregulated_dges_msig_ora"), ranked_dges)
generateAllPlots(downregulated_msig_ora, "MSigDB Downregulated Genes (ORA)", file.path(root_path, "MSIG/graphs/downregulated_dges_msig_ora"), ranked_dges)
generateAllPlots(msig_gsea, "MSigDB All Genes (GSEA)", file.path(root_path, "MSIG/graphs/dges_msig_gsea"), ranked_dges)

