
library(SQMtools)
library(vegan)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggrepel)
library(clusterProfiler)
library(pheatmap)
library(UpSetR)
library(enrichplot)
library(KEGGREST)
library(cowplot)
library(httr)


outpath <- "/Users/sarapour/Desktop/squeeze/deseq/plots_stream1"

stream1 <- loadSQM("/Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream1")

# at least 20 counts
stream1_abundances <- stream1$functions$KEGG$abund[rowSums(stream1$functions$KEGG$abund) >= 20, ]

# metadata
md_stream1 <- data.frame(id = colnames(stream1_abundances))
md_stream1$condition <- ifelse(grepl("_DNA", md_stream1$id),
                               sub("_rep[0-9]+(_DNA)?$", "_DNA", md_stream1$id),
                               sub("_rep[0-9]+$", "", md_stream1$id))



# DEseq2 needs them to be the same to work, but doesnt auto check itself !!!!
stopifnot(identical(md_stream1$id, colnames(stream1_abundances)))


dds_stream1 <- DESeq(DESeqDataSetFromMatrix(countData = stream1_abundances, 
                    colData = md_stream1, 
                    design = ~ condition, 
                    ))

pdf(file = file.path(outpath, "dispersion_estimates_stream1.pdf"), width = 8, height = 6)
plotDispEsts(dds_stream1) 
dev.off()  

padj_threshold <- 0.05
log2FC_threshold <- 0.58  # FC of 1.5

# define comparisons
comparisons <- data.frame(
  comparison_name = c(
    "stream1_BLICH_T2_Control_vs_stream1_3BM_T2_Control",
    "stream1_BLICH_T4_Control_vs_stream1_3BM_T4_Control",
    "stream1_3BM_T2_Control_vs_stream1_3P_T2_Control",
    "stream1_3BM_T4_Control_vs_stream1_3P_T4_Control",
    "stream1_BLICH_T2_Control_vs_stream1_3P_T2_Control",
    "stream1_BLICH_T4_Control_vs_stream1_3P_T4_Control",
    "stream1_BSUB_T2_Control_vs_stream1_3BM_T2_Control",
    "stream1_BSUB_T4_Control_vs_stream1_3BM_T4_Control",
    "stream1_BSUB_T2_Control_vs_stream1_3P_T2_Control",
    "stream1_BSUB_T4_Control_vs_stream1_3P_T4_Control",
    "stream1_PA_T2_Control_vs_stream1_3BM_T2_Control",
    "stream1_PA_T4_Control_vs_stream1_3BM_T4_Control",
    "stream1_PA_T2_Control_vs_stream1_3P_T2_Control",
    "stream1_PA_T4_Control_vs_stream1_3P_T4_Control"
  ),
  denominator = c(
    "stream1_BLICH_T2_Control",
    "stream1_BLICH_T4_Control",
    "stream1_3BM_T2_Control",
    "stream1_3BM_T4_Control",
    "stream1_BLICH_T2_Control",
    "stream1_BLICH_T4_Control",
    "stream1_BSUB_T2_Control",
    "stream1_BSUB_T4_Control",
    "stream1_BSUB_T2_Control",
    "stream1_BSUB_T4_Control",
    "stream1_PA_T2_Control",
    "stream1_PA_T4_Control",
    "stream1_PA_T2_Control",
    "stream1_PA_T4_Control"
  ),
  numerator = c(
    "stream1_3BM_T2_Control",
    "stream1_3BM_T4_Control",
    "stream1_3P_T2_Control",
    "stream1_3P_T4_Control",
    "stream1_3P_T2_Control",
    "stream1_3P_T4_Control",
    "stream1_3BM_T2_Control",
    "stream1_3BM_T4_Control",
    "stream1_3P_T2_Control",
    "stream1_3P_T4_Control",
    "stream1_3BM_T2_Control",
    "stream1_3BM_T4_Control",
    "stream1_3P_T2_Control",
    "stream1_3P_T4_Control"
  ),
  stringsAsFactors = FALSE
)

perform_comparisons_explicit <- function(comparisons_df, dds) {
  group_results <- list()
  for (i in 1:nrow(comparisons_df)) {
    comp_name <- comparisons_df$comparison_name[i]
    numerator <- comparisons_df$numerator[i]
    denominator <- comparisons_df$denominator[i]

    if (!(numerator %in% levels(dds$condition))) {
      stop(paste("Condition", numerator, "not found in the dataset."))
    }
    if (!(denominator %in% levels(dds$condition))) {
      stop(paste("Condition", denominator, "not found in the dataset."))
    }
    
    res <- results(dds, contrast = c("condition", numerator, denominator),
                   independentFiltering = TRUE, alpha = padj_threshold, pAdjustMethod = "BH", parallel = FALSE)
    
    res <- lfcShrink(dds, contrast = c("condition", numerator, denominator), res = res, type = "ashr")
    
    group_results[[comp_name]] <- res
  }
  return(group_results)
}

results_list_stream1 <- perform_comparisons_explicit(comparisons, dds_stream1)

# filter genes based on if present in both numerator and denominator conditions
filter_genes_for_comparison <- function(genes, dds, numerator_condition, denominator_condition, min_count=3) {
  norm_counts <- counts(dds, normalized=TRUE)

  numerator_samples <- which(colData(dds)$condition == numerator_condition)
  denominator_samples <- which(colData(dds)$condition == denominator_condition)
  
  #keep genes that have at least min count in numberator and denominator condiitons, so we're not comparing apples to oranges
  genes[genes %in% rownames(norm_counts) & 
          rowMeans(norm_counts[genes, numerator_samples, drop=FALSE]) >= min_count &
          rowMeans(norm_counts[genes, denominator_samples, drop=FALSE]) >= min_count]
}


# ko ids to good names
get_kegg_descriptions <- function(kegg_ids) {
  sapply(kegg_ids, function(id) {
    tryCatch({
      kegg_info <- keggGet(id)
      if (length(kegg_info) > 0) {
        kegg_info[[1]]$NAME[1]  
      } else {
        NA
      }
    }, error = function(e) NA)
  })
}

ko_definitions <- keggList("ko")
#loop over to do comparisons
for (i in 1:nrow(comparisons)) {
  comparison_name <- comparisons$comparison_name[i]
  numerator <- comparisons$numerator[i]
  denominator <- comparisons$denominator[i]
  print(paste("Comparing:", comparison_name))
  
  comparison_outdir <- file.path(outpath, comparison_name)
  if (!dir.exists(comparison_outdir)) dir.create(comparison_outdir, recursive = TRUE)
  
  res <- results_list_stream1[[comparison_name]]
  res_df <- as_tibble(res, rownames = "ko_id") %>%
    filter(!is.na(log2FoldChange) & !is.na(padj)) %>%
    mutate(
      regulation = case_when(
        padj < padj_threshold & log2FoldChange >= log2FC_threshold ~ "Upregulated",
        padj < padj_threshold & log2FoldChange <= -log2FC_threshold ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )

  ko_definitions <- keggList("ko")
  ko_map <- ko_definitions
  names(ko_map) <- gsub("ko:", "", names(ko_map))
  ko_map_clean <- gsub(";.*", "", ko_map)  
  
  # gene column with better names
  res_df <- res_df %>%
    mutate(gene = ifelse(ko_id %in% names(ko_map_clean), ko_map_clean[ko_id], ko_id))
  top_up_genes <- res_df %>%
    filter(regulation == "Upregulated") %>%
    arrange(padj) %>%
    slice(1:20) %>%
    pull(gene)
  
  top_down_genes <- res_df %>%
    filter(regulation == "Downregulated") %>%
    arrange(padj) %>%
    slice(1:20) %>%
    pull(gene)
  
  res_df <- res_df %>%
    mutate(genelabels = ifelse(gene %in% c(top_up_genes, top_down_genes), gene, ""))

  pvolcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(colour = regulation)) +
    geom_text_repel(aes(label = genelabels)) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    ggtitle(paste(numerator, "vs", denominator)) +
    xlab("log2 Fold Change") +
    ylab("-log10 Adjusted p-value") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, hjust = 0.5), axis.title = element_text(size = 12))
  
  ggsave(file.path(comparison_outdir, paste0("volcano_labels_", comparison_name, ".pdf")), plot = pvolcano, width = 10, height = 8)
  
  # ko ideas for gsekegg
  geneList <- setNames(res_df$log2FoldChange, res_df$ko_id)
  geneList <- sort(geneList, decreasing = TRUE)
  
  # gsea
  gsea_results <- gseKEGG(
    geneList = geneList,
    organism = "ko",
    keyType = "kegg",
    pvalueCutoff = padj_threshold,
    pAdjustMethod = "none",
    nPermSimple = 10000
  )
  
  #dotplot
  gsea_df <- as.data.frame(gsea_results)
  top_gsea <- gsea_df %>%
    arrange(p.adjust) %>%
    head(20)
  
  p_dotplot <- ggplot(top_gsea, aes(x = reorder(Description, NES), y = NES)) +
    geom_point(aes(size = setSize, color = p.adjust)) +
    scale_color_gradient(low = "blue", high = "red") +
    coord_flip() +
    xlab("Pathway") +
    ylab("Normalized Enrichment Score (NES)") +
    ggtitle(paste0("Top GSEA KEGG Pathways: ", numerator, " vs ", denominator)) +
    theme_minimal() +
    theme(plot.title = element_text(size = 14), axis.title = element_text(size = 12))
  
  ggsave(file.path(comparison_outdir, paste0("gsea_dotplot_", comparison_name, ".pdf")), plot = p_dotplot, width = 11, height = 9)
  
  #kegg path enrichment
  up_genes <- res_df %>%
    filter(regulation == "Upregulated") %>%
    pull(ko_id)
  
  down_genes <- res_df %>%
    filter(regulation == "Downregulated") %>%
    pull(ko_id)
  
  if (length(up_genes) > 0) {
    kegg_up_results <- enrichKEGG(
      gene = up_genes,
      organism = "ko",
      keyType = "kegg",
      pvalueCutoff = padj_threshold,
      pAdjustMethod = "BH"
    )
    kegg_up_df <- as.data.frame(kegg_up_results) %>%
      arrange(p.adjust) %>%
      head(20)
    
    p_up <- ggplot(kegg_up_df, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "red") +
      coord_flip() +
      xlab("Pathway") +
      ylab("-log10 Adjusted p-value") +
      ggtitle("Top Enriched KEGG Pathways (Upregulated)") +
      theme_minimal()
    
    ggsave(file.path(comparison_outdir, paste0("kegg_up_barplot_", comparison_name, ".pdf")), plot = p_up, width = 11, height = 9)
  }
  
  if (length(down_genes) > 0) {
    kegg_down_results <- enrichKEGG(
      gene = down_genes,
      organism = "ko",
      keyType = "kegg",
      pvalueCutoff = padj_threshold,
      pAdjustMethod = "BH"
    )
    kegg_down_df <- as.data.frame(kegg_down_results) %>%
      arrange(p.adjust) %>%
      head(20)
    
    p_down <- ggplot(kegg_down_df, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
      geom_bar(stat = "identity", fill = "blue") +
      coord_flip() +
      xlab("Pathway") +
      ylab("-log10 Adjusted p-value") +
      ggtitle("Top Enriched KEGG Pathways (Downregulated)") +
      theme_minimal()
    
    ggsave(file.path(comparison_outdir, paste0("kegg_down_barplot_", comparison_name, ".pdf")), plot = p_down, width = 11, height = 9)
  }
  
  gsea_results_df <- gsea_results@result
  
  # top 5 leading edge genes, with good names
  for (row_idx in seq_len(nrow(gsea_results_df))) {
    les <- gsea_results_df$core_enrichment[row_idx]
    ko_ids <- strsplit(les, "/")[[1]]
    
    # match KO IDs back to res_df to get descriptive names & FCs
    matched <- res_df %>% filter(ko_id %in% ko_ids)
    
    if (nrow(matched) > 0) {
      matched_top5 <- matched %>%
        arrange(desc(log2FoldChange)) %>%
        slice(1:5)
      
      new_genes_str <- paste(matched_top5$gene, collapse = "/")
      
      gsea_results_df$core_enrichment[row_idx] <- new_genes_str
      gsea_results_df$geneID[row_idx] <- new_genes_str
    }
  }
  
  gsea_results@result <- gsea_results_df
  
  p_cnet <- cnetplot(gsea_results, circular = TRUE, colorEdge = TRUE)

  ggsave(
    filename = file.path(comparison_outdir, paste0("cnetplot_", comparison_name, ".pdf")),
    plot = p,
    width = 20, height = 8
  )

  p_total <- plot_grid(
    pvolcano,                                           
    p_dotplot,                                          
    plot_grid(p_up, p_down, p_cnet, ncol = 3,          
              rel_widths = c(1, 1, 1.5),                
              labels = c("C", "D", "E"),                
              label_size = 14),                         
    ncol = 1,                                          
    rel_heights = c(2, 2, 2.5)                          
  )
  ggsave(
    filename = file.path(comparison_outdir, paste0("summary_plot_", comparison_name, ".pdf")),
    plot = p_total,
    width = 24,  
    height = 16   
  )

}


