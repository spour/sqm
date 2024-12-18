
# FOR STREAM2

library(SQMtools)
library(vegan)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(KEGGREST)
library(ggrepel)
library(clusterProfiler)
library(pheatmap)
library(UpSetR)
library(stringr)
library(cowplot)


stream2 <- loadSQM("/Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream2")
outpath <- "/Users/sarapour/Desktop/squeeze/deseq/plots_stream2"
stream2_abundances <- stream2$functions$KEGG$abund[rowSums(stream2$functions$KEGG$abund) >= 20, ]
comparison_outdir <- "/Users/sarapour/Desktop/squeeze/deseq/plots_stream2"
padj_threshold <- 0.05
log2FC_threshold <- 0.58


md_stream2 <- data.frame(id = colnames(stream2_abundances))
md_stream2$condition <- ifelse(grepl("_DNA", md_stream2$id),
                               sub("_rep[0-9]+(_DNA)?$", "_DNA", md_stream2$id),
                               sub("_rep[0-9]+$", "", md_stream2$id))

md_stream2 <- md_stream2 %>%
  mutate(
    condition = as.character(condition)
  ) %>%
  mutate(
    group = ifelse(
      grepl("^Deliverable[0-9]+", condition),
      str_extract(condition, "^Deliverable[0-9]+"),
      NA_character_
    ),
    sample = ifelse(
      grepl("^Deliverable[0-9]+", condition),
      sub("^Deliverable[0-9]+_", "", condition),
      NA_character_
    )
  ) %>%
  mutate(
    strain = ifelse(
      grepl("^stream2_", condition),
      str_split_fixed(condition, "_", 4)[, 2],
      ifelse(
        grepl("^Deliverable6", condition),
        str_split_fixed(sample, "_", 4)[, 1],
        ifelse(
          grepl("^Deliverable4", condition),
          str_split_fixed(sample, "_", 2)[, 1],
          NA_character_
        )
      )
    ),
    timepoint = ifelse(
      grepl("^stream2_", condition),
      str_split_fixed(condition, "_", 4)[, 3],
      ifelse(
        grepl("^Deliverable6", condition),
        str_split_fixed(sample, "_", 4)[, 2],
        NA_character_
      )
    ),
    treatment = ifelse(
      grepl("^stream2_", condition),
      str_split_fixed(condition, "_", 4)[, 4],
      ifelse(
        grepl("^Deliverable6", condition),
        str_split_fixed(sample, "_", 4)[, 3],
        NA_character_
      )
    )
  )

md_stream2$condition <- factor(md_stream2$condition)
stopifnot(identical(md_stream2$id, colnames(stream2_abundances)))

dds_stream2 <- DESeq(DESeqDataSetFromMatrix(countData = stream2_abundances, colData = md_stream2, design = ~ condition))

pdf(file = file.path(outpath, "dispersion_estimates_stream2.pdf"), width = 8, height = 6)
plotDispEsts(dds_stream2) 
dev.off() 


unique_strains <- unique(md_stream2$strain[!is.na(md_stream2$strain)])
unique_treatments <- unique(md_stream2$treatment[!is.na(md_stream2$treatment)])
unique_conditions <- levels(md_stream2$condition)


comparisons <- data.frame(
  comparison_name = character(),
  numerator = character(),
  denominator = character(),
  stringsAsFactors = FALSE
)


for (strain in unique_strains) {
  control_condition <- paste0("stream2_", strain, "_T4_Control")
  if (!(control_condition %in% unique_conditions)) {
    next
  }
  for (treatment in unique_treatments) {
    if (treatment != "Control") {
      treatment_condition <- paste0("stream2_", strain, "_T4_", treatment)
      if (treatment_condition %in% unique_conditions) {
        comparison_name <- paste0(treatment_condition, "_vs_", control_condition)
        comparisons <- rbind(
          comparisons,
          data.frame(
            comparison_name = comparison_name,
            numerator = treatment_condition,
            denominator = control_condition,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
}


filter_genes_for_comparison <- function(genes, dds, numerator_condition, denominator_condition, min_count=3) {
  norm_counts <- counts(dds, normalized=TRUE)
  
  numerator_samples <- which(colData(dds)$condition == numerator_condition)
  denominator_samples <- which(colData(dds)$condition == denominator_condition)

  genes[genes %in% rownames(norm_counts) & 
          rowMeans(norm_counts[genes, numerator_samples, drop=FALSE]) >= min_count &
          rowMeans(norm_counts[genes, denominator_samples, drop=FALSE]) >= min_count]
}

perform_comparisons_explicit <- function(comparisons_df, dds, min_count=10) {
  group_results <- list()
  for (i in 1:nrow(comparisons_df)) {
    comp_name <- comparisons_df$comparison_name[i]
    numerator <- comparisons_df$numerator[i]
    denominator <- comparisons_df$denominator[i]

    if (!(numerator %in% levels(dds$condition))) {
      message(paste("condition", numerator, "not found. Skipping."))
      next
    }
    if (!(denominator %in% levels(dds$condition))) {
      message(paste("condition", denominator, "not found. Skipping."))
      next
    }
    filtered_genes <- filter_genes_for_comparison(
      rownames(dds),
      dds,
      numerator_condition = numerator,
      denominator_condition = denominator,
      min_count = min_count
    )

    res <- results(
      dds,
      contrast = c("condition", numerator, denominator),
      independentFiltering = TRUE,
      alpha = padj_threshold,
      pAdjustMethod = "BH"
    )

    res <- lfcShrink(
      dds,
      contrast = c("condition", numerator, denominator),
      res = res,
      type = "ashr"
    )

    group_results[[comp_name]] <- res
  }
  return(group_results)
}

results_list_stream2 <- perform_comparisons_explicit(comparisons, dds_stream2)
saveRDS(results_list_stream2, file = file.path(outpath, "results_list_stream2.rds"))

results_list_stream2 <- readRDS(file = file.path(outpath, "results_list_stream2.rds"))

ko2pathway <- keggLink("pathway", "ko")  

ko2pathway_df <- data.frame(
  KO = sub("ko:", "", names(ko2pathway)),  
  PathwayID = sub("path:", "", ko2pathway), 
  stringsAsFactors = FALSE
)

pathway_names <- keggList("pathway") 
pathway_df <- data.frame(
  PathwayID = sub("path:", "", names(pathway_names)), 
  NAME = unname(pathway_names) 
)

ko2pathway_df <- ko2pathway_df %>%
  mutate(PathwayID = sub("^ko", "map", PathwayID))

ko2pathway_df <- merge(ko2pathway_df, pathway_df, by = "PathwayID", all.x = TRUE) %>% distinct()

TERM2GENE <- ko2pathway_df %>%
  select(PathwayID, KO, NAME) %>%
  rename(TERM = PathwayID, GENE = KO, DESCRIPTION = NAME)

TERM2GENE$TERM <- sub("^map", "ko", TERM2GENE$TERM)
TERM2NAME <- TERM2GENE %>%
  select(TERM, DESCRIPTION) %>%
  distinct() %>%
  rename(NAME = DESCRIPTION)
saveRDS(TERM2GENE, file = file.path(outpath, "TERM2GENE_KO_to_Pathway.rds"))

TERM2GENE <- readRDS(file = file.path(outpath, "TERM2GENE_KO_to_Pathway.rds"))


for (i in 1:nrow(comparisons)) {
  comparison_name <- comparisons$comparison_name[i]
  numerator <- comparisons$numerator[i]
  denominator <- comparisons$denominator[i]
  print(paste("comparing:", comparison_name))

  comparison_outdir <- file.path(outpath, comparison_name)
  if (!dir.exists(comparison_outdir)) {
    dir.create(comparison_outdir, recursive = TRUE)
  }

  res <- results_list_stream2[[comparison_name]]

  if (is.null(res)) {
    message("no results for comparison: ", comparison_name)
    next
  }
  ko_definitions <- keggList("ko")
  ko_map <- ko_definitions
  names(ko_map) <- gsub("ko:", "", names(ko_map))
  ko_map_clean <- gsub(";.*", "", ko_map)  
  
  res_df <- as_tibble(res, rownames = "ko_id") %>%
    filter(!is.na(log2FoldChange) & !is.na(padj)) %>%
    mutate(
      regulation = case_when(
        padj < padj_threshold & log2FoldChange >= log2FC_threshold ~ "upregulated",
        padj < padj_threshold & log2FoldChange <= -log2FC_threshold ~ "downregulated",
        TRUE ~ "NS"
      )
    )
  
  res_df <- res_df %>%
    mutate(gene = ifelse(ko_id %in% names(ko_map_clean), ko_map_clean[ko_id], ko_id))

  # Identify top 20 up/down-regulated genes
  top_up_genes <- res_df %>%
    filter(regulation == "upregulated") %>%
    arrange(padj) %>%
    slice(1:20) %>%
    pull(ko_id)

  top_down_genes <- res_df %>%
    filter(regulation == "downregulated") %>%
    arrange(padj) %>%
    slice(1:20) %>%
    pull(ko_id)

  # Add gene labels for volcano plot (use KO IDs here, or map to names if available)
  res_df <- res_df %>%
    mutate(genelabels = ifelse(ko_id %in% c(top_up_genes, top_down_genes), ko_id, ""))

  # Volcano plot
  pvolcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(colour = regulation)) +
    geom_text_repel(aes(label = genelabels)) +
    scale_color_manual(values = c("upregulated" = "red", "downregulated" = "blue", "NS" = "grey")) +
    ggtitle(paste(numerator, "vs", denominator)) +
    xlab("log2FC") +
    ylab("-log10 p-vaL") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 12)
    )

  ggsave(file.path(comparison_outdir, paste0("volcano_labels_", comparison_name, ".pdf")), plot = pvolcano, width = 10, height = 8)

  up_genes <- res_df %>% filter(regulation == "upregulated") %>% pull(ko_id)
  down_genes <- res_df %>% filter(regulation == "downregulated") %>% pull(ko_id)

  if (length(up_genes) > 0) {
    up_genes_filtered <- intersect(up_genes, TERM2GENE$GENE)
    if (length(up_genes_filtered) > 0) {
      kegg_up_results <- enrichKEGG(
        gene = up_genes_filtered,
        organism = "ko",
        keyType = "kegg",
        pvalueCutoff = padj_threshold,
        pAdjustMethod = "BH"
      )

      if (!is.null(kegg_up_results) && nrow(kegg_up_results@result) > 0) {
        kegg_up_df <- as.data.frame(kegg_up_results) %>%
          arrange(p.adjust) %>%
          head(20)

        p_up <- ggplot(kegg_up_df, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
          geom_bar(stat = "identity", fill = "red") +
          coord_flip() +
          xlab("pathway") +
          ylab("-log10 p-val") +
          ggtitle("KEGG pathways (upregulated)") +
          theme_minimal()

        ggsave(file.path(comparison_outdir, paste0("kegg_up_barplot_", comparison_name, ".pdf")), plot = p_up, width = 11, height = 9)
        saveRDS(kegg_up_results, file = file.path(comparison_outdir, paste0("kegg_up_results_", comparison_name, ".rds")))
      } else {
        message("no significant enriched pathways for upregulated genes in comparison: ", comparison_name)
      }
    } else {
      message("no upregulated genes match TERM2GENE KO IDs in comparison: ", comparison_name)
    }
  } else {
    message("no upregulated genes for comparison: ", comparison_name)
  }

  if (length(down_genes) > 0) {
    down_genes_filtered <- intersect(down_genes, TERM2GENE$GENE)
    if (length(down_genes_filtered) > 0) {
      kegg_down_results <- enrichKEGG(
        gene = down_genes_filtered,
        organism = "ko",
        keyType = "kegg",
        pvalueCutoff = padj_threshold,
        pAdjustMethod = "BH"
      )

      if (!is.null(kegg_down_results) && nrow(kegg_down_results@result) > 0) {
        kegg_down_df <- as.data.frame(kegg_down_results) %>%
          arrange(p.adjust) %>%
          head(20)

        p_down <- ggplot(kegg_down_df, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
          geom_bar(stat = "identity", fill = "blue") +
          coord_flip() +
          xlab("pathway") +
          ylab("-log10 p-val") +
          ggtitle("KEGG pathways (downregulated)") +
          theme_minimal()

        ggsave(file.path(comparison_outdir, paste0("kegg_down_barplot_", comparison_name, ".pdf")), plot = p_down, width = 11, height = 9)
        saveRDS(kegg_down_results, file = file.path(comparison_outdir, paste0("kegg_down_results_", comparison_name, ".rds")))
      } else {
        message("no significant enriched pathways for downregulated genes in comparison: ", comparison_name)
      }
    } else {
      message("no downregulated genes match TERM2GENE KO IDs in comparison: ", comparison_name)
    }
  } else {
    message("no downregulated genes for comparison: ", comparison_name)
  }

  # GSEA analysis
  geneList <- setNames(res_df$log2FoldChange, res_df$ko_id)
  geneList <- sort(geneList, decreasing = TRUE)
  geneList <- geneList[names(geneList) %in% TERM2GENE$GENE]

  if (length(geneList) > 0) {
    gsea_results <- gseKEGG(
      geneList = geneList,
      organism = "ko",
      keyType = "kegg",
      pvalueCutoff = padj_threshold,
      pAdjustMethod = "none",
      nPermSimple = 10000
    )

    if (!is.null(gsea_results) && nrow(gsea_results@result) > 0) {
      gsea_df <- as.data.frame(gsea_results) %>%
        arrange(p.adjust) %>%
        head(20)

      p_dotplot <- ggplot(gsea_df, aes(x = reorder(Description, NES), y = NES)) +
        geom_point(aes(size = setSize, color = p.adjust)) +
        scale_color_gradient(low = "blue", high = "red") +
        coord_flip() +
        xlab("pathway") +
        ylab("NES") +
        ggtitle(paste0("top GSEA KEGG pathways: ", numerator, " vs ", denominator)) +
        theme_minimal() +
        theme(plot.title = element_text(size = 14), axis.title = element_text(size = 12))

      ggsave(file.path(comparison_outdir, paste0("gsea_dotplot_", comparison_name, ".pdf")), plot = p_dotplot, width = 11, height = 9)

      gsea_results_df <- gsea_results@result
      for (row_idx in seq_len(nrow(gsea_results_df))) {
        les <- gsea_results_df$core_enrichment[row_idx]
        ko_ids <- strsplit(les, "/")[[1]]
        
        matched <- res_df %>% filter(ko_id %in% ko_ids)
        if (nrow(matched) > 0) {
          matched_top5 <- matched %>%
            arrange(desc(log2FoldChange)) %>%
            slice(1:5)
          new_genes_str <- paste(matched_top5$ko_id, collapse = "/")
          gsea_results_df$core_enrichment[row_idx] <- new_genes_str
          gsea_results_df$geneID[row_idx] <- new_genes_str
        }
      }
      gsea_results@result <- gsea_results_df

      p_cnet <- cnetplot(gsea_results, circular = TRUE, colorEdge = TRUE)
      ggsave(
        filename = file.path(comparison_outdir, paste0("cnetplot_", comparison_name, ".pdf")),
        plot = p_cnet,
        width = 20, height = 8
      )
    } else {
      message("no significant pathways found in GSEA for comparison: ", comparison_name)
    }
  } else {
    message("gene list is empty after filtering for TERM2GENE in comparison: ", comparison_name)
  }

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
