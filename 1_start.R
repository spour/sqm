# FILES
# /Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream2 
# /Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream1 
# /Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDMBP
# if (!require("BiocManager", quietly = TRUE)) { install.packages("BiocManager")}; BiocManager::install("SQMtools")

# library("SQMtools")
# Hadza = loadSQM("/Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream1")
# plotTaxonomy(Hadza, rank='phylum', count='percent')

# proteo=subsetTax(Hadza, 'phylum',tax='Firmicutes', rescale_copy_number = F)
# plotTaxonomy(proteo, 'genus','percent', N=10, rescale = T, others = T)

# plotFunctions(Hadza, fun_level = 'KEGG', count = 'tpm', N = 15)

library(SQMtools)
library(vegan)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

outpath <- "/Users/sarapour/Desktop/squeeze/SqueezeOutputs/plots/"


# load SQM data
stream1 <- loadSQM("/Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream1")
stream2 <- loadSQM("/Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream2")

# set parameters for top n 
top_n <- 15 

### seperate
plotTaxonomy(
  stream1, 
  rank = "phylum", 
  count = "percent", 
  N = 15
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),  
    legend.position = "right" 
  )

# ggsave to outpath
ggsave(filename = paste0(outpath, "stream1_phylum.png"), width = 10, height = 6, dpi = 300)


plotTaxonomy(
  stream2, 
  rank = "phylum", 
  count = "percent", 
  N = 15
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6), 
    legend.position = "right"  
  )

ggsave(filename = paste0(outpath, "stream2_phylum.png"), width = 15, height = 6, dpi = 300)


# most abundant taxa in Stream1
top_abundant <- mostAbundant(stream1$taxa$phylum$abund, N = top_n, rescale = TRUE)
plotBars(top_abundant, label_y = "Abundance", label_fill = "Phylum") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 
ggsave(filename = paste0(outpath, "stream1_abundant.png"), width = 10, height = 9, dpi = 300)

top_abundant <- mostAbundant(stream2$taxa$phylum$abund, N = top_n, rescale = TRUE)
plotBars(top_abundant, label_y = "Abundance", label_fill = "Phylum") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))  
ggsave(filename = paste0(outpath, "stream2_abundant.png"), width = 10, height = 9, dpi = 300)

## merge replicates (mayne not needed)
merge_all_replicates <- function(sqm_data) {
  # extract sample names
  sample_names <- colnames(sqm_data$taxa$phylum$abund)
  
  #find unique prefix
  unique_prefixes <- unique(gsub("(_rep[0-9]+.*)$", "", sample_names))
  
  merged_abund_list <- list()
  
  for (prefix in unique_prefixes) {
    pattern <- paste0("^", prefix, "_rep[0-9]+.*$")
    replicate_indices <- grep(pattern, sample_names)
    
    if (length(replicate_indices) > 0) {
      replicates_abund <- sqm_data$taxa$phylum$abund[, replicate_indices, drop = FALSE]
      
      #  mean abundance
      merged_abund <- rowMeans(replicates_abund, na.rm = TRUE)
      
      merged_abund_list[[prefix]] <- merged_abund
    } else {
      message("No samples found for prefix '", prefix, "'.")
    }
  }
  
  merged_abund_matrix <- do.call(cbind, merged_abund_list)
  colnames(merged_abund_matrix) <- names(merged_abund_list)
  
  sqm_data$taxa$phylum$abund <- merged_abund_matrix
  
  return(sqm_data)
}


stream1_merged <- merge_all_replicates(stream1)

top_abundant <- mostAbundant(stream1_merged$taxa$phylum$abund, N = top_n)
plotBars(top_abundant, label_y = "Abundance", label_fill = "Phylum") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6)) 
ggsave(filename = paste0(outpath, "stream1_abundant_merged.png"), width = 10, height = 9, dpi = 300)


non_ribosome_ids <- names(stream1$misc$KEGG_names)[
  !grepl("ribosome|ribosomal", stream1$misc$KEGG_names, ignore.case = TRUE)
]

# non-ribosomal functions
non_ribosomal_stream1 <- subsetFun(
  stream1,
  fun = non_ribosome_ids,
  fixed = TRUE,
  trusted_functions_only = FALSE
)

### kegg filtering
merge_kegg_replicates <- function(sqm_data) {
  merge_replicates <- function(data_matrix) {
    sample_names <- colnames(data_matrix)
    unique_prefixes <- unique(gsub("(_rep[0-9]+.*)$", "", sample_names))
    
    merged_data_list <- list()
    
    for (prefix in unique_prefixes) {
      pattern <- paste0("^", prefix, "_rep[0-9]+.*$")
      replicate_indices <- grep(pattern, sample_names)
      
      if (length(replicate_indices) > 0) {
        merged_data <- rowMeans(data_matrix[, replicate_indices, drop = FALSE], na.rm = TRUE)
        merged_data_list[[prefix]] <- merged_data
      } else {
        message("No samples found for prefix '", prefix, "'.")
      }
    }
    
    merged_data_matrix <- do.call(cbind, merged_data_list)
    colnames(merged_data_matrix) <- names(merged_data_list)
    
    return(merged_data_matrix)
  }

  if (!is.null(sqm_data$taxa$phylum$abund)) {
    sqm_data$taxa$phylum$abund <- merge_replicates(sqm_data$taxa$phylum$abund)
  }
  
  if (!is.null(sqm_data$functions$KEGG$tpm)) {
    sqm_data$functions$KEGG$tpm <- merge_replicates(sqm_data$functions$KEGG$tpm)
  }
  
  if (!is.null(sqm_data$functions$KEGG$abund)) {
    sqm_data$functions$KEGG$abund <- merge_replicates(sqm_data$functions$KEGG$abund)
  }
  if (!is.null(sqm_data$functions$KEGG$bases)) {
    sqm_data$functions$KEGG$bases <- merge_replicates(sqm_data$functions$KEGG$bases)
  }
  if (!is.null(sqm_data$functions$KEGG$cov)) {
    sqm_data$functions$KEGG$cov <- merge_replicates(sqm_data$functions$KEGG$cov)
  }
  if (!is.null(sqm_data$functions$KEGG$cpm)) {
    sqm_data$functions$KEGG$cpm <- merge_replicates(sqm_data$functions$KEGG$cpm)
  }
  if (!is.null(sqm_data$functions$KEGG$copy_number)) {
    sqm_data$functions$KEGG$copy_number <- merge_replicates(sqm_data$functions$KEGG$copy_number)
  }
  
  return(sqm_data)
}

filter_and_plot_kegg <- function(stream_data, name, outpath, top_n = 100, plot_n = 15, exclude_pattern = "ribosome|ribosomal") {
  # KEGG TPM data
  kegg_tpm <- stream_data$functions$KEGG$tpm  
  kegg_ids <- rownames(kegg_tpm)             
  sample_names <- colnames(kegg_tpm)          
  
  total_tpm <- rowSums(kegg_tpm)
  
  kegg_df <- data.frame(
    KEGG_ID = kegg_ids,
    Total_TPM = total_tpm,
    Description = stream_data$misc$KEGG_names[kegg_ids],
    stringsAsFactors = FALSE
  )
  
  top_kegg_df <- kegg_df %>%
    arrange(desc(Total_TPM)) %>%
    head(top_n)
  
  filtered_kegg_df <- top_kegg_df %>%
    filter(!grepl(exclude_pattern, Description, ignore.case = TRUE))
  
  filtered_kegg_ids <- filtered_kegg_df$KEGG_ID
  
  filtered_kegg_tpm <- kegg_tpm[filtered_kegg_ids, ]
  
  # copy to  avoid modifying the original data
  filtered_stream_data <- stream_data
  filtered_stream_data$functions$KEGG$tpm <- filtered_kegg_tpm
  filtered_stream_data$functions$KEGG$abund <- stream_data$functions$KEGG$abund[filtered_kegg_ids, ]
  filtered_stream_data$functions$KEGG$bases <- stream_data$functions$KEGG$bases[filtered_kegg_ids, ]
  filtered_stream_data$functions$KEGG$cov <- stream_data$functions$KEGG$cov[filtered_kegg_ids, ]
  filtered_stream_data$functions$KEGG$cpm <- stream_data$functions$KEGG$cpm[filtered_kegg_ids, ]
  filtered_stream_data$functions$KEGG$copy_number <- stream_data$functions$KEGG$copy_number[filtered_kegg_ids, ]
  
  plot <- plotFunctions(
    filtered_stream_data,
    fun_level = "KEGG",
    count = "tpm",
    N = plot_n 
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
    )
  
  ggsave(
    filename = paste0(outpath, name),
    plot = plot,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  return(filtered_stream_data)
}

stream1_kegg_merged <- merge_kegg_replicates(stream1)
stream1_filtered <- filter_and_plot_kegg(
  stream1_kegg_merged,
  name = "stream1_kegg_filtered.png",
  outpath = outpath,
  top_n = 100,
  plot_n = top_n,
  exclude_pattern = "ribosome|ribosomal"
)

