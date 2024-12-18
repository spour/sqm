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
library(reshape2)
library(scales)

outpath <- "/Users/sarapour/Desktop/squeeze/SqueezeOutputs/plots/"

stream1 <- loadSQM("/Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream1")
stream2 <- loadSQM("/Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream2")


top_n <- 15  

topKEGG = names(sort(rowSums(stream1$functions$KEGG$tpm), decreasing=TRUE))[1:11]
topKEGG = topKEGG[topKEGG!="Unclassified"]

stream1$misc$KEGG_names[topKEGG]

generate_heatmap <- function(data, keyword, stream_name, output_path, N = 20, fontsize = 15) {
  
  stream_signal <- subsetFun(data, keyword)
  
  # most abundant pathways based on TPM
    top_stream <- mostVariable(stream_signal$functions$KEGG$tpm, N=N)
  
  rownames(top_stream) <- paste(rownames(top_stream), stream_signal$misc$KEGG_names[rownames(top_stream)], sep = "; ")
  
  top_stream_long <- melt(top_stream)
  colnames(top_stream_long) <- c("Pathway", "Sample_Rep", "TPM")
  
  top_stream_long <- top_stream_long %>%
    mutate(
      Sample = sub("_rep[0-9]+(_DNA)?$", "", Sample_Rep),
      Replicate = factor(sub(".*_(rep[0-9]+(_DNA)?)$", "\\1", Sample_Rep))
    )
  

  plot <- ggplot(top_stream_long, aes(x = Replicate, y = Pathway, fill = TPM)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = median(top_stream_long$TPM)
    ) +
    labs(
      title = paste("Top", N, keyword, "Pathways for", stream_name),
      x = "Replicate",
      y = "Pathway",
      fill = "TPM"
    ) +
    facet_grid(~ Sample, scales = "free_x", space = "free_x") +
    coord_cartesian(clip = "off") + 
    scale_x_discrete(labels = scales::label_wrap(10)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, size = fontsize, hjust = 1),
      axis.text.y = element_text(size = fontsize),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      strip.text.x = element_text(
        angle = 70,
        hjust = 0.5,
        vjust = 0.5,
        size = fontsize
      ),
      plot.margin = margin(t = 20, b = 10)
    )
  
  ggsave(filename = paste0(output_path, stream_name, "_", gsub(" ", "_", keyword), "_most_VAR_heatmap.png"),
         plot = plot, width = 30, height = 25, units = "in")
  
  message("heatmap saved successfully!")
}

KEYWORDS = c("Signal transduction", "Biofilm Formation", 
            "Bacterial Chemotaxis", "Bacterial Invasion of Epithelial Cells", 
            "Bacterial Toxins", "Bacterial Secretion Systems", 
            "Quorum Sensing", "Biofilm Formation - Multiple Species",
            "Oxidative phosphorylation
"
)

for (keyword in KEYWORDS) {
  generate_heatmap(stream1, keyword, "Stream1", outpath, N = top_n)
  generate_heatmap(stream2, keyword, "Stream2", outpath, N = top_n)
}

plotBins(stream1)
