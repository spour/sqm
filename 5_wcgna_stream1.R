
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


outpath <- "/Users/sarapour/Desktop/squeeze/deseq/plots_stream1"

stream1 <- loadSQM("/Users/sarapour/Desktop/squeeze/SqueezeOutputs/SQMToolsTables/SqueezeCoassemblyISEDStream1")

stream1_abundances <- stream1$functions$KEGG$abund[rowSums(stream1$functions$KEGG$abund) >= 20, ]

md_stream1 <- data.frame(id = colnames(stream1_abundances))
md_stream1$condition <- ifelse(grepl("_DNA", md_stream1$id),
                               sub("_rep[0-9]+(_DNA)?$", "_DNA", md_stream1$id),
                               sub("_rep[0-9]+$", "", md_stream1$id))




stopifnot(identical(md_stream1$id, colnames(stream1_abundances)))


dds_stream1 <- DESeq(DESeqDataSetFromMatrix(countData = stream1_abundances, 
                    colData = md_stream1, 
                    design = ~ condition, 
                    ))
#### function WGCNA
library(WGCNA)

vsd <- vst(dds_stream1, blind = FALSE)
expr_data <- assay(vsd)

datExpr <- t(expr_data)

# check for genes with too many missing values or low variance
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
# OUTLIER SAMPLES
sampleTree <- hclust(dist(datExpr), method = "average")

par(cex = 0.6);
par(mar = c(0,4,2,0))

pdf(file = file.path(outpath, "sample_clustering.pdf"), width = 12, height = 15)
plot(sampleTree, main = "sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
dev.off()
cut.sampleTree <- cutreeStatic(sampleTree, cutHeight = 83, minSize = 10) 

powers <- c(1:20)

sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="soft thresh (power)",
     ylab="sft fit,signed R^2",
     type="n")
# there is a red line at 0.9 which is is often the threshold for determining whether the network approximates a scale-free topology. A value of
# 0.9 or higher is considered good, suggesting it follows a scale-free topology well. 
# you can select power = 17 as your soft-thresholding power for constructing the adjacency matrix
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.90, col="red") 

softPower <- 17
adjacency <- adjacency(datExpr, power = softPower)

TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

geneTree <- hclust(as.dist(dissTOM), method = "average")

dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = 30)

dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors,
                    "dynamic tree cut",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05)
# each branch is group of genes, height on the dendrogram is the dissimilarity between the genes in the branch
# dynamic tree cut algorithm has identified distinct modules (clusters of co-expressed genes) and assigned them colors.

sampleTraits <- as.data.frame(colData(dds_stream1))


condition_dummies <- model.matrix(~ 0 + condition, data=sampleTraits)


rownames(condition_dummies) <- rownames(sampleTraits)

sampleTraits <- as.data.frame(condition_dummies)

sampleTraits <- sampleTraits[match(rownames(datExpr), rownames(sampleTraits)), , drop = FALSE]

dim(datExpr)      
dim(sampleTraits)   

MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes

MEs <- orderMEs(MEs)

moduleTraitCor <- cor(MEs, sampleTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf(file = file.path(outpath, "module_trait_relationships.pdf"), width = 16, height = 20)
par(mar = c(16, 15, 4, 2))  

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(sampleTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,  
               cex.text = 0.7,         
               zlim = c(-1,1),
               main = "module-trait relationships")

dev.off()


