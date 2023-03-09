###Amphioxus scRNAseq analysis: load data for control (of DAPT) 10hpf, 18hpf, 26hpf
library(Matrix)
library(stringr)
library(SingleCellExperiment)
library(scater)
library(scran)
library(msigdbr)
library(igraph)
library(viridis)
library(ggsci)
library(ggthemes)
library(paletteer)
library(DropletUtils)

setwd("~/Dropbox (Cambridge University)/PhD/scRNAseq/2021/Embryos")

#BEST:  10h SITTA2
#       18h SITTB12
#       26h SITTE5

sample_name <- "SITTA2_C_10hpf"
data_path <- paste0(getwd(),"/",sample_name,
                    "/filtered_feature_bc_matrix")
sce <- read10xCounts(data_path)

# Normalise
clusts <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusts)
sce <- logNormCounts(sce)

# Compute UMAP
sce <- runUMAP(sce)

# Plot expression
#p1 <- plotReducedDim(sce,"TSNE",colour_by="BL10935") + scale_color_viridis_c(direction=-1) 
p1 <- plotReducedDim(sce,"UMAP",colour_by="BL23554") + scale_color_viridis_c(direction=-1)
p1

#loop
Genes <- read.csv2("InterestingGenes_Names.csv")
rownames(Genes) <- Genes$Code
setwd(paste0("~/Dropbox (Cambridge University)/PhD/scRNAseq/2021/Embryos/",sample_name,"/Plots"))

for (i in rownames(Genes)) {
  p <- plotReducedDim(sce, dimred="UMAP", point_size = 4, colour_by = i) +
    scale_fill_viridis_c(direction = -1) +
    ggtitle(Genes[i,1])
  ggsave(paste(Genes[i,1], "umap.png"), p, width = 8, height = 5, dpi = 300)
}

#loop for chargenes
Chargenes <- read.csv2("characgenes.csv")
rownames(Chargenes) <- Chargenes$Code
setwd(paste0("~/Dropbox (Cambridge University)/PhD/scRNAseq/2021/Embryos/",sample_name,"/Plots"))

for (i in rownames(Chargenes)) {
  p <- plotReducedDim(sce, dimred="UMAP", point_size = 4, colour_by = i) +
    scale_fill_viridis_c(direction=-1) +
    ggtitle(Chargenes[i,1])
  ggsave(paste(Chargenes[i,1], "umap.png"), p, width = 8, height = 8, dpi = 300)
}


###quick analysis
#feature selection using variance
dec.G11 <- modelGeneVar(sce)
fit.G11 <- metadata(dec.G11)
plot(fit.G11$mean, fit.G11$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.G11$trend(x), col="dodgerblue", add=TRUE, lwd=2)
hvg.G11.var <- getTopHVGs(dec.G11, n=round(dim(sce)[1]*10/100, digits=0))
G11.hvg <- sce[hvg.G11.var,]
altExp(G11.hvg, "original") <- sce
#PCA
set.seed(100) # See below.
G11.hvg <- runPCA(G11.hvg, ncomponents = 40)
plotReducedDim(G11.hvg, dimred="PCA", colour_by = "BL18685")
#UMAP
set.seed(0010100110)
G11.hvg <- runUMAP(G11.hvg, dimred="PCA")
plotReducedDim(G11.hvg, dimred="UMAP", colour_by = "BL19810")
#tsne
set.seed(00101001)
G11.hvg <- runTSNE(G11.hvg, dimred="PCA")
plotReducedDim(G11.hvg, dimred="TSNE", colour_by = "BL13340")


for (i in rownames(Genes)) {
  p <- plotReducedDim(G11.hvg, dimred="UMAP", point_size = 4, colour_by = i) +
    scale_fill_viridis_c(direction = -1) +
    ggtitle(Genes[i,1])
  ggsave(paste(Genes[i,1], "umap.png"), p, width = 8, height = 5, dpi = 300)
}


#cluster
g <- buildSNNGraph(G11.hvg, k=6, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
##add to metadata and plot the clusters
G11.hvg$cluster <- factor(clust)
plott <- plotReducedDim(G11.hvg, "UMAP", colour_by="cluster")
plott

ggsave("clusters_umap_SITTA2_10hpf.png", plott, width = 8, height = 5, dpi = 300)

#search for hvg for each cluster
r_markers_some <- scran::findMarkers(G11.hvg, groups=G11.hvg$cluster, pval.type="all")
data.frame(head(r_markers_some[[5]][1:2], 50))



###select specific clusters
#colData(G11.hvg)["neural"] <- Merged_lognorm_20220122$lineage %in% "neural ectoderm"
#colnames(colData(Merged_lognorm_20220122))

##neural
neural <- G11.hvg[, G11.hvg$cluster %in% c("6","7")]
set.seed(1000)
neural <- runPCA(neural)
neural <- runTSNE(neural)
set.seed(10011)
neural <- runUMAP(neural)
plotReducedDim(neural, dimred="UMAP", colour_by="BL24510") + scale_fill_viridis_c(direction = -1)

#recluster
g <- buildSNNGraph(neural, k=5, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
#add to metadata and plot the clusters
neural$cluster <- factor(clust)
nplott <- plotReducedDim(neural, "UMAP", point_size = 4, colour_by="cluster")
nplott

ggsave("neuralclusters_umap_SITTA2_10hpf.png", nplott, width = 6, height = 5, dpi = 300)

r_markers_some <- scran::findMarkers(neural, groups=neural$cluster, pval.type="all")
data.frame(head(r_markers_some[[1]][1:2], 50))

#candidate gene expression
ngene <- plotReducedDim(neural, dimred="UMAP", point_size = 4, colour_by="BL11117") + scale_fill_viridis_c(direction = -1)
ngene
ggsave("Elav_neural_umap_SITTA2_10hpf.png", ngene, width = 6, height = 5, dpi = 300)


##epidermis
skin <- G11.hvg[, G11.hvg$cluster %in% c("4")]
set.seed(1000)
skin <- runPCA(skin)
skin <- runTSNE(skin)
set.seed(1005)
skin <- runUMAP(skin)

sgene <- plotReducedDim(skin, dimred="UMAP", point_size = 4, colour_by="BL19810") + scale_fill_viridis_c(direction = -1)
sgene
ggsave("Frz58_epidermis_umap_SITTA2_10hpf.png", sgene, width = 6, height = 5, dpi = 300)


#cluster
g <- buildSNNGraph(skin, k=8, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
##add to metadata and plot the clusters
skin$cluster <- factor(clust)
eplott <- plotReducedDim(skin, "UMAP", colour_by="cluster")
eplott

ggsave("epidermisclusters_umap_SITTA2_10hpf.png", eplott, width = 6, height = 5, dpi = 300)

r_markers_some <- scran::findMarkers(skin, groups=skin$cluster, pval.type="all")
data.frame(head(r_markers_some[[4]][1:2], 50))

