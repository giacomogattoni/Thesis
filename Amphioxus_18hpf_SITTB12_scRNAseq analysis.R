#Amphioxus scRNAseq control (of DAPT) 18hpf
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

setwd("~/Dropbox (Cambridge University)/PhD/scRNAseq/2021/Embryos/SITTB12_C_18hpf")

#load
neu2021 <- readRDS("sce.rds")

c <- plotReducedDim(neu2021, dimred="UMAP",
                    point_size = 2) +
  scale_fill_paletteer_d("ggsci::default_ucscgb")
c

#single genes
plotReducedDim(neu2021, dimred="UMAP", colour_by = "BL23931") +
  scale_fill_viridis()

#loop
Genes <- read.csv2("InterestingGenes_Names.csv")
rownames(Genes) <- Genes$Code
setwd("~/Dropbox (Cambridge University)/PhD/scRNAseq/2021/Embryos/SITTB12_C_18hpf/Plots")

for (i in rownames(Genes)) {
  p <- plotReducedDim(neu2021, dimred="UMAP", point_size = 4, colour_by = i) +
    scale_fill_viridis() +
    ggtitle(Genes[i,1])
  ggsave(paste(Genes[i,1], "umap.png"), p, width = 8, height = 8, dpi = 300)
}


###quick analysis
#feature selection using variance
dec.G11 <- modelGeneVar(neu2021)
fit.G11 <- metadata(dec.G11)
plot(fit.G11$mean, fit.G11$var, xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(fit.G11$trend(x), col="dodgerblue", add=TRUE, lwd=2)
hvg.G11.var <- getTopHVGs(dec.G11, n=round(dim(neu2021)[1]*10/100, digits=0))
G11.hvg <- neu2021[hvg.G11.var,]
altExp(G11.hvg, "original") <- neu2021
#PCA
set.seed(100) # See below.
G11.hvg <- runPCA(G11.hvg, ncomponents = 40)
plotReducedDim(G11.hvg, dimred="PCA", colour_by = "BL18685")
#tsne
set.seed(00101001101)
G11.hvg <- runTSNE(G11.hvg, dimred="PCA")
plotReducedDim(G11.hvg, dimred="TSNE", colour_by = "BL13340")
#UMAP
set.seed(00101001101)
G11.hvg <- runUMAP(G11.hvg, dimred="PCA")
plotReducedDim(G11.hvg, dimred="UMAP", colour_by = "BL13340")


#cluster
g <- buildSNNGraph(G11.hvg, k=8, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
##add to metadata and plot the clusters
G11.hvg$cluster <- factor(clust)
plott <- plotReducedDim(G11.hvg, "UMAP", point_size = 4, colour_by="cluster")
plott

ggsave("clusters_new_umap.png", plott, width = 8, height = 5, dpi = 300)

#search for hvg for each cluster
r_markers_some <- scran::findMarkers(G11.hvg, groups=G11.hvg$cluster, pval.type="all")
data.frame(head(r_markers_some[[1]][1:2], 50))

#loop
for (i in rownames(Genes)) {
  p <- plotReducedDim(G11.hvg, dimred="UMAP", point_size = 4, colour_by = i) +
    scale_fill_viridis_c(direction = -1) +
    ggtitle(Genes[i,1])
  ggsave(paste(Genes[i,1], "umap.png"), p, width = 8, height = 5, dpi = 300)
}



###select specific clusters
#colData(G11.hvg)["neural"] <- Merged_lognorm_20220122$lineage %in% "neural ectoderm"
#colnames(colData(Merged_lognorm_20220122))

##neural
neural <- G11.hvg[, G11.hvg$cluster %in% c("4","6","9")]
set.seed(1000)
neural <- runPCA(neural)
neural <- runTSNE(neural)
set.seed(10011)
neural <- runUMAP(neural)
plotReducedDim(neural, dimred="UMAP", colour_by="BL13340") + scale_fill_viridis_c(direction = -1)

#recluster
g <- buildSNNGraph(neural, k=12, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
#add to metadata and plot the clusters
neural$ncluster <- factor(clust)
nplott <- plotReducedDim(neural, "UMAP", point_size = 4, colour_by="cluster")
nplott

ggsave("neuralclusters_umap_SITTB12_10hpf.png", nplott, width = 6, height = 5, dpi = 300)


Ngenes <- read.csv2("neural.csv")
rownames(Ngenes) <- Ngenes$Code

for (i in rownames(Ngenes)) {
  p <- plotReducedDim(neural, dimred="UMAP", point_size = 4, colour_by = i) +
    scale_fill_viridis_c(direction = -1) +
    ggtitle(Ngenes[i,1])
  ggsave(paste(Ngenes[i,1], "umap.png"), p, width = 5, height = 5, dpi = 300)
}




#candidate gene expression
ngene <- plotReducedDim(neural, dimred="UMAP", point_size = 4, colour_by="BL04895") + scale_fill_viridis_c(direction = -1)
ngene
ggsave("Noto_neural_umap_SITTB12_18hpf.png", ngene, width = 6, height = 5, dpi = 300)




r_markers_some <- scran::findMarkers(neural, groups=neural$cluster, pval.type="all")
data.frame(head(r_markers_some[[1]][1:2], 50))









#subset cluster

cl1 <- G11.hvg[, G11.hvg$cluster=="1"]

set.seed(100) 
cl1 <- runPCA(cl1, ncomponents = 5)

set.seed(00101001101)
cl1 <- runTSNE(cl1, dimred="PCA")
cl1 <- runUMAP(cl1, dimred="PCA")
plotReducedDim(cl1, dimred="TSNE", colour_by = "BL11117")

g <- buildSNNGraph(cl1, k=8, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)
##add to metadata and plot the clusters
cl1$cluster <- factor(clust)
plott <- plotReducedDim(cl1, dimred = "TSNE", colour_by="cluster")
plott

r_markers_some <- scran::findMarkers(cl1, groups=cl1$cluster, pval.type="all")
data.frame(head(r_markers_some[[2]][1:2], 20))


colData(cl1)


#select hvg
dec.G11.2 <- modelGeneVar(G11.hvg)
fit.G11.2 <- metadata(dec.G11.2)
dec.G11.2[order(dec.G11.2$bio, decreasing=TRUE),] 
hvg.G11.2.var <- getTopHVGs(dec.G11.2, n=round(dim(G11.hvg)[1]*10/100, digits=0))
str(hvg.G11.2.var)



