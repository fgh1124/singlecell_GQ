library(Seurat)
library(dplyr)
totalcount1 <- read.table("SC_count_E35_E45.txt", sep="\t", header=TRUE, as.is = TRUE)
totalcount1 <- t(totalcount1)
colnames(totalcount1) <- totalcount1[1,]
totalcount1 <- totalcount1[-1,]
embryo_pub <- CreateSeuratObject(counts = totalcount1, project = "SeuratProject", min.cells = 10, min.features = 1000)
embryo_pub$tech <- "pub"
totalcount2 <- read.table("gene_count2.txt",sep="\t",header=TRUE, as.is = TRUE, row.names = 1)
embryo_our <- CreateSeuratObject(counts = totalcount2, project = "SeuratProject", min.cells = 10, min.features = 1000)
embryo_our$tech <- "our"

embryo <- merge(embryo_pub, y = embryo_our, add.cell.ids = c("pub", "our"), project = "embryobig")

data.list <- SplitObject(embryo, split.by = "tech")

for (i in 1:length(data.list)) {
    data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
    data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
}

reference.list <- data.list[c("our", "pub")]
data.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)
data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:50)

library(ggplot2)
library(cowplot)
DefaultAssay(data.integrated) <- "integrated"

all.genes <- rownames(data.integrated)

data.integrated <- ScaleData(data.integrated, features = all.genes)
embryo <- data.integrated
embryo <- RunPCA(embryo, features = VariableFeatures(object = embryo))
print(embryo[["pca"]], dims = 1:5, nfeatures = 5)
pdf ("pca.pdf")
VizDimLoadings(embryo, dims = 1:2, reduction = "pca")
dev.off()
pdf ("pca2.pdf")
DimPlot(embryo, reduction = "pca")
dev.off()

embryo <- JackStraw(embryo, num.replicate = 100)
embryo <- ScoreJackStraw(embryo, dims = 1:20)
pdf ("JackStrawPlot.pdf")
JackStrawPlot(embryo, dims = 1:15)
dev.off()

pdf ("ElbowPlot.pdf")
ElbowPlot(embryo)
dev.off()

embryo <- FindNeighbors(embryo, dims = 1:10)
embryo <- FindClusters(embryo, resolution = 0.05)

embryo <- RunUMAP(embryo, dims = 1:10)
#write.table (embryo$seurat_clusters, "cluster_info.txt")

pdf ("umap1.pdf")
DimPlot(embryo, reduction = "umap", pt.size = 1, label = T)
dev.off()
pdf ("umap2_N.pdf")
DimPlot(embryo, reduction = "umap", group.by = "orig.ident", pt.size = 1, label = F)
dev.off()
pdf ("umap3_N.pdf")
DimPlot(embryo, reduction = "umap", group.by = "tech", pt.size = 1, label = F)
dev.off()

new.cluster.ids <- c("ICM", "TE", "PrE", "epiblast")
names(new.cluster.ids) <- levels(embryo)
embryo <- RenameIdents(embryo, new.cluster.ids)
pdf ("umap4.pdf")
DimPlot(embryo, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
dev.off()

