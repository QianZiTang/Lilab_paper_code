#We use the pipeline of one sample ID : YDJ1filter as an example.

# Load package
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)


# load ST data into Seurat using the Load10X_Spatial function. 
#This input data are directly from the output of the Spaceranger pipeline. 
YDJ1filter <- Load10X_Spatial(data.dir="/WORK/Project/zengbo/spatial/pig_muscle_analysis/count_filter/count-YDJ1filter-final/outs", assay="Spatial", slice="YDJ1filter", filter.matrix=TRUE)
samInfo <- rep("YDJ1filter",nrow(YDJ1filter@meta.data))  
YDJ1filter <- AddMetaData(object = YDJ1filter, metadata = samInfo, col.name = "samInfo") #Add sample name information to metadata


# Data quality control
#Statistics of each mitochondrial gene
YDJ1filter[["ND6"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ND6")
YDJ1filter[["ND5"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ND5")
YDJ1filter[["ND4L"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ND4L")
YDJ1filter[["ND4"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ND4")
YDJ1filter[["ND3"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ND3")
YDJ1filter[["ND2"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ND2")
YDJ1filter[["ND1"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ND1")
YDJ1filter[["CYTB"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "CYTB")
YDJ1filter[["COX3"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "COX3")
YDJ1filter[["COX2"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "COX2")
YDJ1filter[["COX1"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "COX1")
YDJ1filter[["ATP8"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ATP8")
YDJ1filter[["ATP6"]] <- PercentageFeatureSet(YDJ1filter, assay = "Spatial", feature = "ATP6")
YDJ1filter[["percent.mt"]] <- apply(YDJ1filter@meta.data[,5:17], 1, sum) #Sum to get the total proportion of mitochondrial genes
# Export metadate with the counts, gene number and mitochondrial gene proportion information of each spot in the sample
write.csv(YDJ1filter@meta.data, file = "YDJ1filter-spots-metadata-raw.csv")
#Make a violin chart to study data quality
p1 <- VlnPlot(YDJ1filter, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"), ncol = 3, pt.size = 0.1) + NoLegend()
ggsave("YDJ1filter-rawdata-Violin-plot.pdf", plot=p1, width = 8, height =8)
#
#Quantification of spots and genes were carried out by using subset function. 
#We filtered out ~2% of low-quality spots. Spots with over 25% mitochondrial gene expression were also discarded
feature_nums <- colSums(as.matrix(YDJ1filter@assays$Spatial@counts) >0)
mycut_feature <- as.numeric(quantile(feature_nums, 0.02))
count_nums <- colSums(as.matrix(YDJ1filter@assays$Spatial@counts))
mycut_count <- as.numeric(quantile(count_nums, 0.02))
YDJ1filter_f1 <- subset(YDJ1filter, subset = nFeature_Spatial>=mycut_feature | nCount_Spatial>=mycut_count & percent.mt < 25 )
#Genes expressed in fewer than 15 spots were excluded.
min.spots <- 15
num.spots <- rowSums(as.matrix(YDJ1filter_f1@assays$Spatial@counts) > 0)
genes.use <- names(num.spots[which(num.spots >= min.spots)])
mykeepgene <- c(1:nrow(YDJ1filter_f1))[rownames(YDJ1filter_f1)%in%as.character(genes.use)]
YDJ1filter_f2 <- subset(YDJ1filter_f1,features=mykeepgene)
YDJ1filter_f2
# Filter out contaminated genes with specified names.
#Genes related to hemoglobin (considerable variation from blood contents) and Y-chromosome linked genes were removed.
filter_genes <-c("MALAT1", "SLC4A1", "KDM5D", "ANK1", "DDX3Y", "EIF2AK1", "HBQ1", "FTL", "GATA1", "KLF1", "USP9Y", "NFE2", "MT1G", "RPS4Y1", "HBZ", "GYPC", "HEMGN", "SLC25A37", "ALAS2", "EPB41", "AHSP", "GYPA", "UTY", "HBA2", "HBG2", "EIF1AY", "HBA1", "HBM", "HBE1", "HBG1", "MTRNR2L4", "HBB", "MTRNR2L5", "MTRNR2L8", "MTRNR2L10", "MTRNR2L3", "MTRNR2L1", "MTRNR2L7", "MTRNR2L12", "MTRNR2L11", "MTRNR2L13", "MTRNR2L6")
keepgenes <- c(1:nrow(YDJ1filter_f2))[!(rownames(YDJ1filter_f2)%in%as.character(filter_genes))]
YDJ1filter_f3 <- subset(YDJ1filter_f2,features=keepgenes)
write.csv(YDJ1filter_f3@meta.data, file = "YDJ1filter-spots-metadata-clean.csv") # Export the clean spots metadata


# SCT normalization
#The clean expression matrix data were normalized using regularized negative binomial regression.
YDJ1filter_f3  <- SCTransform(YDJ1filter_f3 , assay = "Spatial", return.only.var.genes = FALSE)
DefaultAssay(YDJ1filter_f3) <- "SCT"
write.csv(YDJ1filter_f3@meta.data, file = "YDJ1filter-spots-metadata-clean-SCT.csv")


# Dimensionality reduction (PCA)
#PCA was performed and the 10 most significant components was determined by the DimHeatmap and ElbowPlot function. 
YDJ1filter_f3 <- RunPCA(YDJ1filter_f3)
DimHeatmap(YDJ1filter_f3, dims = 1:15, cells = 2000, balanced = TRUE)
p1 <- ElbowPlot(YDJ1filter_f3)
ggsave("YDJ1filter-PCA-ElbowPlot.pdf", plot=p1, width = 5, height =5)


#Spots Clustering (First round to identify Type I and Type II Myofiber)
YDJ1filter_f3 <- FindNeighbors(YDJ1filter_f3, dims = 1:10)
YDJ1filter_f3 <- FindClusters(YDJ1filter_f3, resolution = 0.14)  #Adjust resolution to change the number of clusters
YDJ1filter_f3 <- RunUMAP(YDJ1filter_f3, dims = 1:10)
#plot UMAP
p1 <- DimPlot(YDJ1filter_f3, reduction = "umap")
ggsave("YDJ1filter-first-cluster-r0.14-dims10-UMAP.pdf", plot=p2, width = 6, height =5)
#
#In order to evaluate the reliability of clustering, we plot the expression abundance heatmap of MYH7, MYH4 and MYH2 genes 
# (3 representative genes of muscle fiber type I ,type 2B and type 2A) on UMAP, and 
# preliminary characterization of the muscle fiber type of cluster were performed by study the expression levels of these 3 genes in each cluster.
g7 <- FeaturePlot(YDJ1filter_f3, features = c("MYH7"))
g4 <- FeaturePlot(YDJ1filter_f3, features = c("MYH4"))
g2 <- FeaturePlot(YDJ1filter_f3, features = c("MYH2"))
ggsave("YDJ1filter-all-spots-MYH7-heatmap.pdf", plot=g7, width = 6, height =5)
ggsave("YDJ1filter-all-spots-MYH4-heatmap.pdf", plot=g4, width = 6, height =5)
ggsave("YDJ1filter-all-spots-MYH-heatmap.pdf", plot=g2, width = 6, height =5)
#
#Visualize the spatial distribution of each cluster to confirm the rationality of the clusters.
p2 <- SpatialDimPlot(YDJ1filter_f3)
ggsave("YDJ1filter-first-cluster-r0.14-dims10-SpatialPlot.pdf", plot=p2, width = 6, height =5)


# Identified differentially expressed genes between TypeI and Type II clusters.
difgene1 <- FindAllMarkers(object=YDJ1filter_f3, assay ="SCT", slot= "counts", min.pct=0.1, logfc.threshold=0.1)
write.csv(difgene1, file="YDJ1filter-first-cluster-diff-genes-LnFC0.1.csv")


# Second round of Clustering (To identify Tpye IIB and IIA)
# In order to better distinguish between the two types of IIB and IIA, we filter out the clusters of type I clusters, and then cluster the remaining spots again.
#Filter out type I spots
keep0 <-c(1:nrow(YDJ1filter_f3@meta.data))[YDJ1filter_f3@meta.data$SCT_snn_res.0.14==0]
keep1 <-c(1:nrow(YDJ1filter_f3@meta.data))[YDJ1filter_f3@meta.data$SCT_snn_res.0.14==1]
YDJ1filter_f3_filter1 <-subset(YDJ1filter_f3, cells= c(keep0))
# Second round of Clustering 
YDJ1filter_f3_filter1 <- RunPCA(YDJ1filter_f3_filter1, verbose = FALSE)
YDJ1filter_f3_filter1 <- RunUMAP(YDJ1filter_f3_filter1, dims = 1:10)
YDJ1filter_f3_filter1 <- FindNeighbors(YDJ1filter_f3_filter1, dims = 1:10)
YDJ1filter_f3_filter1 <- FindClusters(YDJ1filter_f3_filter1, verbose = FALSE, resolution = 0.13)
#MYH genes heatmap
g7 <- FeaturePlot(YDJ1filter_f3_filter1, features = c("MYH7"))
g4 <- FeaturePlot(YDJ1filter_f3_filter1, features = c("MYH4"))
g2 <- FeaturePlot(YDJ1filter_f3_filter1, features = c("MYH2"))
ggsave("YDJ1filter-second-cluster-MYH7-heatmap.pdf", plot=g7, width = 6, height =5)
ggsave("YDJ1filter-second-cluster-MYH4-heatmap.pdf", plot=g4, width = 6, height =5)
ggsave("YDJ1filter-second-cluster-MYH2-heatmap.pdf", plot=g2, width = 6, height =5)
#Umap and Spatial Plot
p1 <- DimPlot(YDJ1filter_f3_filter1, reduction = "umap")
p2 <- SpatialDimPlot(YDJ1filter_f3_filter1)
ggsave("YDJ1filter-second-cluster-r0.13-dims10-UMAP.pdf", plot=p1, width = 6, height =5)
ggsave("YDJ1filter-second-cluster-r0.13-dims10-SpatialPlot.pdf", plot=p2, width = 6, height =5)
#Differential genes
difgene2 <- FindAllMarkers(object= YDJ1filter_f3_filter1, assay ="SCT", slot= "counts", min.pct=0.1, logfc.threshold=0.1)
write.csv(difgene2, file=" YDJ1filter-second-cluster-diff-genes-LnFC0.1.csv")


# Merge these 2 rounds of cluster results
#Rename the cluster and export spots ID 
type_1 <- rownames(YDJ1filter_f3@meta.data)[YDJ1filter_f3@meta.data$SCT_snn_res.0.14==1]
type_2B <- rownames(YDJ1filter_f3_filter1@meta.data)[YDJ1filter_f3_filter1@meta.data$SCT_snn_res.0.13==0]
type_2A <- rownames(YDJ1filter_f3_filter1@meta.data)[YDJ1filter_f3_filter1@meta.data$SCT_snn_res.0.13==1]
#Replace the cluster name to the samInfo column in meta.data
YDJ1filter_f3@meta.data$samInfo[rownames(YDJ1filter_f3@meta.data)%in%as.character(type_1)] <- "I"
YDJ1filter_f3@meta.data$samInfo[rownames(YDJ1filter_f3@meta.data)%in%as.character(type_2A)] <- "IIA"
YDJ1filter_f3@meta.data$samInfo[rownames(YDJ1filter_f3@meta.data)%in%as.character(type_2B)] <- "IIB"
#Plot the final UMAP and SpatialPlot
p1 <- DimPlot(YDJ1filter_f3, reduction = "umap", group.by = "samInfo")
p2 <- SpatialDimPlot(YDJ1filter_f3, group.by = "samInfo")
ggsave("YDJ1filter-final-cluster-UMAP.pdf", plot=p1, width = 6, height =5)
ggsave("YDJ1filter-final-cluster-SpatialPlot.pdf", plot=p2, width = 6, height =5)
# Differential expressed genes of the three clusters after the final merge:
Idents(YDJ1filter_f3) <- YDJ1filter_f3@meta.data$samInfo  #Modify inents according to the final cluster group
diffgene1 <- FindAllMarkers(object=YDJ1filter_f3, assay ="SCT", slot= "counts", min.pct=0.10, logfc.threshold=0.05)
write.csv(diffgene1, file="YDJ1filter-final-cluster-diff-genes.csv")


#Study of Representative spots in each cluster
#Export UMAP coordinates of all spots and grouping information of 3 clusters 
write.table(YDJ1filter_f3@reductions$umap@cell.embeddings, file='3cluster-UMAP-distance.txt', sep='\t', quote=F, row.names=T, col.names=T)
write.table(YDJ1filter_f3@meta.data$samInfo, file='3cluster-group.txt', sep='\t ',quote=F, row.names=T, col.names=T)
#we use a python (v2.7) script programmed by ourselves to obtain representative spots. 
#See the attachment for the script file: get_representative_spots.py .
#
#load the spots list back into R and build a Seurat sample object only content the representative spots.
# Import representative spot information and filter out other spots from the sample
data3<-read.table(file='clusterI-rep-spots187.txt',sep='\t',header=F)
data4<-read.table(file='clusterIIB-rep-spots187.txt',sep='\t',header=F)
data5<-read.table(file='clusterIIA-rep-spots187.txt',sep='\t',header=F)
rep_spots_1 <- data3[,1]
rep_spots_2B <- data4[,1]
rep_spots_2A <- data5[,1]
keep1 <- c(1:nrow(YDJ1filter_f3@meta.data)) [rownames(YDJ1filter_f3@meta.data)%in%as.character(rep_spots_1)]
keep2 <- c(1:nrow(YDJ1filter_f3@meta.data)) [rownames(YDJ1filter_f3@meta.data)%in%as.character(rep_spots_2B)]
keep3 <- c(1:nrow(YDJ1filter_f3@meta.data)) [rownames(YDJ1filter_f3@meta.data)%in%as.character(rep_spots_2A)]
rep_spots <- c(keep1, keep2, keep3)
YDJ1filter_3cluster_rep_spots <- subset(YDJ1filter_f3, cells= rep_spots)
#
#Plot UMAP and SpatialPlot again with group names of representative spots added in cluster names.
YDJ1filter_f3@meta.data$samInfo [rownames(YDJ1filter_f3@meta.data)%in%as.character(rep_spots_1)] <- "I_rep"
YDJ1filter_f3@meta.data$samInfo [rownames(YDJ1filter_f3@meta.data)%in%as.character(rep_spots_2B)] <- "IIB_rep"
YDJ1filter_f3@meta.data$samInfo [rownames(YDJ1filter_f3@meta.data)%in%as.character(rep_spots_2A)] <- "IIA_rep" # Add rep spots group id in samInfo
p1 <- DimPlot(YDJ1filter_f3, reduction = "umap", group.by ="samInfo")
p2 <- SpatialDimPlot(YDJ1filter_f3, group.by ="samInfo")
ggsave("YDJ1filter-final-clusters-UMAP-rep-spots.pdf", plot=p1, width = 6, height =5)
ggsave("YDJ1filter-final-clusters-SpatialPlot-rep-spots.pdf", plot=p2, width = 6, height =5)
#
#Use representative spots to screen the Marker genes.
Idents(YDJ1filter_3cluster_rep_spots) <- YDJ1filter_3cluster_rep_spots@meta.data$samInfo  #Change the cluster inents
markergene1 <- FindAllMarkers(object=YDJ1filter_3cluster_rep_spots, assay ="SCT", slot= "counts", min.pct=0.10, logfc.threshold=0.05)
write.csv(markergene1, file="YDJ1filter-rep-spots-all-marker-genes.csv") #All markers (one cluster vs. all others)
#Markers by Pairwise comparison
difgene1 <- FindMarkers(YDJ1filter_3cluster_rep_spots, ident.1 = "I_rep", ident.2 = "IIB_rep", assay ="SCT", slot= "counts", min.pct=0.10, logfc.threshold=0.01)
write.csv(difgene1, file="YDJ1filter-rep-spots-marker-genes-1-vs-2B.csv")
difgene2 <- FindMarkers(YDJ1filter_3cluster_rep_spots, ident.1 = "I_rep", ident.2 = "IIA_rep", assay ="SCT", slot= "counts", min.pct=0.10, logfc.threshold=0.01)
write.csv(difgene2, file="YDJ1filter-rep-spots-marker-genes-1-vs-2A.csv")
difgene3 <- FindMarkers(YDJ1filter_3cluster_rep_spots, ident.1 = "IIB_rep", ident.2 = "IIA_rep", assay ="SCT", slot= "counts", min.pct=0.10, logfc.threshold=0.01)
write.csv(difgene3, file="YDJ1filter-rep-spots-marker-genes-2B-vs-2A.csv")



