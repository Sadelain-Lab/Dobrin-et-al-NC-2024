#Set the current directory as the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

#Load Required Libraries
library(Matrix)
library(Seurat)
library(patchwork)
library(data.table)
library(reshape2)
library(tidyverse)
#library(monocle3)
library(harmony)
library(ComplexHeatmap)
library(RColorBrewer)
library(stringr)
library(venn)
library(pheatmap)
library(SeuratDisk)
library(DoubletFinder)
library(EnhancedVolcano)
library(openxlsx)
library(Nebulosa)
set.seed(12345)


#### read in data and add group, hashtags
aggr_multi_raw = Read10X(data.dir = "data and code/data/filtered_feature_bc_matrix")
aggr_multi_raw_seurat <- CreateSeuratObject(counts = aggr_multi_raw$`Gene Expression`, min.cells = 3) #19374 cells

## add hashtag data
aggr_multi_raw_seurat[["hashtag"]] <- CreateAssayObject(counts = aggr_multi_raw$`Multiplexing Capture`)
aggr_multi_raw_seurat@meta.data$Cell.BC = substr(rownames(aggr_multi_raw_seurat@meta.data), 1, 16)
aggr_multi_raw_seurat@meta.data$group = substr(rownames(aggr_multi_raw_seurat@meta.data), 18, 18)

aggr_multi_raw_seurat@meta.data$group = plyr::mapvalues(aggr_multi_raw_seurat@meta.data$group, from = c("1", "2"), to = c("low", "high"))


##QC: violin plot to check expression distribution
aggr_multi_raw_seurat[["percent.mt"]] <- PercentageFeatureSet(aggr_multi_raw_seurat, pattern = "^MT-")
aggr_multi_raw_seurat[["percent.rb"]] <- PercentageFeatureSet(aggr_multi_raw_seurat, pattern = "^RP[SL]")
VlnPlot(aggr_multi_raw_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 3, pt.size = 0)

#Cutoffs informed by prior analyses
aggr_multi_raw_seurat <- subset(aggr_multi_raw_seurat, subset = nFeature_RNA > 1000 & 
                                      percent.rb < 35 & percent.mt < 12) 

#Demultiplex cells
assay_hashtags <- NormalizeData(aggr_multi_raw_seurat, assay = "hashtag", normalization.method = "CLR")
hashtag_assigned <- HTODemux(assay_hashtags, assay = "hashtag", positive.quantile = 0.99)
hashtag_assigned$Cell.BC <- substr(hashtag_assigned$Cell.BC, 1, 16)

aggr_multi_raw_seurat@meta.data$hashtag_assigned = "NA"
aggr_multi_raw_seurat@meta.data$hashtag_assigned = hashtag_assigned$hashtag_maxID[match(aggr_multi_raw_seurat@meta.data$Cell.BC, hashtag_assigned$Cell.BC)]

#remove certain car-related features
aggr_multi_raw_seurat = aggr_multi_raw_seurat[-grep("80BB",rownames(aggr_multi_raw_seurat)), ]
aggr_multi_raw_seurat = aggr_multi_raw_seurat[-grep("80BB-junc",rownames(aggr_multi_raw_seurat)), ]
aggr_multi_raw_seurat = aggr_multi_raw_seurat[-grep("80BB-SFG",rownames(aggr_multi_raw_seurat)), ]
aggr_multi_raw_seurat = aggr_multi_raw_seurat[-grep("80BB-P2A",rownames(aggr_multi_raw_seurat)), ]
aggr_multi_raw_seurat = aggr_multi_raw_seurat[-grep("HIT-P2A",rownames(aggr_multi_raw_seurat)), ]
aggr_multi_raw_seurat = aggr_multi_raw_seurat[-grep("HIT-junc1",rownames(aggr_multi_raw_seurat)), ]
aggr_multi_raw_seurat = aggr_multi_raw_seurat[-grep("HIT-junc2",rownames(aggr_multi_raw_seurat)), ]
aggr_multi_raw_seurat = aggr_multi_raw_seurat[-grep("HIT-no80BB",rownames(aggr_multi_raw_seurat)), ]

#### regress out cell cycle effect
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
aggr_multi_raw_seurat_sub_cc <- CellCycleScoring(aggr_multi_raw_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
aggr_multi_raw_seurat_sub_cc <- ScaleData(aggr_multi_raw_seurat_sub_cc, vars.to.regress = c("S.Score", "G2M.Score"))

######Run analyses
aggr_multi_raw_seurat_sub_cc<- NormalizeData(aggr_multi_raw_seurat_sub_cc)
aggr_multi_raw_seurat_sub_cc<- FindVariableFeatures(aggr_multi_raw_seurat_sub_cc, selection.method = "vst", nfeatures = 2500)
aggr_multi_raw_seurat_sub_cc<- ScaleData(object = aggr_multi_raw_seurat_sub_cc, vars.to.regress = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

aggr_multi_raw_seurat_sub_cc <- RunPCA(aggr_multi_raw_seurat_sub_cc, verbose = FALSE)
aggr_multi_raw_seurat_sub_cc <- RunUMAP(aggr_multi_raw_seurat_sub_cc, dims = 1:30)
aggr_multi_raw_seurat_sub_cc <- FindNeighbors(aggr_multi_raw_seurat_sub_cc, dims = 1:30)
aggr_multi_raw_seurat_sub_cc <- FindClusters(aggr_multi_raw_seurat_sub_cc, resolution = 1)

######Remove Doublets
nExp <- round(ncol(aggr_multi_raw_seurat_sub_cc) * 0.04)  # expect 4% doublets
aggr_multi_raw_seurat_sub_cc <- doubletFinder_v3(aggr_multi_raw_seurat_sub_cc, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:30)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(aggr_multi_raw_seurat_sub_cc@meta.data)[grepl("DF.classification", colnames(aggr_multi_raw_seurat_sub_cc@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(aggr_multi_raw_seurat_sub_cc, group.by = "orig.ident") + NoAxes(),
                   DimPlot(aggr_multi_raw_seurat_sub_cc, group.by = DF.name) + NoAxes())

VlnPlot(aggr_multi_raw_seurat_sub_cc, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

aggr_multi_raw_seurat_sub_cc = aggr_multi_raw_seurat_sub_cc[, aggr_multi_raw_seurat_sub_cc@meta.data[, DF.name] == "Singlet"]
dim(aggr_multi_raw_seurat_sub_cc)



#SavePoint
SaveH5Seurat(aggr_multi_raw_seurat_sub_cc, filename = "aggr_multi_raw_seurat_sub_cc_v5.h5seurat" , overwrite = TRUE)
#aggr_multi_raw_seurat_sub_cc <- LoadH5Seurat("aggr_multi_raw_seurat_sub_cc_v5.h5seurat")

######
#### subset HIT and HIT_SFG_CD80BB
Idents(aggr_multi_raw_seurat_sub_cc) = "hashtag_assigned"

aggr_multi_raw_seurat_sub_cc$hashtag_assigned <- plyr::mapvalues(aggr_multi_raw_seurat_sub_cc$hashtag_assigned, 
                                                                                     from=c("hashtag3", "hashtag5"), 
                                                                                     to=c("HIT","HIT + CD80BB"))

aggr_multi_raw_seurat_sub_cc = subset(x = aggr_multi_raw_seurat_sub_cc, subset = hashtag_assigned %in% c("HIT", "HIT + CD80BB")) #412+1864 = 2276


#### Remove mitochondrial genes 
aggr_multi_raw_seurat_sub_cc =aggr_multi_raw_seurat_sub_cc[-grep("MT-",rownames(aggr_multi_raw_seurat_sub_cc)), ]
aggr_multi_raw_seurat_sub_cc = aggr_multi_raw_seurat_sub_cc[-grep("HIST",rownames(aggr_multi_raw_seurat_sub_cc)), ]
aggr_multi_raw_seurat_sub_cc = aggr_multi_raw_seurat_sub_cc[-grep("80BB",rownames(aggr_multi_raw_seurat_sub_cc)), ]

car_seurat = aggr_multi_raw_seurat_sub_cc

car_seurat  <- RunPCA(car_seurat, verbose = FALSE)
ElbowPlot(car_seurat,ndims = 50)

car_seurat  <- RunUMAP(car_seurat , dims = 1:18, min.dist=0.01)
car_seurat  <- FindNeighbors(car_seurat , dims = 1:15)
car_seurat  <- FindClusters(car_seurat , resolution = 1)

#Quality Cntrl
features = c("MS4A1", "CD19", "CD79B", "GNLY", "CD3E", "CD14", "FUT4",  "CD3D", "CD3G")
FeaturePlot(car_seurat, features = features, 
            min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3)

DimPlot(car_seurat, 
        reduction = "umap",group.by = "seurat_clusters",label = TRUE,
        pt.size = 3,
        label.size = 8,
        label.color = "black",
        label.box = TRUE)


DimPlot(car_seurat, 
        reduction = "umap",group.by = "hashtag_assigned",
        pt.size = 3)

###### remove clusters associated with high CD79B, MS4A1
car_seurat <- subset(car_seurat, idents = c(1,2,4,8, 10, 12), invert = TRUE)



car_seurat  <- RunPCA(car_seurat, verbose = FALSE)
ElbowPlot(car_seurat,ndims = 50)
car_seurat  <- RunUMAP(car_seurat , dims = 1:14)
car_seurat  <- FindNeighbors(car_seurat , dims = 1:14)
car_seurat  <- FindClusters(car_seurat , resolution = 0.5)



DimPlot(car_seurat, 
        reduction = "umap",group.by = "hashtag_assigned",label = TRUE,
        pt.size = 3,
        label.size = 8,
        label.color = "black",
        label.box = TRUE)

#SavePoint
SaveH5Seurat(car_seurat, filename = "aggr_multi_raw_seurat_sub_cc_noBclusters_HITSFG_noMITO_noHISTOv5.h5seurat" , overwrite = TRUE)
#car_seurat <- LoadH5Seurat("car_seurat.h5seurat")

#######################
###                 ### 
###   Quickstart    ###
###                 ###
#######################

#Run the below only if you're doing a quickstart
car_seurat <- LoadH5Seurat(file = "aggr_multi_raw_seurat_sub_cc_noBclusters_HITSFG_noMITO_noHISTOvAD.h5seurat")
car_seurat_HIT <- subset(car_seurat, subset = hashtag_assigned == 'HIT')
car_seurat_HITCD80BB <- subset(car_seurat, subset = hashtag_assigned == 'HIT + CD80BB')  


#name the current folder for saving files
#create subfolders for saving
save_folder = 'v3'
#The first time you run a new version, you'll have to initialize the folders
dir.create(paste0(save_folder))
dir.create(paste0(save_folder,'/vis',sep=''))
dir.create(paste0(save_folder,'/GSEA',sep=''))
dir.create(paste0(save_folder,'/DGEs',sep=''))

####Basic Stats and Information
car_seurat
#Cells per cluster and hashtag
table(car_seurat$seurat_clusters)
table(car_seurat$hashtag_assigned)

###UMAP of Clusters
colour_list <- ggplotColours(n=8)
colour_list = c("#FF61CC", "#CD9600", "#C77CFF", "#00BE67", "#00BFC4", "#00A9FF")

a <- DimPlot(car_seurat, 
        reduction = "umap",group.by = "seurat_clusters",label = TRUE,
        pt.size = 3,
        label.size = 8,
        label.color = "black",
        label.box = TRUE,
        cols = colour_list)
        
        
a <- a + labs(title = "") +
               theme(legend.position = "none", 
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               axis.text.y=element_blank(),
               axis.ticks.y=element_blank(),
  axis.title=element_text(size=18,face="bold"))
a


###UMAP of Hashtags

b <- FeaturePlot(car_seurat, features = c('MKI67'), 
                 min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3)
b <- b +
  theme(legend.position = c(0.7, 0.9), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(colour="black", size=18, 
                                   face="bold"))

b

c <- FeaturePlot(car_seurat, features = c('CD8A'), 
                 min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3)
c <- c +
  theme(legend.position = c(0.7, 0.9), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(colour="black", size=18, 
                                          face="bold"))

c

d <- FeaturePlot(car_seurat, features = c('CD4'), 
                 min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3)
d <- d +
  theme(legend.position = c(0.7, 0.9), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(colour="black", size=18, 
                                   face="bold"))

d

f <- FeaturePlot(car_seurat, features = c('IL7R'), 
                 min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3)
f <- f +
  theme(legend.position = c(0.7, 0.9), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(colour="black", size=18, 
                                   face="bold"))

f


e <- DimPlot(car_seurat, 
             reduction = "umap",group.by = "hashtag_assigned",
             pt.size = 3)

e <- e + labs(title = "") +
  theme(legend.position = c(0.7, 0.9), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"))
e <- e + theme(legend.text = element_text(colour="black", size=18, 
                                          face="bold"))
e

pdf(paste0(save_folder,"/vis/FourPlots_2022_10_14.pdf",sep=''), width = 10, height = 10)
a + b + c + d + plot_layout(ncol = 2)
a + a + f + plot_layout(ncol = 2)
dev.off()



#Remove CD4 Clusters
car_seurat_for8  <- FindClusters(car_seurat , resolution = 1.5)
car_seurat_for8 <- subset(car_seurat_for8, idents = c(1,2,3,10,9,8,7), invert = TRUE)

b_for8 <- DimPlot(car_seurat_for8, 
                  reduction = "umap",group.by = "hashtag_assigned",
                  pt.size = 3)

b_for8 <- b_for8 +
  theme(legend.position = c(0.7, 0.9), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(colour="black", size=18, 
                                   face="bold"))

b_for8
a_for8 <- DimPlot(car_seurat_for8, 
                  reduction = "umap",group.by = "seurat_clusters",label = TRUE,
                  pt.size = 3,
                  label.size = 8,
                  label.color = "black",
                  label.box = TRUE)
a_for8 <- a_for8 + labs(title = "") +
  theme(legend.position = "none", 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"))
a_for8
c_for8 <- FeaturePlot(car_seurat_for8, features = c('CD8A'), 
                      min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3)
c_for8 <- c_for8 +
  theme(legend.position = c(0.7, 0.9), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(colour="black", size=18, 
                                   face="bold"))

c_for8

d_for8 <- FeaturePlot(car_seurat_for8, features = c('CD4'), 
                      min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3)
d_for8 <- d_for8 +
  theme(legend.position = c(0.7, 0.9), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_text(colour="black", size=18, 
                                   face="bold"))

d_for8
pdf(paste0(save_folder,"/vis/FourPlots_for8.pdf",sep=''), width = 10, height = 10)
a_for8 + b_for8 + c_for8 + d_for8 + plot_layout(ncol = 2)
dev.off()

a_temp = VlnPlot(car_seurat_for8,features = "IL7R",pt.size = 0) & theme(plot.title = element_text(size=10))
b_temp = VlnPlot(car_seurat_for8,features = "IL7R",pt.size = 0, group.by = "hashtag_assigned") & theme(plot.title = element_text(size=10))
  
pdf(paste0(save_folder,"/vis/Vln_for8.pdf",sep=''), width = 10, height = 10)
a_temp + b_temp + plot_layout(ncol = 2)
dev.off()


###Volcano plot for cluster 4
Idents(car_seurat) = car_seurat$seurat_clusters

i = 4
tmp = FindMarkers(car_seurat, ident.1 = i) 
tmp = tmp %>% filter(p_val_adj < 0.05)
#write.csv(tmp, file = paste0("temp/",i,"vs",j,".csv"))
tmp = tmp[order(tmp$avg_log2FC,decreasing = TRUE),]

pdf(paste0(save_folder,"/vis/Cluster_",i,"_volcano.pdf"), width = 15, height = 10)
print(EnhancedVolcano(tmp,
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      lab = rownames(tmp),
                      selectLab = rownames(tmp)[c(1:20,(length(rownames(tmp))-20):length(rownames(tmp)))],
                      xlim = c(-2, 2),
                      #  ylim = c(0,100),
                      xlab = bquote(~Log[2]~ 'fold change'),
                      ylab = bquote(~-Log[10]~ 'p-value'),
                      gridlines.major = TRUE,
                      gridlines.minor = FALSE,
                      pCutoff = 5e-2,
                      FCcutoff = log(2),
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 0.8,
                      hline = c(10e-1,10e-2),
                      hlineCol = c('grey0','grey25'),
                      colAlpha = 0.5,
                      title = paste0("cluster", i),
                      subtitle = 'Differential Gene Expression',
                      legendLabels=c('Not sig.','Log (base 2) FC','q-value',
                                     'p-value & Log (base e) FC'),
                      legendPosition = 'bottom',
                      legendLabSize = 14,
                      legendIconSize =  6.0,
                      titleLabSize = 20,
                      captionLabSize = 14,
                      axisLabSize = 14,
                      pointSize = 4,
                      labSize = 4,
                      drawConnectors = TRUE,
                      widthConnectors = 0.4,
                      colConnectors = 'grey30'
))
dev.off()

#Feature Map of Canonical Markers
features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
             "CD8A", "NKG7", "CD16", "MS4A1", "MPO", "CD68", "CD14", "FUT4",  "CD3D", "CD3G", "CD247", "CD4", "CD8A")

png(paste0(save_folder,"/vis/Features_CellMarker.png",sep=''), width = 36, height = 36, units = 'cm', res = 500)
FeaturePlot(car_seurat, features = features,
              min.cutoff = "q05", max.cutoff = "q95", order = TRUE)
dev.off()

a <- VlnPlot(car_seurat,features = "percent.mt") & theme(plot.title = element_text(size=10))
b <- VlnPlot(car_seurat,features = "percent.rb") & theme(plot.title = element_text(size=10))
c <- VlnPlot(car_seurat,features = "nCount_RNA") & theme(plot.title = element_text(size=10))
d <- VlnPlot(car_seurat,features = "nFeature_RNA") & theme(plot.title = element_text(size=10))

pdf(paste0(save_folder,"/vis/FourQCPlots.pdf",sep=''), width = 15, height = 15)
a + b + c + d + plot_layout(ncol = 2)
dev.off()


###### DEGs by cluster
Idents(car_seurat) = car_seurat$seurat_clusters
aggr_all_markers_byhashtag = FindAllMarkers(car_seurat)

#Identify top markers for each cluster in top markers, then plot
top_markers <- aggr_all_markers_byhashtag %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
features = top_markers %>% group_by(cluster) %>% top_n(10, avg_log2FC) 

pdf(paste0(save_folder,"/vis/DEGsByCluster_FeatureMap_top10.pdf",sep=''), width = 15, height = 15)
FeaturePlot(car_seurat, features = features$gene,
            min.cutoff = "q05", max.cutoff = "q95", order = TRUE)
dev.off()

pdf(paste0(save_folder,"/vis/DEGsByCluster_Dot_top10_10_12.pdf",sep=''), width = 13.5, height = 4.5)
DotPlot(car_seurat, features = unique(features$gene), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))+
  ggtitle("DEGs By Cluster")
dev.off()

write.csv(top_markers, file = paste0(save_folder, "/DGEs/DGEs_by_cluster.csv"))

## Features from paper
feature_from_louie_etal = c("GZMB","GZMH","PRF1","IFNG","GZMK","NKG7",
                            "GNLY","KLRD1","KLRG1","KLRC1","IL18R1","FCGR3A","FCRL6","HLA-E","HLA-C","CCL5","CXCR4","TGFB1","IL10RA",
                            "KLF2","KLF3","ZEB2","ID2","ZNF683","EGR1","FOS","JUND","ZFP36L2",
                            "LCK","CD27","CD37","IFITM2","STAT3", "JAK1",
                            "CX3CR1","TIGIT","HAVCR2","TOX","TBX21","EOMES","PTPN2","PTPN11","PTPN6",
                            "IL15RA","CD38","HLA-DRA","BIRC2","MKI67","HMGB2",
                            "RORA","MALAT1","ITGA4","IRF1","CD160",
                            "BCL2","NR4A2","LTB","SELL","CD69","TCF7","LEF1","IL7R") 

features_func_sig = c("CD4", "CD8A","GZMA", "GZMH","GZMM","GZMK", "LTB", "CCL5","CCL3", "MKI67", "FASLG", "GNLY", "IL7R","CD40LG")
features_exh_sig = c("ID3", "CD200","CTLA4","BTLA","TIGIT","LAG3","NR4A1","NR4A2","FASLG")
f_selected = c("ID3","NR4A2","CD200","CTLA4","TIGIT","LAG3","BTLA","FASLG","GZMA", "GZMH","GZMK", "LTB")

aa <- FeaturePlot(car_seurat, features = c("ID3","NR4A2","CD200","CTLA4"), 
                 min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3) & NoAxes() + NoLegend()


bb <- FeaturePlot(car_seurat, features = c("TIGIT","LAG3","BTLA","FASLG"),
                  min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3) & NoAxes() + NoLegend()


cc <- FeaturePlot(car_seurat, features = c("GZMA", "GZMH","GZMK", "LTB"),
                  min.cutoff = "q10", max.cutoff = "q95", order = TRUE, pt.size = 3)  & NoAxes() + NoLegend()


pdf(paste0(save_folder,"/vis/features_pub.pdf",sep=''), width = 10, height = 10)
aa + plot_layout(ncol = 2)
bb + plot_layout(ncol = 2)
cc + plot_layout(ncol = 2)
dev.off()




pdf(paste0(save_folder,"/vis/features_pub.pdf",sep=''), width = 16.16, height = 15)
FeaturePlot(car_seurat, features = f_selected, pt.size = 3,
            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) & NoAxes() + NoLegend()

dev.off()



pdf(paste0(save_folder,"/vis/DEGsByCluster_func_sig.pdf",sep=''), width = 15, height = 15)
VlnPlot(car_seurat,features = features_func_sig,pt.size = 0,) & theme(plot.title = element_text(size=10))
dev.off()
pdf(paste0(save_folder,"/vis/DEGsByCluster_exh_sig.pdf",sep=''), width = 15, height = 15)
VlnPlot(car_seurat,features = features_exh_sig,pt.size = 0,) & theme(plot.title = element_text(size=10))
dev.off()


pdf(paste0(save_folder,"/vis/DEGsByCluster_Dot_fromLouieEtAl.pdf",sep=''), width = 15, height = 4)
DotPlot(car_seurat, features = feature_from_louie_etal)+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))+
  ggtitle("DEGs By Cluster")
dev.off()

pdf(paste0(save_folder,"/vis/DEGsByVLN4_8.pdf",sep=''), width = 15, height = 4)
VlnPlot(car_seurat,features = c("CD8A","CD8B","CD4"),pt.size = 0,) & theme(plot.title = element_text(size=10))
dev.off()


###### Feature and Violin Plots


pdf(paste0(save_folder,"/vis/Exhaust_By_Cluster.pdf",sep=''), width = 15, height = 4)
DotPlot(car_seurat, features = unique(features), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))+
  ggtitle("T Cell Genes By Cluster")
dev.off()


#features_exh = c("ID3", "CD200","CTLA4","BTLA","TIGIT","LAG3","HAVCR2","NR4A1","NR4A2","NR4A3", "TOX", "TOX2", "SOX4","PDCD1")
#features_func = c("CD4", "CD8A","GZMA", "GZMH","GZMM","GZMB","GZMK", "LTB", "CCL5", "MKI67", "FASLG","FAS","PRF1","GNLY", "IL7R")

#pdf(paste0(save_folder,"/vis/Feature_Func.pdf",sep=''), width = 20, height = 15)
#FeaturePlot(car_seurat, features = features_func, pt.size = 3,
#            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) & NoAxes() & NoLegend()
#dev.off()

#pdf(paste0(save_folder,"/vis/Feature_Exh.pdf",sep=''), width = 20, height = 15)
#FeaturePlot(car_seurat, features = features_exh, pt.size = 3,
#            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) & NoAxes() & NoLegend()
#dev.off()




for(level in levels(aggr_all_markers_byhashtag$cluster)){
  aggr_all_markers_sub = aggr_all_markers_byhashtag %>% dplyr::filter(cluster == level)
  aggr_all_markers_sub$avg_FC = 2^aggr_all_markers_sub$avg_log2FC
  write.xlsx(aggr_all_markers_sub,paste0(save_folder, "/DGEs/cluster_",level,".xlsx"))
}

###Volcano plot for the two hashtags
Idents(car_seurat) = car_seurat$hashtag_assigned

aggr_all_markers = FindAllMarkers(car_seurat)
aggr_all_markers_sub =  aggr_all_markers %>% dplyr::filter(cluster == 'HIT + CD80BB')

pdf(paste0(save_folder,"/vis/volcano_HIT_resize.pdf"), width = 6, height = 7)
EnhancedVolcano(aggr_all_markers_sub,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                lab = rownames(aggr_all_markers_sub),
                selectLab = rownames(aggr_all_markers_sub)[c(1:100,(length(rownames(aggr_all_markers_sub))-40):length(rownames(aggr_all_markers_sub)))],
                xlim = c(min(aggr_all_markers_sub$avg_log2FC)-0.1, max(aggr_all_markers_sub$avg_log2FC)+0.1),
                #  ylim = c(0,100),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~ 'p-value'),
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                pCutoff = 5e-2,
                FCcutoff = log(1),
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                hline = c(10e-1,10e-2),
                hlineCol = c('grey0','grey25'),
                colAlpha = 0.5,
                #title = paste0("cluster ",i),
                subtitle = 'Differential Gene Expression for HIT',
                legendLabels=c('Not sig.','Log (base 2) FC','q-value',
                               'p-value & Log (base e) FC'),
                legendPosition = 'bottom',
                legendLabSize = 14,
                legendIconSize =  6.0,
                titleLabSize = 20,
                captionLabSize = 14,
                axisLabSize = 14,
                pointSize = 4,
                labSize = 4,
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                colConnectors = 'grey30')
dev.off()

###### DEGs by hashtag
Idents(car_seurat) = car_seurat$hashtag_assigned

aggr_all_markers_byhashtag = FindAllMarkers(car_seurat)

top_markers <- aggr_all_markers_byhashtag %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

pdf(paste0(save_folder,"/vis/DEGSbyHashtag_FeatureMap.pdf",sep=''), width = 15, height = 15)
FeaturePlot(car_seurat, features = top_markers$gene, pt.size = 3,
            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) & NoAxes() & NoLegend()
dev.off()

pdf(paste0(save_folder,"/vis/DEGsByHashtag_Dot.pdf",sep=''), width = 15, height = 3)
DotPlot(car_seurat, features = unique(top_markers$gene), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))+
  ggtitle("DEGs By Cluster")
dev.off()

for(level in levels(aggr_all_markers_byhashtag$cluster)){
  aggr_all_markers_sub = aggr_all_markers_byhashtag %>% dplyr::filter(cluster == level)
  aggr_all_markers_sub$avg_FC = 2^aggr_all_markers_sub$avg_log2FC
  write.xlsx(aggr_all_markers_sub,paste0(save_folder,"/DGEs/hashtag_",level,".xlsx"))
}

write.csv(top_markers, file = paste0(save_folder, "/DGEs/DGEs_by_cluster.csv"))

###### Feature and Violin Plots

pdf(paste0(save_folder,"/vis/Func_By_Hashtag.pdf",sep=''), width = 15, height = 3)
DotPlot(car_seurat, features = unique(features), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))+
  ggtitle("T Cell Genes By Cluster")
dev.off()


features_exh_sig = c("ID3", "CD200","CTLA4","BTLA","TIGIT","LAG3","NR4A1","NR4A2","FASLG")
features_exh_notSig = c("HAVCR2","NR4A3", "TOX", "TOX2", "SOX4","PDCD1")

features_func_sig = c("CD4", "CD8A","GZMA", "GZMH","GZMM","GZMK", "LTB", "CCL5","CCL3", "MKI67", "FASLG", "GNLY", "IL7R","CD40LG")
features_func_notSig= c("GZMB", "FAS","PRF1")


pdf(paste0(save_folder,"/vis/exh_sig plots.pdf",sep=''), width = 20, height = 15)
FeaturePlot(car_seurat, features = features_exh_sig, pt.size = 3,
            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) # & NoXAxes() & NoLegend()
VlnPlot(car_seurat,  group.by = "hashtag_assigned", features = features_exh_sig, pt.size = 0, ncol=6) & theme(axis.title.x = element_blank()) & theme(axis.text.x = element_blank())

dev.off()

pdf(paste0(save_folder,"/vis/exh_notSig plots.pdf",sep=''), width = 20, height = 15)
FeaturePlot(car_seurat, features = features_exh_notSig, pt.size = 3,
            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) # & NoXAxes() & NoLegend()
VlnPlot(car_seurat,  group.by = "hashtag_assigned", features = features_exh_notSig, pt.size = 0, ncol=6) & theme(axis.title.x = element_blank()) & theme(axis.text.x = element_blank())

dev.off()

pdf(paste0(save_folder,"/vis/func_sig plots.pdf",sep=''), width = 20, height = 15)
FeaturePlot(car_seurat, features = features_func_sig, pt.size = 3,
            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) # & NoXAxes() & NoLegend()
VlnPlot(car_seurat,  group.by = "hashtag_assigned", features = features_func_sig, pt.size = 0, ncol=6) & theme(axis.title.x = element_blank()) & theme(axis.text.x = element_blank())

dev.off()

pdf(paste0(save_folder,"/vis/func_notSig plots.pdf",sep=''), width = 20, height = 15)
FeaturePlot(car_seurat, features = features_func_notSig, pt.size = 3,
            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) # & NoXAxes() & NoLegend()
VlnPlot(car_seurat,  group.by = "hashtag_assigned", features = features_func_notSig, pt.size = 0, ncol=6) & theme(axis.title.x = element_blank()) & theme(axis.text.x = element_blank())

dev.off()

pdf(paste0(save_folder,"/vis/vln_prolif.pdf",sep=''), width = 20, height = 20)
VlnPlot(car_seurat,  group.by = "hashtag_assigned", features = c("JPT1", "STMN1", "PTTG1", "LMNB1", "MKI67", "HMGN2", "DDX39A", "SMC4", "TPX2", "NUSAP1", "LBR", "TMPO", "CKS1B", "CENPF", "KNL1"), pt.size = 0, ncol=6) & theme(axis.title.x = element_blank()) & theme(axis.text.x = element_blank())
FeaturePlot(car_seurat, features =  c("JPT1", "STMN1", "PTTG1", "LMNB1", "MKI67", "HMGN2", "DDX39A", "SMC4", "TPX2", "NUSAP1", "LBR", "TMPO", "CKS1B", "CENPF", "KNL1"), pt.size = 3,
            min.cutoff = "q05", max.cutoff = "q95", order = TRUE) # & NoXAxes() & NoLegend()
dev.off()



### GSEA using Escape
BIOCARTA_CTLA4_PATHWAY Reactome CD28
BIOCARTA_41BB_PATHWAY4-1 BB pathway
BiocManager::install("dittoSeq")
library(dittoSeq)
library(escape)

pathways.Biocarta <- getGeneSets(species = "Homo sapiens", library = "C2", subcategory = "CP:BIOCARTA" )
pathways.Hallmark <- getGeneSets(species = "Homo sapiens", library = "H")


ES.seurat <- enrichIt(obj = car_seurat, 
                      gene.sets = pathways.Hallmark, #set which pathway to run 
                      groups = 1000, cores = 2, ssGSEA.norm	= TRUE,
                      min.size = 5)
car_seurat <- Seurat::AddMetaData(car_seurat, ES.seurat)

Idents(car_seurat) = car_seurat$hashtag_assigned
ES2 <- data.frame(car_seurat[[]], Idents(car_seurat))

output <- getSignificance(ES2, 
                          group = "hashtag_assigned", 
                          #group = "seurat_clusters",
                          gene.sets = names(ES.seurat),
                          fit = "Wilcoxon") 
colnames(output)[4] <- 'median.HIT80BB'

output <- output %>% dplyr::mutate(NES_Ratio = median.HIT80BB/median.HIT) %>%
  rownames_to_column("pathway")


write.xlsx(output, file = paste0(save_folder,"/GSEA/Hallmark_By_Hashtag.xlsx",sep=''))
colors <- colorRampPalette(c("#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF"))

#pathway_vars <- c("BIOCARTA_41BB_PATHWAY", "BIOCARTA_IL2RB_PATHWAY", "BIOCARTA_CTLA4_PATHWAY" )
pathway_vars <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
pdf(paste0(save_folder,"/GSEA/GSEA_vln_Hallmark_byhashtag.pdf",sep=''), width = 16, height = 10)
multi_dittoPlot(car_seurat, vars = pathway_vars, 
                group.by = "hashtag_assigned",  plots = c("jitter", "vlnplot", "boxplot"), 
                ylab = "Enrichment Scores", 
                theme = theme_classic() + theme(text = element_text(size=14),
                  plot.title = element_text(size = 18)))
dev.off()

?getGeneSets

library("devtools")
#install_github("GSEA-MSigDB/GSEA_R", force = TRUE)
#remotes::install_github("ctlab/phantasus")
#install_github('immunogenomics/presto')
library(presto)
library(GSEA)
#library(phantasus)
library(msigdbr)
library(fgsea)
library(ExperimentHub)
library(GSEABase)
library(forcats)
library(dplyr)
#install.packages('forcats')


#Computes auROC and Wilcoxon p-value
#can either do by cluster or hashtag
#markers.for.gsea <- wilcoxauc(car_seurat, 'hashtag_assigned')
markers.for.gsea <- wilcoxauc(car_seurat, 'seurat_clusters')
#x <- msigdbr_collections()

pathways.Hallmark <- msigdbr(species = "Homo sapiens", category = "H") %>% split(x = .$gene_symbol, f = .$gs_name)
pathways.Immune <- msigdbr(species = "Homo sapiens", category = "C7", subcategory='IMMUNESIGDB') %>%
  dplyr::filter(grepl("TCELL",gs_name)) %>% #select only T cell immune pathways
  dplyr::filter(!grepl("NEUTROPHIL",gs_name)) %>% #select only T cell immune pathways
  dplyr::filter(!grepl("MACROPHAGE",gs_name)) %>% #select only T cell immune pathways
  dplyr::filter(!grepl("BCELL",gs_name)) %>% #select only T cell immune pathways
  split(x = .$gene_symbol, f = .$gs_name)
pathways.KEGG <- msigdbr(species = "Homo sapiens", category = "C2", subcategory='CP:KEGG') %>% split(x = .$gene_symbol, f = .$gs_name)
pathways.Transcript <- msigdbr(species = "Homo sapiens", category = "C3", subcategory='TFT:GTRD') %>% split(x = .$gene_symbol, f = .$gs_name)

pathway_sets <- c('Hallmark', 'Immune', 'KEGG','Transcript')

#If you don't want to do entire loop, you can just pick one variable and skip it
#i='HIT + CD80BB'
#i='HIT'
#i='9'
#j = pathway_sets[1]

for(i in unique(markers.for.gsea$group)){

  tmp = markers.for.gsea %>%
    dplyr::filter(group == i)  %>%
    arrange(desc(auc)) %>% 
    dplyr::select(feature, auc)
  
  ranks<- deframe(tmp)
  
  nam <- paste("ranks_",gsub(" ", "", i, fixed = TRUE), sep = "")
  assign(nam, ranks)
  
  for(j in pathway_sets) { #iterate through all pathways
  fgseaRes <- fgsea(pathways=eval(parse(text = paste('pathways.',j, sep=''))), 
                    stats=ranks, scoreType = "pos", 
                    minSize=15, maxSize=500, eps = 0, nproc=2)
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% 
    dplyr::filter(padj < 0.01)
  
  nam <- paste("GSEA_", j, "_", gsub(" ", "", i, fixed = TRUE), sep = "")
  assign(nam, fgseaResTidy)
  write.xlsx(fgseaResTidy, file = paste(save_folder, "/GSEA/",nam,".xlsx", sep=""))
  
  fgseaResTidy_head <- fgseaResTidy %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% 
    slice_head(n=20)
  
  fgseaResTidy_tail <- fgseaResTidy %>%
    as_tibble() %>%
    arrange(desc(NES)) %>% 
    slice_tail(n=20)
  
  fgseaResPlot = rbind(fgseaResTidy_head, fgseaResTidy_tail) 
  
  tmp_plot = ggplot(fgseaResPlot , aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill=NES<0)) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=paste0(j," gene sets NES from GSEA, Cluster ", i,sep = "")) + 
    theme_minimal(base_size = 15) + 
    theme(plot.background = element_rect(colour = "white"), legend.position="none")
  print(tmp_plot)
  ggsave(tmp_plot, file=paste0(save_folder, "/GSEA/", i ,"_", j, ".png"), width = 11,height = 12)
  
}}

#Find pathways of interest and plot them using either ranks_HIT or ranks_HIT+CD80BB
a <- plotEnrichment(pathways.immune[["GSE23568_ID3_TRANSDUCED_VS_ID3_KO_CD8_TCELL_UP"]],
               ranks_HIT) + labs(title="HIT+CD80BB, ID3 Transduced vs KO")
b <- plotEnrichment(pathways.immune[["GSE32533_WT_VS_MIR17_OVEREXPRESS_ACT_CD4_TCELL_UP"]],
                    ranks_HIT) + labs(title="HIT+CD80BB, WT_vs_MIR17")
c <- plotEnrichment(pathways.transcript[["NFRKB_TARGET_GENES"]],
                    ranks_HIT) + labs(title="HIT+CD80BB, NFRKB Targets")
d <- plotEnrichment(pathways.transcript[["ASH1L_TARGET_GENES"]],
                    ranks_HIT) + labs(title="HIT+CD80BB, ASH1L Targets")

png(paste0(save_folder,"/GSEA/DGE_Pathways_HIT+CD80BB.png",sep=''), width = 36, height = 36, units = 'cm', res = 500)
a + b + c + d + plot_layout(ncol = 2)
dev.off()

save(file = "workspace_9_24.RData", list = ls(all.names = TRUE))
?save
