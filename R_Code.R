setwd("/media/4T/hmwphd/ospbmconlyfinally")

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggsci)
library(pheatmap)
library(ggplot2)
rm(list=ls())
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

#读取OS患者pbmc
##设置文件目录与样本名称
dir <- dir("/media/4T/hmwphd/Data/OSpriPbmc3/")
dir <- paste0("/media/4T/hmwphd/Data/OSpriPbmc3/", dir)
#查看文件顺序
dir                         
#按文件顺序给样本命名，名称不要以数字开头，中间不能有空格 
samples_name = c("ospbmc1", "ospbmc2", "ospbmc3")
##使用循环命令批量创建seurat对象
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  #不设置min.cells过滤基因会导致CellCycleScoring报错
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = samples_name[i], min.cells = 6, min.features = 300)
  #给细胞barcode加个前缀，防止合并后barcode重名
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
}
#给列表命名
names(scRNAlist) <- samples_name
scRNAospbmc <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNAospbmc$proj <- rep("10x", ncol(scRNAospbmc))
head(scRNAospbmc@meta.data)
save(scRNAospbmc, file = "scRNAospbmc_origin.Rdata")

dir.create("QC")
scRNAospbmc
# An object of class Seurat
# 17062 features across 32720 samples within 1 assay
# Active assay: RNA (17062 features, 0 variable features)
table(scRNAospbmc@active.ident)
# ospbmc1 ospbmc2 ospbmc3 
#   17332    7441    7947

scRNA <- scRNAospbmc
#计算线粒体比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

##绘制质控小提琴图
#设置可能用到的主题
theme.set1 = theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
theme.set2 = theme(axis.title.x=element_blank())
#设置绘图元素
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
group = "orig.ident"

#质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0, 
                                                         features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow = 1)    
ggsave("QC/vlnplot_before_qc.pdf", plot = violin, width = 15, height = 5) 
ggsave("QC/vlnplot_before_qc.png", plot = violin, width = 15, height = 5)  

plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- plot1 + plot2 + plot_layout(guides = "collect")
ggsave(file="QC/nCount_RNAvsPercent.mt.pdf", p, width = 15, height = 5)
ggsave(file="QC/nCount_RNAvsPercent.mt.png", p, width = 15, height = 5)
ggsave(file="QC/nCount_RNAvsPercent.mt.eps", p, width = 15, height = 5)

##数据质控并绘制小提琴图
scRNA <- subset(scRNA, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 10)
scRNA
# An object of class Seurat 
# 17062 features across 26947 samples within 1 assay 
# Active assay: RNA (17062 features, 0 variable features)
table(scRNA@active.ident)
# ospbmc1 ospbmc2 ospbmc3 
#   13878    6665    6404
plots = list()
for(i in seq_along(plot.featrures)){plots[[i]] = VlnPlot(scRNA, group.by = group, pt.size = 0, features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow = 1)     
ggsave("QC/vlnplot_after_qc.pdf", plot = violin, width = 15, height = 5) 
ggsave("QC/vlnplot_after_qc.png", plot = violin, width = 15, height = 5)

##===数据标准化===##
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures() %>% ScaleData(features = rownames(scRNA))

##===高变基因===##
top10 <- head(VariableFeatures(scRNA), 10)
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
p = plot1 + plot2 + plot_layout(guides = "collect")
ggsave(file="QC/High_variable_gene.pdf", plot2, width = 7.5, height = 5)
ggsave(file="QC/High_variable_gene.png", plot2, width = 7.5, height = 5)
ggsave(file="QC/High_variable_gene.eps", plot2, width = 7.5, height = 5)

##保存质控后的结果
save(scRNA, file = "scRNA_QC.Rdata")


library(Seurat)
library(tidyverse)
library(SingleR)
library(harmony)
library(patchwork)
dir.create("Overview")
rm(list=ls())

##===数据概览===##
load("scRNA_QC.Rdata")

##降维聚类
scRNA <- RunPCA(scRNA, verbose = F)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30
scRNA <- scRNA %>% RunTSNE(dims = pc.num) %>% RunUMAP(dims = pc.num) %>% FindNeighbors(dims = pc.num) %>% FindClusters(resolution = 0.3)

##查看批次效应
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident")
ggsave("Overview/batch_overview.png", p1, width = 7.5, height = 5)
ggsave("Overview/batch_overview.pdf", p1, width = 7.5, height = 5)
ggsave("Overview/batch_overview.eps", p1, width = 7.5, height = 5)

#使用分面图查看效果
p <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident", ncol = 3) #+ NoLegend()
ggsave("Overview/batch_facet.pdf", p, width = 15, height = 5)
ggsave("Overview/batch_facet.png", p, width = 15, height = 5)
ggsave("Overview/batch_facet.eps", p, width = 15, height = 5)

#绘制umap图
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Phase")
ggsave("Overview/CellCycle_dimplot.png", p1, width = 7.5, height = 5)
ggsave("Overview/CellCycle_dimplot.pdf", p1, width = 7.5, height = 5)
ggsave("Overview/CellCycle_dimplot.eps", p1, width = 7.5, height = 5)


##保存结果
save(scRNA, file = "Overview/scRNA_overview.Rdata")

##===标准化与批次校正===##
dir.create("Harmony")
rm(list=ls())
load("scRNA_QC.Rdata")

set.seed(1)

##===数据标准化===##
scRNA <- NormalizeData(scRNA) %>% FindVariableFeatures() %>% ScaleData(features = rownames(scRNA))

top10 <- head(VariableFeatures(scRNA), 10)
plot1 <- VariableFeaturePlot(scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plotc <- plot2 + plot_layout(guides = "collect")
ggsave("Harmony/VariableFeatures.pdf", plotc, width = 7.5, height = 5)
ggsave("Harmony/VariableFeatures.png", plotc, width = 7.5, height = 5)
ggsave("Harmony/VariableFeatures.eps", plotc, width = 7.5, height = 5)

##使用harmony整合数据
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", assay.use="RNA", max.iter.harmony = 10)
ElbowPlot(scRNA, ndims = 50)
pc.num=1:30

##===降维聚类与可视化===##
##非线性降维与聚类
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.3)
table(scRNA$seurat_clusters)
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 4082 3151 2869 2835 2744 2081 1632 1585 1361 1231 1001  800  535  500  258  143  139

##查看harmony的整合效果
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters")
p2 <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident")
pc = p1 + p2 
ggsave("Harmony/Harmony_integr.pdf", pc, width = 15, height = 5)
ggsave("Harmony/Harmony_integr.png", pc, width = 15, height = 5)
ggsave("Harmony/Harmony_integr.eps", pc, width = 15, height = 5)

#使用分面图查看效果
p <- DimPlot(scRNA, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident", ncol = 3) + NoLegend()
ggsave("Harmony/Harmony_facet.pdf", p, width = 15, height = 5)
ggsave("Harmony/Harmony_facet.png", p, width = 15, height = 5)

##保存结果
save(scRNA, file = "scRNA_harmony.Rdata")

library(Seurat)
library(tidyverse)
library(patchwork)
dir.create("CellType")
rm(list=ls())
source("/media/4T/Resource/function.R")
scRNA <- get(load("scRNA_harmony.Rdata"))
##===Marker基因鉴定===##
##提取各个Cluster的marker genes
ClusterMarker <- FindAllMarkers(scRNA, assay = "RNA", slot = "data", only.pos = T, logfc.threshold = 0.25, min.pct = 0.1)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'CellType/ClusterMarker.csv', row.names=F)
#ClusterMarker <- read.csv('CellType/ClusterMarker.csv')
#提取差异显著的marker genes
top = 100   #可根据需要调整
TopMarkers1 <- ClusterMarker %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers,'CellType/TopMarkers.csv', row.names=F)

##提取没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", ClusterMarker$gene, ignore.case = F),]
top = 100   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers_noRibo <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers_noRibo,'CellType/TopMarkers_noRibo.csv', row.names=F)

##导出一个方便对比的表格(每行一个cluster的marker gene)
TopMarkers2Lines(data = TopMarkers_noRibo, output = "CellType")

##===人工鉴定===##
##导入人工鉴定结果
cell.type <- c("NK_cells", # 0
               "CD8_T_cells", # 1
               "CD4_T_cells",  # 2
               "CD4_T_cells",  # 3
               "Monocytes",   # 4
               "Plasma_cells", # 5
               "Monocytes",    # 6
               "CD8_T_cells", # 7
               "B_cells", # 8
               "NK_cells", # 9
               "Cycling_plasma_cells", # 10
               "MAITs", # 11
               "Monocytes",  # 12
               "Monocytes",    # 13
               "Unkonwn_cells",  # 14
               "pDCs",    # 15
               "Platelets" # 16
)
Idents(scRNA) <- "seurat_clusters"
names(cell.type) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, cell.type)
scRNA$celltype <- Idents(scRNA)
Idents(scRNA) <- "seurat_clusters"

##结果用UMAP图展示
p1 <- DimPlot(scRNA, reduction = "umap", label = T) #+ NoLegend()
p2 <- DimPlot(scRNA, reduction = "umap", label = T, group.by = "celltype") #+ NoLegend()
pc = p1|p2
ggsave("CellType/CellType_Custom.pdf", pc, width = 15, height = 5)
ggsave("CellType/CellType_Custom.png", pc, width = 15, height = 5)
ggsave("CellType/CellType_Custom.eps", pc, width = 15, height = 5)

##保存结果
save(scRNA, file = "scRNA_classify.Rdata")

table(scRNA$celltype)
#  NK_cells          CD8_T_cells          CD4_T_cells            Monocytes         Plasma_cells              B_cells Cycling_plasma_cells 
#      5313                 4736                 5704                 5411                 2081                 1361                 1001 
#     MAITs        Unkonwn_cells                 pDCs            Platelets 
#       800                  258                  143                  139


##===结果可视化===##
##Stackbar细胞丰度柱状图
tmp <- select(scRNA@meta.data, c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}
source("/media/4T/Resource/function.R")
#按样本统计细胞类型
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = col21) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45))
ggsave('CellType/Stackbar_celltype.pdf', p, width = 8, height = 4)
ggsave('CellType/Stackbar_celltype.png', p, width = 8, height = 4)
ggsave('CellType/Stackbar_celltype.eps', p, width = 8, height = 4)

##Marker基因可视化
Idents(scRNA) <- "celltype"
markers <- c("CD3D", "CD3E", "CD4", "CD8A", "CD8B", "LEF1", "CCR7", "LYZ", "S100A8",
             "GNLY", "NKG7", "MZB1", "JCHAIN", "IGHG1", "MKI67", "TOP2A",
             "SLC4A10", "MS4A1", "CD79A", "PPBP", "PF4", "LILRA4", "IL3RA" # "RORC",
)
markers <- CaseMatch(markers, rownames(scRNA))

if(F){
  dir.create("CellType/MultipleGeneMap")
  OS <- scRNA
  #多基因umap图
  #NK_cells
  OS$NK_cells=apply(OS[['RNA']]@counts[c("KLRF1","PRF1","NKG7"),],2,function(x){ifelse(all(x>0),1,0)})
  p5 <- DimPlot(OS,reduction='umap',group.by='NK_cells',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/NK_cells.png", p5, width = 12, height = 6)
  #CD8_T_cells
  OS$CD8_T_cells=apply(OS[['RNA']]@counts[c("CD8A","CD8B"),],2,function(x){ifelse(all(x>0),1,0)})
  p2 <- DimPlot(OS,reduction='umap',group.by='CD8_T_cells',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/CD8_T_cells.png", p2, width = 12, height = 6)
  #CD4_T_cells
  OS$CD4_T_cells=apply(OS[['RNA']]@counts[c("CD4","LEF1"),],2,function(x){ifelse(all(x>0),1,0)})
  p1 <- DimPlot(OS,reduction='umap',group.by='CD4_T_cells',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/CD4_T_cells.png", p1, width = 12, height = 6)
  #Monocytes
  OS$Monocytes=apply(OS[['RNA']]@counts[c("LYZ","S100A8","S100A9"),],2,function(x){ifelse(all(x>0),1,0)})
  p3 <- DimPlot(OS,reduction='umap',group.by='Monocytes',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/Monocytes.png", p3, width = 12, height = 6)
  #Plasma_cells
  OS$Plasma_cells=apply(OS[['RNA']]@counts[c("MZB1","JCHAIN"),],2,function(x){ifelse(all(x>0),1,0)})
  p11 <- DimPlot(OS,reduction='umap',group.by='Plasma_cells',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/Plasma_cells.png", p11, width = 12, height = 6)
  #B_Cells
  OS$B_Cells=apply(OS[['RNA']]@counts[c("CD79A","CD79B","MS4A1"),],2,function(x){ifelse(all(x>0),1,0)})
  p7 <- DimPlot(OS,reduction='umap',group.by='B_Cells',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/B_Cells.png", p7, width = 12, height = 6)
  #Cycling_plasma_cells
  OS$Cycling_plasma_cells=apply(OS[['RNA']]@counts[c("MKI67","TOP2A"),],2,function(x){ifelse(any(x>0),1,0)})
  p15 <- DimPlot(OS,reduction='umap',group.by='Cycling_plasma_cells',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/Cycling_plasma_cells.png", p15, width = 12, height = 6)
  #MAITs
  OS$MAITs=apply(OS[['RNA']]@counts[c("SLC4A10","RORC"),],2,function(x){ifelse(any(x>0),1,0)})
  p8 <- DimPlot(OS,reduction='umap',group.by='MAITs',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/MAITs.png", p8, width = 12, height = 6)
  #pDCs
  OS$pDCs=apply(OS[['RNA']]@counts[c("LILRA4","IL3RA"),],2,function(x){ifelse(any(x>0),1,0)})
  p17 <- DimPlot(OS,reduction='umap',group.by='pDCs',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/pDCs.png", p17, width = 12, height = 6)
  #Platelets
  OS$Platelets=apply(OS[['RNA']]@counts[c("PPBP","PF4"),],2,function(x){ifelse(any(x>0),1,0)})
  p12 <- DimPlot(OS,reduction='umap',group.by='Platelets',cols=c('grey','red')) + NoLegend()
  ggsave("CellType/MultipleGeneMap/Platelets.png", p12, width = 12, height = 6)
  
  p <- (p5|p2|p1|p3|p11)/(p7|p15|p8|p17|p12)
  ggsave("CellType/MultipleGeneMap/Manymarkers.png",plot = p, width = 30, height = 12)
  ggsave("CellType/MultipleGeneMap/Manymarkers.pdf",plot = p, width = 30, height = 12)
  ggsave("CellType/MultipleGeneMap/Manymarkers.eps",plot = p, width = 30, height = 12)
}



##==================提取细胞子集做亚群鉴定========================##
rm(list=ls())
load("scRNA_classify.Rdata")

##===提取细胞子集===##
dir.create("Subset")

if(F){
  names(table(scRNA$celltype))
  Idents(scRNA) <- "celltype"
  scRNAsub <- subset(scRNA, idents = c("CD8_T_cells"))
}
#净化数据
DefaultAssay(scRNAsub) <- "RNA"
scRNAsub@assays$SCT <- NULL
colnames(scRNAsub@meta.data)
scRNAsub@meta.data <- scRNAsub@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "proj",
                                            "S.Score", "G2M.Score", "Phase","celltype")]
scRNAsub$main.celltype <- scRNAsub$celltype
scRNAsub$celltype <- NULL

##提取子集之后需要从头分析吗？
var.all <- VariableFeatures(scRNA, assay = "RNA")
var.sub <- FindVariableFeatures(scRNAsub, assay = "RNA") %>% VariableFeatures(assay = "RNA")
var.share <- intersect(var.all, var.sub)
length(var.share)

##数据标准化
#log标准化
scRNAsub <- NormalizeData(scRNAsub) %>% FindVariableFeatures() %>% ScaleData(features = rownames(scRNAsub))

save(scRNAsub, file = "Subset/scRNAsub_CD8T_log.Rdata")

##使用harmony整合数据
library(harmony)
scRNAsub <- RunPCA(scRNAsub, npcs = 50, verbose = FALSE)
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
scRNAsub <- RunHarmony(scRNAsub, group.by.vars="orig.ident", assay.use = "RNA", max.iter.harmony = 15)
ElbowPlot(scRNAsub, ndims = 50)
pc.num = 1:30

##降维聚类
scRNAsub <- RunTSNE(scRNAsub, reduction="harmony", dims=pc.num) %>% RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% FindClusters(resolution=0.2) 
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "RNA_snn_res.0.2", label = T) + ggtitle("RNA_snn_res.0.2")
ggsave("Subset/Resolution_test.png", p1, width = 10, height = 8)

#查看结果后采用0.2的分辨率
scRNAsub$seurat_clusters <- scRNAsub$RNA_snn_res.0.2
Idents(scRNAsub) <- "RNA_snn_res.0.2"

##查看harmony的整合效果
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "seurat_clusters")
p2 <- DimPlot(scRNAsub, reduction = "umap", group.by = "orig.ident")
pc = p1 + p2
ggsave("Subset/Harmony_integr.pdf", pc, width = 8, height = 4)
ggsave("Subset/Harmony_integr.png", pc, width = 8, height = 4)
ggsave("Subset/Harmony_integr.eps", pc, width = 8, height = 4)

#使用分面图查看效果
p <- DimPlot(scRNAsub, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident", ncol = 4) + NoLegend()
ggsave("Subset/Harmony_facet.pdf", p, width = 10, height = 6)
ggsave("Subset/Harmony_facet.png", p, width = 10, height = 6)
ggsave("Subset/Harmony_facet.eps", p, width = 10, height = 6)

##保存结果
save(scRNAsub, file = "Subset/scRNAsub_harmony.Rdata")     

##===鉴定数据子集的细胞类型===##
##提取各个Cluster的marker genes
ClusterMarker <- FindAllMarkers(scRNAsub, assay = "RNA", slot = "data", only.pos = T, logfc.threshold = 0.25, min.pct = 0.1)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'Subset/ClusterMarker.csv', row.names=F)
#ClusterMarker <- read.csv('Subset/ClusterMarker.csv')
#提取差异显著的marker genes
top = 50   #可根据需要调整
TopMarkers1 <- ClusterMarker %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers,'Subset/TopMarkers.csv', row.names=F)

##提取没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", ClusterMarker$gene, ignore.case = F),]
top = 50   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers_noRibo <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers_noRibo,'Subset/TopMarkers_noRibo.csv', row.names=F)

##导出一个方便对比的表格(每行一个cluster的marker gene)
source("/media/4T/Resource/function.R")
TopMarkers2Lines(data = TopMarkers_noRibo, output = "Subset")


##导入人工鉴定结果
cell.type <- c("CD8+_Naive_T_cells", "CD8+_Effector_T_cells", "CD8+_Exhausted_T_cells", "CD8+_Memory_T_cells")
Idents(scRNAsub) <- "seurat_clusters"
names(cell.type) <- levels(scRNAsub)
scRNAsub <- RenameIdents(scRNAsub, cell.type)
scRNAsub$celltype <- Idents(scRNAsub)
Idents(scRNAsub) <- "seurat_clusters"

##鉴定结果可视化
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T) 
p2 <- DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5, group.by = "celltype") 
pc = p1|p2
ggsave("Subset/CellType_Custom.pdf", pc, width = 10, height = 4)
ggsave("Subset/CellType_Custom.png", pc, width = 10, height = 4)
ggsave("Subset/CellType_Custom.eps", pc, width = 10, height = 4)

p <- DimPlot(scRNAsub, reduction = 'umap', group.by = 'celltype', split.by = 'orig.ident', ncol = 4)
ggsave("Subset/CellType_orig.ident.pdf", p, width = 9.5, height = 4.5)
ggsave("Subset/CellType_orig.ident.png", p, width = 9.5, height = 4.5)
ggsave("Subset/CellType_orig.ident.eps", p, width = 9.5, height = 4.5)

##保存结果
save(scRNAsub, file = "Subset/scRNAsub_classify.Rdata")

##===结果可视化===##
##Stackbar细胞丰度柱状图
tmp <- select(scRNAsub@meta.data, c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}

#按样本统计细胞类型
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = col21) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45))
ggsave('Subset/Stackbar_celltype.pdf', p, width = 8, height = 4.5)
ggsave('Subset/Stackbar_celltype.png', p, width = 8, height = 4.5)
ggsave('Subset/Stackbar_celltype.eps', p, width = 8, height = 4.5)

##Marker基因小提琴图
p <- VlnPlot(scRNAsub, features = c("LEF1", "CCR7", "SELL", "GNLY", "GZMB", "TIGIT"), slot = "counts", log = TRUE, ncol = 6)
ggsave("Subset/Markers_vlnplot.png", p, width = 15, height = 5)


table(scRNAsub$celltype)
# CD8+_Naive_T_cells  CD8+_Effector_T_cells CD8+_Exhausted_T_cells    CD8+_Memory_T_cells 
#               1486                   1458                   1318                    474 


# Run Monocle2
if(F){
  library(monocle)
  library(Seurat)
  library(dplyr)
  load("scRNAsub_classify.Rdata")
  ##monocle分析
  ###Extract data, phenotype data, and feature data from the SeuratObject
  data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data =scRNAsub@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)
  ###Construct monocle cds
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  ###按cluster选取差异基因
  my_cds<-monocle_cds
  my_cds
  slotNames(my_cds)
  my_cds <- estimateSizeFactors(my_cds)
  my_cds <- estimateDispersions(my_cds)
  my_cds <- detectGenes(my_cds, min_expr = 0.1)
  head(fData(my_cds))
  summary(fData(my_cds)$num_cells_expressed)
  #keep only genes expressed in greater than 5% of cells
  fData(my_cds)$use_for_ordering <- fData(my_cds)$num_cells_expressed > 0.05 * ncol(my_cds)
  table(fData(my_cds)$use_for_ordering)
  #PCA
  plot_pc_variance_explained(my_cds, return_all = F)
  #reduceDimension
  my_cds <- reduceDimension(my_cds,max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',verbose = TRUE)
  #按默认条件聚类
  my_cds <- clusterCells(my_cds, verbose = F)
  plot_cell_clusters(my_cds, color_by = 'as.factor(Cluster)')
  ##Trajectory step 1: choose genes that define a cell's progress
  #提取差异基因
  clustering_DEG_genes <- differentialGeneTest(my_cds,fullModelFormulaStr = '~Cluster',cores = 6)#根据个人电脑CPU的core来选择cores的值，服务器可以选到16
  dim(clustering_DEG_genes)
  clustering_DEG_genes %>% arrange(qval) %>% head()
  my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  my_cds <- setOrderingFilter(my_cds, ordering_genes = my_ordering_genes)
  ##Trajectory step 2: reduce data dimensionality
  my_cds <- reduceDimension(my_cds, method = 'DDRTree')
  ##Trajectory step 3: order cells along the trajectory
  my_cds <- orderCells(my_cds)
  saveRDS(my_cds,file = "my_cds.rds")
  plot_cell_trajectory(my_cds, color_by = "seurat_clusters")
  plot_cell_trajectory(my_cds, color_by = "State")
  plot_cell_trajectory(my_cds, color_by = "Pseudotime")
  ##定义细胞名称
  cell.type <- c("CD8+_Naive_T_cells", "CD8+_Effector_T_cells", "CD8+_Exhausted_T_cells", "CD8+_Memory_T_cells")
  names(cell.type) <- levels(scRNAsub)
  scRNAsub <- RenameIdents(scRNAsub,cell.type)
  #寻找差异基因画热图
  diff_test_res <- differentialGeneTest(my_cds,fullModelFormulaStr = "~sm.ns(Pseudotime)")
  diff_test_res[,c("gene_short_name", "pval", "qval")]
  sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
  plot_pseudotime_heatmap(my_cds[sig_gene_names[1:100],],
                          num_clusters = 3,
                          cores = 40,
                          show_rownames = T)
  ggsave("heatmap.png",p)
  ggsave("heatmap.eps",p)
  ggsave("heatmap.pdf",p)
  
  #画拟时态图
  p <- plot_cell_trajectory(my_cds, color_by = "Pseudotime")
  ggsave("Pseudotime.png",p)
  ggsave("Pseudotime.eps",p)
  ggsave("Pseudotime.pdf",p)
  
  p <- plot_cell_trajectory(my_cds, color_by = "seurat_clusters") 
  ggsave("celltype.png",p)
  ggsave("celltype.eps",p)
  ggsave("celltype.pdf",p)
  
  p <- plot_cell_trajectory(my_cds, color_by = "State")
  ggsave("state.png",p)
  ggsave("state.eps",p)
  ggsave("state.pdf",p)
  
  p <- plot_cell_trajectory(my_cds,markers = c(sig_gene_names[1]),use_color_gradient = TRUE)
  ggsave("monocle_KLHL17.png",p)
  ggsave("monocle_KLHL17.eps",p)
  ggsave("monocle_KLHL17.pdf",p)
  
  p <- plot_cell_trajectory(my_cds,markers = c(sig_gene_names[2]),use_color_gradient = TRUE)
  ggsave("monocle_ISG15.png",p)
  ggsave("monocle_ISG15.eps",p)
  ggsave("monocle_ISG15.pdf",p)
  
  p <- plot_cell_trajectory(my_cds,markers = c(sig_gene_names[3]),use_color_gradient = TRUE)
  ggsave("monocle_TNFRSF18.png",p)
  ggsave("monocle_TNFRSF18.eps",p)
  ggsave("monocle_TNFRSF18.pdf",p)
  
  p <- plot_cell_trajectory(my_cds,markers = c(sig_gene_names[4]),use_color_gradient = TRUE)
  ggsave("monocle_TNFRSF4.png",p)
  ggsave("monocle_TNFRSF4.eps",p)
  ggsave("monocle_TNFRSF4.pdf",p)
  
  p <- plot_cell_trajectory(my_cds,markers = c(sig_gene_names[5]),use_color_gradient = TRUE)
  ggsave("monocle_SDF4.png",p)
  ggsave("monocle_SDF4.eps",p)
  ggsave("monocle_SDF4.pdf",p)
  
  p <- plot_cell_trajectory(my_cds,markers = c(sig_gene_names[6]),use_color_gradient = TRUE)
  ggsave("monocle_B3GALT6.png",p)
  ggsave("monocle_B3GALT6.eps",p)
  ggsave("monocle_B3GALT6.pdf",p)
  
}

##==================提取细胞子集做亚群鉴定========================##
rm(list=ls())
load("scRNA_classify.Rdata")

##===提取细胞子集===##
dir.create("Subset2")
if(F){
  names(table(scRNA$celltype))
  Idents(scRNA) <- "celltype"
  scRNAsub <- subset(scRNA, idents = c("Monocytes"))
}
#净化数据
DefaultAssay(scRNAsub) <- "RNA"
scRNAsub@assays$SCT <- NULL
colnames(scRNAsub@meta.data)
scRNAsub@meta.data <- scRNAsub@meta.data[,c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "proj",
                                            "S.Score", "G2M.Score", "Phase","celltype")]
scRNAsub$main.celltype <- scRNAsub$celltype
scRNAsub$celltype <- NULL

##提取子集之后需要从头分析吗？
var.all <- VariableFeatures(scRNA, assay = "RNA")
var.sub <- FindVariableFeatures(scRNAsub, assay = "RNA") %>% VariableFeatures(assay = "RNA")
var.share <- intersect(var.all, var.sub)
length(var.share)

##数据标准化
#log标准化
scRNAsub <- NormalizeData(scRNAsub) %>% FindVariableFeatures() %>% ScaleData(features = rownames(scRNAsub))

save(scRNAsub, file = "Subset2/scRNAsub_MY_log.Rdata")

##使用harmony整合数据
library(harmony)
scRNAsub <- RunPCA(scRNAsub, npcs = 50, verbose = FALSE)
DimPlot(scRNAsub, reduction = "pca", group.by = "orig.ident")
scRNAsub <- RunHarmony(scRNAsub, group.by.vars="orig.ident", assay.use = "RNA", max.iter.harmony = 15)
ElbowPlot(scRNAsub, ndims = 50)
pc.num = 1:30

##降维聚类
scRNAsub <- RunTSNE(scRNAsub, reduction="harmony", dims=pc.num) %>% 
  RunUMAP(reduction="harmony", dims=pc.num) %>%
  FindNeighbors(reduction="harmony", dims=pc.num) %>% 
  FindClusters(resolution=0.1) 
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T) + ggtitle("RNA_snn_res.0.1")
ggsave("Subset2/Resolution_test.png", p1, width = 10, height = 8)

#查看结果后采用0.1的分辨率
scRNAsub$seurat_clusters <- scRNAsub$RNA_snn_res.0.1
Idents(scRNAsub) <- "RNA_snn_res.0.1"

##查看harmony的整合效果
p1 <- DimPlot(scRNAsub, reduction = "umap", group.by = "seurat_clusters")
p2 <- DimPlot(scRNAsub, reduction = "umap", group.by = "orig.ident")
pc = p1 + p2 
ggsave("Subset2/Harmony_integr.pdf", pc, width = 8, height = 4)
ggsave("Subset2/Harmony_integr.png", pc, width = 8, height = 4)
ggsave("Subset2/Harmony_integr.eps", pc, width = 8, height = 4)

#使用分面图查看效果
p <- DimPlot(scRNAsub, reduction = "umap", group.by = "orig.ident", split.by = "orig.ident", ncol = 4) + NoLegend()
ggsave("Subset2/Harmony_facet.pdf", p, width = 10, height = 6)
ggsave("Subset2/Harmony_facet.png", p, width = 10, height = 6)

##保存结果
save(scRNAsub, file = "Subset2/scRNAsub_MY_harmony.Rdata")     

##===鉴定数据子集的细胞类型===##
##提取各个Cluster的marker genes
ClusterMarker <- FindAllMarkers(scRNAsub, assay = "RNA", slot = "data", only.pos = T, logfc.threshold = 0.25, min.pct = 0.1)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
write.csv(ClusterMarker,'Subset2/ClusterMarker.csv', row.names=F)
#ClusterMarker <- read.csv('Subset/ClusterMarker.csv')
#提取差异显著的marker genes
top = 100   #可根据需要调整
TopMarkers1 <- ClusterMarker %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers,'Subset2/TopMarkers.csv', row.names=F)

##提取没有核糖体的Markers
ClusterMarker_noRibo <- ClusterMarker[!grepl("^RP[SL]", ClusterMarker$gene, ignore.case = F),]
top = 100   #可根据需要调整
TopMarkers1 <- ClusterMarker_noRibo %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_noRibo %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers_noRibo <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers_noRibo,'Subset2/TopMarkers_noRibo.csv', row.names=F)

##导出一个方便对比的表格(每行一个cluster的marker gene)
source("/media/4T/Resource/function.R")
TopMarkers2Lines(data = TopMarkers_noRibo, output = "Subset2")

##使用main.celltype作为参考
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T) 
p2 <- DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5, group.by = "main.celltype")
pc = p1|p2
ggsave("Subset2/CellType_Ref.png", pc, width = 10, height = 4)

##导入人工鉴定结果
cell.type <- c("FCGR3A+_Monocytes", "CD14+_Monocytes", "CD14+_Monocytes", "Plasma-like_Monocytes")
Idents(scRNAsub) <- "seurat_clusters"
names(cell.type) <- levels(scRNAsub)
scRNAsub <- RenameIdents(scRNAsub, cell.type)
scRNAsub$celltype <- Idents(scRNAsub)
Idents(scRNAsub) <- "seurat_clusters"

##鉴定结果可视化
p1 <- DimPlot(scRNAsub, reduction = "umap", label = T) 
p2 <- DimPlot(scRNAsub, reduction = "umap", label = T, label.size = 2.5, group.by = "celltype")
pc = p1|p2
ggsave("Subset2/CellType_Custom.pdf", pc, width = 10, height = 4)
ggsave("Subset2/CellType_Custom.png", pc, width = 10, height = 4)
ggsave("Subset2/CellType_Custom.eps", pc, width = 10, height = 4)

p <- DimPlot(scRNAsub, reduction = 'umap', group.by = 'celltype', split.by = 'orig.ident', ncol = 4)
ggsave("Subset2/CellType_orig.ident.pdf", p, width = 9.5, height = 4.5)
ggsave("Subset2/CellType_orig.ident.png", p, width = 9.5, height = 4.5)
ggsave("Subset2/CellType_orig.ident.eps", p, width = 9.5, height = 4.5)

##保存结果
save(scRNAsub, file = "Subset2/scRNAsub_classify.Rdata")

##===结果可视化===##
##Stackbar细胞丰度柱状图
tmp <- select(scRNAsub@meta.data, c("orig.ident", "celltype"))
df <- data.frame()
for(i in unique(tmp$orig.ident)){
  df_i <- subset(tmp, tmp$orig.ident==i) %>% pull(celltype) %>% table()
  df_i <- data.frame(sample=rep(i, length(df_i)), value=df_i)
  names(df_i) <- c("sample", "celltype", "value")
  df <- rbind(df, df_i)
}
source("/media/4T/Resource/function.R")
#按样本统计细胞类型
p <- ggplot(df, aes(x=sample, y=value, fill=celltype)) + 
  geom_bar(stat= "identity", position = "fill") + 
  scale_fill_manual(values = col21) +
  labs(x = 'Samples', y = 'Relative Abundance', title = 'Samples composition') +
  theme(panel.grid = element_blank(), panel.background = element_blank(), 
        axis.text.x = element_text(angle=45))
ggsave('Subset2/Stackbar_celltype.pdf', p, width = 8, height = 4.5)
ggsave('Subset2/Stackbar_celltype.png', p, width = 8, height = 4.5)
ggsave('Subset2/Stackbar_celltype.eps', p, width = 8, height = 4.5)



##Marker基因可视化
Idents(scRNAsub) <- "celltype"
markers <- c("CD14", "FCGR3A", "LYZ", "S100A8", "S100A9", "MZB1", "JCHAIN")
markers <- CaseMatch(markers, rownames(scRNAsub))

##Marker基因小提琴图
p <- VlnPlot(scRNAsub, features = markers, slot = "counts", pt.size = 0 ,log = TRUE, ncol = 7)
ggsave("Subset2/Markers_vlnplot.png", p, width = 20, height = 6)
ggsave("Subset2/Markers_vlnplot.pdf", p, width = 20, height = 6)
ggsave("Subset2/Markers_vlnplot.eps", p, width = 20, height = 6)

