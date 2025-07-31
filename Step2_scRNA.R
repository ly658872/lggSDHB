rm(list = ls());gc()
if(!require(multtest))install.packages('multtest')
if(!require(metap))install.packages('metap')
if(!require(Seurat))install.packages('Seurat')
if(!require(mindr))install.packages('mindr')
if(!require(tidyverse))install.packages('tidyverse')
if(!require(dplyr))install.packages('dplyr')
library(data.table)
#----1.Input Data-----
####仅有一个稀疏矩阵时的读取方法#####

matrix_data <- fread("GSE117891_all_6148.umi.count.matrix.tsv.gz",data.table = F)
gene <- matrix_data$V1
table(str_split_fixed(colnames(matrix_data),pattern = '_',n = 2)[,1])
exp <- matrix_data %>% 
  dplyr::select(dplyr::starts_with(c('S'))) 
rownames(exp) <- gene
seurat_obj <- CreateSeuratObject(counts = exp)

ifnb.list <- CreateSeuratObject(counts = exp,min.cells = 3,
                                min.features = 300)
length(ifnb.list@meta.data$orig.ident) 
#计算线粒体RNA MT=homo mt=mouse 
ifnb.list[['percent.mt']] <- PercentageFeatureSet(ifnb.list,pattern  = '^MT-')
head(ifnb.list@meta.data)
#计算红细胞数量）
HB.genes<-c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(ifnb.list@assays$RNA)) 
HB.genes <- rownames(ifnb.list@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
ifnb.list[["percent.HB"]]<-PercentageFeatureSet(ifnb.list, features=HB.genes) 
#可视化数量 
VlnPlot(ifnb.list,
        features = c('nCount_RNA','nFeature_RNA',
                     'percent.mt','percent.HB'),
        ncol = 4,pt.size = 0)

ifnb.list <- subset(ifnb.list,
                    subset=nFeature_RNA>300 & nFeature_RNA <8000 & percent.mt < 15 & percent.HB< 0.1)

saveRDS(ifnb.list,file = 'GSE117891_seurat.rds')

### 数据预处理 #####
rm(list = ls());gc()
library(ggplot2)
seurat_obj <- readRDS('GSE117891_seurat.rds')
#合并前将count数据标准化
seurat_obj <- NormalizeData(seurat_obj,
                            normalization.method = 'LogNormalize',
                            scale.factor = 10000)
#寻找VariableFeatures
seurat_obj <- FindVariableFeatures(seurat_obj,selection.method = 'vst',
                                   nfeatures = 2000)
# top10 <- head(VariableFeatures(seurat_obj),10)
# plot1 <- VariableFeaturePlot(seurat_obj)
# plot2 <- LabelPoints(plot1,points = top10,repel = T)
# plot1+plot2

#将数据标准化
seurat_obj <- ScaleData(seurat_obj,
                        features = rownames(seurat_obj)) 
#RunPCA
seurat_obj <- seurat_obj %>% 
  RunPCA(.,features = VariableFeatures(seurat_obj)) 

ElbowPlot(seurat_obj)

# #基于Harmony包对不同样本间结合
library(harmony)
library(ggsci)
seurat_harmony <- RunHarmony(seurat_obj,"orig.ident", plot_convergence = T)
seurat_harmony <- SetIdent(seurat_harmony,value = "orig.ident")

  #plot_annotation(title = "")

seurat_harmony <- seurat_harmony %>%
  FindNeighbors(reduction = "harmony", k.param = 20, dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunTSNE(reduction = "harmony", dims = 1:20, verbose = F) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = F) %>%
  identity()
#seurat_harmony$RNA_snn_res.0.2
DimPlot(seurat_harmony,reduction = "umap")

### 注释 ####

celltype_marker=c(
  'MBP','MOG','PLP1','MAG', #Oligodendrocyte
  'CD3D','CD3E','CD8A', #T cellS 
  'GFAP','AQP4','SOX9','CLU', #Astrocyte 
  'TMEM119','CX3CR1','P2RY12' # Migcroglial Cell
)

Idents(seurat_harmony) <- 'seurat_clusters'
VlnPlot(seurat_harmony,features = celltype_marker,pt.size = 0,ncol = 3)
Oligodendrocyte=c(1)
T_cells=c(10,11,15) 
Astrocyte=c(12,7,2,3,4,5,6,13)
Migcroglial_cell=c(0,8,9,14) #21

current.cluster.ids <- c(Oligodendrocyte,
                         T_cells,
                         Astrocyte,
                         Migcroglial_cell)
new.cluster.ids <- c(rep("Oligodendrocyte",length(Oligodendrocyte)),
                     rep("T_cells",length(T_cells)),
                     rep("Astrocyte",length(Astrocyte)),
                     rep("Migcroglial_cell",length(Migcroglial_cell)))
seurat_harmony@meta.data$Celltype <- plyr::mapvalues(x = as.integer(as.character(seurat_harmony@meta.data$seurat_clusters)),
                                                     from = current.cluster.ids, to = new.cluster.ids)
table(seurat_harmony@meta.data$Celltype)
saveRDS(seurat_harmony,file = 'GSE117891_seurat_anno.rds')

#### Umap ######
rm(list = ls());gc()
library(Seurat)
library(ggunchull)
library(magrittr)

pbmc <- readRDS('GSE117891_seurat_anno.rds')
Idents(pbmc) <- 'Celltype'

umap <- pbmc@reductions$umap@cell.embeddings %>% as.data.frame
cellType <- pbmc@active.ident %>%
  as.data.frame %>%
  `colnames<-` ("Celltype")
umap <- cbind(umap, cellType)
# "#00468BFF" "#ED0000FF" "#42B540FF" "#0099B4FF" "#925E9FFF" "#FDAF91FF" "#AD002AFF" "#ADB6B6FF" "#1B1919FF""
ggsci::pal_lancet('lanonc',alpha = 1)(10)
allcolour <- c("#00468BD8", "#ED0000D8", "#42B540D8", "#0099B4D8","#9370DB")

myTheme <- theme(  panel.background = element_blank(),
                   panel.grid = element_blank(),
                   axis.line = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   legend.background = element_blank(),
                   legend.key = element_blank(),
                   legend.key.size = unit(1, "cm"),
                   legend.text = element_text(size = 12),)  
xstart <- min(umap$umap_1); 
xend <- min(umap$umap_1) + 3
ystart <- min(umap$umap_2); 
yend <- min(umap$umap_2) + 3

p <- ggplot(umap,aes(x= umap_1 ,
                     y = umap_2 ,colour = Celltype)) +
  geom_point(size = 1 , alpha = 0.5) +
  scale_colour_manual(values = allcolour)
p2 <- p + myTheme + 
  guides(colour = guide_legend(override.aes = list(size = 6)))
p3 <- p2 +   geom_segment(aes(x = xstart, y = ystart , xend = xend, yend = ystart),
                          colour = "black", linewidth=1,arrow = arrow(length = unit(0.3,"cm")))+
  geom_segment(aes(x = xstart, y = ystart, xend = xstart , yend = yend),
               colour = "black", linewidth=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = xstart +1.5, y = ystart -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) +
  annotate("text", x = xstart -1, y = ystart + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
p3

cellTypeLoc <- umap %>%
  group_by(Celltype) %>%
  summarise(Loc1 = median(umap_1),
            Loc2 = median(umap_2)  ) %>%
  as.data.frame
##### 添加标签
 p4 <- p3 + geom_text(data = cellTypeLoc,
                      mapping = aes(x = Loc1, y = Loc2,
                                    label = Celltype),
                      color="black",size=5)
 p4
 cellTypeLoc <- umap %>%
   group_by(Celltype) %>%
   summarise(Loc1 = median(umap_1),
             Loc2 = median(umap_2)) %>%
   as.data.frame
##### 添加标签
 p4 <- p3 + geom_text(data = cellTypeLoc,
                      mapping = aes(x = Loc1,
                                    y = Loc2,
                                    label = Celltype),
                      color="black",size=5)
 p4
 library(ggunchull)
 p4 + stat_unchull(fill = "white",
                   alpha = 0,
                   show.legend = FALSE,
                   nsm = 20,
                   nbin = 200,
                   sfac = 1.5)
 ggsave2(filename = 'plot/scRNA_umap_anno.pdf',width = 8,height = 6)
 
##### Dotplot ######
 
 pbmc <- readRDS('GSE117891_seurat_anno.rds')
 Idents(pbmc) <- 'Celltype'
 makers <- c(
     'MBP','MOG','PLP1', #Oligodendrocyte
     'CD3D','CD3E','CD8A', #T cellS 
     'GFAP','AQP4','SOX9', #Astrocyte 
     'TMEM119','CX3CR1','P2RY12' # Migcroglial Cell
     #"CD163","MRC1" #macrophage
   )
 
 top3pbmc.markers<-FindAllMarkers(pbmc,
                                only.pos=TRUE,
                                min.pct=0.25)
 top3pbmc.markers <- top3pbmc.markers[top3pbmc.markers$gene %in% makers,]
 top3pbmc.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
 # saveRDS(top3pbmc.markers,file = 'cell_degs_scRNA.rds')
 top3pbmc.markers <- readRDS('cell_degs_scRNA.rds')
 p<-DotPlot(pbmc,
              features=split(makers,top3pbmc.markers$cluster),
              cols=c("#ffffff","#AD002AFF")
 )+
   RotatedAxis()+#来自Seurat
   theme(
     panel.border=element_rect(color="black"),
     panel.spacing=unit(1,"mm"),
     axis.title=element_blank(),
     axis.text.y=element_blank(),
   )
 p
 p$data$feature.groups2<-factor(p$data$feature.groups,
                                  levels=c("Oligodendrocyte","T_cells","Astrocyte",
                                             "Migcroglial_cell"))
 library(ggh4x)
 allcolour <- c("#00468BD8", "#ED0000D8", "#42B540D8", "#0099B4D8")
 
 strip<-strip_themed(
   background_x=elem_list_rect(fill=allcolour))
 p$data %>%
   ggplot(aes(x=features.plot,
              y=id))+
   geom_point(aes(size=pct.exp,
                  color=avg.exp.scaled))+
   facet_wrap2(~feature.groups2,
               scales="free_x",
               strip=strip,
               nrow=1)+
   theme_classic()+
   theme(axis.text.x=element_text(angle=90,
                                    hjust=0.5,
                                    vjust=0.3,
                                    color="black"),
         axis.title=element_blank(),
         strip.background=element_rect(color="white"),
         axis.text.y=element_blank())+
   scale_color_gradient(low="#ffffff",
                        high="#AD002AFF",
                        name="avg.exp")-> p
 
 p
 df<-data.frame(x=0,y=levels(pbmc),stringsAsFactors=F)
 df$y<-factor(df$y,levels=df$y)
 
 p1<-ggplot(df,aes(x,y,color=factor(y)))+
   geom_point(size=6,shape=15,show.legend=F)+
   scale_color_manual(values=rev(allcolour))+
   theme_classic()+
   scale_x_continuous(expand=c(0,0))+
   theme(
     plot.margin=margin(r=0),
     axis.title=element_blank(),
     axis.text.x=element_blank(),
     axis.text.y=element_text(size=9),
     axis.ticks=element_blank(),
     axis.line=element_blank()
   )
 
 p1
 cowplot::plot_grid(p1,p,align="h",axis="bt",rel_widths=c(1.5,9))
 ggsave2(filename = 'plot/Dot_marker.pdf',width = 8,height = 6)
 
 #### Vlnplot ####
 rm(list = ls());gc()
 library(Seurat)
 library(irGSEA)
 scRNA_Ams <- readRDS('GSE117891_seurat_anno.rds')
 Idents(scRNA_Ams) <- 'Celltype'
 genes <- clusterProfiler::read.gmt('WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA.v2024.1.Hs.gmt')
 genes_Ams=genes$gene
 genes_Ams=genes_Ams[genes_Ams%in%rownames(scRNA_Ams)]
 genes_Ams=list(genes_Ams)
 names(genes_Ams) <- 'METCGs'
 
 scRNA_Ams<-SeuratObject::UpdateSeuratObject(scRNA_Ams)
 pbmc3k.final2<-CreateSeuratObject(counts=CreateAssay5Object(GetAssayData(scRNA_Ams,assay="RNA",slot="counts")),
                                     meta.data=scRNA_Ams[[]])
 pbmc3k.final2<-NormalizeData(pbmc3k.final2)
 pbmc3k.final2<-irGSEA.score(object=pbmc3k.final2,assay="RNA",
                               slot="data",seeds=123,ncores=2,
                               min.cells=3,min.feature=0,
                               custom=T,geneset=genes_Ams,
                               species="Homo sapiens",
                               geneid="symbol",
                               method=c("AUCell","UCell","singscore",
                                        "ssgsea","JASMINE","viper"),
                               aucell.MaxRank=NULL,ucell.MaxRank=NULL,
                               kcdf='Gaussian')
 Assays(pbmc3k.final2)
 pbmc3k.final2@assays$AUCell
 pbmc3k.final2$Celltype
 result.dge<-irGSEA.integrate(object=pbmc3k.final2,
                                group.by="Celltype",
                                method=c("AUCell","UCell","singscore",
                                           "ssgsea","JASMINE","viper"))
 geneset.show<- result.dge$RRA %>%
   dplyr::filter(pvalue<=0.05) %>%
   dplyr::pull(Name) %>%
   unique(.)
 
saveRDS(pbmc3k.final2,file = 'seurat_score.rds')
pbmc3k.final2@assays$AUCell$scale.data
Idents(pbmc3k.final2) <- 'Celltype'
allcolour <- c( "#ED0000D8", "#00468BD8", "#0099B4D8","#42B540D8")

vlnplot<-irGSEA.vlnplot(object=pbmc3k.final2,color.cluster = allcolour,
                          method=c("AUCell","UCell","singscore","ssgsea",
                                     "JASMINE","viper"),show.geneset =1)
vlnplot
cowplot::ggsave2(filename = 'plot/Score_vln.pdf',width = 8,height = 8)
### 统计 #####
library(Seurat)
library(dplyr)

target_celltypes <- c("Astrocyte", "Migcroglial_cell",
                      "Oligodendrocyte", "T_cells")  # 示例名称
# 获取细胞类型标识
celltypes <- as.character(Idents(pbmc3k.final2))
scores <- pbmc3k.final2@assays[['viper']]$data
# 创建统计检验函数
compare_astrocytes <- function(algorithm, geneset_name = NULL) {
  # 提取评分数据（处理向量格式）
  scores <- pbmc3k.final2@assays[[algorithm]]$scale.data
  
  # 创建数据框
  df <- data.frame(
    score = scores,
    celltype = celltypes,
    stringsAsFactors = FALSE
  ) %>% 
    filter(celltype %in% target_celltypes)
  
  # 两两Wilcoxon检验
  results <- list()
  ref_type <- "Astrocyte"
  other_types <- setdiff(target_celltypes, ref_type)
  
  for (type in other_types) {
    group1 <- df$score[df$celltype == ref_type]
    group2 <- df$score[df$celltype == type]
    
    test_res <- wilcox.test(group1, group2, alternative = "two.sided")
    
    results[[paste(ref_type, "vs", type)]] <- data.frame(
      Algorithm = algorithm,
      Comparison = paste(ref_type, "vs", type),
      W_statistic = test_res$statistic,
      P_value = test_res$p.value,
      n_ref = length(group1),
      n_other = length(group2)
    )
  }
  
  do.call(rbind, results)
}

# 执行所有算法检验
algorithms <- c("AUCell", "UCell", "singscore", "ssgsea", "JASMINE", "viper")
all_results <- lapply(algorithms, function(alg) {
  if (alg %in% names(pbmc3k.final2@assays)) {
    tryCatch(compare_astrocytes(alg),
             error = function(e) message("Error in ", alg, ": ", e$message))
  } else {
    message("Skipping ", alg, ": assay not found")
    NULL
  }
}) %>% bind_rows()
# 校正P值（处理可能的NA值）
all_results$Adj_P_value <- p.adjust(all_results$P_value,
                                    method = "bonferroni")
data.table::fwrite(all_results,file = 'Stats_revised.csv',row.names = T)
##### Degs ####
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
obj.markers <- readRDS('cell_degs_scRNA.rds')
Symbol<-mapIds(get("org.Hs.eg.db"),keys=obj.markers$gene,
               keytype="SYMBOL",column="ENTREZID")
ids<-bitr(obj.markers$gene,"SYMBOL","ENTREZID","org.Hs.eg.db")

#ENTREZID to obj.markers
data<-merge(obj.markers,ids,by.x="gene",by.y="SYMBOL")
gcSample<-split(data$ENTREZID,data$cluster)
gcSample
xx<-compareCluster(gcSample,fun="enrichKEGG",
                   organism="hsa",pvalueCutoff=1,qvalueCutoff=1)
res<-xx@compareClusterResult

for(i in 1:dim(res)[1]){
  arr=unlist(strsplit(as.character(res[i,"geneID"]),split="/"))
  gene_names=paste(unique(names(Symbol[Symbol%in%arr])),collapse="/")
  res[i,"geneID"]=gene_names
}

head(res)
data.table::fwrite(res,file = 'enrichKEGG_scRNA.csv',row.names = T)
##pathway iden
table(res$Cluster)
enrich<-res%>%
  group_by(Cluster)%>%
  dplyr::top_n(n=5,wt=-pvalue)%>%
  filter(Cluster %in% c("Astrocyte",'Migcroglial_cell',
                        'Oligodendrocyte','T_cells'))

dt<-enrich
dt<-dt[order(dt$Cluster),]
dt$Description<-factor(dt$Description,levels=unique(dt$Description))
colnames(dt)
#
library(ggplot2)
mytheme<- theme(
  axis.title=element_text(size=13),
  axis.text=element_text(size=11),
  axis.text.y=element_blank(),#在自定义主题中去掉y轴通路标签:
  axis.ticks.length.y=unit(0,"cm"),
  plot.title=element_text(size=13,hjust=0.5,face="bold"),
  legend.title=element_text(size=13),
  legend.text=element_text(size=11),
  #plot.margin=margin(t=5.5,r=10,l=5.5,b=5.5)
)
p<-ggplot(data=dt,aes(x=-log10(pvalue),
                      y=rev(Description),fill=Cluster))+
  scale_fill_manual(values=allcolour)+
  geom_bar(stat="identity",width=0.5,alpha=0.8)+
  scale_x_continuous(expand=c(0,0))+#调整柱子底部与y轴紧贴 "Astrocyte",'Migcroglial_cell','Oligodendrocyte','T_cells'
  labs(x="-Log10(pvalue)",y="",
       title="KEGG Pathway enrichment")+
  #x=0.61用数值向量控制文本标签起始位置
  geom_text(size=3.8,aes(x=0.05,label=Description),hjust=0)+#hjust=0,左对齐
  geom_text(size=3,aes(x=0.05,label=geneID),hjust=0,vjust=2.5,
            color=rep(allcolour,each=5))+#hjust=0,左对齐
  theme_classic()+
  mytheme+
  NoLegend()

p
ggsave(filename="pathway_celltypes.pdf",width=5.2,height=5.1,plot=p)

#### plot #####
library(ggplot2)
library(cowplot)
cols = c("gray", "red")
plot9 <- FeaturePlot(pbmc, features = 'SDHB',cols = cols, pt.size = 0.8)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 
cowplot::ggsave2(plot9,filename = 'SDHB_Feature_dis.pdf',width = 5,height = 4)

plot10 <- DotPlot(pbmc, features = 'SDHB',group.by = 'Celltype')+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
cowplot::ggsave2(plot10,filename = 'SDHB_Dot_dis.pdf',width = 5,height = 4)


