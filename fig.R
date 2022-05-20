setwd("D:/R/Rstudy/data/mouse/TAC/fb")

#R packge
library(Seurat)
library(cowplot)
library(tidyverse)
library(patchwork)
library(monocle)
library(Seurat)
library(cowplot)
library(tidyverse)
library(patchwork)
library(monocle)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)
library(clusterProfiler)
library(paletteer)
library(viridis)
library(lessR)
require(qdapTools)
require(REdaS)
library(ggrepel)
library(harmony)
library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ComplexHeatmap)
library(dorothea)
library(dplyr)
library(Seurat)
library(tibble)
library(pheatmap)
library(tidyr)
library(viper)
library(reshape2)
library(ggpubr)
library(stringr)
library(tidyr)
library(dplyr)
library(hash)
library(monocle)
library(scales)
library(enrichplot)
library(enrichR)
library(pheatmap)
library(clusterProfiler)
library(msigdbr)
library(fgsea)
library(DOSE)
library(RColorBrewer)
library(circlize)
library(viridis)
library(patchwork)
library(hrbrthemes)
library(circlize)
setwd('D:/R/Rstudy/data/mouse/TAC/rds')

#Sham1 data processing
Sham1 <-Read10X('D:/R/Rstudy/data/mouse/TAC/GSE155882_RAW/sham1/')
Sham1 <-CreateSeuratObject(counts = Sham1, project = "Sham1", min.cells = 3, min.features = 200)
Sham1[["percent.mt"]] <- PercentageFeatureSet(Sham1, pattern = "^mt-")
Sham1[["percent.rp"]]  = PercentageFeatureSet(Sham1, pattern = "^Rp[sl][[:digit:]]")
Sham1<- subset(Sham1, subset = nFeature_RNA >300 & nFeature_RNA < 3000 & percent.mt<12& percent.rp<20& nCount_RNA < 6000)
VlnPlot(Sham1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4,pt.size =0.1)
mm=Sham1
library(DoubletFinder)
mm=NormalizeData(mm, verbose =T) 
mm <- FindVariableFeatures(mm, selection.method = "vst")
mm <- ScaleData(mm, verbose = FALSE)
mm <- RunPCA(mm, verbose = T)

pct <- mm [["pca"]]@stdev / sum( mm [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

mm=RunTSNE(mm,dims = 1:20, verbose = T)
mm=FindNeighbors(mm, dims = 1:20,verbose = T)
mm=FindClusters(mm,resolution =0.3)

sweep.res.list <- paramSweep_v3(mm, PCs = 1:20, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)  
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() 

DoubletRate = ncol(mm)*8*1e-6

homotypic.prop <- modelHomotypic(mm$seurat_clusters)
nExp_poi <- round(DoubletRate*ncol(mm)) 
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

mm <- doubletFinder_v3(mm, PCs = 1:20, pN = 0.25, pK = pK_bcmvn, 
                       nExp = nExp_poi.adj, reuse.pANN = F, sct = F)

p=DimPlot(mm, group.by = "DF.classifications_0.25_0.2_447")

dpi=300
png(file="Sham1_feature.png", width = dpi*6, height = dpi*5, units = "px",res = dpi,type='cairo')
p
dev.off()


p=VlnPlot(mm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 2,pt.size =0.1,group.by = "DF.classifications_0.25_0.2_447")
dpi=300
png(file="Sham1_vln.png", width = dpi*10, height = dpi*12, units = "px",res = dpi,type='cairo')
p
dev.off()
mi.list <- SplitObject(mm, split.by = "DF.classifications_0.25_0.2_447")
mm=mi.list$'Singlet'
saveRDS(mm,"Sham1.rds")



----------------------------------------------------------------------------------------------------------
#merge and quality control,clustering etc
  mm <- merge(sham1, c(sham2,TAC1,TAC2,TAC_JQ1,TAC_JQ2,TAC_JQ_withdrawn1,TAC_JQ_withdrawn2), add.cell.ids=c("sham1","sham2","TAC1","TAC2","TJ1","TJ2","TJW1","TJW2"))

mm <- NormalizeData(mm, verbose = T)
mm <-FindVariableFeatures(mm, selection.method = "vst")
mm <- ScaleData(mm, verbose = T)#,vars.to.regress = "percent.mt"
mm <- RunPCA(mm, verbose = T)


pct <- mm [["pca"]]@stdev / sum( mm [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()



library(harmony)
mm=RunHarmony(mm,"orig.ident", plot_convergence = TRUE)
mm=RunUMAP(mm,reduction = "harmony", dims = 1:16) 
mm=FindNeighbors(mm,reduction = "harmony", dims = 1:16,verbose = T)
mm=FindClusters(mm,resolution =0.5)

----------------------------------------------------------------------------------------------------
#load data
mm <- readRDS("D:/R/Rstudy/data/mouse/TAC/rds/all_anno.rds")


#fig 1b
p=DimPlot(mm,cols =pal,label =F,repel = T)+umap_theme+NoLegend()

dpi=300
png(file="1b.png", width = dpi*8, height = dpi*8, units = "px",res = dpi,type='cairo')
p
dev.off()

#fig 1c
mm=subset(mm,idents=c("Ly6a+ CF","Activated CF","Recruited MAC","Resident MAC","G0s2+ CF"))


#fig 1c
mm=subset(mm,idents=c("Ly6a+ CF","Activated CF","G0s2+ CF"))
VlnPlot(mm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0.1)
mm<- subset(mm, subset = nFeature_RNA >200 & nFeature_RNA < 3000 & percent.mt<8)#& percent.rp<25
VlnPlot(mm, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size =0.1)
mm <- NormalizeData(mm, verbose = T)
mm <-FindVariableFeatures(mm, selection.method = "vst")
mm <- ScaleData(mm, verbose = T)#,vars.to.regress = "nFeature_RNA"
mm <- RunPCA(mm, verbose = T)
pct <- mm [["pca"]]@stdev / sum( mm [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

mm=RunHarmony(mm,"condition", plot_convergence = TRUE)
mm=RunUMAP(mm,reduction = "harmony", dims = 1:12) 
mm=FindNeighbors(mm,reduction = "harmony", dims = 1:12,verbose = T)
mm=FindClusters(mm,resolution =0.8)

p1=VlnPlot(mm, features = c("nFeature_RNA"))
p2=VlnPlot(mm, features = c("percent.mt"))
p1 | p2
png(file="zhikong.png", width = dpi*16, height = dpi*6, units = "px",res = dpi,type='cairo')
p1 | p2
dev.off()
DimPlot(mm,label = T,split.by = "condition")
DimPlot(mm, label = TRUE)
DimPlot(mm, label = TRUE,split.by = "orig.ident",pt.size = 1)

mm.markers <- FindAllMarkers(mm, only.pos =TRUE, min.pct =0.25, logfc.threshold =0.25)
top5=mm.markers %>% group_by(cluster) %>% top_n(n =20, wt = avg_log2FC)
DoHeatmap(mm, features = top5$gene)+theme(plot.title = element_text(hjust = 0.5,size=30))+scale_fill_gradientn(colors = c("blue", "white", "red"))
mm=subset(mm,idents=c("1","2","3","4","5","6","7","8","9"))

mm <- NormalizeData(mm, verbose = T)
mm <-FindVariableFeatures(mm, selection.method = "vst")
mm <- ScaleData(mm, verbose = T)#,vars.to.regress = "nFeature_RNA"
mm <- RunPCA(mm, verbose = T)
pct <- mm [["pca"]]@stdev / sum( mm [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

mm=RunHarmony(mm,"condition", plot_convergence = TRUE)
mm=RunUMAP(mm,reduction = "harmony", dims = 1:11) 
mm=FindNeighbors(mm,reduction = "harmony", dims = 1:11,verbose = T)
mm=FindClusters(mm,resolution =0.7)
DimPlot(mm, label = TRUE)
p1=VlnPlot(mm, features = c("nFeature_RNA"))
p2=VlnPlot(mm, features = c("percent.mt"))
p1 | p2

mm.markers <- FindAllMarkers(mm, only.pos =TRUE, min.pct =0.25, logfc.threshold =0.25)
top5=mm.markers %>% group_by(cluster) %>% top_n(n =30, wt = avg_log2FC)
DoHeatmap(mm, features = top5$gene)+theme(plot.title = element_text(hjust = 0.5,size=30))+scale_fill_gradientn(colors = c("blue", "white", "red"))


p=DimPlot(mm,label =T,repel = T,
          label.size = 8,pt.size = 0.8)+umap_theme+NoLegend()
dpi=300
png(file="1c.png", width = dpi*7, height = dpi*6, units = "px",res = dpi,type='cairo')
p
dev.off()

#fig 1d
DimPlot(mm,label = T)
mm$celltype=Idents(mm) 
#table(mm$celltype)
av <-AverageExpression(mm,
                       group.by = "celltype",
                       assays = "RNA")
av=av[[1]]
head(av)

cg=names(tail(sort(apply(av, 1, sd)),2000))
#View(av[cg,])
#View(cor(av[cg,],method = 'spearman'))
p1=pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),color = colorRampPalette(c( "blue","white", "red"))(256))
dpi=300
png(file="1d.png", width = dpi*4, height = dpi*4, units = "px",res = dpi,type='cairo')
p1
dev.off()


#fig 1e
features=c("Spry2","Actb","Marf1","Nfia","Bgn","Cd9","Serpinh1","Col8a1","Ly6a","Cd248","Ly6c1","Pi16","Apoe","Lpl","G0s2","Socs3","Cilp","Postn","Ctgf","Meox1")
mm$celltype=Idents(mm)
features=features
mm@meta.data$CB=rownames(mm@meta.data) 
bubble.df=as.matrix(mm[["RNA"]]@data[features,])
bubble.df=t(bubble.df)
bubble.df=as.data.frame(scale(bubble.df))
bubble.df$CB=rownames(bubble.df)
bubble.df=merge(bubble.df,mm@meta.data[,c("CB","celltype")],by = "CB")
bubble.df$CB=NULL
celltype_v=c()
gene_v=c()
mean_v=c()
ratio_v=c()
for (i in unique(bubble.df$celltype)) {
  bubble.df_small=bubble.df%>%filter(celltype==i)
  for (j in features) {
    exp_mean=mean(bubble.df_small[,j])
    exp_ratio=sum(bubble.df_small[,j] > min(bubble.df_small[,j])) / length(bubble.df_small[,j])
    celltype_v=append(celltype_v,i)
    gene_v=append(gene_v,j)
    mean_v=append(mean_v,exp_mean)
    ratio_v=append(ratio_v,exp_ratio)
  }
}
plotdf=data.frame(
  celltype=celltype_v,
  gene=gene_v,
  exp=mean_v,
  ratio=ratio_v
)
plotdf$celltype=factor(plotdf$celltype,levels = sort(unique(plotdf$celltype)))
plotdf$gene=factor(plotdf$gene,levels = rev(as.character(features)))
plotdf$exp=ifelse(plotdf$exp>1.5,1.5,plotdf$exp)
p=plotdf%>%ggplot(aes(x=celltype,y=gene,size=ratio,color=exp))+geom_point()+
  scale_x_discrete("")+scale_y_discrete("")+
  scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[8:11])))+
  scale_size_continuous(limits = c(0,1))+theme_bw()+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )
dpi=300
png(file="1e.png", width = dpi*5, height = dpi*4, units = "px",res = dpi,type='cairo')
p
dev.off()

#fig 1f
mm=RenameIdents(mm,"2"="C1")#.......
p=DimPlot(mm,cols =pal,label =F,repel = T,
          label.size = 8,split.by = 'condition',pt.size = 0.8)+umap_theme+NoLegend()

png(file="1f.png", width = dpi*15, height = dpi*5, units = "px",res = dpi,type='cairo')
p
dev.off()
p=DimPlot(mm,split.by = 'condition')
png(file="1f.png", width = dpi*5, height = dpi*4, units = "px",res = dpi,type='cairo')
p
dev.off()

#fig 1g
p=FeaturePlot(mm,"Apoe",cols = viridis(256))+umap_theme+
  theme(plot.title = element_text(hjust = 0.5,size=30))+NoLegend()
dpi=300
png(file="1g-6.png", width = dpi*6, height = dpi*6, units = "px",res = dpi,type='cairo')
p
dev.off()


#fig 1h
mm$celltype=Idents(mm)
phe =mm@meta.data
metadata <- phe
meta.temp <- metadata[,c("celltype", "condition")]
prop.table.error <- list()
for(i in 1:length(unique(meta.temp$condition))){
  vec.temp <- meta.temp[meta.temp$condition==unique(meta.temp$condition)[i],"celltype"]
  table.temp <- freqCI(vec.temp, level = c(.95))
  prop.table.error[[i]] <- print(table.temp, percent = TRUE, digits = 3)
}
names(prop.table.error) <- unique(meta.temp$condition)
tab.1 <- as.data.frame.array(do.call(rbind, prop.table.error))
b <- c()
a <- c()
for(i in names(prop.table.error)){
  a <- rep(i,nrow(prop.table.error[[1]]))
  b <- c(b,a)
}
tab.1$condition <- b
aa <- gsub(x = row.names(tab.1), ".1", "")
aa <- gsub(x = aa, ".2", "")
aa <- gsub(x = aa, ".3", "")
aa <- gsub(x = aa, ".4", "")
tab.1$cell <- aa  
colnames(tab.1)[1] <- "lower"
colnames(tab.1)[3] <- "upper"
p<- ggplot(tab.1, aes(x=condition, y=Estimate, group=cell)) +
  geom_bar(stat = "identity", aes(fill=cell)) + facet_grid(cols =  vars(cell)) + 
  theme_bw() + scale_fill_manual(values =pal)+ 
  theme(axis.text = element_text(angle = 45, hjust=1, vjust=0.5,size = 20), legend.position= "none") + 
  xlab("") + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,position=position_dodge(0.05)) 

png(file="1e.png", width = dpi*10, height = dpi*8, units = "px",res = dpi,type='cairo')
p
dev.off()




---------------------------------------------------------------------------------------------------------------------------------------
#fig 2a-d
object.markers <- FindMarkers(mm, ident.1 = 'Sham',logfc.threshold = 0.2,min.pct = 0,pseudocount.use=0.01)
object.markers$names <- rownames(object.markers)
object.markers$group=0
for (i in 1:nrow(object.markers)) {
  if (object.markers$avg_log2FC[i] >= 0.5 & object.markers$Difference[i] >= 0.1 & object.markers$pct.2[i]<=0.5){
    object.markers$group[i]='up'
  }
  else if(object.markers$avg_log2FC[i] <= -0.5 & object.markers$Difference[i] <= -0.1 & object.markers$pct.1[i]<=0.5){
    object.markers$group[i]='down'
  }
  else {
    object.markers$group[i]='no'
  }
}
p=ggplot(object.markers, aes(x=Difference, y=avg_log2FC)) +
  geom_point(size=2,aes(color=group)) +
  scale_color_manual(values=c('blue','grey','red'))+
  geom_label_repel(data=subset(object.markers, group !='no'), aes(label=names), segment.size =0.25,size=10)+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()+umap_theme+NoLegend()
dpi=300
png(file="2a-1.png", width = dpi*10, height = dpi*10, units = "px",res = dpi,type='cairo')
p 
dev.off()

#fig 2e
deg1=FindMarkers(object = mm, ident.1 = 'Sham',
                 min.pct = 0.01, logfc.threshold = 0)

deg2=FindMarkers(object = mm, ident.1 = 'TAC',
                 min.pct = 0.01, logfc.threshold = 0) 

deg3=FindMarkers(object = cc, ident.1 = 'TAC+JQ1',
                 min.pct = 0.01, logfc.threshold = 0)
deg4=FindMarkers(object = cc, ident.1 = 'TAC+JQ1 withdrawn',
                 min.pct = 0.01, logfc.threshold = 0)

geneList= deg4$avg_log2FC 
names(geneList)= toupper(rownames(deg4))
geneList=sort(geneList,decreasing = T)
head(geneList)
gmtfile ='c2.cp.kegg.v7.5.1.symbols.gmt'
gmtfile ='c5.go.bp.v7.5.1.symbols.gmt'
library(GSEABase) 
geneset <- read.gmt( gmtfile ) 
geneset$term <- gsub('KEGG_','',geneset$term)
geneset$term <- tolower(geneset$term)
geneset$term <- gsub('_',' ',geneset$term)


length(unique(geneset$term))
egmt <- GSEA(geneList, TERM2GENE=geneset, 
             minGSSize = 1,
             pvalueCutoff = 0.99,
             verbose=FALSE)
head(egmt)
egmt@result 
gsea_results_df <- egmt@result 
rownames(gsea_results_df)
write.csv(gsea_results_df,file = 'kegg_deg4.csv')

library(fmsb)
data=read.table('go-merge.txt')
colnames(data)=c('oxidative phosphorylation',	'supramolecular fiber organization','atp metabolic process',
                 'tissue remodeling',	'negative regulation of response to oxidative stress','cell cell adhesion','myeloid leukocyte activation')
rownames(data) <- c('1','2','Sham','TAC','TAC+JQ1',"TAC+JQ1 withdrawn")
colors_border=c(rgb(0.8549,0.64706,0.12549,1),rgb(1,0.49804,0.31373,1), rgb(0.81569,0.12549,0.56471,1),rgb(0.56723,0.198745,0.3,1))
colors_in=c(rgb(0.8549,0.64706,0.12549,0.5),rgb(1,0.49804,0.31373,0.5), rgb(0.81569,0.12549,0.56471,0.5),rgb(0.56723,0.198745,0.3,0.5))

radarchart( data , axistype=0,
            pcol=colors_border , pfcol=colors_in , plwd=3 , plty=1,pty=32,
            cglcol="black", cglty=2, cglwd=0.6)

dpi=300
png(file="2b-2.png", width = dpi*8, height = dpi*8, units = "px",res = dpi,type='cairo')
radarchart( data , axistype=0,
            pcol=colors_border , pfcol=colors_in , plwd=3 , plty=1,pty=32,
            cglcol="black", cglty=2, cglwd=0.6)
dev.off()

#fig 2f
GO_DATA <- get_GO_data("org.Mm.eg.db", "ALL", "ENTREZID")

findGO <- function(pattern, method = "key"){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  if(method == "key"){
    pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
  } else if(method == "gene"){
    pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
  }
  
  colnames(pathways) = "pathway"
  
  if(length(pathways) == 0){
    cat("No results!\n")
  } else{
    return(pathways)
  }
}

getGO <- function(ID){
  
  if(!exists("GO_DATA"))
    load("GO_DATA.RData")
  allNAME = names(GO_DATA$PATHID2EXTID)
  if(ID %in% allNAME){
    geneSet = GO_DATA$PATHID2EXTID[ID]
    names(geneSet) = GO_DATA$PATHID2NAME[ID]
    return(geneSet)     
  } else{
    cat("No results!\n")
  }
}


findGO("fibrosis")
getGO("GO:0030198")
a=getGO("GO:0030198")
a=as.data.frame(a)
a=bitr(a$extracellular.matrix.organization,fromType = "ENTREZID",toType = c("SYMBOL"),OrgDb = "org.Mm.eg.db")

features=a$SYMBOL


mat<- GetAssayData(mm, slot = "counts")
mat=log10(mat+1)
gene_features <-features
cluster_info <- sort(mm$condition)
mat <- as.matrix(mat[gene_features, names(cluster_info)])

devtools::install_github("caleblareau/BuenColors")
library("BuenColors")
col <- jdb_color_maps[1:9]
names(col) <- levels(cluster_info)

top_anno <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = col),
                       labels = levels(cluster_info), 
                       labels_gp = gpar(cex = 0.5, col = "black"))) 

mark_gene <- c("Hif1a","Hmox1","Sod1","Nfe2l2","Nfe2l1","Nos3","Sod2","Sod3","Gsta3","Gpx3","Gstm1")
gene_pos <- which(rownames(mat) %in% mark_gene)
row_anno <-  rowAnnotation(mark_gene = anno_mark(at = gene_pos,labels = mark_gene))
p=Heatmap(mat,
          cluster_rows = T,
          cluster_columns = F,
          show_column_names = F,
          show_row_names = F,
          column_split = cluster_info,
          top_annotation = top_anno,
          right_annotation = row_anno,
          column_title = NULL)

dpi=300
png(file="2f.png", width = dpi*6, height = dpi*4, units = "px",res = dpi,type='cairo')
p
dev.off()


#fig 2g
a=mm$condition
dorothea_regulon_mouse <- get(data("dorothea_mm", package = "dorothea"))
regulon <- dorothea_regulon_mouse %>%
  dplyr::filter(confidence %in% c("A","B","C"))
mm <- run_viper(mm, regulon,
                options = list(method = "scale", minsize = 4, 
                               eset.filter = FALSE, cores = 1, 
                               verbose = T))
DefaultAssay(mm) <- "dorothea"
mm <- ScaleData(mm)
mm <- RunPCA(mm, features = rownames(mm), verbose = FALSE)

pct <- mm [["pca"]]@stdev / sum( mm [["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
mm<- RunUMAP(mm, dims = 1:8)
mm <- FindNeighbors(mm, dims = 1:8, verbose = FALSE)
mm <- FindClusters(mm, resolution = 0.5, verbose = FALSE)

Idents(mm)=a

viper_scores_df <- GetAssayData(mm, slot = "scale.data", 
                                assay = "dorothea") %>%
  data.frame(check.names = F) %>%
  t()


CellsClusters <- data.frame(cell = names(Idents(mm)), 
                            cell_type = as.character(Idents(mm)),
                            check.names = F) #也可以使用其他的分类信息

viper_scores_clusters <- viper_scores_df  %>%
  data.frame() %>% 
  rownames_to_column("cell") %>%
  gather(tf, activity, -cell) %>%
  inner_join(CellsClusters)


summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))


highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(100, var) %>%
  distinct(tf)

## We prepare the data for the plot
summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE) 
palette_length = 100
my_color = colorRampPalette(c("#00FFFF", "white","red"))(palette_length)

my_breaks <- c(seq(min(summarized_viper_scores_df), 0, 
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(summarized_viper_scores_df)/palette_length, 
                   max(summarized_viper_scores_df), 
                   length.out=floor(palette_length/2)))

viper_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=15, 
                       fontsize_row = 15, 
                       color=my_color, breaks = my_breaks, 
                       main = "TF activity", 
                       treceheight_col = 0,  border_color = "#000000",angle_col = 45) 

dpi=300
png(file="2f.png", width = dpi*4, height = dpi*8, units = "px",res = dpi,type='cairo')
viper_hmap
dev.off()

---------------------------------------------------------------------------------------------------------------------------------
#fig 3a-d
mm$celltype=Idents(mm)
table(mm$celltype)
data <- as(as.matrix(mm@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = mm@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
HSMM<-monocle_cds
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.10 )
expressed_genes <- row.names(subset(fData(HSMM),
                                    num_cells_expressed >= 100))
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
ordering_genes=disp.genes
diff <- differentialGeneTest(HSMM[expressed_genes,],fullModelFormulaStr="~celltype",cores=2) 
ordering_genes <- row.names (subset(diff, qval <1e-10))
ordering_genes <-VariableFeatures(mm)
saveRDS(ordering_genes,"ordering_genes.rds")
HSMM <- setOrderingFilter(HSMM,ordering_genes)
plot_ordering_genes(HSMM)
HSMM<- reduceDimension(HSMM, max_components = 2,method = 'DDRTree')
HSMM <- orderCells(HSMM)

a=plot_cell_trajectory(HSMM, color_by = "celltype") + 
  scale_color_manual(breaks = c("X", "Y", "Z"), values=pal)+umap_theme
b=plot_cell_trajectory(HSMM, color_by = "Pseudotime")+umap_theme+NoLegend()
c=DimPlot(mm,label=T,cols=pal,repel = T,label.size = 8)+umap_theme+NoLegend()
d=plot_cell_trajectory(HSMM, color_by = "State", show_state_number = F, show_tree = T) + 
  facet_wrap(~condition,nrow = 1) + scale_color_manual(breaks = c("X", "Y", "Z"), values=pal1) +
  theme(axis.line = element_line(size = 6), axis.text.x = element_text(size=16, face = "italic"),
        axis.text.y = element_text(size=18),
        axis.title = element_text(size=18), legend.position = "top", strip.text = element_text(size=18))
plot_grid(a,b,c,d)

dpi=300
png(file="3d.png", width = dpi*20, height = dpi*6, units = "px",res = dpi,type='cairo')
d
dev.off()

#fig 3d(propotion)
anno <- c("Sham","TAC","TAC+JQ1","TAC+JQ1 withdrawn")
Data=mm
Data.Monocle=HSMM
States <- Data.Monocle$State
names(States) <- rownames(Data.Monocle@phenoData@data)
Data <- AddMetaData(Data, States, "State")
A <- data.frame(Timepoint = Data@meta.data$condition, State = Data@meta.data$State)
a <- as.numeric(table(A$Timepoint))
x <- aggregate(A, by = list(A$Timepoint, A$State), FUN = length)
x <- A%>%count(Timepoint,State)%>% complete(Timepoint, State, fill = list(n = 0))
prop =NULL
for (j in anno){
  y <- filter(x, Timepoint == j)
  s <- sum(y$n)
  pobs<- c(y$n/s)
  prop = c(prop, pobs)
}
x <- data.frame(x, Prop= prop)
p1=ggplot(x, aes(x=factor(1), y= prop, fill = State))+ 
  geom_bar(color = "black", stat = "identity", width = 1)+
  coord_polar(theta = "y", clip = "on")+scale_fill_manual(values =pal1 )+
  scale_x_discrete(breaks = "", labels = "")+
  facet_grid(~Timepoint)+theme_bw()+
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank())
p1
dpi=300
png(file="3d-prop.png", width = dpi*12, height = dpi*4, units = "px",res = dpi,type='cairo')
p1
dev.off()







#fig 3e

BEAM_res <- BEAM(HSMM[ordering_genes,], branch_point = 2, cores = 2)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

tmp1=plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                                       qval < 1e-50)),],
                                 branch_point =2,
                                 num_clusters =3, 
                                 cores = 1,
                                 branch_labels = c("Cell fate 1", "Cell fate 2"),
                                 hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                 branch_colors = c("#979797", "#F05662", "#7990C8"), 
                                 use_gene_short_name = F,
                                 show_rownames = F,
                                 return_heatmap = T 

png(file="3e.png", width = dpi*4, height = dpi*5, units = "px",res = dpi,type='cairo')
tmp1$ph_res
dev.off()
gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)
write.csv(gene_group,"gene_group_branch.csv")

allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}

write.csv(allcluster_go,"go_branch.csv")

#fig 3h-i

States <- HSMM@phenoData@data$State
names(States) <- rownames(HSMM@phenoData@data)
mm <- AddMetaData(mm, States, "State")

X <- data.frame(State = mm@meta.data$State, Postn = expm1(mm@assays$RNA@data["Postn", ]), Sparc = expm1(mm@assays$RNA@data["Sparc", ]),  Cilp= expm1(mm@assays$RNA@data["Cilp", ]))
X <- data.frame(State = mm@meta.data$State, Nrf2 = expm1(mm@assays$RNA@data["Nfe2l2", ]), Hmox1 = expm1(mm@assays$RNA@data["Hmox1", ]), Sod2 = expm1(mm@assays$RNA@data["Sod2", ]))
X.m <- melt(X)
X.m$State=factor(X.m$State,
                 levels=c( "branch 1", "branch 2","branch 3"))

p=ggplot(X.m, aes(x= State, y= value, fill = State))+
  geom_jitter(height = 0, color = "grey80")+
  geom_violin(scale = "width")+ 
  facet_wrap(~variable, scales = "free") +scale_fill_manual(values =pal1 )+
  scale_y_continuous(name = "normalized UMI counts (log10)", trans = "log1p")+
  theme(axis.line = element_line(size = 2), axis.text.x = element_text(size=20, angle=45, hjust=1, vjust=1),
        axis.text.y = element_text(size=12), 
        axis.title = element_text(size=20), legend.position = "top", strip.text = element_text(size=30))+NoLegend()

png(file="3h.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
p
dev.off()

-----------------------------------------------------------------------------------------------------
#fig 4b

s.genes <-c( "Gstm1","Gpx3","Ggt5","Gsta3")
s.genes <-c( "Actb","Actg1","Fn1","Anxa1")
s.genes <-c( "Ctgf","Col3a1","Postn","Sparc")


p=plot_genes_branched_pseudotime(HSMM[s.genes,],
                                 branch_point = 2,
                                 color_by = "State",
                                 ncol = 2)+ scale_color_manual(breaks = c("X", "Y", "Z"), values=pal1)+NoLegend()+theme(plot.title = element_text(size = 40))
p
png(file="4b-1.png", width = dpi*6, height = dpi*6, units = "px",res = dpi,type='cairo')
p
dev.off()


#fig 4c
Gene_list1 = c("Slc2a1", "Hk1", "Hk2", "Gpi1", "Pfkm", "Pfkp", "Aldoa", "Gapdh", "Eno1", "Eno3", "Pkm","Pgk1", 
              "Cs", "Aco2", "Idh1", "Idh2", "Idh3a", "Idh3b","Pdha1", "Fh1", "Mdh1", "Mdh2", "Ldha", "Ldhb","Sdha", "Sdhc", "Sucla2", "Pcx","Ogdh", "Acly",
              "Uqcrc1", "Ndufb8", "Ndufa11", "Ndufb4", "Ndufb7", "Uqcrb","Uqcrq", "Uqcrh", "Uqcrfs1", "Cox7a2", "Cox7a2l", "Cox7c","Cox8a",
              "Acat1","Acat2","Hmgcr","Acaca","Apoe","Lpl","Ldlr","Vldlr","Acsl4","Cpt1a","Acadm","Acadsb","Acadl","Lipe")

Gene_list2 = c("G6pdx","Pgls","Rpia","Rpe","Tkt","Taldo1","Aldoa",
              "Ugp2","Gys1","Gbe1","Gk",
              "Pygm","Pygb","Pygl","Agl","Dbr1","G6pc3",
              "Pcx","Pck2","Mdh1","Mdh2","Got1","Got2")

for (i in Gene_list){
  if(length(Data@assays$RNA@data[i,] == length(rownames(Data@meta.data)))){print(i)}
}
EndMT <-c("branch 1","branch 2","branch 3")
FC <- NULL
p.val <- NULL
Ti <- NULL
for (i in 1:length(EndMT)){
  print(EndMT[i])
  if(i == 1){
    p.val <- rep(NA, length(Gene_list))}
  else{
    Tester <- FindMarkers(Data, features = Gene_list, test.use = "bimod", logfc.threshold = 0, min.pct = 0, return.thresh = 1,
                          ident.2 ="branch 2",ident.1 = EndMT[i])
    Tester <- Tester[Gene_list, ]
    p.val <- cbind(p.val, Tester$p_val_adj)}
}
G <- NULL
Means <- NULL
for (i in EndMT){
  for(j in Gene_list){
    x <- mean(expm1(Data@assays$RNA@data[j, WhichCells(Data, cells = which(Data@meta.data$EndMT.counter == i))]))
    G <- c(G,x)
  }
  Means= cbind(Means, G)
  G <- NULL
}
rownames(Means) <- Gene_list
rownames(p.val) <- Gene_list
colnames(Means) <- EndMT
colnames(p.val) <- EndMT
#Means. <- Means
#for (i in 1:length(p.val)){
#if(p.val[i] > 0.05 & is.na(p.val[i]) == F){
#  Means[i] <- NA
# }
#}
colo <- c("blue", "grey80", "red")
my_palette <- colorRampPalette(colo)(255)
Pathway <- c(rep("Glycolysis", times = 12), rep("TCA cycle", times = 18),rep("Mitochondrial respiratory chain", times = 13),rep("Fat acid metabolism", times = 14))
J_dataframe <- data.frame(Pathway = Pathway)
rownames(J_dataframe) <- rownames(Means)
p=pheatmap::pheatmap(Means, scale = "row", cluster_rows = F, cluster_cols = F, 
                     color = my_palette,angle_col = 45, na_col = "white", border_color = "white", gaps_row = c(12,30,43), annotation_row = J_dataframe)

png(file="4c.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
p
dev.off()

---------------------------------------------------------------------------------------------------------
#fig 5(L-R)
  
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))

ligand_target_matrix = readRDS('ligand_target_matrix.rds')
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")

mm$celltype=Idents(mm)
mm@meta.data %>% head()
mm@meta.data$celltype %>% table()

sender_celltypes =c("branch 1" ,"branch 1" ,"branch 3")



nichenet_output = nichenet_seuratobj_aggregate(
  seurat_obj = mm, 
  receiver =c("Arterial EC","venous EC","Capillary EC") ,
  sender = sender_celltypes, condition_colname = "condition", condition_oi = "TAC+JQ1", condition_reference = "TAC",
  ligand_target_matrix = ligand_target_matrix, lr_network = lr_network, 
  weighted_networks = weighted_networks, organism = "mouse",assay_oi='RNA',
  top_n_ligands =20,filter_top_ligands=T )


p=DotPlot(mm, features = nichenet_output$top_receptors %>% rev(), cols = "RdYlBu") + RotatedAxis()
png(file="5e.png", width = dpi*12, height = dpi*4, units = "px",res = dpi,type='cairo')
p
dev.off()


p=DotPlot(mm, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
png(file="5f.png", width = dpi*8, height = dpi*6, units = "px",res = dpi,type='cairo')
p
dev.off()

p=nichenet_output$ligand_receptor_heatmap
png(file="5g.png", width = dpi*12,height = dpi*4, units = "px",res = dpi,type='cairo')
p
dev.off()


-------------------------------------------------------------------------------------------------------
#fig 6(a)

features=c("Jun","Fos","Egr1","Atf3","Junb","Fosb","Cebpb","Jund")
p=VlnPlot(mm,features,cols = pal1,ncol = 4)
png(file="6a.png", width = dpi*10, height = dpi*10, units = "px",res = dpi,type='cairo')
p
dev.off()

#fig 6b
p=VlnPlot(cc,features,cols = col,ncol = 4,group.by = 'condition')
png(file="6b.png", width = dpi*12, height = dpi*10, units = "px",res = dpi,type='cairo')
p
dev.off()


#fig 6c
p=plot_genes_branched_pseudotime(HSMM[features,],
                                 branch_point = 2,
                                 color_by = "State",
                                 ncol = 4)+ scale_color_manual(breaks = c("X", "Y", "Z"), values=pal1)+NoLegend()

png(file="6c.png", width = dpi*12, height = dpi*8, units = "px",res = dpi,type='cairo')
p
dev.off()








