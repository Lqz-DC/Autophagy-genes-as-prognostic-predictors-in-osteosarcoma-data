#Set the working directory
setwd("F:\\scRNA\\GSE162454_RAW")
####Obtain the list of folders. There are 6 sample folders placed under this path.####
folder_list <- list.dirs(path = "F:\\scRNA\\GSE162454_RAW", full.names = TRUE, recursive = FALSE)
####Traverse the folder####
for (folder in folder_list) {  file_list <- list.files(path = folder, full.names = TRUE, recursive = FALSE)     
for (file in file_list) {    
  extension <- tools::file_ext(file)        
  # Compressed file  
  if (extension %in% c("tsv","mtx")) { gz_file <- paste0(file,".gz")
  cmd <- paste("gzip", file)            
  # Run the system command to compress the file    
  system(cmd)            
  # Rename the compressed files     
  file.rename(paste0(file,".gz"),gz_file)    
  cat("Compressed file:",file,"\n")    
  } 
}
}
####01 Read the file and create a Seurat object####
library(Seurat)
library(harmony)
library(tidyverse)
library(dplyr)
library(patchwork)
library(tidydr)
library(ggplot2)
library(cowplot)
library(ggunchull)
library(ggrepel)
library(presto)
# Create an empty list of Seurat objects
seurat_list <- list()
# Loop through each sample folder
for (folder in folder_list) {  

  data <- Read10X(data.dir = folder)   
  # Convert to Seurat object
  seurat <- CreateSeuratObject(counts = data, project = basename(folder),min.cells = 3,min.features = 300)    
  # Add the Seurat object to the list 
  seurat_list[[basename(folder)]] <- seurat
}
seurat_list

for(i in 1:length(seurat_list)){
  sc <- seurat_list[[i]]  # Obtain the i-th Seurat object from the seurat_list and calculate the mitochondrial proportion.
  sc[["mt_percent"]] <- PercentageFeatureSet(sc, pattern = "^MT-")  # Calculate the proportion of genes starting with "MT-".
  
  # Assign the value of sc to seurat_list[[i]] 
  seurat_list[[i]] <- sc
  
  # Delete "sc"
  rm(sc)
}

####02 Batch generation of pre-quality control violin plots#### 
violin_before <- list()
for(i in 1:length(seurat_list)){
  violin_before[[i]] <- VlnPlot(seurat_list[[i]],
                                layer ="counts",
                                features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),
                                pt.size = 0.01, 
                                ncol=4)
}
violin_before## Output the violin plot before quality control
violin_before[[2]]#Click on the number to view the data of a certain sample.

####03 Batch filtration of cells and MT####
seurat_list <- lapply(seurat_list, function(x) {
  x <- subset(x,
              subset = nFeature_RNA > 500 & 
                nFeature_RNA < 6500 & 
                mt_percent < 10 & 
                nCount_RNA <30000& 
                nCount_RNA > 1000)
  return(x)
})
view(seurat_list[[1]]@meta.data)

##### 04 Merge multiple Seurat objects####
seurat_merge<-merge(x=seurat_list[[1]],y=seurat_list[-1], add.cell.ids = names(seurat_list))

#The data after quality control (violin pick)
VlnPlot(seurat_merge,	
        features = c("nFeature_RNA","nCount_RNA","mt_percent","HB_percent"),	
        split.by = "orig.ident",
        layer = "counts",
        pt.size = 0.01, 
        ncol=4)
## Count the number of cells
table(seurat_merge[[]]$orig.ident)

#### 05 Standardization & Variable-variant Genes & Normalization & PCA####
seurat_merge <- NormalizeData(seurat_merge)
seurat_merge <- FindVariableFeatures(seurat_merge)
seurat_merge <- ScaleData(seurat_merge) #Eliminate the influence of mitochondria
seurat_merge <- RunPCA(seurat_merge,verbose=F)

###Integrated with harmony
scRNA_harmony <- IntegrateLayers(object =seurat_merge,
                                 method = HarmonyIntegration, 
                                 orig.reduction ="pca",
                                 new.reduction ="harmony", 
                                 verbose =FALSE)

#### 06 JoinLevers combines sample counts and data####
scRNA_harmony[["RNA"]] <- JoinLayers(scRNA_harmony[["RNA"]])
save(scRNA_harmony,file = "scRNA_harmony.Rdata")
load("scRNA_harmony.Rdata")	

#The cluster tree determines the appropriate number of clusters
library(clustree)	
clustree(scRNA_harmony)

####Dimensionality Reduction Clustering Umap Plot####	
ElbowPlot(scRNA_harmony,ndims = 50)	
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30)	

scRNA_harmony <- FindClusters(scRNA_harmony,resolution = 0.2)
scRNA_harmony <- RunUMAP(scRNA_harmony,dims = 1:30,reduction = "harmony")
scRNA_harmony <- RunTSNE(scRNA_harmony,dims = 1:30,reduction ="harmony")	

#Save the data after dimensionality reduction and clustering
save(scRNA_harmony,file = "scRNA_harmony_resulition0.2.Rdata")

#Load data
scRNA_harmony <- load("scRNA_harmony_resulition0.2.Rdata")	


#Look at the UMAP graph to see the integration status of Harmony.	
DimPlot(scRNA_harmony, reduction = "umap", group.by = "orig.ident")+
  ggtitle("Harmony")

scRNA_harmony$RNA_snn_res.0.2

Idents(scRNA_harmony) <-"RNA_snn_res.0.2"#Set the default resolution to 0.2

DimPlot(scRNA_harmony, reduction = "umap",label =TRUE)

# Drawing
umap_integrated1 <- DimPlot(scRNA_harmony, reduction = "umap",group.by = "orig.ident")
tsne_integrated2<- DimPlot(scRNA_harmony, reduction ="tsne",label =TRUE) 
umap_integrated3 <- DimPlot(scRNA_harmony, reduction = "umap",label=TRUE)

save(scRNA_harmony,file = "scRNA_harmony.Rdata")
load("scRNA_harmony.Rdata")

#### Differential gene annotation####

options(stringsASFactors=F)
setwd("F:\\scRNA\\GSE162454_RAW")
load("scRNA_harmony.Rdata")

rm(list = ls())
load("celltype.scRNA.Rdata")
####Part One: Difference Analysis of FindAllMarkers
markers <- FindAllMarkers(object = scRNA_harmony, test.use="wilcox",
                          only.pos =TRUE,
                          logfc.threshold=0.25)

#Screen the maker genes of each calculated cluster
all.markers=markers%>%dplyr::select(gene,everything())%>%subset(p_val<0.05)
top10 <- all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
view(top10)
write.csv(top10,"14cluster_top10.csv",row.names = T)
write.csv(all.markers,"allmarkers.csv",row.names = T)

DefaultAssay(scRNA_harmony)="RNA"
scRNA_harmony <- FindClusters(scRNA_harmony,resolution = 0.2)
#Create a marker set based on the literature
marker <- c('ALPL','RUNX2', 'IBSP',#osteoblastic OS cells 
            'LYZ','CD68',#myeloid cells
              'ACP5', 'CTSK',# osteoclastic cells 
              'FBLN1','ACTA2','TAGLN','COL3A1', 'COL6A1',#fibroblasts (CAFs)
            'CD2', 'CD3D','CD3E', 'CD3G','GNLY', 'NKG7', 'KLRD1', 'KLRB1',#NK/T cells
            'EGFL7', 'PLVAP',#endothelial cells
              'MS4A1', 'CD79A',#B cells 
            'IGHG1', 'MZB1')#plasma cells

par(mar = c(5, 4, 8, 3))
p<- DotPlot(scRNA_harmony,features = marker,group.by = "celltype")+RotatedAxis()+ 
  coord_flip()
#Add comments to the cluster
scRNA_harmony$celltype <- recode(scRNA_harmony@meta.data$seurat_clusters,
                                 "0"="myeloid cells","1"="NK/T cells","2"="osteoblastic OS cells",
                                 "3"="myeloid cells","4"="myeloid cells","5"="osteoblastic OS cells",
                                 "6"="fibroblasts (CAFs)","7"="osteoblastic OS cells","8"="plasma cells",
                                 "9"="endothelial cells","10"="osteoclastic cells","11"="B cells","12"="myeloid cells","13"="osteoblastic OS cells")
                              
table(scRNA_harmony@meta.data$celltype)


# Define the sequence of cell types (keeping it consistent with the color vector)
celltype_levels <- c("osteoblastic OS cells","myeloid cells","NK/T cells",
                     "fibroblasts (CAFs)","plasma cells","endothelial cells",
                     "osteoclastic cells","B cells")

# Convert the "celltype" column to a factor and sort it.
scRNA_harmony$celltype <- factor(scRNA_harmony$celltype, 
                                 levels = celltype_levels,
                                 ordered = TRUE)

# Use the custom colors you provided (the number of colors should match the number of cell types)
custom_colors <- c("#76b7b2","#e15759","#ff9da7","#59a14f","#edc949","#af7aa1","#ff9f40","#9c755f")
names(custom_colors) <- celltype_levels

meta <- cbind(scRNA_harmony@meta.data, 
              scRNA_harmony@reductions$umap@cell.embeddings)

# Calculate the position of the label
celltype_med <- meta %>% 
  group_by(celltype) %>% 
  summarise(umap_1 = median(umap_1),
            umap_2 = median(umap_2))

# Create  UMAP diagram
p_umap <- ggplot(meta, aes(x = umap_1, y = umap_2)) +
  stat_unchull(aes(fill = celltype,color=celltype), 
               alpha = 0.05, size =1, lty = 2, delta = 0.25,show.legend = F) +
  geom_point(aes(color = celltype), size = 0.3, alpha = 0.8) +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  geom_text_repel(
    aes(label = celltype, color = "black"),
    data = celltype_med,
    fontface = "bold", 
    size = 4,
    box.padding = 0.8,
    max.overlaps = Inf
  ) +
  theme_classic(base_size = 14) +
  theme(
    aspect.ratio = 1,
    axis.line = element_line(arrow = arrow(type = "closed", length = unit(2, "mm"))),
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(x = "UMAP1", y = "UMAP2")
p_umap
####DotPlot####
p_dot <- DotPlot(scRNA_harmony, 
                 features = marker, 
                 group.by = "celltype",
                 cols = c("lightskyblue", "#b2182b")) + # 双色渐变
  RotatedAxis() +
  scale_y_discrete(limits = rev(levels(scRNA_harmony$celltype))) + # 保持顺序一致
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_dot

####Model gene distribution expression map####
modelgene <- c("PEA15","MYC","SAR1A","BNIP3")
p5 <- FeaturePlot(scRNA_harmony,features =modelgene,ncol =2,cols = c("lightskyblue", "#b2182b"))
p5

p6 <- VlnPlot(scRNA_harmony,
        layer ="counts",
        features = modelgene,
        #pt.size = 0, 
        group.by = "celltype",
        ncol=4)+theme(axis.text.x = element_text(angle = 45, hjust = 1))

