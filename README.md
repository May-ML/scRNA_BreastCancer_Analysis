Name: Meiheng Liang
Programming Language: [Python]
Date: [Dec 2023]
Description:
This script is a demonstration the application of scRNA analysis using seurat pipeline in analyzing cellular composition of invasive ductal breast cancer tissue. The data which available from https://10xgenomics.com

############################################################################
Required files:

filtered_feature_bc_matrix

############################################################################
Required packages:
Tidyverse
patchwork
Seurat

############################################################################

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)


##Load libraries
  library(dplyr)
  library(Seurat)
  library(patchwork)
```

## R Markdown
```{r dataset prep}
  #Check if the file is in your working directory
  dir()
  #Change your working directory to location of files
  setwd("path/to/your/directory")
  #Load in 10x Genomics data
  bre_c.data <- Read10X("filtered_feature_bc_matrix")
  #Create RNA object, where the counts
  bc <- CreateSeuratObject(counts = bre_c.data, project = "nameyourproject", min.cells = 3, min.features = 200) #only genes expressed within minimal of 3 cells and cells with minimal 200 gene expression
  #brain # Features, samples, 1 assay
  head(bc) # We can look at different column data
```

```{r Data  preprocessing}
head(bc) # We can look at different column data, notice theres now percent.mt column
  VlnPlot(bc, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3) # Display these column types 
  
  # The mitochondria dna lead to a lot of noise
  bc <- subset(bc, subset = nFeature_RNA < 3000 & nCount_RNA < 30000 & percent.MT < 15)
  VlnPlot(bc, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), ncol = 3)

## Data normalization
  # Highly variable gene shown pre-normalization
  VlnPlot(bc,features="Hbb-bs") # Try Inpp5d as well, a marker for microglia
  # Normalization
  bc <- NormalizeData(bc, normalization.method = "LogNormalize", scale.factor = 10000) # Scale factor is the magnitude of your normalization. 10000 is recommend by Satija
  # Highly variable gene shown pre-normalization
  VlnPlot(bc,features="Hbb-bs") # Look at the difference after normalization
  
## Identification key variable 
  # adjust thereshold as necessary
  bc <- FindVariableFeatures(bc, selection.method = "vst", nfeatures = 2000) # Identifies features that are outliers on a 'mean variability plot'.
  
  # The code below is just for fun, not really required for analysis.
  top10 <- head(VariableFeatures(bc), 10) # Find the top ten variable genes, display them
  top10
  plot1 <- VariableFeaturePlot(bc)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = FALSE)
  plot1 + plot2
  plot1
  plot2
  
## Scaling the data  #scaling for 10,000 here
  all.genes <- rownames(bc)
  bc <- ScaleData(bc, features = all.genes)
```


```{r Perform linear dimensional reduction}
  # principal component 1,2,3 decenting importance in impacting data
  bc <- RunPCA(bc, features = VariableFeatures(object = bc))
  VizDimLoadings(bc, dims = 1:2, reduction = "pca")
  DimPlot(bc, reduction = "pca")
  
  # There is clear association though in the first 7 or 8 PC
  # Question: how do we make sure that everything is going correctly
  DimHeatmap(bc, dims = 1, cells = 500, balanced = TRUE)
  DimHeatmap(bc, dims = 1:15, cells = 500, balanced = TRUE) ##principal component determines good separation between upregulated or downregulated genes, if the separation of gene expression in heatmap is not clear, wrong with principal component
  ##ideally, heat map contain top expression gene and least expression gene in control versus disease, determine a separation point/cutoff point for PCA based on graph visualization
  
```


```{r Determine dimensionality}
  bc <- JackStraw(bc, num.replicate = 100)
  bc <- ScoreJackStraw(bc, dims = 1:20)
  ElbowPlot(bc,ndims = 50)   ##becaz looking for highly variated genes, when the elbow tilting toward 0 means less vairate and less meaningful

## Cells Clustering
  # Use elbow plot to determine how many dimensions to look at, 
  # Here we can go to 1:25 but even 1:40 works
  bc <- FindNeighbors(bc, dims = 1:25) #
  bc <- FindClusters(bc, resolution = 0.5) #resolution adjust according to cell group, higher the resolution smaller number of cell in scope, more detail for each
                                                  # Think of resolutions like microscopes. You can use different resolutions
                                                 # to change how "zoomed in" you want to be

  #bc@active.ident # Looking at the current amount of clusters we are actually using
  #Idents(bc) <- "RNA_snn_res.0.6" # Can switch between any of the resolutions with this command
  ##returns 0-23, 0: biggest cluster, and 23 the smallest
## Non-linear dimensional reduction (UMAP/tSNE)
  bc <- RunUMAP(bc, dims = 1:25)  # 
  DimPlot(bc, reduction = "umap", label = TRUE)
  # Example below on looking at one specific gene, Hbb-bs
  FeaturePlot(bc, features = c("Hbb-bs")) # Let's look at Inpp5d
  
##Save your R file to save computation time in the future
## saveRDS(bc,file = "file_name.rds")

```

```{r Identify differentially expressed features and Annotation}
  # Get all genes that are positive + assign markers
  bc.markers <- FindAllMarkers(bc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  head(bc.markers)
  # How to see the 1st, second, third cluster?
  head(subset(bc.markers, cluster =="0"))
  # Then go to https://panglaodb.se/samples.html and manually look up genes, 
  # Look at each cluster, significant genes (p-value), 
  # make 25 different names (0:24)
  # Rename clusters ID  - step for annotation
  new.cluster.ids <- c("", "", "", "", "", " ",
                       "", "", "", "", "",
                       "", "", "", "", "", "",
                       "", "", "", "", "", "",
                       "")
  names(new.cluster.ids) <- levels(bc)
  bc <- RenameIdents(bc, new.cluster.ids)
  # Plot the Umap but with the annotated clusters
  DimPlot(bc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
  


############################################################################
# Output files: 
.html
