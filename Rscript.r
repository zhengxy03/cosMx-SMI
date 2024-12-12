library(Seurat)
setwd("E:/project/ESCC/data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images")

data <- ReadNanostring(data.dir = ".", type = "centroids")
cents <- CreateCentroids(data$centroids)
coords <- CreateFOV(coords = list("centroids" = cents), type = "centroids")
nano.obj <- CreateSeuratObject(counts = data$matrix)
nano.obj[["fov"]]<-subset(coords$boundaries, cell=Cells(nano.obj))
setwd("E:/project/ESCC/data")
azimuth.data <- readRDS("../../ref/nanostring_data.Rds")
nano.obj <- AddMetaData(nano.obj, metadata = azimuth.data$annotations)
nano.obj[["proj.umap"]] <- azimuth.data$umap
Idents(nano.obj) <- nano.obj$predicted.annotation.l1
options(future.globals.maxSize = 8000 * 1024^2)
nano.obj <- SCTransform(nano.obj, assay = "Nanostring", clip.range = c(-10, 10), verbose = FALSE)

# text display of annotations and prediction scores
head(slot(object = nano.obj, name = "meta.data")[2:5])
ImageDimPlot(nano.obj, fov = "lung5.rep1", axes = TRUE, cols = "glasbey")