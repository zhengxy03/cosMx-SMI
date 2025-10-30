# cosMx-SMI使用手册
CosMx SMI是首个可对福尔马林固定石蜡包埋（FFPE）和新鲜冷冻（FF）组织等样本进行细胞和亚细胞分辨率空间多组学分析的高通量原位分析平台。SMI可对多达6000种RNA和64种经过验证的蛋白质实现快速定量和可视化分析。CosMx SMI是一个灵活的空间单细胞成像平台，能推动对细胞图谱、组织表型、细胞间相互作用和生物标志物发现的进一步发展。
# 1 介绍
## 1.1 技术原理
CosMx SMI是一种基于杂交的单分子条形码检测的、无酶、无核酸扩增的分析方法。
* ISH Target Probe: CosMx panel中的ISH 靶向探针（包含两个部分）：  
    * ①Target Binding Domain：通过杂交与内源mRNA结合；  
    * ②Reporter Readout Domain：为4个连续的reporter结合区域，每个区域仅与唯一对应的reporter结合；<br>
* Reporter Probe：与ISH靶向探针中reporter结合区域杂交配对的荧光报告探针；<br>
* Reporter Set：每轮杂交添加4种荧光探针。<br>

![技术原理](./pic/技术原理.png "技术原理")<br>
## 1.2 工作流程
* 固定并透化5μm FFPE或新鲜冷冻组织切片；
* RNA探针与组织样本中的内源RNA杂交；
* 清洗组织样本，然后与寡核苷酸标记的抗体孵育，进行形态学标志物染色；
* 上样至SMI仪器进行形态学标志物成像；
* 选择组织上所需的成像区域FOV(最大300mm2)；
* 多轮报告子结合和荧光成像，读取每个成像的RNA探针或蛋白质抗体的荧光条形码信息。<br>

![工作流程](./pic/工作流程.png "工作流程")<br>
# 2 上游实验与处理
## 2.1 样本准备
## 2.2 CosMx SMI 实验流程
* 探针杂交
* RNA 读取与成像
* 细胞分割
* 质量控制
# 3 下游数据处理
## 3.1 数据下载
这里使用NSCLS公开数据集中的一个重复进行示例：
```
curl -LO https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/Lung5_Rep1/Lung5_Rep1+SMI+Flat+data.tar.gz

tar xvfz Lung5_Rep1+SMI+Flat+data.tar.gz
```
解压之后可以看到这些文件：<br>
![数据下载](./pic/数据下载.png "数据下载")<br>
文件里的具体内容可以看[这里](https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/SMI-ReadMe.html)<br>

使用Seurat创建对象进行后续分析：
```
library(Seurat)
setwd("E:/project/ESCC/data")
nano.obj <- LoadNanostring(data.dir = "Lung5_Rep1+SMI+Flat+data/Lung5_Rep1/Lung5_Rep1-Flat_files_and_images", fov = "lung5.rep1")
```
* fov 是指 "field of view" 的缩写，它在单细胞基因表达数据的上下文中通常用来指代特定的视野或区域。在空间转录组学或类似的技术中，组织样本被分割成多个视野或区域，每个区域可以单独分析，以获得该区域的基因表达数据。
## 3.2 细胞注释
对于这个数据集，Seurat v5官网上并没有进行无监督分析，而是将Nanostring的分析结果与Azimuth健康人类肺脏参考数据库进行对比，这个数据库是通过单细胞RNA测序（scRNA-seq）技术建立的。使用的是Azimuth软件的0.4.3版本以及人类肺脏参考数据库的1.0.0版本。可以从[这个链接](https://seurat.nygenome.org/vignette_data/spatial_vignette_2/nanostring_data.Rds)下载预先计算好的分析结果，这些结果包括了注释信息、预测分数以及UMAP的可视化图。每个细胞平均检测到的转录本数量是249。<br>

但是这个数据库是来自健康样本，这里分析的组织是来自肿瘤样本，因此这种注释可能不准确，我后面又进行了常规无监督聚类分析及手动细胞注释。
```
azimuth.data <- readRDS("../nanostring_data.Rds")
nano.obj <- AddMetaData(nano.obj, metadata = azimuth.data$annotations)
nano.obj[["proj.umap"]] <- azimuth.data$umap
Idents(nano.obj) <- nano.obj$predicted.annotation.l1
```
## 3.3 预处理与标准化
SCTransform函数是一种数据标准化和变量选择的方法，它通过正则化和缩放数据来减少技术变异，并提高数据的可比性。
```
# set to avoid error exceeding max allowed size of globals
options(future.globals.maxSize = 8000 * 1024^2)
nano.obj <- SCTransform(nano.obj, assay = "Nanostring", clip.range = c(-10, 10), verbose = FALSE)

# text display of annotations and prediction scores
head(slot(object = nano.obj, name = "meta.data")[2:5])
```
* 在处理大型单细胞RNA-seq数据时，内存和计算资源是一个重要的考虑因素。SCTransform 函数尤其需要大量的内存，因为它涉及到复杂的数值计算和矩阵操作。确保在继续之前解决这些潜在的问题。
## 3.4 UMAP降维
我们可以可视化 Nanostring 细胞和注释，并使用UMAP降维。请注意，对于此 NSCLC 样本，肿瘤样本被注释为“基础”，这是健康参考中最接近的细胞类型匹配。
```
library(ggsci)
npg_pal <- pal_npg()(10)
npg_extended <- colorRampPalette(npg_pal)(38)
p1 <- DimPlot(nano.obj, cols = npg_extended)
ggsave("nano_cluster.pdf",plot=p1)
```
![umap](./pic/umap.png "umap")<br>
## 3.5 常规方法--无监督聚类
这里`标准化`也可以使用单细胞测序的常规方法：
```
nano.obj <- NormalizeData(nano.obj)
nano.obj <- ScaleData(nano.obj)
nano.obj <- FindVariableFeatures(nano.obj)
```
然后`PCA降维`，这里的维度Nanostring使用手册上默认50，因此我这里也设置50：
```
nano.obj <- RunPCA(nano.obj, features = VariableFeatures(object = nano.obj))
nano.obj <- RunUMAP(nano.obj, dims = 1:50)
```
`聚类`：
```
nano.obj <- FindNeighbors(nano.obj, reduction = "pca", dims = 1:50)
nano.obj <- FindClusters(nano.obj, resolution = 0.3)

DimPlot(nano.obj, raster = FALSE, label = TRUE)
```
`细胞注释`，这里的手动注释是查看Cellmarker2.0：
```
cluster_markers <- FindAllMarkers(nano.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
list_marker <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
df_marker=data.frame(p_val = list_marker$p_val,
                     avg_log2FC = list_marker$avg_log2FC,
                     pct.1 = list_marker$pct.1,
                     pct.2 = list_marker$pct.2,
                     p_val_adj = list_marker$p_val_adj,
                     cluster = list_marker$cluster,
                     gene = list_marker$gene)
write.csv(df_marker,"marker.csv")

nano.obj <- RenameIdents(nano.obj, '0' = "Tumor", '1' = "Macro", '2' = "naive T cell", '3' = "Fibroblast", '4' = "Neutrophil", '5' = "Endothelial cell", '6' = "B", '7' = "Plasma", '8' = "Macro", '9' = "Proliferative", '10' = "MPC", '11' = "CTL")

DimPlot(nano.obj, reduction = "umap", label = TRUE, pt.size = 0.5)
```
![无监督聚类](./pic/无监督聚类.png "无监督聚类")<br>
# 4 细胞类型和表达定位模式的可视化
ImageDimPlot() 这个函数会根据细胞在空间上的分布位置来绘制它们，并依据细胞被指定的类型来对它们进行颜色标记。可以观察到，基底细胞群在空间上的排列非常紧凑有序，这与我们的预期是一致的。
```
ImageDimPlot(nano.obj, fov = "lung5.rep1", axes = TRUE, cols = "glasbey")
```
 * cols = "glasbey" 这样的参数设置通常出现在绘图函数中，用于指定颜色方案。<br>
![细胞分布](./pic/细胞分布.png "细胞分布")<br>
可以突出显示一些选定组的定位：
```
ImageDimPlot(nano.obj, fov = "lung5.rep1", cells = WhichCells(nano.obj, idents = c("Basal", "Macrophage", "Smooth Muscle", "CD4 T")), cols = c("red", "green", "blue", "orange"), size = 0.6)
```
![部分细胞](./pic/部分细胞.png "部分细胞分布")<br>
# 5 可视化基因表达标记
KRT17是上皮基底细胞的标志物，下面以它为例进行可视化。
## 5.1 VlnPlot
```
VlnPlot(nano.obj, features = "KRT17", assay = "Nanostring", layer = "counts", pt.size = 0.1, y.max = 30) + NoLegend()
```
* 在Seurat中，一个assay是一组相关的数据，可以是基因表达矩阵、蛋白质表达矩阵等。<br>
![VlnPlot](./pic/VlnPlot.png "VlnPlot")<br>
## 5.2 FeaturePlot
FeaturePlot 函数用于生成特定基因或特征的表达水平图。这种图通常以点图的形式展示，可以直观地看到不同细胞中特定基因的表达水平。
```
FeaturePlot(nano.obj, features = "KRT17", max.cutoff = "q95")
```
* max.cutoff = "q95" 是一个参数，用于指定表达水平的上限。超过这个值的数据将不被显示。
![FeaturePlot](./pic/FeaturePlot.png "FeaturePlot")
## 5.3 ImageFeaturePlot
```
ImageFeaturePlot(nano.obj, fov = "lung5.rep1", features = "KRT17", max.cutoff = "q95")
```
![ImageFeaturePlot](./pic/ImageFeaturePlot.png "ImageFeaturePlot")
## 5.4 ImageDimPlot
```
ImageDimPlot(nano.obj, fov = "lung5.rep1", alpha = 0.3, molecules = "KRT17", nmols = 10000) + NoLegend()
```
* alpha = 0.3 设置了图中点的透明度。
* nmols = 10000 是一个参数，用于指定在图中显示的最大分子数量。nmols 参数限制了图中显示的分子点的数量，这有助于控制图像的复杂性和可读性。
![ImageDimPlot](./pic/ImageDimPlot.png "ImageDimPlot")<br>
* 还可以共同可视化多个标记物的表达，包括 KRT17（基底细胞）、C1QA（巨噬细胞）、IL7R（T 细胞）和 TAGLN（平滑肌细胞）。
```
ImageDimPlot(nano.obj, fov = "lung5.rep1", group.by = NA, alpha = 0.3, molecules = c("KRT17", "C1QA", "IL7R", "TAGLN"), nmols = 20000)
```
![多个标记物分布](./pic/多个标记物分布.png "多个标记物分布")<br>
* 还可以使用 Crop() 函数放大一个富含基底的区域。放大后，可以看见单个细胞边界。
```
basal.crop <- Crop(nano.obj[["lung5.rep1"]], x = c(159500, 164000), y = c(8700, 10500))
nano.obj[["zoom1"]] <- basal.crop
DefaultBoundary(nano.obj[["zoom1"]]) <- "segmentation"
ImageDimPlot(nano.obj, fov = "zoom1", cols = "polychrome", coord.fixed = FALSE)
```
![zoom-in](./pic/zoom-in.png "zoom-in")<br>
* 注释细胞和标志物
```
ImageDimPlot(nano.obj, fov = "zoom1", cols = "polychrome", alpha = 0.3, molecules = c("KRT17", "IL7R", "TPSAB1"), mols.size = 0.3, nmols = 20000, border.color = "black", coord.fixed = FALSE)
```
![marker-cell](./pic/marker-cell.png "marker-cell")<br>

---
下面使用淋巴结的数据进行查看：
```
curl -LO https://smi-public.objects.liquidweb.services/LN28_6k/seurat.zip
unzip seurat.zip
```
这里提供了Seurat格式的数据直接下载转化为seurat对象即可：
```
library(Seurat)
lymph.obj <- readRDS("seurat.rds")
```
数据不全，没有空间信息。。。<br>

下载原始数据后想绘制umap图但是内存不够。。。

# 互作分析
```
# 获取坐标信息
coords <- nano.obj@images$lung5.rep1@boundaries$centroids@coords
head(coords)
# 从 centroids 对象获取正确的 cells/barcodes
centroid_cells <- nano.obj@images$lung5.rep1@boundaries$centroids@cells
rownames(coords) <- centroid_cells

# 获取细胞类型标签 - 使用 predicted.annotation.l1
labels <- as.character(nano.obj@meta.data$predicted.annotation.l1)
names(labels) <- rownames(nano.obj@meta.data)

# 进行匹配
common <- intersect(rownames(coords), names(labels))
print(paste("共同barcodes数量:", length(common)))

# 检查细胞类型分布
print("细胞类型分布:")
print(table(labels))

# 处理数据
coords <- coords[common, , drop = FALSE]
labels <- labels[common]
barcodes <- common

# 移除NA值的细胞
valid_cells <- !is.na(labels)
coords <- coords[valid_cells, , drop = FALSE]
labels <- labels[valid_cells]
barcodes <- barcodes[valid_cells]

print(paste("有效细胞数量:", length(labels)))
print("最终细胞类型分布:")
print(table(labels))

# 继续处理...
xy <- coords
colnames(xy) <- c("x","y")

print("开始计算Delaunay三角剖分...")
library(deldir)
x <- xy[,1]
y <- xy[,2]
rw <- c(min(x), max(x), min(y), max(y))

deld <- deldir(x, y, rw = rw)
segs <- deld$delsgs

print(paste("生成三角边数量:", nrow(segs)))

# 使用索引
edges <- cbind(segs$ind1, segs$ind2)

# 去除自环与重复
edges <- edges[edges[,1] != edges[,2], , drop=FALSE]
edges <- t(apply(edges, 1, function(x) sort(x)))
edges <- unique(edges)
edges_df <- data.frame(from = edges[,1], to = edges[,2])

print(paste("最终边数量:", nrow(edges_df)))

# 获取所有细胞类型
types <- sort(unique(labels))
K <- length(types)

print(paste("细胞类型数量:", K))
print("细胞类型:")
print(types)

# 将索引映射到类型
type_by_index <- labels
index_to_type <- type_by_index

# For edges_df, get types:
t1 <- index_to_type[edges_df$from]
t2 <- index_to_type[edges_df$to]

# 构建矩阵
mat_obs <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))

for(i in seq_along(t1)){
    a <- t1[i]
    b <- t2[i]
    mat_obs[a,b] <- mat_obs[a,b] + 1
    mat_obs[b,a] <- mat_obs[b,a] + 1
}

print("观察到的细胞-细胞接触矩阵：")
print(mat_obs)

# 排列检验
library(doParallel)
library(foreach)

nperm <- 1000
ncores <- parallel::detectCores() - 1
ncores <- max(1, ncores)

print(paste("使用", ncores, "个核心进行排列检验"))

cl <- makeCluster(ncores)
registerDoParallel(cl)

from_idx <- edges_df$from
to_idx <- edges_df$to

perm_counts <- foreach(p = 1:nperm, .packages = c(), .combine = rbind) %dopar% {
    set.seed(p + 12345)
    perm_labels <- sample(index_to_type, length(index_to_type), replace = FALSE)
    
    pt1 <- perm_labels[from_idx]
    pt2 <- perm_labels[to_idx]
    
    mat <- matrix(0, nrow = K, ncol = K, dimnames = list(types, types))
    
    for(i in seq_along(pt1)){
        a <- pt1[i]
        b <- pt2[i]
        mat[a,b] <- mat[a,b] + 1
        mat[b,a] <- mat[b,a] + 1
    }
    as.vector(mat)
}

stopCluster(cl)

print("排列检验完成")

# 继续剩余的分析代码...
obs_vec <- as.vector(mat_obs)
mu_rand <- colMeans(perm_counts)
sd_rand <- apply(perm_counts, 2, sd)

# z score
z_vec <- (obs_vec - mu_rand) / (sd_rand + 1e-8)

# 经验 p 值
p_emp <- sapply(seq_along(obs_vec), function(i){
  perm_i <- perm_counts[, i]
  obs_i <- obs_vec[i]
  mu_i <- mu_rand[i]
  p_val = (sum(abs(perm_i - mu_i) >= abs(obs_i - mu_i)) + 1) / (nperm + 1)
  p_val
})

# 转回矩阵形式
mat_mu <- matrix(mu_rand, nrow = K, ncol = K, dimnames = list(types, types))
mat_sd <- matrix(sd_rand, nrow = K, ncol = K, dimnames = list(types, types))
mat_z <- matrix(z_vec, nrow = K, ncol = K, dimnames = list(types, types))
mat_p <- matrix(p_emp, nrow = K, ncol = K, dimnames = list(types, types))

# z-score 热图
library(pheatmap)
pheatmap(mat_z,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Z-score of contact enrichment (Obs vs Random)",
         fontsize = 10)

print("分析完成！")

# 保存结果
save(mat_obs, mat_z, mat_p, mat_mu, file = "spatial_contact_analysis_results.RData")

#显著性
# 创建显著性标注矩阵
signif_symbols <- matrix("", nrow = K, ncol = K, dimnames = list(types, types))

# 设置显著性阈值
signif_symbols[mat_p < 0.001] <- "***"
signif_symbols[mat_p < 0.01 & mat_p >= 0.001] <- "**"
signif_symbols[mat_p < 0.05 & mat_p >= 0.01] <- "*"

# 绘制带显著性标注的热图
p4 <- pheatmap(mat_z,
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         main = "Z-score of contact enrichment (Obs vs Random)",
         fontsize = 10,
         display_numbers = signif_symbols,  # 显示显著性星号
         number_color = "black",           # 星号颜色
         fontsize_number = 12)             # 星号大小
ggsave("nano_zscore.pdf",plot=p4,width=10,height=8)


## 修复同配性分析
library(igraph)

# 创建图对象
g <- graph_from_edgelist(as.matrix(edges_df), directed = FALSE)

# 将细胞类型转换为数值因子
cell_type_factor <- as.factor(index_to_type)
V(g)$cell_type_numeric <- as.numeric(cell_type_factor)
V(g)$cell_type_original <- index_to_type

# 计算同配性系数
assortativity_global <- assortativity_nominal(g, types = cell_type_factor)
print(paste("全局同配性系数:", round(assortativity_global, 4)))

# 解释同配性系数
if(assortativity_global > 0.3) {
  print("强同配性: 同类型细胞倾向于相互接触")
} else if(assortativity_global > 0.1) {
  print("中等同配性: 一定程度同类型接触偏好")
} else if(assortativity_global > -0.1) {
  print("随机混合: 细胞类型间接触接近随机")
} else if(assortativity_global > -0.3) {
  print("中等异配性: 一定程度不同类型接触偏好")
} else {
  print("强异配性: 不同类型细胞倾向于相互接触")
}

# 计算每对细胞类型的同配性
pairwise_assortativity <- function(g, type1, type2) {
  # 创建二分图：type1和type2 vs 其他所有类型
  binary_types <- ifelse(V(g)$cell_type_original %in% c(type1, type2), "pair", "other")
  binary_factor <- as.factor(binary_types)
  assortativity_nominal(g, types = binary_factor)
}

# 对所有类型对进行计算
pairwise_results <- matrix(NA, nrow = K, ncol = K, dimnames = list(types, types))

for(i in 1:K) {
  for(j in 1:K) {
    if(i == j) {
      # 对角线：同类型接触的同配性
      binary_types <- ifelse(V(g)$cell_type_original == types[i], "same", "other")
      binary_factor <- as.factor(binary_types)
      pairwise_results[i,j] <- assortativity_nominal(g, types = binary_factor)
    } else {
      pairwise_results[i,j] <- pairwise_assortativity(g, types[i], types[j])
    }
  }
}

print("成对同配性矩阵:")
print(round(pairwise_results, 4))

# 连通区域分析
comp_analysis <- components(g)
print(paste("连通区域总数:", comp_analysis$no))
print("连通区域大小分布:")
print(table(comp_analysis$csize))

# 分析每个细胞类型在连通区域中的分布
print("各细胞类型在连通区域中的分布:")
for(type in types) {
  type_vertices <- which(V(g)$cell_type_original == type)
  if(length(type_vertices) > 0) {
    type_components <- comp_analysis$membership[type_vertices]
    unique_comps <- length(unique(type_components))
    avg_component_size <- mean(table(type_components))
    print(paste(type, ": 分布在", unique_comps, "个连通区域中，平均区域大小:", round(avg_component_size, 2)))
  }
}

# 可视化连通区域
if(comp_analysis$no > 1) {
  # 给每个连通区域分配颜色
  comp_colors <- rainbow(comp_analysis$no)
  V(g)$color <- comp_colors[comp_analysis$membership]
  
  # 简化绘图（如果图太大）
  if(vcount(g) > 1000) {
    # 随机抽样或只绘制最大的几个连通区域
    largest_comp <- which.max(comp_analysis$csize)
    sub_vertices <- which(comp_analysis$membership == largest_comp)
    sub_g <- induced_subgraph(g, sub_vertices)
    plot(sub_g, vertex.size = 3, vertex.label = NA, 
         main = paste("最大连通区域 (大小:", comp_analysis$csize[largest_comp], ")"))
  } else {
    plot(g, vertex.size = 3, vertex.label = NA, 
         main = paste("空间接触网络 (", comp_analysis$no, "个连通区域)"))
  }
}


#可视化连通性
# 在空间图像上叠加连通图
library(ggplot2)

# 准备数据
spatial_df <- data.frame(
  x = xy[,1],
  y = xy[,2],
  cell_type = labels,
  component = comp_analysis$membership
)

# 获取边数据
edge_coords <- do.call(rbind, lapply(1:nrow(edges_df), function(i) {
  from_idx <- edges_df$from[i]
  to_idx <- edges_df$to[i]
  data.frame(
    x = spatial_df$x[from_idx],
    y = spatial_df$y[from_idx],
    xend = spatial_df$x[to_idx],
    yend = spatial_df$y[to_idx],
    component = spatial_df$component[from_idx]
  )
}))

# 绘制空间连通图
p5 <- ggplot() +
  # 绘制细胞点
  geom_point(data = spatial_df, aes(x = x, y = y, color = cell_type), 
             size = 1, alpha = 0.7) +
  # 绘制连通边
  geom_segment(data = edge_coords, 
               aes(x = x, y = y, xend = xend, yend = yend, group = component),
               color = "gray40", alpha = 0.3, linewidth = 0.2) +
  # 按连通区域绘制轮廓
  geom_density_2d_filled(data = spatial_df, 
                        aes(x = x, y = y, group = component), 
                        alpha = 0.1) +
  theme_minimal() +
  coord_fixed() +
  labs(title = "空间细胞连通网络",
       subtitle = paste("包含", comp_analysis$no, "个连通区域"))

ggsave("liant.pdf",plot=p5)

#zoom的可视化
# 正确的方法：使用与之前一致的坐标获取方式
coords <- nano.obj@images$lung5.rep1@boundaries$centroids@coords
centroid_cells <- nano.obj@images$lung5.rep1@boundaries$centroids@cells
rownames(coords) <- centroid_cells

# 获取zoom区域的坐标信息
zoom_coords <- GetTissueCoordinates(nano.obj[["zoom1"]])
zoom_cells <- rownames(zoom_coords)

print("Zoom区域细胞barcodes:")
print(head(zoom_cells))
print("总坐标细胞barcodes:")
print(head(rownames(coords)))

# 提取zoom区域内的细胞
zoom_idx <- which(rownames(coords) %in% zoom_cells)
zoom_xy <- coords[zoom_idx, ]
zoom_labels <- labels[zoom_idx]  # 使用之前定义好的labels
zoom_barcodes <- rownames(coords)[zoom_idx]

print(paste("Zoom区域细胞数量:", length(zoom_labels)))
print("Zoom区域细胞类型分布:")
print(table(zoom_labels))

# 为zoom区域重新计算Delaunay三角网
library(deldir)
x_zoom <- zoom_xy[,1]
y_zoom <- zoom_xy[,2]
rw_zoom <- c(min(x_zoom), max(x_zoom), min(y_zoom), max(y_zoom))

deld_zoom <- deldir(x_zoom, y_zoom, rw = rw_zoom)
segs_zoom <- deld_zoom$delsgs

# 创建zoom区域的边
edges_zoom <- cbind(segs_zoom$ind1, segs_zoom$ind2)
edges_zoom <- edges_zoom[edges_zoom[,1] != edges_zoom[,2], , drop=FALSE]
edges_zoom <- t(apply(edges_zoom, 1, function(x) sort(x)))
edges_zoom <- unique(edges_zoom)
edges_df_zoom <- data.frame(from = edges_zoom[,1], to = edges_zoom[,2])

print(paste("Zoom区域边数量:", nrow(edges_df_zoom)))

# 继续之前的可视化代码...
library(igraph)
g_zoom <- graph_from_edgelist(as.matrix(edges_df_zoom), directed = FALSE)
V(g_zoom)$cell_type <- zoom_labels

comp_zoom <- components(g_zoom)
print(paste("Zoom区域连通区域数量:", comp_zoom$no))

# 准备绘图数据
spatial_df_zoom <- data.frame(
  x = zoom_xy[,1],
  y = zoom_xy[,2],
  cell_type = zoom_labels,
  component = comp_zoom$membership,
  barcode = zoom_barcodes
)

# 获取边坐标
edge_coords_zoom <- do.call(rbind, lapply(1:nrow(edges_df_zoom), function(i) {
  from_idx <- edges_df_zoom$from[i]
  to_idx <- edges_df_zoom$to[i]
  data.frame(
    x = spatial_df_zoom$x[from_idx],
    y = spatial_df_zoom$y[from_idx],
    xend = spatial_df_zoom$x[to_idx],
    yend = spatial_df_zoom$y[to_idx],
    component = spatial_df_zoom$component[from_idx]
  )
}))

zoom_cells_modified <- paste0(zoom_cells, "_1")
print("修改后的zoom barcodes:")
print(head(zoom_cells_modified))

# 提取zoom区域内的细胞
zoom_idx <- which(rownames(coords) %in% zoom_cells_modified)

if(length(zoom_idx) == 0) {
  # 方法2：如果还是匹配不上，检查是否有其他模式
  print("尝试其他匹配模式...")
  
  # 查看zoom坐标的完整结构
  print("Zoom坐标结构:")
  str(zoom_coords)
  
  # 检查是否有cells信息在zoom图像中
  if("cells" %in% slotNames(nano.obj[["zoom1"]])) {
    zoom_cells_correct <- nano.obj[["zoom1"]]@cells
    print("从zoom图像获取cells:")
    print(head(zoom_cells_correct))
    zoom_idx <- which(rownames(coords) %in% zoom_cells_correct)
  }
}

if(length(zoom_idx) == 0) {
  # 方法3：使用空间坐标范围来匹配
  print("使用空间坐标范围匹配...")
  x_range <- range(zoom_coords$x)
  y_range <- range(zoom_coords$y)
  
  print(paste("Zoom区域坐标范围: x(", x_range[1], "-", x_range[2], 
              "), y(", y_range[1], "-", y_range[2], ")"))
  
  # 在总坐标中查找在zoom区域范围内的细胞
  in_zoom_x <- coords[,1] >= x_range[1] & coords[,1] <= x_range[2]
  in_zoom_y <- coords[,2] >= y_range[1] & coords[,2] <= y_range[2]
  zoom_idx <- which(in_zoom_x & in_zoom_y)
}

print(paste("找到的zoom区域细胞数量:", length(zoom_idx)))

if(length(zoom_idx) > 0) {
  zoom_xy <- coords[zoom_idx, ]
  zoom_labels <- labels[zoom_idx]
  zoom_barcodes <- rownames(coords)[zoom_idx]
  
  print(paste("Zoom区域细胞数量:", length(zoom_labels)))
  print("Zoom区域细胞类型分布:")
  print(table(zoom_labels))
  
  # 为zoom区域重新计算Delaunay三角网
  x_zoom <- zoom_xy[,1]
  y_zoom <- zoom_xy[,2]
  rw_zoom <- c(min(x_zoom), max(x_zoom), min(y_zoom), max(y_zoom))
  
  deld_zoom <- deldir(x_zoom, y_zoom, rw = rw_zoom)
  segs_zoom <- deld_zoom$delsgs
  
  # 创建zoom区域的边
  edges_zoom <- cbind(segs_zoom$ind1, segs_zoom$ind2)
  edges_zoom <- edges_zoom[edges_zoom[,1] != edges_zoom[,2], , drop=FALSE]
  edges_zoom <- t(apply(edges_zoom, 1, function(x) sort(x)))
  edges_zoom <- unique(edges_zoom)
  edges_df_zoom <- data.frame(from = edges_zoom[,1], to = edges_zoom[,2])
  
  print(paste("Zoom区域边数量:", nrow(edges_df_zoom)))
  
  # 继续可视化...
  library(igraph)
  g_zoom <- graph_from_edgelist(as.matrix(edges_df_zoom), directed = FALSE)
  V(g_zoom)$cell_type <- zoom_labels
  
  comp_zoom <- components(g_zoom)
  print(paste("Zoom区域连通区域数量:", comp_zoom$no))
  
  # 准备绘图数据
  spatial_df_zoom <- data.frame(
    x = zoom_xy[,1],
    y = zoom_xy[,2],
    cell_type = zoom_labels,
    component = comp_zoom$membership,
    barcode = zoom_barcodes
  )
  
  # 获取边坐标
  edge_coords_zoom <- do.call(rbind, lapply(1:nrow(edges_df_zoom), function(i) {
    from_idx <- edges_df_zoom$from[i]
    to_idx <- edges_df_zoom$to[i]
    data.frame(
      x = spatial_df_zoom$x[from_idx],
      y = spatial_df_zoom$y[from_idx],
      xend = spatial_df_zoom$x[to_idx],
      yend = spatial_df_zoom$y[to_idx],
      component = spatial_df_zoom$component[from_idx]
    )
  }))
  
  # 可视化
  library(ggplot2)
    p_zoom <- ggplot() +
  geom_segment(data = edge_coords_zoom, 
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "darkgray", alpha = 0.6, linewidth = 0.3) +
  geom_point(data = spatial_df_zoom, 
             aes(x = x, y = y, color = cell_type), 
             size = 2.5, alpha = 0.8) +
  geom_text(data = spatial_df_zoom %>% 
             dplyr::group_by(component) %>% 
             dplyr::summarise(x = mean(x), y = mean(y), count = dplyr::n()),
            aes(x = x, y = y, label = paste("区域", component, "\n(", count, "细胞)")),
            size = 2.5, fontface = "bold", color = "darkred") +
    scale_color_manual(values = npg_extended) +
  theme_minimal() +
  coord_fixed() +
  labs(title = "Zoom区域细胞连通网络",
       subtitle = paste(comp_zoom$no, "个连通区域,", nrow(edges_df_zoom), "条连接边"),
       color = "细胞类型")
}

ggsave("nano_zoom_liant.pdf",plot=p_zoom,width = 20,height = 15)


#同配性热图+zscore显著性
library(pheatmap)

# 创建显著性标注
#signif_symbols <- matrix("", nrow = K, ncol = K, dimnames = list(types, types))

# 设置显著性阈值（基于经验）
#signif_symbols[abs(pairwise_results) > 0.3] <- "***"
#signif_symbols[abs(pairwise_results) > 0.2 & abs(pairwise_results) <= 0.3] <- "**"
#signif_symbols[abs(pairwise_results) > 0.1 & abs(pairwise_results) <= 0.2] <- "*"

p6 <- pheatmap(pairwise_results,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = signif_symbols,
         number_color = "black",
         fontsize_number = 12,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         main = "成对同配性矩阵\n(值越高表示越倾向于共同定位)",
         fontsize = 10)
ggsave("assort.pdf",plot=p6,width=12,height=8)
```
# 互作机制分析
* 识别强互作对
```
# 基于z-score和同配性筛选强互作对
strong_interactions <- data.frame(
  cell_type1 = rep(types, each = K),
  cell_type2 = rep(types, times = K),
  z_score = as.vector(mat_z),
  p_value = as.vector(mat_p),
  assortativity = as.vector(pairwise_results)
) %>%
  filter(cell_type1 != cell_type2) %>%
  filter(p_value < 0.05 & abs(z_score) > 2) %>%  # 显著且强互作
  arrange(desc(abs(z_score)))

print("强互作细胞对:")
print(strong_interactions)
```
* 自选互作对分析
```
# 检查Basal和Myofibroblast在矩阵中的位置
basal_idx <- which(types == "Basal")
myofibro_idx <- which(types == "Myofibroblast")

print(paste("Basal位置:", basal_idx))
print(paste("Myofibroblast位置:", myofibro_idx))

# 提取这对细胞的互作指标
basal_myofibro_stats <- data.frame(
  z_score = mat_z[basal_idx, myofibro_idx],
  p_value = mat_p[basal_idx, myofibro_idx],
  assortativity = pairwise_results[basal_idx, myofibro_idx],
  observed_contacts = mat_obs[basal_idx, myofibro_idx],
  expected_contacts = mat_mu[basal_idx, myofibro_idx]
)

print("Basal-Myofibroblast互作统计:")
print(basal_myofibro_stats)


identify_basal_myofibro_edges <- function(edges_df, labels, basal_type = "Basal", myofibro_type = "Myofibroblast") {
  interacting_pairs <- c()
  
  for(i in 1:nrow(edges_df)) {
    cell1_type <- labels[edges_df$from[i]]
    cell2_type <- labels[edges_df$to[i]]
    
    # 检查是否是Basal-Myofibroblast互作
    if((cell1_type == basal_type & cell2_type == myofibro_type) |
       (cell1_type == myofibro_type & cell2_type == basal_type)) {
      interacting_pairs <- c(interacting_pairs, i)
    }
  }
  
  return(edges_df[interacting_pairs, ])
}

# 获取互作边
bm_edges <- identify_basal_myofibro_edges(edges_df, labels)
print(paste("Basal-Myofibroblast直接接触对数量:", nrow(bm_edges)))

# 获取互作细胞
bm_interacting_cells <- unique(c(
  bm_edges$from, bm_edges$to
))

print(paste("参与互作的细胞数量:", length(bm_interacting_cells)))
```
* 提取互作区域的基因表达特征
```
#基底

# 首先设置细胞标识
nano.obj$is_direct_contact <- ifelse(colnames(nano.obj) %in% basal_direct_barcodes, 
                                    "Direct_Contact", "Non_Contact")

# 设置活跃标识
Idents(nano.obj) <- nano.obj$is_direct_contact

print(paste("直接接触Basal细胞:", sum(nano.obj$is_direct_contact == "Direct_Contact")))
print(paste("非接触Basal细胞:", sum(nano.obj$is_direct_contact == "Non_Contact")))

# 现在使用正确的FindMarkers语法
if(sum(nano.obj$is_direct_contact == "Direct_Contact") >= 10 & 
   sum(nano.obj$is_direct_contact == "Non_Contact") >= 10) {
  
  basal_direct_de <- FindMarkers(nano.obj,
                                ident.1 = "Direct_Contact",
                                ident.2 = "Non_Contact",
                                logfc.threshold = 0.25,
                                min.pct = 0.1)
  
  basal_direct_sig <- basal_direct_de %>%
    mutate(gene = rownames(.)) %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
    arrange(desc(avg_log2FC))
  
  print(paste("直接接触Basal细胞显著基因数:", nrow(basal_direct_sig)))
  print(head(basal_direct_sig, 10))
  
  # 生成直接接触Basal火山图数据（处理NA和Inf）
  basal_direct_volcano_data <- basal_direct_de %>%
    mutate(
      gene = rownames(.),
      p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj),  # 处理0值避免Inf
      y_value = -log10(p_val_adj),
      y_value = ifelse(is.na(y_value), 0, y_value),  # 替换NA
      significance = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, "Significant", "Not significant"),
      label = ifelse(gene %in% head(basal_direct_sig$gene, 8), gene, "")
    )
  
  # 计算y轴上限
  max_y <- max(basal_direct_volcano_data$y_value[is.finite(basal_direct_volcano_data$y_value)], na.rm = TRUE)
  y_upper <- ifelse(is.finite(max_y), max_y * 1.1, 30)
  
  # 绘制火山图（移除网格线，保留坐标轴）
  p_basal_direct <- ggplot(basal_direct_volcano_data, aes(x = avg_log2FC, y = y_value)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 2) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, box.padding = 0.5) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
    ylim(0, y_upper) +
    # 主题设置：移除网格线，保留坐标轴
    theme_bw() +  # 使用带坐标轴的主题为基础
    theme(
      panel.grid.major = element_blank(),  # 移除主网格线
      panel.grid.minor = element_blank(),  # 移除次网格线
      axis.line = element_line(color = "black")  # 确保坐标轴线显示
    ) +
    labs(title = "Basal细胞: 直接接触 vs 非接触",
         subtitle = paste("显著基因:", nrow(basal_direct_sig), "个"),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)")
  
  print(p_basal_direct)
  
  # 保存为PDF（可调整宽度和高度）
  ggsave(
    filename = "basal_direct_contact_volcano.pdf",
    plot = p_basal_direct,
    device = "pdf",
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  print("直接接触Basal细胞数量不足")
  basal_direct_sig <- data.frame()
}

#成纤维
print("=== Myofibroblast直接接触分析 ===")

# 设置Myofibroblast的接触标识
nano.obj$myofibro_contact <- ifelse(colnames(nano.obj) %in% myofibro_direct_barcodes, 
                                   "Myofibro_Contact", "Myofibro_Non_Contact")

# 只对Myofibroblast细胞进行分析
myofibro_cells <- colnames(nano.obj)[nano.obj$cell_type == "Myofibroblast"]
nano.obj_myofibro <- subset(nano.obj, cells = myofibro_cells)

# 设置标识
Idents(nano.obj_myofibro) <- nano.obj_myofibro$myofibro_contact

print(paste("直接接触Myofibroblast细胞:", sum(nano.obj_myofibro$myofibro_contact == "Myofibro_Contact")))
print(paste("非接触Myofibroblast细胞:", sum(nano.obj_myofibro$myofibro_contact == "Myofibro_Non_Contact")))

if(sum(nano.obj_myofibro$myofibro_contact == "Myofibro_Contact") >= 10 & 
   sum(nano.obj_myofibro$myofibro_contact == "Myofibro_Non_Contact") >= 10) {
  
  myofibro_direct_de <- FindMarkers(nano.obj_myofibro,
                                   ident.1 = "Myofibro_Contact",
                                   ident.2 = "Myofibro_Non_Contact",
                                   logfc.threshold = 0.25,
                                   min.pct = 0.1)
  
  myofibro_direct_sig <- myofibro_direct_de %>%
    mutate(gene = rownames(.)) %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>%
    arrange(desc(avg_log2FC))
  
  print(paste("直接接触Myofibroblast细胞显著基因数:", nrow(myofibro_direct_sig)))
  print(head(myofibro_direct_sig, 10))
  
  # 生成直接接触Myofibroblast火山图数据（处理NA和Inf）
  myofibro_direct_volcano_data <- myofibro_direct_de %>%
    mutate(
      gene = rownames(.),
      p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj),  # 处理0值避免Inf
      y_value = -log10(p_val_adj),
      y_value = ifelse(is.na(y_value), 0, y_value),  # 替换NA值
      significance = ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5, "Significant", "Not significant"),
      label = ifelse(gene %in% head(myofibro_direct_sig$gene, 8), gene, "")
    )
  
  # 计算y轴上限（带缓冲）
  max_y <- max(myofibro_direct_volcano_data$y_value[is.finite(myofibro_direct_volcano_data$y_value)], na.rm = TRUE)
  y_upper <- ifelse(is.finite(max_y), max_y * 1.1, 30)  # 设置默认上限以防极端情况
  
  # 绘制火山图（移除网格线，保留坐标轴）
  p_myofibro_direct <- ggplot(myofibro_direct_volcano_data, aes(x = avg_log2FC, y = y_value)) +
    geom_point(aes(color = significance), alpha = 0.7, size = 2) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, box.padding = 0.5) +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
    ylim(0, y_upper) +  # 设置y轴范围避免顶部截断
    # 主题设置：移除网格线，保留坐标轴线
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),  # 移除主网格线
      panel.grid.minor = element_blank(),  # 移除次网格线
      axis.line = element_line(color = "black")  # 确保坐标轴线显示
    ) +
    labs(title = "Myofibroblast细胞: 直接接触 vs 非接触",
         subtitle = paste("显著基因:", nrow(myofibro_direct_sig), "个"),
         x = "log2 Fold Change", y = "-log10(adjusted p-value)")
  
  print(p_myofibro_direct)
  
  # 保存为PDF
  ggsave(
    filename = "myofibroblast_direct_contact_volcano.pdf",
    plot = p_myofibro_direct,
    device = "pdf",
    width = 10,
    height = 8,
    dpi = 300
  )
  
} else {
  print("直接接触Myofibroblast细胞数量不足")
  myofibro_direct_sig <- data.frame()
}

#共同热图
# 分析共同变化的基因
if(exists("basal_direct_sig") && nrow(basal_direct_sig) > 0 && 
   exists("myofibro_direct_sig") && nrow(myofibro_direct_sig) > 0) {
  
  common_direct_genes <- intersect(basal_direct_sig$gene, myofibro_direct_sig$gene)
  print(paste("直接接触共同显著基因:", length(common_direct_genes)))
  
  if(length(common_direct_genes) > 0) {
    print("共同基因列表:")
    print(common_direct_genes)
    
    # 可视化共同基因
    common_genes_expr <- GetAssayData(nano.obj, slot = "scale.data")[common_direct_genes, ]
    
    # 创建注释
    annotation_common <- data.frame(
      Cell_Type = nano.obj$cell_type,
      Contact_Status = nano.obj$is_direct_contact,
      row.names = colnames(nano.obj)
    )
    
    # 筛选需要展示的细胞（与之前分析一致）
    selected_cells <- colnames(nano.obj) %in% c(basal_direct_barcodes, myofibro_direct_barcodes)
    
    # 热图显示共同基因并保存为PDF
    pheatmap(
      common_genes_expr[, selected_cells],
      annotation_col = annotation_common,
      show_colnames = FALSE,
      main = "直接接触共同显著基因表达",
      fontsize_row = 8,
      # 保存为PDF参数
      filename = "common_direct_genes_heatmap.pdf",  # 文件名
      device = "pdf",                               # 设备类型
      width = 10,                                   # 宽度
      height = 8 + length(common_direct_genes)*0.1  # 高度（根据基因数量动态调整）
    )
  }
}
```
# 受体配体分析
```
library(CellChat)

# 重新创建Basal-Myofibroblast子集
bm_cells <- colnames(nano.obj)[nano.obj$cell_type %in% c("Basal", "Myofibroblast")]
nano.obj_bm <- subset(nano.obj, cells = bm_cells)

print(paste("子集细胞数:", ncol(nano.obj_bm)))
print(table(nano.obj_bm$cell_type))

# 正确的CellChat分析流程
cellchat_bm <- createCellChat(object = nano.obj_bm, group.by = "cell_type")

# 1. 首先设置数据库
CellChatDB <- CellChatDB.human
cellchat_bm@DB <- CellChatDB

# 2. 子集化信号基因数据
cellchat_bm <- subsetData(cellchat_bm)  # 这一步必须在identifyOverExpressedGenes之前

# 3. 识别过表达基因
cellchat_bm <- identifyOverExpressedGenes(cellchat_bm)

# 4. 识别过表达相互作用
cellchat_bm <- identifyOverExpressedInteractions(cellchat_bm)

# 5. 计算通信概率（使用默认参数）
cellchat_bm <- computeCommunProb(cellchat_bm)

# 6. 如果细胞数量较少，可能需要调整参数
# cellchat_bm <- computeCommunProb(cellchat_bm, type = "truncatedMean", trim = 0.1)

# 7. 计算通路水平的通信概率
cellchat_bm <- computeCommunProbPathway(cellchat_bm)

# 8. 整合通信网络
cellchat_bm <- aggregateNet(cellchat_bm)

print("CellChat分析成功完成！")

#提取显著受体配体对
# 提取详细的受体-配体对信息
extract_detailed_LR_pairs <- function(cellchat_obj, source, target, prob_threshold = 0.01) {
  # 获取通信概率
  prob_matrix <- cellchat_obj@net$prob
  lr_probs <- prob_matrix[source, target, ]
  
  # 筛选显著的对
  significant_idx <- which(lr_probs > prob_threshold)
  
  if(length(significant_idx) > 0) {
    significant_probs <- lr_probs[significant_idx]
    lr_names <- names(significant_probs)
    
    # 从数据库获取详细信息
    interaction_db <- cellchat_obj@DB$interaction
    pair_details <- interaction_db[interaction_db$interaction_name %in% lr_names, ]
    
    # 合并信息
    result <- data.frame(
      interaction_name = lr_names,
      communication_prob = significant_probs,
      ligand = pair_details$ligand,
      receptor = pair_details$receptor,
      pathway = pair_details$pathway_name
    ) %>%
      arrange(desc(communication_prob))
    
    return(result)
  } else {
    warning(paste("No significant LR pairs found between", source, "and", target, "with threshold", prob_threshold))
    return(data.frame())
  }
}

# 提取双向通信对
bm_LR <- extract_detailed_LR_pairs(cellchat_bm, "Basal", "Myofibroblast", 0.01)
mb_LR <- extract_detailed_LR_pairs(cellchat_bm, "Myofibroblast", "Basal", 0.01)

print("Basal -> Myofibroblast 显著受体-配体对:")
if(nrow(bm_LR) > 0) {
  print(head(bm_LR, 15))
} else {
  print("无显著对，尝试降低阈值...")
  bm_LR <- extract_detailed_LR_pairs(cellchat_bm, "Basal", "Myofibroblast", 0.001)
  print(head(bm_LR, 10))
}

print("Myofibroblast -> Basal 显著受体-配体对:")
if(nrow(mb_LR) > 0) {
  print(head(mb_LR, 15))
} else {
  print("无显著对，尝试降低阈值...")
  mb_LR <- extract_detailed_LR_pairs(cellchat_bm, "Myofibroblast", "Basal", 0.001)
  print(head(mb_LR, 10))
}


#可视化气泡图
library(ggplot2)
library(dplyr)
library(tidyr)

# 可靠方法：手动提取数据创建自定义双向通信图
create_custom_bidirectional_plot <- function(cellchat_obj) {
  # 提取通信概率矩阵
  prob_matrix <- cellchat_obj@net$prob
  
  # 提取两个方向的通信概率
  bm_probs <- prob_matrix["Basal", "Myofibroblast", ]
  mb_probs <- prob_matrix["Myofibroblast", "Basal", ]
  
  # 获取所有相互作用名称
  all_interactions <- unique(c(names(bm_probs), names(mb_probs)))
  
  # 创建数据框
  plot_data <- data.frame(
    interaction_name = all_interactions,
    Basal_to_Myofibro = ifelse(all_interactions %in% names(bm_probs), 
                              bm_probs[all_interactions], 0),
    Myofibro_to_Basal = ifelse(all_interactions %in% names(mb_probs), 
                              mb_probs[all_interactions], 0)
  ) %>%
    # 筛选显著相互作用（概率 > 0.01）
    filter(Basal_to_Myofibro > 0.01 | Myofibro_to_Basal > 0.01) %>%
    # 计算总概率用于排序
    mutate(total_prob = Basal_to_Myofibro + Myofibro_to_Basal) %>%
    arrange(desc(total_prob)) %>%
    # 选择前30个最重要的相互作用
    head(30) %>%
    # 转换为长格式
    pivot_longer(cols = c("Basal_to_Myofibro", "Myofibro_to_Basal"),
                names_to = "direction",
                values_to = "probability") %>%
    filter(probability > 0)  # 只保留概率大于0的
  
  # 创建气泡图
  p <- ggplot(plot_data, aes(x = direction, y = reorder(interaction_name, total_prob))) +
    geom_point(aes(size = probability, color = probability), alpha = 0.8) +
    scale_size_continuous(range = c(2, 10), 
                         breaks = c(0.02, 0.05, 0.1, 0.2),
                         name = "Communication\nProbability") +
    scale_color_viridis_c(name = "Communication\nProbability",
                         breaks = c(0.02, 0.05, 0.1, 0.2)) +
    scale_x_discrete(labels = c("Basal_to_Myofibro" = "Basal → Myofibro", 
                               "Myofibro_to_Basal" = "Myofibro → Basal")) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 9),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    labs(title = "Basal-Myofibroblast双向配体-受体通信",
         x = "Communication Direction", 
         y = "Ligand-Receptor Pair")
  
  return(p)
}

# 生成并保存图形
custom_bidirectional_plot <- create_custom_bidirectional_plot(cellchat_bm)
print(custom_bidirectional_plot)

# 保存高质量图片
ggsave("custom_bidirectional_communication.pdf", 
       plot = custom_bidirectional_plot, 
       width = 10, height = 12, dpi = 300, bg = "white")
