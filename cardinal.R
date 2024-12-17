library("Cardinal")
library("BiocParallel")
register(SerialParam()) 
register(SnowParam(workers = 24))
nedc_cmc=readMSIData("D://SKLxXia//Spatal_Metabolomics//20240608_best_plant msi//Treatment//leave_treatment_1-total ion count.imzML")
cmc_pre <- process(nedc_cmc)
cmc_peaks <- peakPick(cmc_pre,method="filter",SNR=1)
cmc_all_peaks <- peakAlign(cmc_peaks,tolerance=200, units="ppm")
cmc_peaksf=fData(cmc_all_peaks)
peaks_cmcfilt<- as.data.frame(cmc_peaksf)
write.table(peaks_cmcfilt, file = "D://SKLxXia//Spatal_Metabolomics//20240608_best_plant msi//Treatment//cmc_pick.csv", sep = "\t", row.names = FALSE, col.names = TRUE)
cmc_pca=PCA(cmc_all_peaks,ncomp=10)
cmc_nmf=NMF(cmc_all_peaks,ncomp=10)
# 提取 model slot 以获取 activation 矩阵
model_data <- NMF_test@model
activation_matrix <- model_data$activation

# 提取 featureData slot
feature_data <- NMF_test@featureData

# 提取峰大小（例如 mz 列）
peak_sizes <- feature_data@listData$mz

# 检查提取的峰大小数组
head(peak_sizes)

# 创建最终的数据框，其中包括峰大小和激活矩阵
final_data <- data.frame(mz = peak_sizes, activation_matrix)

# 检查合并结果
head(final_data)

# 导出最终数据到 CSV 文件
write.csv(final_data, "activation_with_peak_sizes.csv", row.names = FALSE)



png("nmf9.png", width = 10, height = 8, units = "in", res = 300)

plot(NMF_test,
      style="dark",
      grid=FALSE,
      scale=TRUE,
      smooth="gaussian",
     select=c(9),
      layout=c(1,1))
dev.off()



plot(cmc_all_peaks,linewidth=1, annPeaks=10)



plot(leaf_m,xlim=c(316.5,317.1),ylim=c(0,1))
plot(leaf_pre,xlim=c(338,340),ylim=c(0,1))

png("other_peaks_wt.png", width = 18, height = 5, units = "in", res = 300)

image(leaf_pre, mz=c(243.0582

), tolerance=0.15,layout=c(1,6),
      units=c("mz"), style="dark",grid=FALSE,scale=TRUE,smooth="gaussian")


dev.off()

image(leaf_m, mz=c(339.0277

), tolerance=0.15,layout=c(1,1),
units=c("mz"), style="dark",grid=FALSE,scale=TRUE,smooth="gaussian")


image(leaf_pre, mz=c(191.019734

), tolerance=0.15,layout=c(1,1),
units=c("mz"), style="dark",grid=FALSE,scale=TRUE,smooth="gaussian")




plot(cmc_pca,smooth="adaptive", superpose=FALSE,enhance="histogram",type="rotation", layout=c(3,4), linewidth=0.5)
png("mutant_nmf.png", width = 20, height = 16, units = "in", res = 300)
image(NMF_test,
      style="dark",
      grid=FALSE,
      scale=TRUE,
      smooth="gaussian", 
      superpose=FALSE,
      
      layout=c(3,4))
dev.off()


png("mutant_nmfmerge2,4,5,7.png", width = 3, height = 4, units = "in", res = 300)
image(NMF_test,
      
      grid=FALSE,
      scale=TRUE,
      smooth="gaussian", 
      select=c(2,4,5,7),
      layout=c(1,1))
dev.off()



png("mutant_nmf.png", width = 16, height = 12, units = "in", res = 600)

image(wt_nmf,
      style="dark",
      grid=FALSE,
      scale=TRUE,
      smooth="gaussian", 
      superpose=FALSE,
      layout=c(3,4))
dev.off()

a <- 286.0482871

b <- paste(a, "png",sep = ".")
png(b, width = 6,height = 8,units = "in",res = 600)

image(leafmutant_peaks, mz=a

                             
  ,tolerance=0.15,
  layout=c(1,1),
  units=c("mz"), 
  style="dark",grid=FALSE,scale=TRUE,smooth="gaussian")

dev.off()


a <- 164.071702594
b <- paste(a, "png",sep = ".")
png(b, width = 6,height = 8,units = "in",res = 600)

image(leaf_peaks, mz=a
      ,tolerance=0.15,
      layout=c(1,1),
      units=c("mz"), 
      style="dark",grid=FALSE,scale=TRUE,smooth="gaussian")

dev.off()













png("nmf.png", width = 16, height = 12, units = "in", res = 600)


plot(wt_nmf,style="dark",grid=FALSE)

dev.off()


png("mutant_nmf.png", width = 6, height = 8, units = "in", res = 300)

image(mutant_nmf,
      style="dark",
      grid=FALSE,
      scale=TRUE,
      smooth="gaussian",select=c(2,4,5,7,10))

dev.off()






k <- 10  # 您可以根据需要调整簇的数量
kmeans_result <- kmeans(umap_data, centers = k)

# 2. 将聚类结果添加到数据框中
umap_data$Cluster <- as.factor(kmeans_result$cluster)

# 3. 为每个簇计算凸包
hulls <- lapply(1:k, function(cluster) {
  cluster_points <- umap_data[umap_data$Cluster == cluster, c("UMAP1", "UMAP2", "UMAP3")]
  if(nrow(cluster_points) < 4) return(NULL)  # 如果点太少，跳过
  ch <- chull(cluster_points)
  return(cluster_points[c(ch, ch[1]), ])  # 闭合凸包
})

# 4. 创建3D UMAP图，包括簇和凸包
umap_wt <- plot_ly()

# 添加散点
umap_wt <- umap_wt %>% add_trace(
  data = umap_data,
  x = ~UMAP1, y = ~UMAP2, z = ~UMAP3,
  color = ~Cluster,
  colors = c("purple","#912C2C","#F2BB6B","#C2ABC8","#329845","#AED185","#276C9E","#A3C9D5"),
  type = "scatter3d", mode = "markers",
  marker = list(size = 0.5, opacity = 0.8)
)

# 添加凸包
for(i in 1:length(hulls)) {
  if(!is.null(hulls[[i]])) {
    umap_wt <- umap_wt %>% add_trace(
      x = hulls[[i]]$UMAP1, 
      y = hulls[[i]]$UMAP2, 
      z = hulls[[i]]$UMAP3,
      type = "mesh3d",
      opacity = 0.2,
      color = i,
      hoverinfo = "none"
    )
    
    # 添加簇的数字标注
    cluster_center <- colMeans(hulls[[i]][, c("UMAP1", "UMAP2", "UMAP3")])
    umap_wt <- umap_wt %>% add_trace(
      x = cluster_center[1],
      y = cluster_center[2],
      z = cluster_center[3],
      type = "scatter3d",
      mode = "text",
      text = as.character(i),
      textposition = "middle center",
      textfont = list(size = 20, color = "black"),
      showlegend = FALSE
    )
  }
}

# 设置布局
umap_wt <- umap_wt %>% layout(
  scene = list(
    xaxis = list(title = 'UMAP1'),
    yaxis = list(title = 'UMAP2'),
    zaxis = list(title = 'UMAP3')
  ),
  title = "UMAP 3D Projection with Clusters",
  legend = list(
    font = list(size = 30),  # 大幅增加字体大小
    itemsizing = "constant",
    itemwidth = 50,  # 增加图例项的宽度
    itemheight = 30  # 增加图例项的高度
)
)


# 保存为HTML文件
htmlwidgets::saveWidget(umap_wt, "umap_3d_with_clustersgai.html")





selected_clusters <- c(1, 5, 6,7)

# 2. 创建一个新的数据框，只包含选定的簇
selected_data <- umap_data[kmeans_result$cluster %in% selected_clusters, ]

# 3. 获取原始图像的坐标（如果还没有这个信息）
coord_data <- coord(leafmutant_peaks)

# 4. 创建包含原始空间坐标和簇信息的数据框
image_data <- data.frame(
  x = coord_data$x,
  y = coord_data$y,
  cluster = factor(kmeans_result$cluster)
)

# 5. 只保留选定的簇
image_data_selected <- image_data[image_data$cluster %in% selected_clusters, ]

# 6. 创建颜色映射
color_mapping <- c("1" = "purple", "5" = "#329845", "6" = "#276C9E", "7" = "#A3C9D5")

# 7. 创建UMAP映射图，只显示选定的簇
r

复制
# 假设 kmeans_result 已经计算完成，umap_data 包含了UMAP的结果

# 1. 选择指定的簇
selected_clusters <- c(2)

# 2. 创建一个新的数据框，只包含选定的簇
selected_data <- umap_data[kmeans_result$cluster %in% selected_clusters, ]

# 3. 获取原始图像的坐标（如果还没有这个信息）
coord_data <- coord(leafmutant_peaks)

# 4. 创建包含原始空间坐标和簇信息的数据框
image_data <- data.frame(
  x = coord_data$x,
  y = coord_data$y,
  cluster = factor(kmeans_result$cluster)
)

# 5. 只保留选定的簇
image_data_selected <- image_data[image_data$cluster %in% selected_clusters, ]

# 6. 创建颜色映射
color_mapping <- c("1" = "purple", "5" = "#4DAF4A", "7" = "#377EB8")

# 7. 创建UMAP映射图，只显示选定的簇，并减小点的大小
umap_map <- ggplot(image_data_selected, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.01, alpha = 0.7) +  # 减小点的大小到0.1，并添加一些透明度
  scale_color_manual(values = color_mapping) +
  theme_minimal() +
  ggtitle("Selected UMAP Clusters") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_fixed(ratio = 1) +
  labs(color = "Cluster", x = "X", y = "Y")

# 8. 显示图形
print(umap_map)

# 9. 保存图形
ggsave("umap_selected.png", umap_map, width = 10, height = 8, dpi = 600)




# 9. 保存图形
ggsave("umap_selected.png", umap_map, width = 10, height = 8, dpi = 600)


"purple","#912C2C","#F2BB6B","#C2ABC8","#329845","#AED185","#276C9E","#A3C9D5"
,"#276C9E","#AED185","#912C2C","#F2BB6B","purple"









# 1. 选择指定的簇
selected_clusters <- c(10)

# 2. 创建一个新的数据框，只包含选定的簇
selected_data <- umap_data[kmeans_result$cluster %in% selected_clusters, ]

# 3. 获取原始图像的坐标（如果还没有这个信息）
coord_data <- coord(leafmutant_peaks)

# 4. 创建包含原始空间坐标和簇信息的数据框
image_data <- data.frame(
  x = coord_data$x,
  y = coord_data$y,
  cluster = factor(kmeans_result$cluster)
)

# 5. 只保留选定的簇
image_data_selected <- image_data[image_data$cluster %in% selected_clusters, ]

# 6. 创建颜色映射
color_mapping <- c("10" = "purple", "5" = "#4DAF4A", "7" = "#377EB8")

# 7. 创建UMAP映射图，只显示选定的簇，并减小点的大小
umap_map <- ggplot(image_data_selected, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.01, alpha = 0.7) +  # 减小点的大小到0.1，并添加一些透明度
  scale_color_manual(values = color_mapping) +
  theme_minimal() +
  ggtitle("Selected UMAP Clusters") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  coord_fixed(ratio = 1) +
  labs(color = "Cluster", x = "X", y = "Y")

# 8. 显示图形
print(umap_map)

library(ggplot2)
library(viridis)
library(scales)
library(plotly)
library(Cardinal)  # 假设你使用 Cardinal 包处理 MS imaging 数据
library(uwot)  # 确保已加载 uwot 包
library(Rtsne) 

load("/home/bingxu/spatial_metabolomics/data/wt_h/wt_h.RData")
setwd("/home/bingxu/spatial_metabolomics/data/wt_h/umap_tsne/")


intensity_data <- spectra(leafwt_peaks)
intensity_matrix <- t(as.matrix(intensity_data))  # 转置矩阵，使每行代表一个像素
# 2. 预处理数据
# 标准化每个特征
intensity_matrix_scaled <- scale(intensity_matrix)

umap_result <- umap(intensity_matrix_scaled, 
                    n_neighbors = 50, 
                    n_components = 3, 
                    metric = "euclidean", 
                    n_epochs = 200, 
                    learning_rate = 1, 
                    min_dist = 0.1, 
                    spread = 1, 
                    n_threads = 40)

# 只提取前 3 个维度用于三维绘图
umap_data <- as.data.frame(umap_result[, 1:3])
colnames(umap_data) <- c("UMAP1", "UMAP2", "UMAP3")

# 创建三维 UMAP 图，并调整点的大小
umap_wt <- plot_ly(umap_data, 
             x = ~UMAP1, 
             y = ~UMAP2, 
             z = ~UMAP3, 
             color = ~UMAP3, 
             colors = c("purple","#912C2C","#F2BB6B","#C2ABC8","#329845","#AED185","#276C9E","#A3C9D5"), 
             marker = list(size = 0.5)) %>%  # 调整点的大小
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'UMAP1'),
                      yaxis = list(title = 'UMAP2'),
                      zaxis = list(title = 'UMAP3')),
         title = "UMAP 3D Projection")

# 保存为 HTML 文件
htmlwidgets::saveWidget(umap_wt, "umap_3d_projectiongai.html")



# 1. 创建颜色映射函数
color_map <- function(values, colors) {
  scaled <- scales::rescale(values)
  colorRamp(colors)(scaled)
}

# 2. 为 UMAP 创建颜色 (使用第3维)
umap_colors <- color_map(umap_result[,3], c("purple","#912C2C","#F2BB6B","#C2ABC8","#329845","#AED185","#276C9E","#A3C9D5"))
umap_colors_hex <- rgb(umap_colors, maxColorValue = 255)


# 4. 获取原始图像的坐标
coord_data <- coord(leafwt_peaks)

# 5. 创建包含原始空间坐标和颜色的数据框
image_data <- data.frame(
  x = coord_data$x,
  y = coord_data$y,
  umap_color = umap_colors_hex
)

# 6. 创建 UMAP 映射图
umap_map <- ggplot(image_data, aes(x = x, y = y, color = umap_color)) +
  geom_point(size = 0.5) +
  scale_color_identity() +
  theme_minimal() +
  ggtitle("UMAP Mapped to Original Image Space") +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1)



png("umap.png", width = 8, height = 12, units = "in", res = 600)

ggplot(image_data, aes(x = x, y = y, color = umap_color)) +
  geom_point(size = 0.5) +
  scale_color_identity() +
  theme_minimal() +
  ggtitle("UMAP Mapped to Original Image Space") +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1)
  
dev.off()











tsne_result <- Rtsne(intensity_matrix_scaled, dims = 3, perplexity = 50, verbose = TRUE, max_iter = 300)

# 将结果转换为数据框
tsne_data <- as.data.frame(tsne_result$Y)
colnames(tsne_data) <- c("TSNE1", "TSNE2", "TSNE3")

# 创建三维 t-SNE 图，并调整点的大小
tsne_wt <- plot_ly(tsne_data, 
             x = ~TSNE1, 
             y = ~TSNE2, 
             z = ~TSNE3, 
             color = ~TSNE3, 
             colors = c("#912C2C","#F2BB6B","#C2ABC8","#329845","#AED185","#276C9E","#A3C9D5","purple"), 
             marker = list(size = 1)) %>%  # 调整点的大小
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'TSNE1'),
                      yaxis = list(title = 'TSNE2'),
                      zaxis = list(title = 'TSNE3')),
         title = "t-SNE 3D Projection")

# 保存为 HTML 文件
htmlwidgets::saveWidget(tsne_wt, "tsne_3dgai_projection.html")



# 假设 umap_result 和 tsne_result 已经计算完成

# 1. 创建颜色映射函数
color_map <- function(values, colors) {
  scaled <- scales::rescale(values)
  colorRamp(colors)(scaled)
}


# 3. 为 t-SNE 创建颜色 (使用第3维)
tsne_colors <- color_map(tsne_result$Y[,3], c("#912C2C","#F2BB6B","#C2ABC8","#329845","#AED185","#276C9E","#A3C9D5","purple"))
tsne_colors_hex <- rgb(tsne_colors, maxColorValue = 255)

# 4. 获取原始图像的坐标
coord_data <- coord(leafwt_peaks)

# 5. 创建包含原始空间坐标和颜色的数据框
image_data <- data.frame(
  x = coord_data$x,
  y = coord_data$y,
  tsne_color = tsne_colors_hex
)


# 7. 创建 t-SNE 映射图
tsne_map <- ggplot(image_data, aes(x = x, y = y, color = tsne_color)) +
  geom_point(size = 0.01,alpha = 0.7) +
  scale_color_identity() +
  theme_minimal() +
  ggtitle("t-SNE Mapped to Original Image Space") +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1)

save.image("wt_.RData")

png("sne.png", width = 8, height = 12, units = "in", res = 600)

ggplot(image_data, aes(x = x, y = y, color = tsne_color)) +
  geom_point(size = 0.01,alpha = 0.7) +
  scale_color_identity() +
  theme_minimal() +
  ggtitle("t-SNE Mapped to Original Image Space") +
  theme(legend.position = "none") +
  coord_fixed(ratio = 1)
  
dev.off()




