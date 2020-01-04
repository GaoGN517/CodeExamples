library(Seurat)

## load the cluster four cell distance matrix
C4 <- readRDS("Liu_Wu_scVI_C4.rds")
C4.dist <- readRDS("Euclidean_C4_Liu_Wu.rds")
C4.dist <- as.matrix(C4.dist)
dim(C4.dist)
## 1:992 are Wu cells, 993:1863 are Liu cells. 
## we can divide the table into 1:992 rows and 993:1863 rows. 
Wu_C4_to_all.dist <- C4.dist[1:992, ]
Liu_C4_to_all.dist <- C4.dist[993:1863, ]
dim(Wu_C4_to_all.dist)
Wu_C4_cluster_change <- apply(Wu_C4_to_all.dist, 1, function(x){
  x_ordered <- order(x, decreasing = F)
  x_first <- x_ordered[1]
  x_second <- x_ordered[2]
  ifelse((x_first <= 992) & (x_second <= 992), 1, 0)
})
table(Wu_C4_cluster_change)
# 0   1 
# 72 920 
920/992*100 # 92.74 %

Liu_C4_cluster_change <- apply(Liu_C4_to_all.dist, 1, function(x){
  x_ordered <- order(x, decreasing = F)
  x_first <- x_ordered[1]
  x_second <- x_ordered[2]
  ifelse((x_first <= 992) | (x_second <= 992), 0, 1)
})
table(Liu_C4_cluster_change)
# 0   1 
# 71 800 
800/887*100 # 90.19 %

#apply(Liu_C4_to_all.dist[1:10,], 1, which.min)
