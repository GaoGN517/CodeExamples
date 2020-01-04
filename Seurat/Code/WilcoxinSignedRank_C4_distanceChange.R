library(Seurat)

## Data preparation
C4 <- readRDS("Liu_Wu_scVI_C4.rds")
C4_cells <- C4@assays$wu_liu_scVI_imputation@data@Dimnames[[2]]

Wu_orig <- readRDS("Wu_obj.RDS")
Liu_orig <- readRDS("Liu_obj.RDS")

Wu_orig_C4 <- subset(x = Wu_orig, cells = C4_cells)
Liu_orig_C4 <- subset(x = Liu_orig, cells = C4_cells)

## From Seurat objects to get the data matrices
C4.data <- GetAssayData(C4, slot = "data")
Wu_orig_C4.data <- GetAssayData(Wu_orig_C4, slot = "data")
Liu_orig_C4.data <- GetAssayData(Liu_orig_C4, slot = "data")

## Calculate the distance and store for future use.
## Wu_orig_C4.dist <- dist(t(as.matrix(Wu_orig_C4.data)))


## load the previously stored distance matrices. 
C4.dist <- readRDS("Euclidean_C4_Liu_Wu.rds")
Wu_orig_C4.dist <- readRDS("Euclidean_C4_Wu_Orig.rds")
Liu_orig_C4.dist <- readRDS("Euclidean_C4_Liu_Orig.rds")

## Make sure that the distance are stored in matrices. 
C4.dist.matrix <- as.matrix(C4.dist)

Wu_orig_C4.dist.matrix <- as.matrix(Wu_orig_C4.dist)
Liu_orig_C4.dist.matrix <- as.matrix(Liu_orig_C4.dist)

Wu_C4.dist.matrix <- C4.dist.matrix[1:992, 1:992]
Liu_C4.dist.matrix <- C4.dist.matrix[993:1863, 993:1863]

Liu_orig_C4.dist.matrix.scale <- Liu_orig_C4.dist.matrix / mean(rowMeans(Liu_orig_C4.dist.matrix)) 
Liu_C4.dist.matrix.scale <- Liu_C4.dist.matrix / mean(rowMeans(Liu_C4.dist.matrix))

## Wilcoxin Signed Rank Tests
## Liu
Liu_C4.matrix <- rbind(Liu_orig_C4.dist.matrix, Liu_C4.dist.matrix)
system.time(
  Liu_C4.WilcoxonSignedRank <- apply(Liu_C4.matrix, 2, function(X) {
    n <- length(X)
    X_orig <- X[1:(n/2)]
    X_Sim <- X[(n/2+1):n]
    wilcox.test(X_orig, X_Sim, paired = T)$p.value
  })
)
plot(Liu_C4.WilcoxonSignedRank)
plot(p.adjust(Liu_C4.WilcoxonSignedRank))


## Wu
Wu_C4.matrix <- rbind(Wu_orig_C4.dist.matrix, Wu_C4.dist.matrix)

system.time(
  Wu_C4.WilcoxonSignedRank <- apply(Wu_C4.matrix, 2, function(X) {
    n <- length(X)
    X_orig <- X[1:(n/2)]
    X_Sim <- X[(n/2+1):n]
    wilcox.test(X_orig, X_Sim, paired = T)$p.value
  })
)
plot(Wu_C4.WilcoxonSignedRank)
plot(p.adjust(Wu_C4.WilcoxonSignedRank))

## log transform the p-values. 
Liu_C4.WilcoxonSignedRank.log <- log10(Liu_C4.WilcoxonSignedRank)
Wu_C4.WilcoxonSignedRank.log <- log10(Wu_C4.WilcoxonSignedRank)

## Get the cell total expression values
Liu_orig_C4.data.cellExpression <- colSums(as.matrix(Liu_orig_C4.data))
Wu_orig_C4.data.cellExpression <- colSums(as.matrix(Wu_orig_C4.data))
C4.data.cellExpression <- colSums(as.matrix(C4.data))


## Plot the -log_10(pval) vs cellExpression 
plot(log10(Wu_orig_C4.data.cellExpression), -Wu_C4.WilcoxonSignedRank.log,
     xlab = "log(Expression Level)", ylab = "-log(Wilcoxon Signed Rank Test p-value)",
     main = "")
plot(log10(Liu_orig_C4.data.cellExpression), -Liu_C4.WilcoxonSignedRank.log)
