library(limma)

## Data
EXP <- readRDS("Expression.rds")
GRP <- readRDS("GroupInfor.rds")
head(EXP)
head(GRP)
## There are many genes that had no counts in any of the samples. Remove them 
## from the analysis.
EXP <- EXP[rowSums(EXP) > 0, ]
rowSums(EXP)
colSums(EXP)


## Design matrix.
#limma_dge <- DGEList(DTA, lib.size = colSums(DTA),group=GRP[,2])
limma_X<-model.matrix(~ Age + Group + Sex + Age * Group + RIN,data = GRP)

## Voom transformation and differential expression analysis.
limma_voom <- voom(EXP, design = limma_X)
limma_fit <- lmFit(EXP, design = limma_X)
limma_eb <- eBayes(limma_fit)
## For some reason, no estimated FDRs fall below 0.95.
limma_FDRs <- topTable(limma_eb, n = nrow(EXP), coef = c(1,5))
# Display rownames with FDR (condition)
rownames(limma_FDRs[limma_FDRs$adj.P.Value<0.1,])

hist(limma_FDRs$P.Value)
quantile(limma_FDRs$adj.P.Val)

