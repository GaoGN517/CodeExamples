####
#### Validation with EdgeR
####
library(edgeR)
####
#### Data 
####

EXP <- readRDS("Expression.rds")
GRP <- readRDS("GroupInfor.rds")
head(EXP)
head(GRP)
## There are many genes that had no counts in any of the samples. Remove them 
## from the analysis.
EXP <- EXP[rowSums(EXP) > 0, ]
rowSums(EXP)
colSums(EXP)

## Create a DGElist object for holding the counts along with the library sizes (the total 
## number of counts) and the group membership information for each sample. 
edgeR_dge <- DGEList(EXP, lib.size = colSums(EXP), group = GRP[, 2])
edgeR_dge

mydesign1<-model.matrix(~Age+Group+Sex+Age*Group+Sex*Age + RIN, data=GRP)
D1 <- estimateGLMCommonDisp(edgeR_dge, mydesign1) 

# D is the DGEList of the normalized counts
#D1 <- estimateGLMTrendedDisp(D1, mydesign1)
fit1 <- glmQLFit(D1, mydesign1)
results<-glmQLFTest(fit1,coef = 5)
edgeR_out <- topTags(results, n = nrow(EXP), sort.by = "PValue")
edgeR_out_sig_01 <- subset(edgeR_out$table, FDR < 0.1)
edgeR_out_sig_005 <- subset(edgeR_out$table, FDR < 0.05)

nms_edgeR_sig <- rownames(edgeR_out_sig_01)
rownames(EXP[edgeR_out_sig_01, ])

nms_edgeR_sig_n <- rownames(edgeR_out_sig_005)
rownames(EXP[edgeR_out_sig_005, ])

fit2 <- glmFit(D1, mydesign1)
