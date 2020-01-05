library(DESeq2)
####
#### Data cleaning
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

## We tested 6 different models containing different different covariates. 
## Some of the combinations are not tested here because the related matrices
## does not reach full rank.
## The models we would test:
deseq2_deseqds1 <- DESeqDataSetFromMatrix(countData = EXP, colData = GRP, design = ~ Group + Age + Sex + Group * Age)
deseq2_deseqds2 <- DESeqDataSetFromMatrix(countData = EXP, colData = GRP, design = ~ Group + Age + Sex + Sex * Age + Group * Age)

deseq2_deseqds1r <- DESeqDataSetFromMatrix(countData = EXP, colData = GRP, design = ~ Group + Age + Sex + Group * Age + RIN)
deseq2_deseqds2r <- DESeqDataSetFromMatrix(countData = EXP, colData = GRP, design = ~ Group + Age + Sex + Sex * Age + Group * Age + RIN)

## Apply rlog transformation.
deseq2_rld1 <- rlog(deseq2_deseqds1)
deseq2_rld2 <- rlog(deseq2_deseqds2)

deseq2_rld1r <- rlog(deseq2_deseqds1r)
deseq2_rld2r <- rlog(deseq2_deseqds2r)

## Apply variance stabilizing transformation (VST) transformation:
deseq2_vst1 <- varianceStabilizingTransformation(deseq2_deseqds1, blind = TRUE, fitType = "parametric")
deseq2_vst2 <- varianceStabilizingTransformation(deseq2_deseqds2, blind = TRUE, fitType = "parametric")
deseq2_vst1r <- varianceStabilizingTransformation(deseq2_deseqds1r, blind = TRUE, fitType = "parametric")
deseq2_vst2r <- varianceStabilizingTransformation(deseq2_deseqds2r, blind = TRUE, fitType = "parametric")


## Cluster the tested dogs. All models generate the same clusters. 
## There appear to be three clusters. The third female clusters by itself, and 
## the first male clusters by itself. Referring to Table S3 from the supplementary 
## information of this publication, there are some aspects of these two samples 
## that are distinct from the others.
deseq2_hc1 <- hclust(dist(t(assay(deseq2_rld1))))

deseq2_hc1_vst <- hclust(dist(t(assay(deseq2_vst1))))

par(mfrow = c(1, 2))
plot(deseq2_hc1)
plot(deseq2_hc1$height)


par(mfrow = c(1, 2))
plot(deseq2_hc1_vst)
plot(deseq2_hc1_vst$height)

## show the dogs in each of the clusters:
GRP$Cluster <- cluster0 <- cutree(deseq2_hc1, k=3) 
c10 <- EXP[, cluster0 == 1] 
c20 <- EXP[, cluster0 == 2] 
c30 <- EXP[, cluster0 == 3] 
colnames(c10)
colnames(c20)
colnames(c30)

GRP1 <- GRP[cluster0 == 1,] 
GRP2 <- GRP[cluster0 == 2,] 
GRP3 <- GRP[cluster0 == 3,] 

## PCA plot and analysis. Using DESeq2's plotPCA, it is clear that two samples are 
## different from the others. Running the PCA ourselves and looking at the loadings, we 
## can see, as expected, that these are the two samples mentioned above.
plotPCA(deseq2_rld1, ntop = nrow(assay(deseq2_rld1)), intgroup = "Age")
plotPCA(deseq2_rld1, ntop = nrow(assay(deseq2_rld1)), intgroup = "Sex")
PC1 <- plotPCA(deseq2_rld1, ntop = nrow(assay(deseq2_rld1)), intgroup = "Group", return = TRUE)

plotPCA(deseq2_vst1, ntop = nrow(assay(deseq2_vst1)), intgroup = "Age")
plotPCA(deseq2_vst1, ntop = nrow(assay(deseq2_vst1)), intgroup = "Sex")
PC1vst <- plotPCA(deseq2_rld1, ntop = nrow(assay(deseq2_vst1)), intgroup = "Group", return = TRUE)

## Show the percentage of each of the principle components:
deseq2_svd1 <- svd(t(scale(t(assay(deseq2_rld1)), center = TRUE, scale = FALSE)))
round(deseq2_svd1$d ^ 2 / sum(deseq2_svd1$d ^ 2), 2)

## Differential expression analysis. 
deseq2_de1 <- DESeq(deseq2_deseqds1)
deseq2_de_out1 <- results(deseq2_de1)
deseq2_de2 <- DESeq(deseq2_deseqds2)
deseq2_de_out2 <- results(deseq2_de2)

deseq2_de1r <- DESeq(deseq2_deseqds1r)
deseq2_de_out1r <- results(deseq2_de1r)
deseq2_de2r <- DESeq(deseq2_deseqds2r)
deseq2_de_out2r <- results(deseq2_de2r)

summary(deseq2_de_out1)

## Show the miRNAs with an estimated FDR < 0.05 and FDR < 0.1.
deseq2_de_out_sig_0051 <- subset(deseq2_de_out1, padj < 0.05)
(nms_deseq2_sig_0051 <- rownames(deseq2_de_out_sig_0051))

deseq2_de_out_sig_011 <- subset(deseq2_de_out1, padj < 0.1)
(nms_deseq2_sig_011 <- rownames(deseq2_de_out_sig_011))

deseq2_de_out_sig_0052 <- subset(deseq2_de_out2, padj < 0.05)
(nms_deseq2_sig_0052 <- rownames(deseq2_de_out_sig_0052))

deseq2_de_out_sig_012 <- subset(deseq2_de_out2, padj < 0.1)
(nms_deseq2_sig_012 <- rownames(deseq2_de_out_sig_012))


deseq2_de_out_sig_0051r <- subset(deseq2_de_out1r, padj < 0.05)
(nms_deseq2_sig_0051r <- rownames(deseq2_de_out_sig_0051r))

deseq2_de_out_sig_011r <- subset(deseq2_de_out1r, padj < 0.1)
(nms_deseq2_sig_011r <- rownames(deseq2_de_out_sig_011r))

deseq2_de_out_sig_0052r <- subset(deseq2_de_out2r, padj < 0.05)
(nms_deseq2_sig_0052r <- rownames(deseq2_de_out_sig_0052r))

deseq2_de_out_sig_012r <- subset(deseq2_de_out2r, padj < 0.1)
(nms_deseq2_sig_012r <- rownames(deseq2_de_out_sig_012r))

###
### test stattistic and FDR of selected miRNAs
###
deseq2_de_out_sig_011
