
# Load expression data
exprmatrx1 <- load(Filepath)

## Creates a DGEList object from a table of counts (rows=features, columns=samples), group indicator for each column,
## library size (optional) and a table of feature annotation (optional).
dge <- edgeR::DGEList(counts = assays(exprMatrix)$counts, samples = colData(exprMatrix))
## Calculate normalization factors to scale the raw library sizes. 
dge <- edgeR::calcNormFactors(dge)

## Generate design matrix
#colnames(dge$samples)[length(colnames(dge$samples))] <- "Group"
formula1 <- ~ 0 + Group
design  <- model.matrix(formula1, data = dge$samples)

# Fit the model
## Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to 
## compute appropriate observational-level weights. The data are then ready for linear modelling.
voom    <- limma::voom(dge, design, plot = FALSE)
fit     <- limma::lmFit(voom, design)


# Fit contrast
contrast.matrix <- limma::makeContrasts(
  "Treatment1_vs_Control" = Treatment1-Control,
  "Treatment2_vs_Treatment1" = Treatment2-Treatment1
  )

fit <- limma::contrasts.fit(fit, contrast.matrix)

# Bayesian correction
fit <- limma::eBayes(fit)
