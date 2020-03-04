if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")

## Load the library
library("biomaRt")

## Look at the list of Marts
listMarts()

## use the mart "ensembl"
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)

## Only look at human genes
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
datasets <- listDatasets(ensembl)
head(datasets)

## In the dataset there are only 1:22 chromatin genes, so we only look at genes on these chromosomes.
## This restriction can be anything. 
affyids=as.character(1:22)

## Select gene id, gene name, chromatin, starting site, ending site from the data source. 
## Using chromasome name as a filter. And the filtered values are stored in affyids.
Genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name','chromosome_name',
                              'start_position', 'end_position'),
               filters = 'chromosome_name', 
               values = affyids, 
               mart = ensembl)

## Check whether the start and end points could be reverse
sum((Genes$End_10000 - Genes$Start_10000) < 0) # no

## We usually look at a larger area than just the gene area. TSS, TES
Genes$Start_10000 <- Genes$start_position - 10000
Genes$End_10000 <- Genes$end_position - 10000

## Store the gene dataframe. 
saveRDS(Genes, "EnsembleHumanGenes1_22Chr.rds")
