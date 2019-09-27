#BiocManager::install("biomaRt")
library(biomaRt)
#listMarts() # check the list of current available marts (the filters that can be used)

traget_genes <- as.character(unique(read.csv("traget_genes.csv")$gene_name)) #129
#hyper$external_gene_name <- toupper(hyper$external_gene_name)
hypo <- as.character(unique(read.csv("Osr1.hypo20p_5kb_p.05.csv")$external_gene_name)) #101
#################################################################################
## Using getBM to download Go base data. 
## Problem, the downloaded dataset is much smaller than expected.
mart<- useMart(biomart="ensembl", dataset = "mmusculus_gene_ensembl")

inflammation <- unique(getBM(attributes=c('ensembl_gene_id', 'external_gene_name','go_id'),
                   filters = 'go', values = "GO:0006954", mart = mart)$external_gene_name)#321

## http://www.informatics.jax.org/vocab/gene_ontology/GO:0008152


inflammation_target <- intersect(traget_genes, inflammation)


merge_target <- unique(Reduce(append, list(inflammation_target, category2_target, 
                                          category3_target, ...)))
other_target <- setdiff(traget_genes, merge_target) 

################################################################################################
## Instead, directly download the dataset. (Find the uppermost target in the category you want. )
inflammatory_table <- unique(read.table("GO_inflammatory_20190913_173834.txt", header = F, fill = T)[, 1:2]$V2) 
inflammatory_target_table <- intersect(traget_genes, inflammatory_table)


merge_target_table <- unique(Reduce(append, list(inflammation_target_table, category2_target_table, 
                                           category3_target_table, ...)))
other_target_table <- setdiff(traget_genes, merge_target_table) 

##################################################################################
## Add to Workbooks
library(xlsx)
wb = createWorkbook()

sheet = createSheet(wb, "target")
addDataFrame(inflammation_target, sheet=sheet, startColumn=1, row.names=FALSE)
addDataFrame(category2_target, sheet=sheet, startColumn=3, row.names=FALSE)

sheet = createSheet(wb, "target_table_download")
addDataFrame(inflammation_target_table, sheet=sheet, startColumn=1, row.names=FALSE)
addDataFrame(category2_target_table, sheet=sheet, startColumn=3, row.names=FALSE)

saveWorkbook(wb, "GenesFoundInPathways.xlsx")
