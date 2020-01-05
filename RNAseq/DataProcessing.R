################
## Expression
################
## Read Raw Data.
filenames <- list.files("RNAseq/RawData/")
filenames_len <- length(filenames)
files <- list()
for(i in 1:filenames_len) {
  files[[i]] <- read.table(paste0("RNAseq/RawData/", filenames[i]))
}
head(files[[27]])
## Now we have a list with 26 elements, each element is the result for a sample. 
## The sample name is in the first row, second column. The rest rows are for each
## gene (mRNA).
#sample_names <- unlist(lapply(files, '[', 2, ))
sample_names <- as.character(do.call(rbind, lapply(files, head, 1))[, 2])

## Then merge the genes from each samples. 
#RNAs <- do.call(function(x, y) merge(x, y, by = V1, all = T), lapply(files, function(x) x[-1, ]))
RNAs <- lapply(files, function(x) x[-1, ])

## Merge each element in the list
Expression <- Reduce(function(x, y) merge(x, y, by = "V1", all = T), RNAs)

## Change the column names to sample names
colnames(Expression) <- c("Gene", sample_names)
rowname_Exp <- Expression[, 1]
Expression <- Expression[, -1]
head(Expression)

Expression <- apply(Expression, 2, as.numeric)
rownames(Expression) <- rowname_Exp
## Save data.
saveRDS(Expression, file = "Expression.rds")

################
## Sample Information
################
GRP <- read.csv("RNAseq/GroupInfor.csv")
RIN <- read.csv("miRNA_RIN.csv")
RIN$Sample <- gsub("/", "_", RIN$Sample)
RIN$Sample <- gsub("\\.", "_", RIN$Sample)
RIN$Sample
GRP$Sample <- as.character(GRP$Sample)
GRP$Sample
GRP_RIN <- merge(GRP, RIN, by = "Sample", all.x = T)
saveRDS(GRP_RIN, file = "GroupInfor.rds")
