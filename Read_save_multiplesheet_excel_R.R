###################################################################
## Save to excel file with multiple sheets.
library(xlsx)
wb = createWorkbook()

sheet = createSheet(wb, "sheetname1")
addDataFrame(dataframe1, sheet=sheet, startColumn=1, row.names=T)
sheet = createSheet(wb, "sheetname2")
addDataFrame(dataframe2, sheet=sheet, startColumn=1, row.names=T)

saveWorkbook(wb, "filename.xlsx")

###################################################################
## read excel file with multiple sheets.
library(tidyverse)
library(readxl)
## check which sheets are there in the file.
excel_sheets("filename.xlsx")
## read the sheets we need. 
dataframe1 <- read_excel("filename.xlsx", sheet = "sheetname1")
dataframe2 <- read_excel("filename.xlsx", sheet = "sheetname2")