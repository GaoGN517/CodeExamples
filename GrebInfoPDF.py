## Load libraries
import PyPDF2 
import re
import os
import glob
import pandas as pd
import csv

## for one file extract all contents
pdfFileObj = open('Filename.pdf', 'rb') 
pdfReader = PyPDF2.PdfFileReader(pdfFileObj) 
print(pdfReader.numPages) 
# creating a page object 
pageObj = pdfReader.getPage(0) 
# extracting text from page 
Text = pageObj.extractText()
print(Text) 

## Function to extract based on a pattern
def findInfo(pattern, inputText):
    pattern_Group= re.search(pattern, Text)
    if (pattern_Group != None):
        Output = pattern_Group.group(1)
    else:
        Output = None
    return(Output)

## Here set the partten to be the string between two required strings. 
patterns = ["Beforestrin1(.*?)Afterstring1", "Beforestrin2(.*?)Afterstring2","Beforestrin3(.*?)Afterstring3",...]
len(patterns)

## For each extraction, we give it a name key
Keys = ["Key1","Key2","Key3",...]
len(Keys)
# Check whether length of patterns matches length of keys. 

## Function to greb information from file with patterns
def GetInfoFromFile(keys, patterns, FileName):
    Info = {}
    i = 0
    for pattern in patterns:
        Info[keys[i]] = findInfo(pattern, Text)
        i +=1
    return(Info)
   
## Test the function using Text we extracted above
GetInfoFromFile(Keys, patterns, Text)

## list all the files we want to extract
Filenames = glob.glob("*.pdf")
#print(Filenames)
len(Filenames)

## Extract information by patterns
FileInfo = {}
for filename in Filenames:
    pdfFileObj = open(filename, 'rb') 
    pdfReader = PyPDF2.PdfFileReader(pdfFileObj) 
    pageObj = pdfReader.getPage(0) 
    Text = pageObj.extractText()
    FileInfo[filename] = GetInfoFromFile(Keys, patterns, Text)
    # closing the pdf file object 
    pdfFileObj.close() 

## We want the output to be a table so we use pd.DataFrame function to store it as a dataframe. 
## And we want the rows to be each file, column to be the variables read. So we need to transpose the dataframe. 
df = pd.DataFrame(FileInfo).T

## Store the result to a csv file. 
df.to_csv("output.csv")





