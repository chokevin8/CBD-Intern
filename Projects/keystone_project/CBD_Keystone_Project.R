#this is the keystone project
#credits: 

#remove all objects
rm(list = ls())

#set working directory to desktop
setwd("C:/Users/choke/OneDrive/Desktop")

#call in all relevent packages:
library(ROntoTools)
library(graph)
library(tidyverse)
library(limma)
library(edgeR)
library(dtplyr)
library(data.table)
library(psych)

#function to read rna seq data file & cleaning up
read_rnaseq <- function(rnaseq_file){
    require(tidyverse)
    require(data.table)
    rnaseq <- fread(rnaseq_file)
    
    #change the column name "V1" in rnaseq data to "gene"
    rnaseq <- setnames(rnaseq, "V1", "gene")
    
    #remove all strings "\\?" and "\\|.*$" from all elements of rnaseq$gene
    
    rnaseq$gene <- str_remove_all(rnaseq$gene, "\\?")
    rnaseq$gene <- str_remove_all(rnaseq$gene, "\\|.*$")
    
    #remove all empty spaces
    rnaseq[gene == ""] <- NA
    rnaseq <- drop_na(rnaseq, gene)
    
    #remove any duplicates
    rnaseq <- distinct(rnaseq, gene, .keep_all = TRUE)
    
    #from rnaseq$gene, take its column names and make them the rownames of rnaseq 
    rnaseq <- column_to_rownames(rnaseq, var = "gene")
    
    #This function is dependent to tidyverse/data.table package.
    #Therefore, you need to install those, if it is not installed yet.
    #If there is changes in those packages make sure everything is OK.
    return(rnaseq)
}

#use read_rnaseq function to read & clean up the read data
rnaseq_patients <- read_rnaseq("RNASeq_PAAD")


#make a formal DGEList (a different class required for voom)
y <- DGEList(rnaseq_patients, genes = rownames(rnaseq_patients))
