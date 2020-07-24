#this is the keystone project
#credits: 

#Remove all objects
rm(list = ls())

#Set working directory to desktop
setwd("C:/Users/choke/OneDrive/Desktop")

#Call in all relevant packages:
library(graph)
library(tidyverse)
library(limma)
library(edgeR)
library(dtplyr)
library(data.table)
library(psych)

#Function to read rna-seq data file & clean up
#This function is dependent to tidyverse/data.table package
#Therefore, you need to install those, if it is not installed yet
#If there is changes in those packages make sure everything is ok
vi
read_rnaseq <- function(rnaseq_file){
    
    rnaseq <- fread(rnaseq_file)
   
    
    #change the column name "V1" in rna-seq data to "gene"
    rnaseq <- setnames(rnaseq, "V1", "gene")
    
    #remove all strings "\\?" and "\\|.*$" from all elements of rnaseq$gene
    
    rnaseq$gene <- str_remove_all(rnaseq$gene, "\\?")
    rnaseq$gene <- str_remove_all(rnaseq$gene, "\\|.*$")
    
    #remove all empty spaces for genes
    rnaseq[gene == ""] <- NA
    rnaseq <- drop_na(rnaseq, gene)
    
    #remove any duplicates
    rnaseq <- distinct(rnaseq, gene, .keep_all = TRUE)
    
    #from rnaseq$gene, take its column names and make them the rownames of rnaseq 
    rnaseq <- column_to_rownames(rnaseq, var = "gene")
    
    return(rnaseq)
}

#use read_rnaseq function to read & clean up the cancer and normal patient's data
rnaseq_cancer <- read_rnaseq("RNASeq_PAAD")
rnaseq_normal <- read_rnaseq("RNASeq_PAAD_Normal")

#rnaseq_patients <- cbind(rnaseq_cancer,rnaseq_normal) 

#make a formal DGEList (a different class required for voom) for cancer and normal:
dge_cancer <- DGEList(rnaseq_cancer, genes = rownames(rnaseq_cancer))
dge_normal <- DGEList(rnaseq_normal, genes = rownames(rnaseq_normal))

#calculate normalization factors for voom
dge_cancer <- calcNormFactors(dge_cancer)
dge_normal <- calcNormFactors(dge_normal)

#perform TMM normalization & log2CPM, voom does this all by itself, store log2CPM
#normalized values
log2CPM_cancer <- voom(dge_cancer, plot = FALSE)$E
log2CPM_normal <- voom(dge_normal, plot = FALSE)$E

#from here we need to calculate log2FC for each gene: rowSums...?
#after we get log2FC value for each gene, we are done with edgeR/limma + voom.



#now we download and parse homo sapiens keggPathways using kpg:
library(ROntoTools)
kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)

#now we store pathwaynames in object kpn so we can select
kpn <- keggPathwayNames("hsa", updateCache = TRUE, verbose = TRUE)

#now perform pDis
pDisRes <- pDis()
#x= fold change vector


#these are to map kegg pathway to entrez id
library(org.Hs.eg.db)
#somehow use org.Hs.egPATH
#look at: https://davetang.org/muse/2013/12/16/bioconductor-annotation-packages/
#maybe also: http://www.imsbio.co.jp/RGM/R_rdfile?f=org.Hs.eg.db/man/org.Hs.egBASE.Rd&d=R_BC

    #select vector of fc values?)



