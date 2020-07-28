#this is the keystone project
#credits: 

#Remove all objects
rm(list = ls())

#Set working directory to desktop
#setwd("C:/Users/choke/OneDrive/Desktop")

#Call in all relevant packages:
library(graph)
library(tidyverse)
library(limma)
library(edgeR)
library(dtplyr)
library(data.table)

#Function to read rna-seq data file & clean up
#This function is dependent to tidyverse/data.table package
#Therefore, you need to install those, if it is not installed yet
#If there is changes in those packages make sure everything is ok

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


#make a formal DGEList (a different class required for voom) for cancer and normal:
dge_cancer <- DGEList(rnaseq_cancer, genes = rownames(rnaseq_cancer))
dge_normal <- DGEList(rnaseq_normal, genes = rownames(rnaseq_normal))

#calculate normalization factors for voom
dge_cancer <- calcNormFactors(dge_cancer)
dge_normal <- calcNormFactors(dge_normal)

#need group matrix for DEG analysis:
#need to combine rnaseq_cancer and rnaseq_normal to make dimensions equal!!
#cancer_list <- rep("Cancer",ncol(rnaseq_cancer))
#normal_list <- rep("Normal",ncol(rnaseq_normal))
list <- c(cancer_list,normal_list)
Group <- factor(list, levels = c("Cancer","Normal"))
design <- model.matrix(~0 + Group)

#perform TMM normalization & log2CPM, voom does this all by itself, store log2CPM
#normalized values
voom_cancer <- voom(dge_cancer,design, plot = FALSE)$E
voom_normal <- voom(dge_normal,design, plot = FALSE)$E

#now perform contrasts and perform DEG:





#from here we need to calculate log2FC for each gene
#log2CPM_normal_value <- rowSums(log2CPM_normal, na.rm = FALSE, dims = 1)/ncol(log2CPM_normal)
#log2CPM_cancer_value <- rowSums(log2CPM_cancer, na.rm = FALSE, dims = 1)/ncol(log2CPM_cancer)

#log2CPM_normal <- log2CPM_normal %>% cbind(log2CPM_normal_value)
#log2CPM_cancer <- log2CPM_cancer %>% cbind(log2CPM_cancer_value)

#log2FC_value <- log2CPM_cancer_value / log2CPM_normal_value

#we got log2FC value for each gene, so we are done with preprocessing with edgeR/limma + voom.

#we need to map kegg pathway to entrez id, but tcga given data is in HGNC gene symbols
#first convert HGNC gene symbol to entrez id using gprofiler2 package
library(gprofiler2)
entrez_id <- as.vector(gconvert(query = (as.vector(rownames(rnaseq_cancer))), organism = "hsapiens", target = "ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)$target)

#can transform entrez to kegg genes by putting on a "hsa:"
KEGG_gene_id <- paste0("hsa:",entrez_id)

#then need to map entrez id to KEGG hsa (homo sapiens) pathways using KEGGREST package
library(KEGGREST)

#when supplying KEGG_gene_id at once, it cannot handle the amount of inputs. In fact,
#it handles around 400 inputs at once. KEGG_gene_id has 17260 outputs. We need to create a
#for loop...
i = 100
KEGG_pathway_id = c()
for (i in length(KEGG_gene_id)) {
    temp <- KEGG_gene_id[i - 99:i]
    temp_list <- keggLink("pathway", temp)
    KEGG_pathway_id <- c(KEGG_pathway_id,temp_list)
    #if (i > length(KEGG_gene_id)) {
       # break
    #}
    i = i + 100
}
#when you get KEGG_pathway_id, merge them all together. 

a <- keggLink("pathway",KEGG_gene_id)


#now we download and parse homo sapiens keggPathways using kpg:
library(ROntoTools)
kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)

#now we store pathwaynames in object kpn so we can select
kpn <- keggPathwayNames("hsa", updateCache = TRUE, verbose = TRUE)

#now finally perform pDis
pDisRes <- pDis(x = log2FC_value, graphs = kpg, nboot = 2000, verbose = FALSE )
