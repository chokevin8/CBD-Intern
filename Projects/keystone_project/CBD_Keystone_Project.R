#this is the keystone project
#credits: 

#Remove all objects
rm(list = ls())

#Set working directory to project folder
setwd("C:/Users/choke/OneDrive/Desktop/CBD_Intern/Projects/keystone_project")

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

#combine these two into one
rnaseq_all <- cbind(rnaseq_cancer,rnaseq_normal)

#make a formal DGEList (a different class required for voom) for all:
dge_all <- DGEList(rnaseq_all, genes = rownames(rnaseq_all))

#calculate normalization factors for voom
dge_all <- calcNormFactors(dge_all)

#need group matrix for DEG analysis:
cancer_all <- rep("Cancer",ncol(rnaseq_cancer))
normal_all <- rep("Normal",ncol(rnaseq_normal))
vec_all <- c(cancer_all,normal_all)
Group <- factor(vec_all, levels = c("Cancer","Normal"))
design <- model.matrix(~0 + Group)

#perform TMM normalization & log2CPM, voom does this all by itself
voom_all <- voom(dge_all,design, plot = FALSE)$E

#now perform contrasts and perform DEG:
fit <- lmFit(voom_all,design)
contr <- makeContrasts(GroupCancer - GroupNormal, levels = colnames(coef(fit)))
temp <- contrasts.fit(fit, contr)
temp <- eBayes(temp)
DEGTable <- topTable(temp, sort.by = "P", number = Inf)

#we got log2FC value for each gene, so we are done with preprocessing with edgeR/limma + voom.

#more preprocessing: we need to map kegg pathway to entrez id, but tcga given data is in HGNC gene symbols
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


#now we download and parse homo sapiens keggPathways using kpg:
library(ROntoTools)
kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)

#now we store pathwaynames in object kpn so we can select
kpn <- keggPathwayNames("hsa", updateCache = TRUE, verbose = TRUE)

#now finally perform pDis
pDisRes <- pDis(x = DEGTable$logFC, graphs = kpg, ref = rownames(DEGTable), nboot = 2000, verbose = FALSE )
