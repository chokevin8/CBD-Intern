#This is the keystone project
#credits: 

#Remove all objects.
rm(list = ls())

#Set working directory to project folder.
#setwd("C:/Users/choke/OneDrive/Desktop/CBD_Intern/Projects/keystone_project")

#Call in all relevant packages.
library(graph)
library(tidyverse)
library(limma)
library(edgeR)
library(dtplyr)
library(data.table)

#Function to read rna-seq data file & clean up.
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

#Use read_rnaseq function to read & clean up the cancer and normal patient's data
rnaseq_cancer <- read_rnaseq("/Data/RNASeq_PAAD")
rnaseq_normal <- read_rnaseq("/Data/RNASeq_PAAD_Normal")

#Combine cancer and normal data into one. 
rnaseq_all <- cbind(rnaseq_cancer,rnaseq_normal)

#since duplicate names will happen for patients with normal data, add string "-NORMAL" to them. 
#Then assign those names back to column of rnaseq_all

temp <- as.vector(colnames(rnaseq_all))
colnames(rnaseq_all) <- make.unique(temp, sep = "_Normal")
colnames(rnaseq_all) <- str_replace_all(colnames(rnaseq_all),"_Normal1","-Normal")

#Make a formal DGEList (a different class required for voom) for all.
dge_all <- DGEList(rnaseq_all, genes = rownames(rnaseq_all))

#Calculate normalization factors for voom.
dge_all <- calcNormFactors(dge_all)

#Create a group matrix for DEG analysis:
cancer_all <- rep("Cancer",ncol(rnaseq_cancer))
normal_all <- rep("Normal",ncol(rnaseq_normal))
vec_all <- c(cancer_all,normal_all)
Group <- factor(vec_all, levels = c("Cancer","Normal"))
design <- model.matrix(~0 + Group)

#Perform TMM normalization & log2CPM, voom does this all by itself.
voom_all <- voom(dge_all,design, plot = FALSE)$E

#Now perform contrasts and finally perform DGE:
fit <- lmFit(voom_all,design)
contr <- makeContrasts(GroupCancer - GroupNormal, levels = colnames(coef(fit)))
temp <- contrasts.fit(fit, contr)
temp <- eBayes(temp)
DEGTable <- topTable(temp, sort.by = "P", number = Inf)

#We are done with preprocessing with edgeR/limma + voom since we have a DEGTable. 

#More preprocessing: We need to map KEGG pathway to Entrez id, but input TCGA data is in HGNC gene symbols.
#First, convert HGNC gene symbol to Entrez id using gprofiler2 package.
library(gprofiler2)
Entrez_id <- as.vector(gconvert(query = (as.vector(rownames(rnaseq_cancer))), organism = "hsapiens", target = "ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)$target)

#Transform Entrez to KEGG genes by putting on a "hsa:".
KEGG_gene_id <- paste0("hsa:",Entrez_id)

#Then, map Entrez id to KEGG hsa (homo sapiens) pathways using KEGGREST package.
library(KEGGREST)

#When supplying KEGG_gene_id at once, it cannot handle the amount of inputs. In fact,
#it handles around 400 inputs at once. KEGG_gene_id has 17260 outputs. Therefore, we need to create a for loop. 

i = 100
KEGG_pathway_id = c()
for (i in length(KEGG_gene_id)) {
    temp <- KEGG_gene_id[i - 99:i]
    temp_list <- keggLink("pathway", temp)
    KEGG_pathway_id <- c(KEGG_pathway_id,temp_list)
    if (i > length(KEGG_gene_id)) {
       break
    }
    i = i + 100
}

#Now we download and parse homo sapiens keggPathways using kpg.

library(ROntoTools)
kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)

#Store pathway names in kpn.
kpn <- keggPathwayNames("hsa", updateCache = TRUE, verbose = TRUE)

#Finally perform primary disregulation (pDis).
pDisRes <- pDis(x = DEGTable$logFC, graphs = kpg, ref = rownames(DEGTable), nboot = 2000, verbose = FALSE )

#We have the result of most perturbed pathways ordered by pValue. Stored in "Result". 
Result <- Summary(pDisRes, pathNames = kpn, totalpDis = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "ppDis")

