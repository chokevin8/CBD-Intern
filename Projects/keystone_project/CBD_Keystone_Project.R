#This is the keystone project.

#Remove all objects.
rm(list = ls())

#Call in all relevant packages.
library(graph)
library(ROntoTools)
library(GeneBook)
library(tidyverse)
library(limma)
library(edgeR)
library(dtplyr)
library(plyr)
library(data.table)

#Function to read rna-seq data file & clean up.
read_rnaseq <- function(rnaseq_file){
    
    rnaseq <- fread(rnaseq_file)
   
    #Change the column name "V1" in rna-seq data to "gene".
    rnaseq <- setnames(rnaseq, "V1", "gene")
    
    #Remove all strings "?" and zero or more instances of "|." from all elements of rnaseq$gene.
    rnaseq$gene <- str_remove_all(rnaseq$gene, "\\?")
    rnaseq$gene <- str_remove_all(rnaseq$gene, "\\|.*$")
    
    #Remove all empty spaces for genes.
    rnaseq[gene == ""] <- NA
    rnaseq <- drop_na(rnaseq, gene)
    
    #Remove any duplicates.
    rnaseq <- distinct(rnaseq, gene, .keep_all = TRUE)
    
    #From rnaseq$gene, take its column names and make them the rownames of rnaseq. 
    rnaseq <- column_to_rownames(rnaseq, var = "gene")
    
    return(rnaseq)
}

#Use read_rnaseq function to read & clean up the cancer and normal patient's data.
rnaseq_cancer <- read_rnaseq(str_c(as.character(getwd()),"/Data/RNASeq_LIHC"))
rnaseq_normal <- read_rnaseq(str_c(as.character(getwd()),"/Data/RNASeq_LIHC_Normal"))

#Combine cancer and normal data into one. 
rnaseq_all <- cbind(rnaseq_cancer,rnaseq_normal)

#Since duplicate names will happen for patients with normal data, add string "-NORMAL" to them. 
#Then assign those names back to column of rnaseq_all.
temp <- as.vector(colnames(rnaseq_all))
colnames(rnaseq_all) <- make.unique(temp, sep = "_Normal")
colnames(rnaseq_all) <- str_replace_all(colnames(rnaseq_all),"_Normal1","-Normal")
rm(temp)

#Make a formal DGEList (a different class required for voom) for all.
dge_all <- DGEList(rnaseq_all, genes = rownames(rnaseq_all))

#Calculate normalization factors for voom.
dge_all <- calcNormFactors(dge_all)

#Create a design matrix for DEG analysis.
cancer_all <- rep("Cancer",ncol(rnaseq_cancer))
normal_all <- rep("Normal",ncol(rnaseq_normal))
vec_all <- c(cancer_all,normal_all)
group <- factor(vec_all, levels = c("Cancer","Normal"))
design <- model.matrix(~0 + group)

#Perform TMM normalization & log2CPM, voom does this all by itself.
voom_all <- voom(dge_all,design, plot = FALSE)$E

#Now perform contrasts and finally perform DGE.
fit <- lmFit(voom_all,design)
contr <- makeContrasts(groupCancer - groupNormal, levels = colnames(coef(fit)))
temp <- contrasts.fit(fit, contr)
temp <- eBayes(temp)
degTable <- topTable(temp, sort.by = "P", number = Inf)

#Take out genes that contain "LOC" in degTable.
degTable <- dplyr::filter(degTable, !grepl("LOC", rownames(degTable)))
rm(temp)

#Change input gene ID's to official gene ID's that are recognized by the GeneCard package that will subsequently be used. 
base::rownames(degTable) <- make.names(alias2SymbolTable(rownames(degTable), species = "Hs"), unique = TRUE)

#We are done with edgeR/limma + voom since we have a DEGTable.
#Now, convert input TCGA gene symbol to Entrez id using GeneCard package.
#Clean up genecard_id data.
gc_list <- as.data.frame((GeneBook::genecard_id))
gc_list <- dplyr::filter(gc_list, grepl("Entrez", gc_list$subname))
gc_list$subname <- gsub('"','',gc_list$subname)
gc_list$subname <- str_replace_all(gc_list$subname, "Entrez Gene: ", "hsa:")
base::colnames(gc_list) <- c("gene_symbol", "entrez")

#Convert to Entrez id by using inner_join.
degTable <- degTable %>% select("logFC","P.Value","adj.P.Val") %>% rownames_to_column(var = "gene_symbol") %>%
    inner_join(gc_list, by = "gene_symbol") %>% column_to_rownames(var = "gene_symbol")

#Alternatively use gprofiler2 package to find entrez id's.
#library(gprofiler2)
#entrez <- as.vector(gconvert(query = temp, organism = "hsapiens", target = "ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)$target)
    
#Now we download and parse homo sapiens keggPathways using kpg.
kpg <- keggPathwayGraphs("hsa", updateCache = FALSE, verbose = FALSE)

#Alternatively:
#kpg <- keggPathwayGraphs("hsa", targRelTypes = c("GErel", "PCrel", "PPrel"), relPercThres = 0.9, nodeOnlyGraphs = FALSE, updateCache = FALSE, verbose = TRUE)
kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype", edgeWeightByType = list(activation = 1, inhibition = -1, expression = 1, repression = -1), defaultWeight = 0)
pv <- degTable$P.Value
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)

#Alternatively:
#kpg <- setNodeWeights(kpg, weights = alpha1MR(pv), defaultWeight = 1)

#Store pathway names in kpn.
kpn <- keggPathwayNames("hsa", updateCache = FALSE, verbose = TRUE)

#This is to perform cut-off-free analysis.
x <- degTable$logFC
names(x) <- degTable$entrez
ref <- degTable$entrez

#Finally perform primary disregulation (pDis).
pDisRes <- pDis(x = x, graphs = kpg, ref = ref, nboot = 1000000, verbose = FALSE )

#We have the result of most perturbed pathways ordered by pValue of pDis. 
result_ppDis <- Summary(pDisRes, pathNames = kpn, totalpDis = TRUE, pORA = FALSE, comb.pv = NULL, order.by = "ppDis")
View(result_ppDis)

#We have the result of most perturbed pathways ordered by total value of pDis.
result_totalpDis <- Summary(pDisRes, pathNames = kpn, totalpDis = TRUE, pORA = FALSE, comb.pv = NULL, order.by = "totalpDis")
View(result_totalpDis)

#Save the results in to a csv file. We are done now!
write_csv(result_ppDis, path = "/home/jjkim/Documents/r-project/CBD_Intern/Results/PAAD_ppDis_1000000.csv")
write_csv(result_totalpDis, path = "/home/jjkim/Documents/r-project/CBD_Intern/Results/PAAD_totalpDis_1000000.csv")

#Optional: 
#This is if you want to set TCGA gene name to Affy Probe ID's. 
#Note: didn't filter for sensitivity or specificity.
#gene_list <- rownames(degTable)
#gene_list <- gene_list %>% str_remove_all("[:blank:]") %>% str_remove_all("\\-AS[:digit:]") %>% str_remove_all("\\-AS") %>% str_remove_all("\\-") %>% str_to_upper()

#read probe data (annotation):
#mapped_probes <- read_tsv(str_c(as.character(getwd()),"/Data/GeneAnnot_191128_unmerge.tsv"))

#mapped_probes <- mapped_probes %>% 
#dplyr::select(c(AffyID, GeneSymbol, Sensitivity, Specificity)) %>%
#dplyr::rename(gene_name = GeneSymbol) %>% 
#drop_na(c(Sensitivity, Specificity, gene_name)) 
#mapped_probes$gene_name <- mapped_probes$gene_name %>% str_remove_all("[:blank:]") %>% str_remove_all("\\-AS[:digit:]") %>% str_remove_all("\\-AS") %>% str_remove_all("\\-") %>% str_to_upper()

#Now set rownames of degTable to AffyID.
#temp <- cbind(mapped_probes$AffyID,mapped_probes$gene_name) %>% as.data.frame()
#row_names <- temp[which(temp$V2 %in% gene_list),]
#rm(temp)
#row_names <- row_names[match(gene_list,row_names$V2),]$V1
#row_names <- make.names(row_names, unique = TRUE)
#row.names(degTable) <- row_names 

#Clean up row names of degTable.
#base::rownames(degTable) <- str_remove_all(base::rownames(degTable),"^X")
#degTable <- degTable[base::grep("NA.", base::rownames(degTable), invert = TRUE),] 
