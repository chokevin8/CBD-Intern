#This is the keystone project.

#Remove all objects.
rm(list = ls())

#Call in all relevant packages.
library(graph)
library(ROntoTools)
library(KEGGREST)
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
rnaseq_cancer <- read_rnaseq(str_c(as.character(getwd()),"/Data/RNASeq_PAAD"))
rnaseq_normal <- read_rnaseq(str_c(as.character(getwd()),"/Data/RNASeq_PAAD_Normal"))

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

#Create a group matrix for DEG analysis:
cancer_all <- rep("Cancer",ncol(rnaseq_cancer))
normal_all <- rep("Normal",ncol(rnaseq_normal))
vec_all <- c(cancer_all,normal_all)
group <- factor(vec_all, levels = c("Cancer","Normal"))
design <- model.matrix(~0 + group)

#Perform TMM normalization & log2CPM, voom does this all by itself.
voom_all <- voom(dge_all,design, plot = FALSE)$E

#Now perform contrasts and finally perform DGE:
fit <- lmFit(voom_all,design)
contr <- makeContrasts(groupCancer - groupNormal, levels = colnames(coef(fit)))
temp <- contrasts.fit(fit, contr)
temp <- eBayes(temp)
degTable <- topTable(temp, sort.by = "P", number = Inf)

#Take out genes that contain "LOC" in degTable: 
degTable <- dplyr::filter(degTable, !grepl("LOC", rownames(degTable)))

#We are done with edgeR/limma + voom since we have a DEGTable.
#Now, convert input TCGA gene symbol to Entrez id using GeneCard package. 
    
gc_list <- as.data.frame((GeneBook::genecard_id))
gc_list <- dplyr::filter(gc_list, grepl("Entrez", gc_list$subname))
gc_list$subname <- gsub('"','',gc_list$subname)
gc_list$subname <- str_replace_all(gc_list$subname, "Entrez Gene: ", "hsa:")
base::colnames(gc_list) <- c("gene_symbol", "entrez")

degTable <- degTable %>% select("logFC","P.Value","adj.P.Val") %>% rownames_to_column(var = "gene_symbol") %>%
    inner_join(gc_list, by = "gene_symbol") %>% column_to_rownames(var = "gene_symbol")

#a <- alias2Symbol(rownames(degTable, species = "Hs"))


#.SD = subsets all of gc_list except gene (so subname) specified in the "by" argument, and we use lapply to paste the characters together for 
#subnames and group it by the gene names. 
#gc_list <- as.data.table(gc_list)
#collapsed_gc_list <- gc_list[, base::lapply(.SD, paste0, collapse=" "), by=gc_list$gene]
#collapsed_gc_list <- collapsed_gc_list %>% as.matrix() %>% str_remove_all("[:blank:]") %>% str_to_upper()

#Now match gc_list to rownames(degTable)
#temp <- rownames(degTable)
#temp <- str_to_upper(temp)
#temp <- as.data.frame(temp)
#colnames(temp) <- "V1"
#a <- temp[match(temp$V1,collapsed_gc_list),] 
#temp <- c(temp)

#entrez <- as.vector(gconvert(query = temp, organism = "hsapiens", target = "ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)$target)
#entrez_a <- as.vector(gconvert(query = a, organism = "hsapiens", target = "ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)$target)

#Now set TCGA gene name to Affy Probe ID's. 
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

#Then, map Entrez id to KEGG hsa (homo sapiens) pathways using KEGGREST package.
#When supplying KEGG_gene_id at once, it cannot handle the amount of inputs. In fact,
#it handles around 400 inputs at once. KEGG_gene_id has 17260 outputs. Therefore, we need to create a for loop. 

#i = 2
#KEGG_pathway_id = c()
#for (i in length(entrez)) {
#    temp <- entrez[i - 1:i]
 #   temp_list <- keggLink("pathway", temp)
#    KEGG_pathway_id <- c(entrez,temp_list)
 #   i = i + 1
#}
#temp_list <- keggLink("pathway", entrez)
    
#Now we download and parse homo sapiens keggPathways using kpg.
kpg <- keggPathwayGraphs("hsa", updateCache = FALSE, verbose = TRUE)
kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype", edgeWeightByType = list(activation = 1, inhibition = -1, expression = 1, repression = -1), defaultWeight = 0)
pv <- degTable$P.Value
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv))

#Store pathway names in kpn.
kpn <- keggPathwayNames("hsa", updateCache = FALSE, verbose = TRUE)

#Finally perform primary disregulation (pDis).
pDisRes <- pDis(x = degTable$logFC, graphs = kpg, ref = degTable$entrez, nboot = 2000, verbose = FALSE )

#We have the result of most perturbed pathways ordered by pValue. Stored in "Result". 
result <- Summary(pDisRes, pathNames = kpn, totalpDis = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "ppDis")


