#This is the keystone project
#credits: 

#Remove all objects.
rm(list = ls())

#Call in all relevant packages.
library(graph)
library(ROntoTools)
library(KEGGREST)
library(gprofiler2)
library(annotate)
library(hgu133a2.db)
library(tidyverse)
library(limma)
library(edgeR)
library(dtplyr)
library(plyr)
library(data.table)

#Function to read rna-seq data file & clean up.
read_rnaseq <- function(rnaseq_file){
    
    rnaseq <- fread(rnaseq_file)
   
    #change the column name "V1" in rna-seq data to "gene"
    rnaseq <- setnames(rnaseq, "V1", "gene")
    
    #remove all strings "?" and zero or more instances of "|." from all elements of rnaseq$gene
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
rnaseq_cancer <- read_rnaseq(str_c(as.character(getwd()),"/Data/RNASeq_PAAD"))
rnaseq_normal <- read_rnaseq(str_c(as.character(getwd()),"/Data/RNASeq_PAAD_Normal"))

#Combine cancer and normal data into one. 
rnaseq_all <- cbind(rnaseq_cancer,rnaseq_normal)

#since duplicate names will happen for patients with normal data, add string "-NORMAL" to them. 
#Then assign those names back to column of rnaseq_all
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

#We are done with preprocessing with edgeR/limma + voom since we have a DEGTable. 


#More preprocessing: We need to map KEGG pathway to Entrez id, but input TCGA data is in HGNC gene symbols.
#First, convert HGNC gene symbol to Entrez id using gprofiler2 package.
entrez_id <- as.vector(gconvert(query = (as.vector(rownames(rnaseq_all))), organism = "hsapiens", target = "ENTREZGENE_ACC", mthreshold = Inf, filter_na = TRUE)$target)

#Transform Entrez to KEGG genes by putting on a "hsa:".
kegg_gene_id <- paste0("hsa:",entrez_id)

#work in progress: bind KEGG_gene_id to DEGTable
#DEGTable <- rbind.fill(as.data.frame(t(DEGTable)),as.data.frame(t(KEGG_gene_id)))
#set column name to "entrez"

#Now we make the gene names of DEGTable in to probe names of affymetrix. 
#gene_list <- rownames(degTable)
#mapped_probes <- mappedkeys(hgu133a2SYMBOL)
#probe_n_sym <- as.data.frame(hgu133a2SYMBOL[mapped_probes])
#result <- probe_n_sym[which(probe_n_sym$symbol %in% gene_list),]
#result_1 <- result[match(gene_list,result$symbol),]

gene_list <- rownames(degTable)
gene_list <- gene_list %>% str_remove_all("[:blank:]") %>% str_remove_all("\\-AS[:digit:]") %>% str_remove_all("\\-AS") %>% str_remove_all("\\-") %>% str_to_upper()
mapped_probes <- read_tsv(str_c(as.character(getwd()),"/Data/GeneAnnot_191128_unmerge.tsv"))
mapped_probes <- mapped_probes %>% 
dplyr::select(c(AffyID, GeneSymbol, Sensitivity, Specificity)) %>%
dplyr::rename(Gene_Name = GeneSymbol) %>% 
drop_na(c(Sensitivity, Specificity, Gene_Name)) 
#didn't filter for sensitivity or specificity

mapped_probes$Gene_Name <- mapped_probes$Gene_Name %>% str_remove_all("[:blank:]") %>% str_remove_all("\\-AS[:digit:]") %>% str_remove_all("\\-AS") %>% str_remove_all("\\-") %>% str_to_upper()
probe_n_sym <- cbind(mapped_probes$AffyID,mapped_probes$Gene_Name)
probe_n_sym <- as.data.frame(probe_n_sym)
result <- probe_n_sym[which(probe_n_sym$V2 %in% gene_list),]
result_1 <- result[match(gene_list,result$V2),]

#do separate_rows(불러온데이터, GeneSymbol, sep = "[:blank:]" or sep = "") in above?

#Then, map Entrez id to KEGG hsa (homo sapiens) pathways using KEGGREST package.
#When supplying KEGG_gene_id at once, it cannot handle the amount of inputs. In fact,
#it handles around 400 inputs at once. KEGG_gene_id has 17260 outputs. Therefore, we need to create a for loop. 

#i = 2
#KEGG_pathway_id = c()
#for (i in length(KEGG_gene_id)) {
#    temp <- KEGG_gene_id[i - 1:i]
 #   temp_list <- keggLink("pathway", temp)
#    KEGG_pathway_id <- c(KEGG_pathway_id,temp_list)
 #   i = i + 1
#}

#temp_list <- keggLink("pathway", kegg_gene_id)
    
# Both gives an error:"Error in curl::curl_fetch_memory(url, handle = handle) : 
# Send failure: Connection reset by peer"".
# But trying this will not give an error:

# a <- kegg_gene_id[1:100]
# b <- keggLink("pathway", a)
# head(b) matches gene to pathway.
# So i am wondering what the error is about. 


#Now we download and parse homo sapiens keggPathways using kpg.
kpg <- keggPathwayGraphs("hsa", updateCache = TRUE, verbose = TRUE)

#Store pathway names in kpn.
kpn <- keggPathwayNames("hsa", updateCache = TRUE, verbose = TRUE)

#Finally perform primary disregulation (pDis).
pDisRes <- pDis(x = degTable$logFC, graphs = kpg, ref = rownames(degTable), nboot = 2000, verbose = FALSE )

#We have the result of most perturbed pathways ordered by pValue. Stored in "Result". 
result <- Summary(pDisRes, pathNames = kpn, totalpDis = FALSE, pORA = FALSE, comb.pv = NULL, order.by = "ppDis")


