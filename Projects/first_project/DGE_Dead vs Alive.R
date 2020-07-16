#credits: https://github.com/jjnote/NeXT-VcoS/blob/master/SAHARA/Panc/SAHARA_Index_Calculator_Panc.R
#https://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
#https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf

rm(list = ls())
setwd("C:/Users/choke/OneDrive/Desktop/")

#call in all relevant library
library(survival)
library(limma) 
library(edgeR)
library(survminer)
library(tidyverse)
library(dtplyr)
library(data.table)
library(psych)
library(corrplot)

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

#cpm(y) calculates counts per million for each gene and patient, and adds all the patients' cpm value,
#and if that is greater than ncol(rnaseq_patients)/2 we include the genes in matrix y, and leave out the others.
#ncol(rnaseq_patients)/2 = top 50%
isexpr <- rowSums( cpm(y) > 1 ) >= ncol(rnaseq_patients)/2 

#apply logical isexpr matrix to y, only keep genes that we want
y <- y[isexpr, ][[1]]

#scale normalization
y <- apply(y,1:2,calcNormFactors)


#read in clinical factors and put in barcodes as column names of Clinical
Clinical <- read_tsv("ClinicalFactor_PAAD")
Barcodes <- intersect(colnames(rnaseq_patients), Clinical$Barcode)
Clinical <- Clinical[match(Barcodes, Clinical$Barcode),]

#create design matrix for a specific phenotype. Let's choose Dead vs Alive:
#create factors with two levels dead and alive, and then design a matrix with 1's and 
# 0's for dead/alive
Group <- factor(Clinical$Status, levels = c("Dead","Alive"))
design <- model.matrix(~0 + Group)
colnames(design) <- c("Dead","Alive")

#use voom to fit a linear model to log2 CPM and calculate residuals. $E gives 
#normalized log2 counts
v <- voom(y, design, plot = TRUE)$E
rnaseq_patients_v <- v

#use limma function (lmfit and eBayes) to perform differential expression testing, objects in fit are related to statistical testing

fit <- lmFit(v, design)
cont.matrix <- makeContrasts(DeadvsAlive = Dead - Alive, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

#print results
DEG_results <- topTable(fit2, adjust.method = "BH", sort.by = "p")
View(DEG_results)
print(rownames(DEG_results))

