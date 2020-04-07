#extract the gene ID from ensembl database

# rm(list=ls())
# invisible(dev.off())

# source("http://bioconductor.org/biocLite.R")
# biocLite("BiocInstaller")
# biocLite("Biostrings")
# biocLite("biomaRt")
# biocLite("S4Vectors")
# biocLite("BiocGenerics")
# install.packages("readr")
# install.packages("dplyr")
# sessionInfo()                // show the version of packages & R


#setwd("~/Dropbox (ASU)/poneglyph/Script")

library(BiocInstaller)
library(BiocGenerics)
library(S4Vectors)
library(Biostrings)
library(biomaRt)
library(readr)
library(dplyr)

main=function(ouFile){
  human = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  # x = listDatasets()                  ## 190 datasets + version
  # y = listAttributes(human)$name
  # z = listFilters(human)$name
  
  #attribute name may change sometimes. 
  attributes = c("ensembl_gene_id",
                 "ensembl_gene_id_version",
                 "ptroglodytes_homolog_ensembl_gene",
                 "ptroglodytes_homolog_orthology_type"
  )
  
  
  Filt = "with_ccds"
  
  orth = getBM(attributes,filters=Filt, values=list(TRUE), mart = human, uniqueRows=TRUE)
  
  #Element-wise comparison single &. 
  orth_filtered = orth %>% filter(`ptroglodytes_homolog_orthology_type`=="ortholog_one2one")
  
  #Gene with mane transcript
  orth_gene = orth_filtered[,c(1,3)]
  
  # ouFile = "../Data/geneId.txt"
  
  write.table(orth_gene,ouFile,quote=FALSE,sep='\t',row.names=F,col.names = F)
}



args=commandArgs(trailingOnly=TRUE)
main(args[1])

