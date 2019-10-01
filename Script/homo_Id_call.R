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

#setwd("~/Dropbox (ASU)/poneglyph/Script")

library(BiocInstaller)
library(BiocGenerics)
library(S4Vectors)
library(Biostrings)
library(biomaRt)
library(readr)
library(dplyr)

main=function(output1,output2,output3){
  ensembl=useMart("ensembl")
  human = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  
  # x = listDatasets()                  ## 190 datasets + version
  # y = listAttributes(human)$name
  # z = listFilters(human)$name
  
  #attribute name may change sometimes. 
  attributes = c("ensembl_gene_id",
                 "ensembl_gene_id_version",
                 "ensembl_transcript_id",
                 "ensembl_transcript_id_version",
                 "cgchok1gshd_homolog_ensembl_gene",
                 "cgchok1gshd_homolog_orthology_type",
                 "mmusculus_homolog_ensembl_gene",
                 "mmusculus_homolog_orthology_type" ,
                 "rnorvegicus_homolog_ensembl_gene",
                 "rnorvegicus_homolog_orthology_type"
  )
  
  
  Filt      = c("with_ccds","mane_select","transcript_appris")
  Filt.less = c("with_ccds","transcript_appris")
  
  orth      = getBM(attributes,filters=Filt, values=list(TRUE,TRUE,TRUE), mart = human, uniqueRows=TRUE)
  orth.more = getBM(attributes,filters=Filt.less, values=list(TRUE,TRUE), mart = human, uniqueRows=TRUE)
  
  #Element-wise comparison single &. 
  orth_filtered = orth %>% filter(`mmusculus_homolog_orthology_type`=="ortholog_one2one" & 
                                    `rnorvegicus_homolog_orthology_type`=="ortholog_one2one" &
                                    `cgchok1gshd_homolog_orthology_type`=="ortholog_one2one")
  
  orth_filtered.more = orth.more %>% filter(`mmusculus_homolog_orthology_type`=="ortholog_one2one" & 
                                              `rnorvegicus_homolog_orthology_type`=="ortholog_one2one" &
                                              `cgchok1gshd_homolog_orthology_type`=="ortholog_one2one")
  
  
  #Gene with mane transcript
  orth_gene = orth_filtered[,c(1,5,7,9)]
  ref_transcript = orth_filtered[,c(1,2,3,4)]
  
  
  #Gene without mane transcript
  orth_gene_more = orth_filtered.more[!duplicated(orth_filtered.more[,1]),]
  orth_gene_left = orth_gene_more[!(orth_gene_more[,1] %in% orth_filtered[,1]),]
  
  orth_gene_left = orth_gene_left[,c(1,5,7,9)]
  
  
  # output1 = "../Raw_data/geneId_with_MANE.txt"
  # output2 = "../Raw_data/transId_with_MANE.txt"
  # output3 = "../Raw_data/geneId_no_MANE.txt"
  
  write.table(orth_gene,output1,quote=FALSE,sep='\t',row.names=F)
  write.table(ref_transcript,output2,quote=FALSE,sep='\t',row.names=F)
  write.table(orth_gene_left,output3,quote=FALSE,sep='\t',row.names=F)
  
}





args=commandArgs(trailingOnly=TRUE)
main(args[1],args[2],args[3])

