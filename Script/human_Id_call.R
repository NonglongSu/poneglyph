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

main=function(output1,output2){
  ensembl=useMart("ensembl")
  human = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  
  # x = listDatasets()                  ## 190 datasets + version
  # y = listAttributes(human)$name
  # z = listFilters(human)$name
  
  #attribute name may change sometimes. 
  attributes = c("ensembl_gene_id",
                 "ensembl_gene_id_version",
                 "ensembl_transcript_id",
                 "ensembl_transcript_id_version"
  )
  
  
  Filt      = c("with_ccds","mane_select","transcript_appris")
  Filt.less = c("with_ccds","transcript_appris")
  
  Id.MANE      = getBM(attributes,filters=Filt, values=list(TRUE,TRUE,TRUE), mart = human, uniqueRows=TRUE)
  Id.Appr      = getBM(attributes,filters=Filt.less, values=list(TRUE,TRUE), mart = human, uniqueRows=TRUE)
  

  #Genes without MANE transcript
  Id.no.MANE   = Id.Appr[!(Id.Appr[,1] %in% Id.MANE[,1]),]
  
  # output1 = "../Raw_data/humanId_with_MANE.txt"
  # output2 = "../Raw_data/humanId_no_MANE.txt"

  write.table(Id.MANE,output1,quote=FALSE,sep='\t',row.names=F)
  write.table(Id.no.MANE,output2,quote=FALSE,sep='\t',row.names=F)
  
}





args=commandArgs(trailingOnly=TRUE)
main(args[1],args[2],args[3])

