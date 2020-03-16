library(tidyverse)
library(seqinr)
library(Biostrings)

# setwd("~/Dropbox (ASU)/poneglyph/Script")

# inFile = "../Data/geneId_cut.txt"
# inDir  = "../Data/Test/"
# ouFile = "../Data/Results/iso_parsi.txt"
# i=1971

main = function(rscript,inFile,inDir,ouFile){
  
  #introduced with script
  #source("fhi_gap_sub.R")
  source(rscript)
  
  #read from the homo-geneId file
  nList = read_delim(inFile,"\t", col_names = FALSE)
  
  DF = c()
  for(i in 1:length(nList[[1]])){
    gene.stem  = nList[[1]][i]
    trans.set  = list.files(path=inDir,pattern=paste0(gene.stem))
    #Extract all transcripts from that specific gene
    trans.set.sorted = trans.set[order(regmatches(trans.set,regexpr("[0-9]_",trans.set)))]
    dup = rle(gsub("(_).*","",trans.set.sorted))[1][[1]][1] 
    
    #Model
    Best_iso = model(inDir, trans.set.sorted, dup)
    #Final dataframe
    df = data.frame("geneId"=gene.stem,"human_trans"=Best_iso[1,],"chimp_trans"=Best_iso[2,])
    DF = rbind(DF,df)
  }
  write.table(DF,file = ouFile,sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
}

args= commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3],args[4])
