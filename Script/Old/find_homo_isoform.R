library(tidyverse)
library(seqinr)
library(Biostrings)


# setwd("~/Dropbox (ASU)/poneglyph/Script")

# inFile = "../Data/geneId_cut.txt"
# inDir  = "../Data/Align_mafft/"
# ouFile = "../Data/Results/iso_align_tb.txt"

main = function(inFile,inDir,ouFile){
  #read from the homo-geneId file
  nList = read_delim(inFile,"\t", col_names = FALSE)
  #Create a final dataframe
  DF = c()
  for(i in 1:length(nList[[1]])){
    gene.stem  = nList[[1]][i]
    trans.set  = list.files(path=inDir,pattern=paste0(gene.stem))
    #Extract all transcripts from that specific gene
    trans.set.sorted = trans.set[order(regmatches(trans.set,regexpr("[0-9]_",trans.set)))]
    dup = rle(gsub("(_).*","",trans.set.sorted))[1][[1]][1] 
    #Cal. the score of msa and subset each transcript. 
    Score = c()
    Misma = c()
    for(j in 1:length(trans.set.sorted)){
      dna   = readDNAStringSet(paste0(inDir,trans.set.sorted[j]),format = "fasta")
      geneA = str_split(dna[[1]],"")[[1]]
      geneB = str_split(dna[[2]],"")[[1]]
      #constant,linear,affine,convex
      score = length(which(geneA=='-'))+length(which(geneB=='-'))
      Score = c(Score,score)
      #hamming distance
      misMatch = length(which(geneA[geneA != geneB] != '-'))
      Misma    = c(Misma,misMatch)
    }
    trans.set.sorted.split = split(trans.set.sorted,ceiling(seq_along(trans.set.sorted)/dup))
     
    Score.split  = split(Score,ceiling(seq_along(Score)/dup))
    Misma.split  = split(Misma,ceiling(seq_along(Misma)/dup))
    index.sc.min = lapply(Score.split, function(x){which(min(x)==x)})
    
    for(k in 1:length(index.sc.min)){
      if(length(index.sc.min[[k]])>1){#amount of gaps are the same!
        mis.score = Misma.split[[k]][index.sc.min[[k]]]
        index.sc.min[[k]] = index.sc.min[[k]][which(min(mis.score)==mis.score)]
      }
    }
    #Create a data.frame with Gene Id, homo_transcript. 
    best_hit   = Map(`[`, trans.set.sorted.split, index.sc.min)
    best_hits  = sapply(best_hit,function(x){names(readDNAStringSet(paste0(inDir,x),format = "fasta"))})   
    best_hits1 = lapply(best_hits,function(x){gsub("\\..*","",x)})
    best_hits2 = matrix(unlist(best_hits1), nrow = 2, byrow = FALSE)
    
    df = data.frame("geneId"=gene.stem,"human_trans"=best_hits2[1,],"chimp_trans"=best_hits2[2,])
    DF = rbind(DF,df)
  }
  write.table(DF,file = ouFile,sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
}

args= commandArgs(trailingOnly = TRUE)
main(args[1],args[2],arg[3])
main