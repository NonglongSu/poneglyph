#Use parsimony method to find the match. 
#Consecutive indels are treated as a single event. 

library(seqinr)
library(stringr)
library(Biostrings)

model = function(mRNA, Loc.set){
  
  #Cal. the score of multiple pairwise alignment and select the winner. 
  flag  =  rep(1,length(mRNA)) #accumulation flag
  Index =  c()
  for(i in 1:length(Loc.set[[1]])){
    loc.i = unlist(lapply(Loc.set,`[[`,i), use.names = FALSE)
    score = c()
    for (j in 1:length(mRNA)) {
      dna   = readDNAStringSet(paste0(Dir,mRNA[j]),format = "fasta")
      #Extract the exon part
      geneA = str_split(dna[[1]],"")[[1]][flag[j]:loc.i[j]]
      geneB = str_split(dna[[2]],"")[[1]][flag[j]:loc.i[j]]
      
      #constant,linear,affine,convex
      geneA.1 = DNAStringSet(paste0(geneA,collapse = ""),use.names = FALSE)
      geneB.1 = DNAStringSet(paste0(geneB,collapse = ""),use.names = FALSE)
      
      gaps    = lapply(str_split(c(geneA.1,geneB.1),''),function(x){IRanges(x=='-')}) 
      gapCost = sum(unlist(lapply(IRangesList(gaps),function(x){length(x)})))
      #Hamming distance
      misMatch = length(which(geneA!='-')) - length(which(geneB=='-')) - length(which(geneA==geneB))
      score    = c(score,sum(gapCost,misMatch))
    }
    #print(score[1])
    #print(score[4])
    #Find where the best exon from.
    index = which(min(Score)==Score)[1]   #Randomly pick one hit for now. (avoid the duplicated gene assumption)
    Index = c(Index,index)
    #update the flag
    flag  = unlist(lapply(Loc.set,`[[`,i), use.names = FALSE) + 1 
  }
  #trans.set.sorted.split = split(mRNA,ceiling(seq_along(mRNA)/dup))
  
  #Score.split  = split(Score,ceiling(seq_along(Score)/dup))
  # index.sc.min = lapply(Score.split, function(x){which(min(x)==x)})
  # 
  # best_hit   = Map(`[`, trans.set.sorted.split, index.sc.min)
  # best_hits  = sapply(best_hit,function(x){names(readDNAStringSet(paste0(inDir,x),format = "fasta"))})   
  # best_hits1 = lapply(best_hits,function(x){gsub("\\..*","",x)})
  # best_hits2 = matrix(unlist(best_hits1), nrow = 2, byrow = FALSE)
  
  return(Index)
  
}
