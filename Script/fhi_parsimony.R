#Use parsimony method to find the match. 
#Consecutive indels are treated as a single event. 

library(seqinr)
library(Biostrings)

model = function(inDir, trans.set.sorted, dup){
  Score = c()
  #Cal. the score of msa and subset each transcript. 
  for(j in 1:length(trans.set.sorted)){
    dna   = readDNAStringSet(paste0(inDir,trans.set.sorted[j]),format = "fasta")
    geneA = str_split(dna[[1]],"")[[1]]
    geneB = str_split(dna[[2]],"")[[1]]
    #constant,linear,affine,convex
    gaps     = lapply(str_split(dna,''),function(x){IRanges(x=='-')}) 
    gap.cost = length(IRangesList(gaps)[[1]])
    #hamming distance
    misMatch = length(which(geneA[geneA != geneB] != '-'))
    Score = c(sum(gap.cost,misMatch),Score)
  }
  trans.set.sorted.split = split(trans.set.sorted,ceiling(seq_along(trans.set.sorted)/dup))
  
  Score.split  = split(Score,ceiling(seq_along(Score)/dup))
  index.sc.min = lapply(Score.split, function(x){which(min(x)==x)})
  
  best_hit   = Map(`[`, trans.set.sorted.split, index.sc.min)
  best_hits  = sapply(best_hit,function(x){names(readDNAStringSet(paste0(inDir,x),format = "fasta"))})   
  best_hits1 = lapply(best_hits,function(x){gsub("\\..*","",x)})
  best_hits2 = matrix(unlist(best_hits1), nrow = 2, byrow = FALSE)
  
  return(best_hits2)
  
}
