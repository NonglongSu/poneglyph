library(seqinr)
library(stringr)
library(Biostrings)

model = function(inDir, trans.set.sorted, dup){
   Score = c()
   Misma = c()
   #Cal. the score of msa and subset each transcript. 
    for(j in 1:length(trans.set.sorted)){
      dna   = readDNAStringSet(paste0(inDir,trans.set.sorted[j]),format = "fasta")
      geneA = str_split(dna[[1]],"")[[1]]
      geneB = str_split(dna[[2]],"")[[1]]
      #constant,linear,affine,convex
      score = length(which(geneA=='-'))+length(which(geneB=='-'))
      Score = c(Score,score)
      #hamming distance
      mismatch = length(which(geneA!='-')) - length(which(geneB=='-')) - length(which(geneA==geneB))
      Misma    = c(Misma,mismatch)
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

    best_hit   = Map(`[`, trans.set.sorted.split, index.sc.min)
    best_hits  = sapply(best_hit,function(x){names(readDNAStringSet(paste0(inDir,x),format = "fasta"))})   
    best_hits1 = lapply(best_hits,function(x){gsub("\\..*","",x)})
    best_hits2 = matrix(unlist(best_hits1), nrow = 2, byrow = FALSE)
    
    return(best_hits2)
        
}
  
