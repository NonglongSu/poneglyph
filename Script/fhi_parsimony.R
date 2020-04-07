# Consecutive indels are treated as a single event originally. 
# MP method
library(stringr)
<<<<<<< HEAD
# library(seqinr)
# library(Biostrings)
=======
library(Biostrings)
library(foreach)
library(doParallel)
>>>>>>> f6574de58f14f532f4e68afbde9796b9f37a0e93

s_map = function(gene) {
  H_state = c()
  for (k in 1:length(gene)) {
    if (grepl("-", gene[k])) { 
      #indels
      h_state = 0
    } else {
      #subs
      h_state = 1
    }
    H_state = c(H_state, h_state)
  }
  val = length(rle(H_state)$values)
  return(val)
}

p_score = function(geneA, geneB) {
  geneA.1 = DNAStringSet(paste0(geneA,collapse = ""), use.names = FALSE)
  geneB.1 = DNAStringSet(paste0(geneB,collapse = ""), use.names = FALSE)
  
  indel  = lapply(str_split(c(geneA.1, geneB.1), ''), function(x){IRanges(x == '-')}) 
  indels = sum(unlist(lapply(IRangesList(indel), function(x){length(x)})))
  # Hamming distance
  subs  = length(which(geneA != '-')) - length(which(geneB == '-')) - length(which(geneA == geneB))
  score = sum(indels, subs)
  return(score)
}

model = function(mRNA, Loc.set) {
  
  # Cal. the score of multiple pairwise alignment and select the winner. 
  flag  =  rep(1, length(mRNA)) #accumulation flag
  Index =  c()
  
  # Start Parallelization
  numCores <- detectCores()
  registerDoParallel(numCores)

  foreach (i = 1:length(Loc.set[[1]])) %dopar% {
    loc.i = unlist(lapply(Loc.set,`[[`,i), use.names = FALSE)
    score = c()
    for (j in 1:length(mRNA)) {
      dna   = readDNAStringSet(paste0(Dir, mRNA[j]), format = "fasta")
      #Extract the exon 
      geneA = str_split(dna[[1]],"")[[1]][flag[j]:loc.i[j]]
      geneB = str_split(dna[[2]],"")[[1]][flag[j]:loc.i[j]]
      
      if (('-' %in% geneA) || ('-' %in% geneB) ) {
        if (all(geneA == '-') || all(geneB == '-') ) {
          #alternative splicing NNNNNN/------ or ------/NNNNNN
          score = Inf
        } else {
          valA = s_map(geneA)
          valB = s_map(geneB)
          #Rule 1 -- gaps exist before(after) the exon in ref.
          if (valA == 2){
            #Only one consecutive gaps
            geneA.upd = geneA[which(geneA != '-')]
            geneB.upd = geneB[which(geneA != '-')]
            if (all(geneB.upd == '-')){#alternative splicing ---NNN/NNN---
              score = Inf
            } else {
              #partial overlap ---NNN/NNNNNN
              score = p_score(geneA.upd, geneB.upd)
            }
          } else {
              # Rule 2 -- gaps exist before (after) the exon in target (treat them as individual indel events)
              if (valA == 1 && valB == 2){
                #NNNNNN/-----N or NNNNNN/NNNNN-
                pScore = p_score(geneA, geneB)
                gaps   = length(which(geneB == '-'))
                score  = pScore + gaps - 1
              } else {
                score = p_score(geneA, geneB)
            }
          }
        }
      #Only substitutions occur.
      } else {
        score = p_score(geneA, geneB)
      }
      #Record each score.
      Score = c(Score, score)
    }
    #Find where the best exon from.
    index = which(min(Score) == Score)[1]   #Randomly pick one hit for now. (ignore the 'duplicated gene' assumption)
    Index = c(Index, index)
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
