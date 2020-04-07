<<<<<<< HEAD
# Use jukes cantor  + HMM
# ML method
=======
#Use jukes cantor + HMM
>>>>>>> f6574de58f14f532f4e68afbde9796b9f37a0e93

# library(Biostrings)
library(Matrix)
<<<<<<< HEAD
library(stringr)
=======
library(foreach)
library(doParallel)
>>>>>>> f6574de58f14f532f4e68afbde9796b9f37a0e93

obs_map = function(seq1, seq2, B) {
  Obs_prob = c()
  for (k in 1:length(seq1)) {
<<<<<<< HEAD
    hit = paste0(c(seq1[k],seq2[k]),collapse = "")
    if(grepl("-",hit)){# indel state
      obs_prob = B[5,1]
    }else{# substitution state
      if(grepl("AA|GG|CC|TT",hit)){
=======
    hit = paste0(c(seq1[k], seq2[k]), collapse = "")
    if (grepl("-", hit)) { 
      #indel state
      obs_prob = B[5,1]
    } else { 
      #substitution state
      if (grepl("AA|GG|CC|TT", hit)) {
>>>>>>> f6574de58f14f532f4e68afbde9796b9f37a0e93
        obs_prob = B[1,1]
      } else {
        obs_prob = B[1,2]
      }
    }
    Obs_prob = c(Obs_prob, obs_prob)
  }
  return(Obs_prob)
}


# mRNA--- 9 alignment files.
# Loc.set--- exon locus across 9 files.

<<<<<<< HEAD
model = function(mRNA, Loc.set){
  
  # Building the Q matrix.
  q_matrix = c(-1,1/3,1/3,1/3,
               1/3,-1,1/3,1/3,
               1/3,1/3,-1,1/3,
               1/3,1/3,1/3,-1)
  nuc_freq = 1/4
  Q = matrix (nuc_freq*q_matrix,4,4,byrow = T)
  
  # Cal. the score of msa and subset each transcript. 
  flag  = rep(1,length(mRNA)) # accumulation flag.
  Index = c()
  Score = c()
  for (i in 1:length(Loc.set[[1]])) {
    loc.i = unlist(lapply(Loc.set, `[[`, i), use.names = FALSE)
    score = c()
    for (j in 1:length(mRNA)){
    dna   = readDNAStringSet(paste0(Dir, mRNA[j]), format = "fasta")
    # Extract the exon.
    geneA = str_split(dna[[1]], "")[[1]][flag[j]:loc.i[j]]
    geneB = str_split(dna[[2]], "")[[1]][flag[j]:loc.i[j]]
    
    # number of mismatches
    subs     = length(which(geneA != '-')) - length(which(geneB == '-')) 
    mismatch = length(which(geneA != '-')) - length(which(geneB == '-')) - length(which(geneA == geneB))
    # number of indels
    gaps     = lapply(str_split(dna, ''), function(x){IRanges(x == '-')}) 
    indels   = length(IRangesList(gaps)[[1]])
    
    # number of sub->sub transitions
    non.gaps     = lapply(str_split(dna,''),function(x){IRanges(x!='-')})
=======
model = function(inDir, trans.set.sorted, dup) {
  
  #Building Q matrix
  q_matrix = c(-1, 1/3, 1/3, 1/3,
               1/3, -1, 1/3, 1/3,
               1/3, 1/3, -1, 1/3,
               1/3, 1/3, 1/3, -1)
  nuc_freq = 1/4
  Q = matrix(nuc_freq * q_matrix, 4, 4, byrow = T)
  
  numCores <- detectCores()
  registerDoParallel(numCores)

  #Calculate the score of msa and subset each transcript. 
  Score = c()
  foreach (j = 1:length(trans.set.sorted)) %dopar% {
    dna   = readDNAStringSet(paste0(inDir, trans.set.sorted[j]), format = "fasta")
    geneA = str_split(dna[[1]],"")[[1]]
    geneB = str_split(dna[[2]],"")[[1]]
    
    #number of mismatches
    subs = length(which(geneA != '-')) - length(which(geneB == '-')) 
    mismatch = length(which(geneA != '-')) - length(which(geneB == '-')) - length(which(geneA == geneB))
    #number of indels
    gaps     = lapply(str_split(dna, ''), function(x){IRanges(x == '-')}) 
    indels   = length(IRangesList(gaps)[[1]])
    
    #number of sub->sub transitions
    non.gaps     = lapply(str_split(dna, ''), function(x){IRanges(x!='-')})
>>>>>>> f6574de58f14f532f4e68afbde9796b9f37a0e93
    num.subs     = length(findOverlaps(non.gaps[[1]],non.gaps[[2]], maxgap=-1L))
    sub.sub      = subs - num.subs
    # number of sub->indel transitions
    sub.gap = sum(unlist(lapply(dna, function(x){length(which((gregexpr('[^-]-',x)[[1]] > -1)))})))
    # number of indel->indel transitions
    gap.gap = sum(unlist(lapply(gaps, function(x){width(x)})) - 1)
    # number of indel->sub transitions
    gap.sub = sum(unlist(lapply(dna, function(x){length(which((gregexpr('-[^-]',x)[[1]] > -1)))})))
   
<<<<<<< HEAD
    # cal. the initiation prob (pai).
    init_prob  = c(subs/(subs+indels),indels/(subs+indels)) 
    init_state = paste0(c(geneA[1],geneB[1]),collapse = "")
    if(grepl("-",init_state)){init_p = init_prob[2]}else{init_p = init_prob[1]}
    
    # cal. the transition matrix
    transition_matrix = matrix(c(sub.sub/(sub.sub+sub.gap),sub.gap/(sub.sub+sub.gap)
                                 ,gap.sub/(gap.sub+gap.gap),gap.gap/(gap.sub+gap.gap)),2,2,byrow=T)
    # cal. the  emission matrix
    Ds = mismatch / subs
    t  = -(3/4)*log((1-4/3*Ds),exp(1)) 
    P  = expm(Q*t)                                      # cal the branch length and the P matrix (substitution)
    emission_matrix = rbind(P,c(0.25,0.25,0.25,0.25))   # combine the indel weight with P matrix
    
    # Create a sequence of state and cal. the likelihood. 
    obs.state = obs_map(geneA,geneB,emission_matrix)
    hidden.state = transition_matrix^(c(sub.sub,gap.sub,sub.gap,gap.gap))
    # Cal. the likelihood
    likelihood = log(prod(c(init_p,obs.state,hidden.state)),exp(1))
    
    # hamming distance
    Score = c(Score,likelihood)
=======
    #calculate the initiation prob (pai).
    init_prob  = c(subs / (subs + indels), indels / (subs + indels)) 
    init_state = paste0(c(geneA[1], geneB[1]), collapse = "")
    if (grepl("-", init_state)) {
      init_p = init_prob[2]
    } else {
      init_p = init_prob[1]
    }
    
    #calculate the transition matrix
    transition_matrix = matrix(c(sub.sub / (sub.sub + sub.gap), sub.gap / (sub.sub + sub.gap),
                                 gap.sub / (gap.sub + gap.gap), gap.gap / (gap.sub + gap.gap)),
                                 2, 2, byrow = T)
    #calculate the  emission matrix
    Ds = mismatch / subs
    t  = -(3/4) * log(1 - ((4 / 3) * Ds), exp(1)) 
    P  = expm(Q * t)                                        #cal the branch length and the P matrix (substitution)
    emission_matrix = rbind(P, c(0.25, 0.25, 0.25, 0.25))   #combine the indel weight with P matrix
    
    #Create a sequence of state and cal. the likelihood. 
    obs.state = obs_map(geneA, geneB, emission_matrix)
    hidden.state = transition_matrix^(c(sub.sub, gap.sub, sub.gap, gap.gap))
    #Cal. the likelihood
    likelihood = log(prod(c(init_p, obs.state, hidden.state)), exp(1))
    
    #hamming distance
    Score = c(Score, likelihood)
>>>>>>> f6574de58f14f532f4e68afbde9796b9f37a0e93
  }

  trans.set.sorted.split = split(trans.set.sorted, ceiling(seq_along(trans.set.sorted) / dup))
  
  Score.split  = split(Score,ceiling(seq_along(Score) / dup))
  index.sc.max = lapply(Score.split, function(x){which(max(x) == x)})
  
  best_hit   = Map(`[`, trans.set.sorted.split, index.sc.max)
  best_hits  = sapply(best_hit, function(x){names(readDNAStringSet(paste0(inDir, x), format = "fasta"))})   
  best_hits1 = lapply(best_hits, function(x){gsub("\\..*", "", x)})
  best_hits2 = matrix(unlist(best_hits1), nrow = 2, byrow = FALSE)
}
  return(best_hits2)
  
}

