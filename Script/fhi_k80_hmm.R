# Use Kimura_1980 + HMM
# ML method
# Indel emission via geometric method 

library(Matrix)
library(stringr)
library(doParallel)
library(stats4)

s_map = function(gene){
  H_state = c()
  for (k in 1:length(gene)) {
    if(grepl("-", gene[k])){#indels
      h_state = 0
    }else{#subs
      h_state = 1
    }
    H_state = c(H_state, h_state)
  }
  val = length(rle(H_state)$values)
  return(val)
}

obs_map = function(seq1, seq2, B){
  obs_prob = 1
  for (k in 1:length(seq1)) {
    hit = paste0(c(seq1[k], seq2[k]), collapse = "")
    if(grepl("-", hit)){# indel state
      obs_prob = obs_prob * B[5, 1]
    }else{# substitution state
      if(grepl("AA|GG|CC|TT", hit)){
        obs_prob = obs_prob * B[1, 1]              
      }else if(grepl("AG|GA|CT|TC", hit)){
        obs_prob = obs_prob * B[1, 2]              
      }else{
        obs_prob = obs_prob * B[1, 3]              
      }
    }
  }
  return(obs_prob)
}

#Define alpha and beta. 
para_def = function(seq1, seq2, mismatch){
    num_trans = 0
    for (k in 1:length(seq1)) {
      hit = paste0(c(seq1[k], seq2[k]), collapse = "")
      if(grepl("AG|GA|CT|TC", hit)){
         num_trans = num_trans + 1           
        }
      }
    num_transver = mismatch - num_trans
    res = (c(num_trans, num_transver))
    return(res)
  }

k80_hmm = function(geneA, geneB, I){
  # number of subs.
  subs   = length(which(geneA != '-')) - length(which(geneB == '-')) 
  # number of indels
  exon.dna = DNAStringSet(c(paste0(geneA, collapse = ""), paste0(geneB, collapse = "")))
  gaps     = lapply(str_split(exon.dna, ''), function(x){IRanges(x == '-')}) 
  indels   = sum(unlist(mclapply(IRangesList(gaps), function(x){length(x)})))
  
  # number of S->S transitions
  non.gaps  = lapply(str_split(exon.dna, ''), function(x){IRanges(x != '-')})
  num.olaps = length(findOverlaps(non.gaps[[1]], non.gaps[[2]], maxgap = -1L))
  sub.sub   = subs - num.olaps
  
  # number of I->D transitions(pretty rare).
  trans.pos = mclapply(exon.dna, function(x){str_locate_all(x, "[^-]-|-[^-]")[[1]][, 1]})
  gap.gap   =  length(which(trans.pos[[1]] %in% trans.pos[[2]]))
  gap.num   = sum(unlist(lapply(gaps, function(x){width(x)})))    # total length of gaps
  
  # number of S->I/D transitions.
  sub.gap  = sum(unlist(lapply(exon.dna, function(x){length(which((gregexpr('[^-]-', x)[[1]] > -1)))}))) - gap.gap
  
  # number of I/D->S transitions.
  gap.sub  = sum(unlist(lapply(exon.dna, function(x){length(which((gregexpr('-[^-]', x)[[1]] > -1)))}))) - gap.gap
  
  # cal. the initiation prob (pi).
  init_prob  = c(subs / (subs + indels), indels / (subs + indels)) 
  init_state = paste0(c(geneA[1], geneB[1]), collapse = "") 
  if (grepl("-", init_state)) {
    init_p = init_prob[2]
  }else{
    init_p = init_prob[1]
  }
  
  # cal. the transition matrix.
  transition_matrix = matrix(c(sub.sub / (sub.sub + sub.gap), sub.gap / (sub.sub + sub.gap)
                               ,gap.sub/ (gap.sub + gap.gap), gap.gap / (gap.sub + gap.gap)), 2, 2, byrow=T)
  
  # Build the Q matrix and find the distance K. 
  mismatch = subs - length(which(geneA == geneB))
  a = para_def(geneA, geneB, mismatch)[1]           # By using the apparent number of transitions/transversions.
  b = para_def(geneA, geneB, mismatch)[2] 
  p = a / subs
  q = b / subs
  
  q_matrix = c( -(a+2*b),       a,        b,         b,
                      a, -(a+2*b),        b,         b,
                      b,        b, -(a+2*b),         a,
                      b,        b,        a,   -(a+2*b))     # A,G,C,T
  nuc_freq = 1/4                                             # assume equal base freqs.                  
  Q        = matrix (nuc_freq * q_matrix, 4, 4, byrow = T)
  
  if(q >= 0.5 || (2*p+q) >= 1){# the limitation of the K2P model.
    return(Inf)
  }else{
    K  = -(0.5) * log((1 - 2 * p - q) * (sqrt(1 - 2 * q)), exp(1))    # K = -(1/2)*ln( (1-2p-q)*(sqrt(1-2q)) )
    P  = expm(Q * K)                                      # cal. the branch length and the P matrix (substitution)
  } 
  # Generate the emission matrix
  emission_matrix = rbind(P, I)                      
    
  # Cal. the product of the observed prob.
  obs.prob = obs_map(geneA, geneB, emission_matrix)
  # Cal. the product of the hidden prob. 
  trans_vec = as.vector(transition_matrix)
  trans_com = c(sub.sub, gap.sub, sub.gap, gap.gap)
  num = which(trans_com != 0)
  hidden.prob = prod(trans_vec[num] ^ (trans_com[num]))      # (ii)^n * (ij)^m * (ji)^p * (jj)^q
  # Cal. the (Pr(x,pi)).
  L = -log(prod(c(init_p, obs.prob, hidden.prob)), exp(1))
  return(L)
}

# mle_ab = function(dna){ R = (-2t / ln(1 - 2n2/n)) - 1}


# mRNA--- 9 alignment files.
# Loc.set--- exon locus across 9 files.

model = function(mRNA, Loc.set){
  
  # Assume the emission prob via the geometric method
  I = c(0.25, 0.25, 0.25, 0.25)                      
  
  # Cal. the L and pick the best exon candidate.
  flag  = rep(1, length(mRNA)) # accumulation flag.
  Index = c()
  for (i in 1:length(Loc.set[[1]])) {
    loc.i = unlist(lapply(Loc.set, `[[`, i), use.names = FALSE)
    Score = c()
    for (j in 1:length(mRNA)) {
      dna   = readDNAStringSet(paste0(Dir, mRNA[j]), format = "fasta")
      # Extract the exon.
      geneA = str_split(dna[[1]], "")[[1]][flag[j]:loc.i[j]]
      geneB = str_split(dna[[2]], "")[[1]][flag[j]:loc.i[j]]
      
      if(('-' %in% geneA) || ('-' %in% geneB) ){
        if(all(geneA == '-') || all(geneB == '-')){# alternative splicing NNNNNN/------ or ------/NNNNNN
          score = Inf
        }else{
          valA = s_map(geneA)
          valB = s_map(geneB)
          # Rule1--gaps exist before(after) the exon in ref.
          if(valA == 2){# only one consecutive gaps
            geneA.upd = geneA[which(geneA != '-')]
            geneB.upd = geneB[which(geneA != '-')]
            if(all(geneB.upd == '-')){# alternative splicing ---NNN/NNN---
               score = Inf
            }else{# partial overlap ---NNN/NNNNNN
               score = k80_hmm(geneA.upd, geneB.upd, I)
            }
          }else{
            score = k80_hmm(geneA, geneB, I)
          }
      }
    }else{ # Only substitutions occur.
      score = k80_hmm(geneA, geneB, I)
    }
      Score = c(Score, score)
   }
   #Find where the best exon from.
   index = which(min(Score) == Score)[1]   #Randomly pick one hit for now. (ignore the 'duplicated gene' assumption)
   Index = c(Index, index)
   #update the flag
   flag  = unlist(lapply(Loc.set, `[[`, i), use.names = FALSE) + 1 
 }
  
 return(Index)
}

