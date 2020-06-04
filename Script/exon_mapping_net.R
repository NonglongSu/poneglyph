library(seqinr)
library(Biostrings)
library(tidyverse)
library(dplyr)        # bind_rows

# Find the best alignment based on parsimony model.
#   best_hit   = Map(`[`, trans.set.sorted.split, index.sc.min)
#   best_hits  = sapply(best_hit,function(x){names(readDNAStringSet(paste0(inDir,x),format = "fasta"))})
#   best_hits1 = lapply(best_hits,function(x){gsub("\\..*","",x)})
#   best_hits2 = matrix(unlist(best_hits1), nrow = 2, byrow = FALSE)
#   res = list(trans.set.sorted.split,best_hit,best_hits2)
#   return(res)
# }

#Find the exon-intron boundary in ref seq.
extendLoc = function(dna, loc){
  geneA   = str_split(dna[[1]], "")[[1]]
  loc.2.0 = c()
  p = 0
  q = 0
  for (k in 1:length(loc)) {
    while (q < loc[k] ) {
      p = p + 1
      if(geneA[p] != '-'){
        q = q + 1
      }
    }
    loc.2.0 = c(loc.2.0, p)  
  }
  return(loc.2.0)
}
##
Uplocate = function(mRNA, loc){
  loc.set = list() 
  for(j in 1:length(mRNA)){
    dna     = readDNAStringSet(paste0(Dir, mRNA[j]), format = "fasta")
    new.loc = extendLoc(dna, loc)
    loc.set[[j]] = new.loc
  }
  returnValue(loc.set)
}

#Recheck the 0/1 state.
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
#Find the chimp exon Id.
RecordTestId = function(p, q, cId, cLoc){
  left.edge  = nchar(str_remove_all(p, '-'))
  right.edge = nchar(str_remove_all(q, '-'))
  
  if(all(left.edge<=cLoc)){
    low = 1
  }else{
    low = tail(which(left.edge >= cLoc), 1) + 1
  }
  high = which(right.edge <= cLoc)[1] 
  return(cId[low:high])
}
#Return the best exon seq/Id from chimp.
RecovermRNA = function(mRNA.Id, End, Start, cId, cLoc){
  dna  = readDNAStringSet(paste0(Dir,mRNA.Id),format = "fasta")
  sub1 = str_split(dna[[1]], "")[[1]][Start:End]
  sub2 = str_split(dna[[2]], "")[[1]][Start:End]
  #Test_state
  State1 = s_map(sub1)
  if(State1 == 2 ){
    subseq1 = paste0(sub1[which(sub1 != '-')], collapse = "")
    subseq2 = paste0(sub2[which(sub1 != '-')], collapse = "")
    #Update the start/end
    Start = which(sub1 != '-')[1]
    End   = tail(which(sub1 != '-'), 1)
  }else{
    subseq1 = toString(substr(dna[[1]], start=Start, stop=End))
    subseq2 = toString(substr(dna[[2]], start=Start, stop=End))
  }
  ###before/after in target seq.
  before.sub = toString(substr(dna[[2]], 1, Start))
  after.sub  = toString(substr(dna[[2]], 1, End))
  TestId     = RecordTestId(before.sub, after.sub, cId, cLoc)
  
  res = list(subseq1, subseq2, TestId)
  return(res)
}


file   = "../Data/geneId_update.txt"
d1     = "../Data/Mega_test/Results/From/Exon_id/"
d2     = "../Data/Mega_test/Results/From/Exon_loc/"
inDir  = "../Data/Align_mafft/"
ouD1   = "../Data/Mega_test/Results/To/Iso_id/"
ouD2   = "../Data/Mega_test/Results/To/Iso_seq/"
seed   = "../Data/Mega_test/Results/From/Exon_id/ENSG00000196547"

# setwd("~/Dropbox (ASU)/poneglyph/Script/Mega_test")
main = function(Rscr, seed, file, inDir, d1, d2, ouD1, ouD2){
  
  # Source script must be in the same folder as main script.
  source(Rscr) 
  #source("fhi_parsimony.R")
  #source("fhi_jc69_hmm.R")
  #source("fhi_k80_hmm.R")
  
  #read from the homo-geneId file
  Dir <<- inDir
  geneId = read_delim(file, "\t", col_names = FALSE)
  
  sd = str_extract(basename(seed), "[^.]+")
  gene.stem  = geneId[which(geneId[[1]] == sd), ]
  
  file = paste0(ouD1,gene.stem[1],".txt")
  if(file.exists(file)){
    quit()
  }
  #Match obj in database
  h.exoId    = read_delim(paste0(d1, gene.stem[1]), "\t", col_names = TRUE)
  h.exoLoc   = readLines(paste0(d2, gene.stem[1]))
  c.exoId    = read_delim(paste0(d1, gene.stem[2]), "\t", col_names = TRUE)
  c.exoLoc   = readLines(paste0(d2, gene.stem[2]))
  
  
  #Extract all transcripts from that specific gene
  mRNA.set       = list.files(path=inDir, pattern=paste0(gene.stem[1]))
  pre.order      = as.numeric(unlist(lapply((regmatches(mRNA.set, regexpr("[0-9]+_", mRNA.set))), 
                                  function(x){str_extract(x, "[^_]+")})))
  mRNA.sorted    = mRNA.set[order(pre.order)]
  dup            = rle(gsub("(_).*","", mRNA.sorted))[1][[1]][1] 
  
  mRNA.sorted.sp = split(mRNA.sorted, ceiling(seq_along(mRNA.sorted)/dup))       #18
  mRNA.sorted.sp = lapply(mRNA.sorted.sp, 
                          function(x){x[order(as.numeric(gsub("\\..*", "", gsub("^.*\\_", "", x))),x)]})
  
  h.Loc   = lapply(h.exoLoc, function(x){as.numeric(str_split(x, " ")[[1]])})    #18
  c.Loc   = lapply(c.exoLoc, function(x){as.numeric(str_split(x, " ")[[1]])})    #9
  #chimp Id
  c.Id.sp = split(c.exoId$Exon, factor(c.exoId$Transcript, levels = unique(c.exoId$Transcript)) )   #9 
  
  Recom.mRNA1 = list()
  Recom.mRNA2 = list()
  Recom.Id    = list()
  for (i in 1:length(mRNA.sorted.sp)) {
    Loc.set = Uplocate(mRNA.sorted.sp[[i]], h.Loc[[i]])   #ref exon loc -> 9
    #parsimony, JC69, ...
    Flag  = 1 
    recom.mRNA1 = c()
    recom.mRNA2 = c()
    recom.Id    = list()
    Index = model(mRNA.sorted.sp[[i]], Loc.set)                                  #9
    for (j in 1:length(Index)) {
      mRNA.Id = mRNA.sorted.sp[[i]][Index[j]] #This transcript owns the best exon
      Loc.Id  = Loc.set[[Index[j]]][j]        #The order of the exon
      mRNAplusId = RecovermRNA(mRNA.Id, Loc.Id, Flag, c.Id.sp[[Index[j]]], c.Loc[[Index[j]]])
      #
      recom.mRNA1    = paste0(recom.mRNA1,mRNAplusId[[1]], sep="", collapse = NULL)
      recom.mRNA2    = paste0(recom.mRNA2,mRNAplusId[[2]], sep="", collapse = NULL)
      recom.Id[[j]]  = mRNAplusId[[3]]
      #
      Flag = Loc.set[[Index[j+1]]][j] + 1     #The flag always has to be the recent one. 
    }
    Recom.mRNA1[[i]] = recom.mRNA1
    Recom.mRNA2[[i]] = recom.mRNA2
    Recom.Id[[i]]    = recom.Id
  }
  
  
  #write out fasta file
  order = 1:length(h.Loc)
  write.fasta(sequences = Recom.mRNA1, names = order, nbchar=80,
              open = "w", as.string = TRUE, file.out = paste0(ouD2,gene.stem[1],".fa"))
  write.fasta(sequences = Recom.mRNA2, names = order, nbchar=80,
              open = "w", as.string = TRUE, file.out = paste0(ouD2,gene.stem[2],".fa"))
  
  #Write out exon table
  Target.ExonId   = lapply(Recom.Id, function(x){as.data.frame(as.matrix(x))})
  Target.ExonId.1 = bind_rows(Target.ExonId, .id = "V1")                       #automatically generate index
  colnames(Target.ExonId.1) = c("Target.Index", "Target.Exon")
  Exon.map = cbind(h.exoId[, 1:2], Target.ExonId.1)
  Exon.map$Target.Exon = vapply(Exon.map$Target.Exon, paste, collapse = ", ", character(1L))
  write.table(Exon.map, file = paste0(ouD1,gene.stem[1],".txt"), 
              sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE, append = FALSE)
  
}

args= commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])