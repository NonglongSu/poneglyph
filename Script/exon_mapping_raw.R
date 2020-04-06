library(tidyverse)
library(seqinr)
library(Biostrings)
library(dplyr)   #bind_rows
#setwd("~/Dropbox (ASU)/poneglyph/Script")

#Find the best alignment based on parsimony model. 
#   best_hit   = Map(`[`, trans.set.sorted.split, index.sc.min)
#   best_hits  = sapply(best_hit,function(x){names(readDNAStringSet(paste0(inDir,x),format = "fasta"))})
#   best_hits1 = lapply(best_hits,function(x){gsub("\\..*","",x)})
#   best_hits2 = matrix(unlist(best_hits1), nrow = 2, byrow = FALSE)
#   
#   res = list(trans.set.sorted.split,best_hit,best_hits2)
#   return(res)
# }

# #Find the hidden sequence structure
# state_map = function(seq1,seq2){
#   Hidden_state = c()
#   for (k in 1:length(seq1)) {
#     hit = paste0(c(seq1[k],seq2[k]),collapse = "")
#     if(grepl("-",hit)){#indel state
#       hidden_state = 0
#     }else{#substitution state
#       hidden_state = 1
#     }
#     Hidden_state = c(Hidden_state,hidden_state)
#   }
#   #rle(Hidden_state)
#   for (k  in 1:(length(Hidden_state)-1) ) {
#     if (Hidden_state[k] == Hidden_state[k+1]) {
#       Hidden_state[k] = 'M'
#     }
#   }
#   Hs = Hidden_state[Hidden_state == 0 | Hidden_state == 1]
#   return(Hs)
# }




#Find the exon-intron boundary in ref alignment.
extendLoc = function(dna,loc){
  geneA   = str_split(dna[[1]],"")[[1]]
  loc.2.0 = c()
  for (k in 1:length(loc)) {
    p = q = 1
    while (q < loc[k] ) {
      if(geneA[p] != '-'){
        q = q + 1
      }
      p = p + 1
    }
    loc.2.0 = c(loc.2.0,p)  
  }
  return(loc.2.0)
}
###
Uplocate = function(mRNA,loc){
  loc.set = list()
  for(j in 1:length(mRNA)){
    dna     = readDNAStringSet(paste0(Dir,mRNA[j]),format = "fasta")
    new.loc = extendLoc(dna,loc)
    loc.set[[j]] = new.loc
  }
  returnValue(loc.set)
}

#Find the corresponding chimp exon Id.
RecordTestId = function(p,q,cId,cLoc){
  left.edge  = nchar(str_remove_all(p,'-'))
  right.edge = nchar(str_remove_all(q,'-'))
  
  if(all(left.edge<=cLoc)){
    low = 1
  }else{
    low = tail(which(left.edge >= cLoc),1) + 1
  }
  high = which(right.edge <= cLoc)[1] 
  return(cId[low:high])
}
#check the 0/1 state.
s_map = function(gene){
  H_state = c()
  for (k in 1:length(gene)) {
    if(grepl("-",gene[k])){#indels
      h_state = 0
    }else{#subs
      h_state = 1
    }
    H_state = c(H_state,h_state)
  }
  val = length(rle(H_state)$values)
  return(val)
}
#Recover the exon-intron boundary from chimp alignment.
RecovermRNA = function(mRNA.Id,End,Start,cId,cLoc){
  dna  = readDNAStringSet(paste0(Dir,mRNA.Id),format = "fasta")
  sub1 = str_split(dna[[1]],"")[[1]][Start:End]
  sub2 = str_split(dna[[2]],"")[[1]][Start:End]
  #Tester
  State1 = s_map(sub1)
  if(State1 == 2 ){
    subseq1 = paste0(sub1[which(sub1!='-')],collapse = "")
    subseq2 = paste0(sub2[which(sub1!='-')],collapse = "")
    #Update the start/end
    Start = which(sub1!='-')[1]
    End   = tail(which(sub1!='-'),1)
  }else{
    subseq1 = toString(substr(dna[[1]],start=Start,stop=End))
    subseq2 = toString(substr(dna[[2]],start=Start,stop=End))
  }
  ###before/after in target seq.
  before.sub = toString(substr(dna[[2]],1,Start))
  after.sub  = toString(substr(dna[[2]],1,End))
  TestId = RecordTestId(before.sub,after.sub,cId,cLoc)
  
  res = list(subseq1,subseq2,TestId)
  return(res)
}


# a=c("A","A","A","-","-","-","A","A","A")
# b=c("A","A","A","A","A","A","-","-","-")


file   = "../Data/Results/From/geneId.txt"
f1     = "../Data/Results/From/human_exonId.txt"
f2     = "../Data/Results/From/human_exonLoc.txt"
f3     = "../Data/Results/From/chimp_exonId.txt"
f4     = "../Data/Results/From/chimp_exonLoc.txt"
inDir  = "../Data/Test/"

ouF1  = "../Data/Results/To/Parsi/exon_map_ref.fa"
ouF2  = "../Data/Results/To/Parsi/exon_map_tes.fa"
ouF3  = "../Data/Results/To/Parsi/exon_map.txt"
#i=1971


main = function(Rscr,file,f1,f2,f3,f4,inDir,ouF1,ouF2,ouF3){
  
  source(Rscr)
  #source("fhi_parsimony.R")
  
  Dir <<- inDir
  #read from the homo-geneId file
  geneId = read_delim(file,"\t", col_names = FALSE)
  gene.stem  = geneId[[1]]
  
  h.exoId    = read_delim(f1,"\t",col_names = TRUE)
  h.exoLoc   = readLines(f2)
  c.exoId    = read_delim(f3,"\t",col_names = TRUE)
  c.exoLoc   = readLines(f4)
  
  #Extract all transcripts from that specific gene
  mRNA.set       = list.files(path=Dir,pattern=paste0(gene.stem))
  mRNA.sorted    =  mRNA.set[order(regmatches(mRNA.set,regexpr("[0-9]_", mRNA.set)))]
  mRNA.sorted    =  mRNA.sorted[order(nchar(mRNA.sorted),mRNA.sorted)]
  dup            = rle(gsub("(_).*","", mRNA.sorted))[1][[1]][1] 
  mRNA.sorted.sp = split(mRNA.sorted,ceiling(seq_along(mRNA.sorted)/dup))       #18
  #
  h.Loc   = lapply(h.exoLoc, function(x){as.numeric(str_split(x," ")[[1]])})    #18
  c.Loc   = lapply(c.exoLoc, function(x){as.numeric(str_split(x," ")[[1]])})    #9
  #
  c.Id.sp = split(c.exoId$Exon,factor(c.exoId$Transcript,levels=unique(c.exoId$Transcript)) )   #9 
  
  #h.Loc.Bset = list() 
  Recom.mRNA1 = list()
  Recom.mRNA2 = list()
  Recom.Id    = list()
  for (i in 1:length(mRNA.sorted.sp)) {
    Loc.set = Uplocate(mRNA.sorted.sp[[i]],h.Loc[[i]])   #ref exon loc
    #parsimony, JC69, ...
    flag  = 1 
    recom.mRNA1 = c()
    recom.mRNA2 = c()
    recom.Id   = list()
    Index = model(mRNA.sorted.sp[[i]],Loc.set)                                  #9
    for (j in 1:length(Index)) {
      mRNA.Id = mRNA.sorted.sp[[i]][Index[j]] #This transcript owns the best exon
      Loc.Id  = Loc.set[[Index[j]]][j]        #The order of the exon
      mRNAplusId = RecovermRNA(mRNA.Id, Loc.Id, flag, c.Id.sp[[Index[j]]], c.Loc[[Index[j]]])
      #
      recom.mRNA1    = paste0(recom.mRNA1,mRNAplusId[[1]],sep="",collapse = NULL)
      recom.mRNA2    = paste0(recom.mRNA2,mRNAplusId[[2]],sep="",collapse = NULL)
      recom.Id[[j]]  = mRNAplusId[[3]]
      #
      flag = Loc.set[[Index[j+1]]][j] + 1     #The flag always has to be the recent one. 
      }
    Recom.mRNA1[[i]] = recom.mRNA1
    Recom.mRNA2[[i]] = recom.mRNA2
    Recom.Id[[i]]    = recom.Id
    }
  

    
  #write out fasta file
  order = 1:length(h.Loc)
  write.fasta(sequences = Recom.mRNA1,names = order, nbchar=80,
              open = "w",as.string = TRUE, file.out = ouF1)
  write.fasta(sequences = Recom.mRNA2,names = order, nbchar=80,
              open = "w",as.string = TRUE, file.out = ouF2)
  
  #Write out exon table
  Target.ExonId   = lapply(Recom.Id,function(x){as.data.frame(as.matrix(x))})
  Target.ExonId.1 = bind_rows(Target.ExonId, .id = "V1")                       #automatically generate index
  colnames(Target.ExonId.1) = c("Target.Index","Target.Exon")
  Exon.map = cbind(h.exoId,Target.ExonId.1)
  Exon.map$Target.Exon = vapply(Exon.map$Target.Exon, paste, collapse = ", ", character(1L))
  write.table(Exon.map,file = ouF3,sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE,append = FALSE)

}

args= commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9],args[10])