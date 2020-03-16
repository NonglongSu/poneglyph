library(tidyverse)
library(seqinr)
library(Biostrings)

#setwd("~/Dropbox (ASU)/poneglyph/Script")

#Find the best alignment based on parsimony model. 
# parsi_model = function(inDir,mRNA.sorted,dup){
#   Score = c()
#   #Cal. the score of msa and subset each transcript. 
#   for(j in 1:length(mRNA.sorted)){
#     dna   = readDNAStringSet(paste0(inDir,mRNA.sorted[j]),format = "fasta")
#     # geneA = str_split(dna[[1]],"")[[1]]
#     # geneB = str_split(dna[[2]],"")[[1]]
#     genes = extract_dna(inDir,mRNA.sorted[j])
#     #constant,linear,affine,convex
#     gaps     = lapply(str_split(dna,''),function(x){IRanges(x=='-')}) 
#     gap.cost = length(IRangesList(gaps)[[1]])
#     #hamming distance
#     mismatch = length(which(genes[[1]]!='-')) - length(which(genes[[2]]=='-')) - length(which(genes[[1]]==genes[[2]]))
#     Score = c(Score,sum(gap.cost,mismatch))
#   }
#   #trans.set.sorted.split = split(trans.set.sorted,ceiling(seq_along(trans.set.sorted)/dup))
#   
#   Score.split  = split(Score,ceiling(seq_along(Score)/dup))
#   index.sc.min = lapply(Score.split, function(x){which(min(x)==x)}[1])  #random pick one
#   
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

# #Extract the gene content
# extract_dna = function(inDir,stem){
#   dna   = readDNAStringSet(paste0(inDir,stem),format = "fasta")
#   geneA = str_split(dna[[1]],"")[[1]]
#   geneB = str_split(dna[[2]],"")[[1]]
#   
#   res = list(geneA,geneB)
#   return(res)
# }


# exon_match = function(inDir,Exon,Set,BestHit){
#   Gene.state = list()
#   for (p in 1:length(Set)) {
#     genes      = extract_dna(inDir,Set[p])
#     gene.state = state_map(genes[[1]],genes[[2]])
#     Gene.state[[p]] = gene.state
#   }
#   BestState = Gene.state[which(BestHit == Set)]
#   
# }


#Find the exon-intron boundary in human alignment.
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
#####################
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
  low  = tail(which(left.edge >= cLoc),1) + 1
  high = which(right.edge <= cLoc)[1]  

  return(cId[low:high])
}
#Recover the exon-intron boundary from chimp alignment.
RecovermRNA = function(mRNA.Id,End,Start,cId,cLoc){
  dna     = readDNAStringSet(paste0(Dir,mRNA.Id),format = "fasta")
  sub.dna = toString(substr(dna[[2]],start=Start,stop=End))
  if(Start>1){
    before.sub = toString(substr(dna[[2]],1,Start))
    after.sub  = toString(substr(dna[[2]],1,End))
    TestId = RecordTestId(before.sub,after.sub,cId,cLoc)
  }else{
    TestId = cId[1]
  }
  res = list(sub.dna,TestId)
  return(res)
}


a=c("A","A","A","-","-","-","A","A","A")
b=c("A","A","A","A","A","A","-","-","-")
q=p=1


file   = "../Data/Results/From/geneId.txt"
f1     = "../Data/Results/From/human_exonId.txt"
f2     = "../Data/Results/From/human_exonLoc.txt"
f3     = "../Data/Results/From/chimp_exonId.txt"
f4     = "../Data/Results/From/chimp_exonLoc.txt"
inDir  = "../Data/Test/"

ouFile  = "../Data/Results/To/exon_map.txt"
ouDir   = "../Data/Results/To/"
#i=1971


main = function(file,f1,f2,f3,f4,inDir,ouFile,ouDir){
  
  source("fhi_parsimony.R")
  
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
  #h.Id    = unique(h.exoId[1])
  h.Loc   = lapply(h.exoLoc, function(x){as.numeric(str_split(x," ")[[1]])})    #18
  c.Loc   = lapply(c.exoLoc, function(x){as.numeric(str_split(x," ")[[1]])})    #9
  #
  c.Id.sp = split(c.exoId$Exon,c.exoId$Transcript)                              #9 
  #
  #h.Loc.Bset = list() 
  Recom.mRNA = list()
  Recom.Id   = list()
  for (i in 1:length(mRNA.sorted.sp)) {
    Loc.set = Uplocate(mRNA.sorted.sp[[i]],h.Loc[[i]])   
    #parsimony, JC69, ...
    flag  = 1 
    recom.mRNA = c()
    recom.Id   = list()
    Index = model(mRNA.sorted.sp[[i]],Loc.set)                                  #9
    for (j in 1:length(Index)) {
      mRNA.Id = mRNA.sorted.sp[[i]][Index[j]] #The transcript owning the best exon
      Loc.Id  = Loc.set[[Index[j]]][j]        #The order of the exon
      mRNAplusId = RecovermRNA(mRNA.Id, Loc.Id, flag, c.Id.sp[[Index[j]]], c.Loc[[Index[j]]])
      #
      recom.mRNA    = paste0(recom.mRNA,mRNAplusId[[1]],sep="",collapse = NULL)
      recom.Id[[j]] = mRNAplusId[[2]]
      #
      flag = Loc.Id + 1
      }
    Recom.mRNA[[i]] = recom.mRNA
    Recom.Id[[i]]   = recom.Id
    }
  
  
    #Exon_matching
    
    
    
    #Final dataframe
    df = data.frame("geneId"=gene.stem,"human_trans"=Best_iso[1,],"chimp_trans"=Best_iso[2,])
    DF = rbind(DF,df)
  
  write.table(DF,file = ouFile,sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
  
}

args= commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])