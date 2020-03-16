#We can start make user-interaction
#Extract the IDs and locations of exons from any random reference gene that we select.

# setwd("~/Dropbox (ASU)/poneglyph/Script")

library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)
library(BiocGenerics) 
library(dplyr)
library(purrr) # accumulate

Find_id_loc = function(ETable,TRANS.ID,EXON.ID){
  Exon.loc = list()   #Relative location  
  Exon.id = list()   #Functional Id 
  server = "http://rest.ensembl.org"
  for(i in 1:length(TRANS.ID)){
    name = gsub('\"', "", TRANS.ID[i], fixed = TRUE)
    ext  = paste("/lookup/id/",name,"?expand=1",sep = "")
    r    = GET(paste(server, ext, sep = ""), content_type("application/json"))
    stop_for_status(r)
    id   = fromJSON(toJSON(content(r)))
    
    #locate translation range
    cds.box = c(id$Translation$start,id$Translation$end)
    
    #extract exon location
    id_exon = id$Exon
    exon.df = do.call(rbind,Map(data.frame,start=id_exon$start,end=id_exon$end))
    #5'3'UTR removal
    UTR.5 = which(exon.df$start <= cds.box[1])
    UTR.3 = which(exon.df$end >= cds.box[2])
    #locate the translation start/end
    if(exon.df$start[head(UTR.3,1)] <= cds.box[1]  && exon.df$end[tail(UTR.5,1)] >= cds.box[2]){#translation start inside the last exon
      cds.start = NULL
      cds.end   = cds.box[2] - cds.box[1] + 1 
    }else{#translation start before the last exon
      cds.start = exon.df$end[tail(UTR.5,1)] - cds.box[1] + 1 
      cds.end   = cds.box[2] - exon.df$start[head(UTR.3,1)] + 1 
    }
    
    #find exon location
    exon.fix.df  = exon.df[-c(UTR.5,UTR.3),]
    middle.len   = exon.fix.df$end - exon.fix.df$start+1
    exon.loc     = c(cds.start,c(middle.len,cds.end)) %>% accumulate(`+`)
    Exon.loc[[i]]= exon.loc  
    #filter the exon id
    exon.id.tmp    = EXON.ID[which(ETable$Transcript==TRANS.ID[i])]
    if(length(UTR.5)+length(UTR.3) <= 2){
      exonId.final = exon.id.tmp
    }else{
      exonId.final   = exon.id.tmp[-c(head(UTR.5,-1),tail(UTR.3,-1))]
    }
    Exon.id[[i]]   = exonId.final
  }
  Exon.ID = do.call(rbind,Map(data.frame,Transcript = TRANS.ID, Exon = Exon.id))
  res = list(Exon.ID,Exon.loc)
  return(res)
}


# inDir  = "../Data/Exon_table/"
# inFile = "../Data/geneId_cut.txt"
# ouDir  = "../Data/Results/From/"

main = function(inDir,inFile,ouDir){
  geneList   = read_delim(inFile,"\t", col_names = FALSE)
  name = c("human","chimp")
  #input: user 
  #gene.stem = "ENSG00000003402"
  #gene.stem = geneList[1971,]
  rand = sample(nrow(geneList),1,replace = FALSE)
  gene.stem  = geneList[rand,]
  
  ID.LIST = list()
  for (i in 1:length(gene.stem)) {
    exon.table = read_delim(paste0(inDir,gene.stem[[i]],".txt"),"\t",col_names = TRUE)
    trans.id = rle(exon.table[[1]])[[2]]
    exon.id  = rle(exon.table[[2]])[[2]]
    
    ID.list = Find_id_loc(exon.table,trans.id,exon.id)
    ID.LIST[[i]] = ID.list
  }
  
  
  for (j in 1:length(ID.LIST)) {
    write.table(ID.LIST[[j]][[1]],paste0(ouDir,name[[j]],"_exonId.txt"),quote = FALSE,
                row.names=FALSE,col.names = TRUE,append = FALSE,sep = "\t")
    lapply(ID.LIST[[j]][[2]], write, paste0(ouDir,name[[j]],"_exonLoc.txt"), append=TRUE, ncolumns=1000)
  }
  cat(unlist(gene.stem),file = paste0(ouDir,"geneId.txt"),sep = "\t")
 
}

args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2],args[3])