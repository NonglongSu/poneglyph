#Extract the exon Id from ensembl API.
#Never call it from the server. 

# setwd("~/Dropbox (ASU)/poneglyph/Script")

library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)
library(BiocGenerics) 

# Name1="ENSG00000003402"
# Name2="ENSPTRG00000012798"

main = function(Name1,Name2){
  
  nList=c(Name1,Name2)
  
  server = "http://rest.ensembl.org"
  for(i in 1:length(nList)){
    name = gsub('\"', "", nList[i], fixed = TRUE)
    ext  = paste("/lookup/id/",name,"?expand=1",sep = "")
    r    = GET(paste(server, ext, sep = ""), content_type("application/json"))
    stop_for_status(r)
    id   = fromJSON(toJSON(content(r)))
    id_trans  = id$Transcript[["id"]]
    id_trans  = id_trans[id$Transcript$biotype=="protein_coding"]
    id_exon   = id$Transcript$Exon[id$Transcript$biotype=="protein_coding"]
    #call the exon & cds
    exon.list=list()
    cds.list =list()
    for (j in 1:length(id_trans)) {
      ext_1 = paste("/sequence/id/",id_trans[[j]],"?type=cds",sep = "")
      r_1   = GET(paste(server, ext_1, sep = ""), content_type("text/x-fasta"))
      if(r_1$status_code==400){
        next
      }else{
        exon.list[[j]] = unlist(id_exon[[j]]$id)
        stop_for_status(r_1)
        cds_     = (content(r_1))
        cds.list = c(cds.list,cds_)
      }
    }
    exon.lt = Filter(Negate(is.null), exon.list)
    df = do.call(rbind,Map(data.frame,Transcript = id_trans, Exon = exon.lt))
   
    output1="../Data/cds/"
    output2="../Data/Exon_table/"
    
    write.table(cds.list,paste0(output1,name,".fa"),quote = FALSE,
                row.names=FALSE,col.names = FALSE,append = FALSE,sep = "\n")
    write.table(df,paste0(output2,name,".txt"),quote = FALSE,
                row.names=FALSE,col.names = TRUE,append = FALSE,sep = "\t")
  }
  
}

args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2])