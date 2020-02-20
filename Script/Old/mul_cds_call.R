#Extract the transcript Id from ensembl API.
#Never call it from the server. 

# setwd("~/Dropbox (ASU)/poneglyph/Script")

library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)


# Name1="ENSG00000114735"
# Name2="ENSMAMG00000016211"


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
    id_trans  = id_trans[id$Transcript$biotype=="protein_coding"]             #protein coding sequence only
    
    cds.list=list()
    for (i in 1:length(id_trans)) {
      ext_1 = paste("/sequence/id/",id_trans[[i]],"?type=cds",sep = "")
      r_1   = GET(paste(server, ext_1, sep = ""), content_type("text/x-fasta"))
      if(r_1$status_code==400){
        next
      }else{
        stop_for_status(r_1)
        cds_  = (content(r_1))
        cds.list = c(cds.list,cds_)
      }
    }
    
    output="../Data/cds/"
    write.table(cds.list,paste0(output,name,".fa"),quote = FALSE,
                row.names=FALSE,col.names = FALSE,append = FALSE,sep = "\n")
  }
  
}

args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2])

