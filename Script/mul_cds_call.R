#Extract the exon Id from ensembl API.

# setwd("~/Dropbox (ASU)/poneglyph/Script")
library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)


#Name="ENSMUSG00000046567"

main = function(Name){
  
  server = "http://rest.ensembl.org"
  
  name = gsub('\"', "", Name, fixed = TRUE)
  ext  = paste("/lookup/id/",name,"?expand=1",sep = "")
  r    = GET(paste(server, ext, sep = ""), content_type("application/json"))
  stop_for_status(r)
  id   = fromJSON(toJSON(content(r)))
  id_  = id$Transcript[["id"]]
  
  cds.list=list()
  for (i in 1:length(id_)) {
    ext_1 = paste("/sequence/id/",id_[[i]],"?type=cds",sep = "")
    r_1   = GET(paste(server, ext_1, sep = ""), content_type("text/x-fasta"))
    if(r_1$status_code==400){
      next
    }else{
      stop_for_status(r_1)
      cds_  = (content(r_1))
      cds.list = c(cds.list,cds_)
    }
  }
  
  output="../Raw_data/focal_cds/"
  write.table(cds.list,paste0(output,name,".fa"),quote = FALSE,
              row.names=FALSE,col.names = FALSE,append = FALSE,sep = "\n")
  
}

args = commandArgs(trailingOnly = TRUE)
main(args[1])

