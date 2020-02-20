#Extract the canonical transcript from ensembl API.

#setwd("~/Dropbox (ASU)/Indel_project/Script")
library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(Biostrings)
library(readr)

server = "http://rest.ensembl.org"

inFile="../Raw_data/geneId_cano.txt"
ouDir="../Raw_data/ref_can_cds/"                        

#name="ENSG00000158796" 

cds_call = function(name){
      name = gsub('\"', "", name, fixed = TRUE)
      ext  = paste("/lookup/id/",name,"?expand=1",sep = "")
      r    = GET(paste(server, ext, sep = ""), content_type("application/json"))
      stop_for_status(r)
      id   = fromJSON(toJSON(content(r)))
      id_  = id$Transcript[["id"]][id$Transcript$is_canonical==1]
      #print(id)
      
      ext_1 = paste("/sequence/id/",id_,"?type=cds",sep = "")
      r_1   = GET(paste(server, ext_1, sep = ""), content_type("text/x-fasta"))
      stop_for_status(r_1)
      cds_  = (content(r_1))
      
      seq.list = c(name,cds_)
      return(seq.list)
           
}

main = function(inFile,ouDir){
  
  cano.Id = read_tsv(inFile,col_names = FALSE)
  for (i in 1:nrow(cano.Id[,1])) {
    seq.out = cds_call(cano.Id[i,1])
    write(seq.out[2],paste0(ouDir,seq.out[1],".fa"),append = FALSE,sep = "\t")
  }
  
}



args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2])


