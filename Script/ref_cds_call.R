#Extract the refer gene's mane cds from ensembl API

# setwd("~/Dropbox (ASU)/poneglyph/Script")
library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)


# file.1 = "../Raw_data/transId.txt"
# ouFile="../Raw_data/transId_cds.txt"

#Name="ENST00000373115"

main = function(Name){

  server = "https://rest.ensembl.org"
 
  name = gsub('\"',"",Name,fixed=TRUE)
  ext  = paste("/sequence/id/",name,"?type=cds",sep = "")
  r    = GET(paste(server, ext, sep = ""), content_type("text/x-fasta"))
  stop_for_status(r)
  cds  = (content(r))
  
  output="../Raw_data/ref_cds/"
  write.table(cds,paste0(output,name,".fa"),quote = FALSE,
              row.names=FALSE,col.names = FALSE,append = FALSE,sep = "\t")
}
    
args=commandArgs(trailingOnly = TRUE)
main(args[1])
  

