# Extract the exon Id from ensembl API.
# Using 'repeat' deal with "too many request"
# Two arguments required: 

library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)
library(BiocGenerics) 
library(stringr)
# setwd("~/Dropbox (ASU)/poneglyph/Script")

# Name1 = "ENSG00000003402"
# Name2 = "ENSPTRG00000012798"
# name = "ENSG00000284690"
# name = "ENSG00000003402"

extractAPI = function(s,x){
  R = tryCatch(
    GET(paste(s, x, sep = ""), content_type("application/json")),
    error = function(err) {
      message(err)
      return(NULL)
    }
  )
  stop_for_status(R)
  return(R)
}

# file = "../Data/geneId_cut.txt"
# output1 = "../Data/cds/"
# output2 = "../Data/Exon_table/"
main = function(n, file, output1, output2) {
    
    geneId = read_delim(file, "\t", col_names = FALSE)  
  
    n          = str_extract(basename(n), "[^.]+")
    gene.stem  = geneId[which(geneId[[1]] == n), ]
    
    
    for (i in 1:length(gene.stem)) {
    
      server = "http://rest.ensembl.org"
      name   = gene.stem[[i]]
      ext    = paste("/lookup/id/", name, "?expand=1", sep = "")
      
      repeat{
        r = extractAPI(server, ext)
        if(!is.null(r)){
          break
        }
      }
      
      id        = fromJSON(toJSON(content(r)))
      id_trans  = id$Transcript[["id"]]
      id_trans  = id_trans[id$Transcript$biotype == "protein_coding"]
      id_exon   = id$Transcript$Exon[id$Transcript$biotype == "protein_coding"]
  
      # call the cds && exons
      if(length(id_trans) == 0){# no trans_id exists.
        if(i == 2){
          file.remove(paste0(output1, gene.stem[[1]], ".fa"))
          file.remove(paste0(output2, gene.stem[[1]], ".txt"))
        }
        quit()
      }else{
        exon.list = list()
        cds.list  = list()
        for (j in 1:length(id_trans)) {
          ext_1 = paste("/sequence/id/", id_trans[[j]], "?type=cds",sep = "")
          r_1   = GET(paste(server, ext_1, sep = ""), content_type("text/x-fasta"))
          stop_for_status(r_1)
          
          if (r_1$status_code == 400) {# keep the exon and transcript synchronized.
            next
          }else {
            exon.list[[j]] = unlist(id_exon[[j]]$id)
            cds_     = (content(r_1))
            cds.list = c(cds.list, cds_)
          }
        }
        exon.lt = Filter(Negate(is.null), exon.list)
        df = do.call(rbind, Map(data.frame, Transcript = id_trans, Exon = exon.lt))
        ##########
        write.table(cds.list, paste0(output1, name, ".fa"), quote = FALSE,
                    row.names=FALSE, col.names = FALSE, append = FALSE, sep = "\n")
        write.table(df, paste0(output2, name,".txt"),quote = FALSE,
                    row.names=FALSE, col.names = TRUE, append = FALSE, sep = "\t")
      }
      
    }
}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3], args[4])