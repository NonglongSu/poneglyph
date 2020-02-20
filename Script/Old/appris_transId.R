#Extract the rest of genes with appris_principal_1

# setwd("~/Dropbox (ASU)/Indel_project/test_human_mouse_rat/Script")

library(readr)
library(Biostrings)


file.1 ="../Raw_data/Appris_tag/Appris_1.txt"
file.2 ="../Raw_data/Appris_tag/Appris_2.txt"
file.3 ="../Raw_data/Appris_tag/Appris_3.txt"
file.4 ="../Raw_data/Appris_tag/Appris_4.txt"
file.5 ="../Raw_data/Appris_tag/Appris_5.txt"

inFile ="../Raw_data/humanId_no_MANE.txt"
ouFile ="../Raw_data/humanId_cano.txt"

ouFi.1 = "../Raw_data/ApprisId_1.txt"
ouFi.2 = "../Raw_data/ApprisId_2.txt"
ouFi.3 = "../Raw_data/ApprisId_3.txt"
ouFi.4 = "../Raw_data/ApprisId_4.txt"
ouFi.5 = "../Raw_data/ApprisId_5.txt"

File.in = c(file.1,file.2,file.3,file.4,file.5,inFile)
File.ou = c(ouFi.1,ouFi.2,ouFi.3,ouFi.4,ouFi.5,ouFile)

#####################################
appris.df = list()
  for (i in 1:length(File.in)) {
    appris.Id      = read_tsv(File.in[i],col_names = FALSE)
    appris.df[[i]] = data.frame(matrix(unlist(appris.Id),nrow=nrow(appris.Id),byrow=FALSE))
  }
  
   geneId.list = appris.df[[6]]
   
   # The results show us there is no coverage of appris_x. 
appr.df = list()
   for (j in 1:(length(appris.df)-1) ) {
     geneId.appr = geneId.list[(geneId.list[,1]  %in% appris.df[[j]][,1]), ]
     geneId.rest = geneId.list[!(geneId.list[,1] %in% appris.df[[j]][,1]), ]
    
     appr.df[[j]] = geneId.appr 
     geneId.list  = geneId.rest
  
  }
  
   for (k in 1:length(appr.df)) {
     write.table(appr.df[[k]],File.ou[[k]],quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = FALSE)
  }
  
     write.table(geneId.list,File.ou[[6]],quote = FALSE, sep="\t",row.names = FALSE,col.names = FALSE,append = FALSE)



  