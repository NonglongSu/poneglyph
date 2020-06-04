#Extract the IDs and locations of exons from all reference genes in the database(2000).

#setwd("~/Dropbox (ASU)/poneglyph/Script")

library(methods)
library(httr)
library(jsonlite)
library(xml2)
library(readr)
library(Biostrings)
library(BiocGenerics) 
library(dplyr) 
library(purrr) # accumulate
library(stringr)

findExonLoc = function(TRANS.ID, Align.set, flag){
  exon.Df = list()
  cds.Box = list()
  j = 0
  T.id = TRANS.ID
  for(i in 1:length(TRANS.ID)){
    name = gsub('\"', "", TRANS.ID[i], fixed = TRUE)
    ext  = paste("/lookup/id/", name, "?expand=1", sep = "")
    r    = GET(paste(server, ext, sep = ""), content_type("application/json"))
    warn_for_status(r)
    id   = fromJSON(toJSON(content(r)))
    if(is.null(id$error) == FALSE){# http 400 error "bad request"
      j = c(j,i)
      next 
    }
    #locate translation range
    cds.box = c(id$Translation$start,id$Translation$end)
    #extract exon locus
    id_exon = id$Exon
    exon.df = do.call(rbind, Map(data.frame, start = id_exon$start, end = id_exon$end))
    exon.df = exon.df[order(exon.df$start), ]
    
    exon.Df[[i]] = exon.df
    cds.Box[[i]] = cds.box
    
  }
  # Quality control
  if(length(j) > 1){
    T.id = TRANS.ID[-j] 
    for (w in 1:length(j)) {
      if(flag == 1){
        mAtch  = as.numeric(unlist(lapply((regmatches(Align.set, regexpr("[0-9]+_", Align.set))), function(x){str_extract(x, "[^_]+")})))
      }else{
        mAtch  = as.numeric(unlist(lapply(Align.set, function(x){str_extract(gsub("^.*\\_", "", x), "[^.]+")})))
      }
      rm.set = Align.set[which(mAtch == w)]
      lapply(rm.set, function(x){file.remove(paste0(Dir, x))})   
    }
  }
  
  E.df = Filter(Negate(is.null), exon.Df)
  C.ds = Filter(Negate(is.null), cds.Box)
  res  = list(E.df, C.ds, T.id)
  return(res)
  
} 

find_id_loc = function(ETable, EXON.ID, E.Loc) {
  
  Exon.df  = list()
  Exon.loc = list()    #Relative location  
  Exon.id  = list()    #Functional Id 
  for(i in 1:length(E.Loc[[1]])){
    
    cds.box = E.Loc[[2]][[i]]
    exon.df = E.Loc[[1]][[i]]
    #5'3'UTR locus.
    UTR.5 = which(cds.box[1] >= exon.df$start)
    UTR.3 = which(cds.box[2] <= exon.df$end)
    #Test if start&end points are on the same/different exons .
    if (exon.df$start[head(UTR.3, 1)] <= cds.box[1] && exon.df$end[tail(UTR.5, 1)] >= cds.box[2]) { #translation start/end on the same exon
      exon.df.1 = data.frame(start = cds.box[1], end = cds.box[2])
      cds.start = NULL
      cds.end   = cds.box[2] - cds.box[1] + 1 
    } else { #translation start/end on different exons.
      exon.start = c(cds.box[1], exon.df[1][exon.df[1] > cds.box[1] & exon.df[1] <= cds.box[2]])
      exon.end   = c(exon.df[2][exon.df[2] >= cds.box[1] & exon.df[2] < cds.box[2]], cds.box[2])
      exon.df.1  = data.frame(start = exon.start, end = exon.end)
      cds.start = exon.df$end[tail(UTR.5, 1)] - cds.box[1] + 1 
      cds.end   = cds.box[2] - exon.df$start[head(UTR.3, 1)] + 1 
    }
    
    #find the middle exon locus.
    exon.fix.df  = exon.df[-c(UTR.5,UTR.3),]
    middle.len   = exon.fix.df$end - exon.fix.df$start + 1
    exon.loc     = c(cds.start, c(middle.len,cds.end)) %>% purrr::accumulate(`+`)
    
    #filter the exon id
    exon.id.tmp  = EXON.ID[which(ETable$Transcript == E.Loc[[3]][i])]
    if (length(UTR.5) + length(UTR.3) <= 2) { # all
      exonId.final = exon.id.tmp
    }else { # partial
      exonId.final = exon.id.tmp[-c(head(UTR.5, -1), tail(UTR.3, -1))]
    }
    
    Exon.df[[i]]  = exon.df.1
    Exon.id[[i]]  = exonId.final
    Exon.loc[[i]] = exon.loc  
    
  }
  E.start = lapply(Exon.df, `[[`, 1)
  E.end   = lapply(Exon.df, `[[`, 2)
  Exon.ID = do.call(rbind, Map(data.frame, Transcript = E.Loc[[3]], Exon = Exon.id, Start = E.start, End = E.end))
  
  res = list(Exon.ID, Exon.loc)
  return(res)
}


# file = "../Data/geneId_update.txt"
# inD1 =  "../Data/Exon_table/"
# inD2 =  "../Data/Align_mafft/"
# ouD1 = "../Data/Mega_test/Results/From/Exon_id/"
# ouD2 = "../Data/Mega_test/Results/From/Exon_loc/"

main = function(n, file, inD1, inD2, ouD1, ouD2){
  
  # input: user 
  # Test
  # name = "ENSG00000173757"
  
  # numCores <- detectCores()
  # registerDoParallel(numCores)
  name   = str_extract(basename(n), "[^.]+")
  geneId = read_delim(file, "\t", col_names = FALSE)  
  gene.stem  = geneId[which(geneId[[1]] == name), ]
  
  # list all alignment of that gene.
  align.set = list.files(path = inD2, pattern = gene.stem[[1]])
  Dir <<- inD2
  
  for (i in 1:length(gene.stem)) {
    f1 = paste0(ouD1, gene.stem[[i]])
    f2 = paste0(ouD2, gene.stem[[i]])
    exon.table = read_delim(paste0(inD1, gene.stem[[i]], ".txt"), "\t", col_names = TRUE)
    
    if(file.exists(f1) || file.exists(f2)){
      quit()
    }else{
      
      trans.id = rle(exon.table[[1]])[[2]]
      exon.id  = rle(exon.table[[2]])[[2]]
      server <<- "http://rest.ensembl.org"
      eLoc.list   = findExonLoc(trans.id, align.set, i)
      id.loc.list = find_id_loc(exon.table, exon.id, eLoc.list)
      
      write.table(id.loc.list[[1]],
                  f1,
                  quote = FALSE,
                  row.names = FALSE, 
                  col.names = TRUE, 
                  append = FALSE, 
                  sep = "\t")
      lapply(id.loc.list[[2]], 
             write,
             f2,
             append=TRUE,
             ncolumns = 1000)
    }
  }
  
}

args = commandArgs(trailingOnly = TRUE)
main(args[1],args[2], args[3], args[4], args[5], args[6])