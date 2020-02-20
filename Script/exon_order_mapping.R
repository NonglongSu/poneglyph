library(seqinr)
library(Biostrings)
library(tidyverse)

setwd("~/Dropbox (ASU)/poneglyph/Script")

inFile1 = "../Data/geneId_cut.txt"
inFile2 = "../Data/trans_align_tb.txt"
inDir   = "../Data/Exon_table/"

main = function(inFile1, inFile2, inDir,ouFile){
  #read from the homo-geneId file
  gene.lst = read_delim(inFile1,"\t", col_names = FALSE)
  tran.lst = read_delim(inFile2,"\t", col_names = FALSE)
}
