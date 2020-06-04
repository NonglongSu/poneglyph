# Reorganize the isoforms for alignment
# Three arguments: 
# Input File, Input Directory of Coding Sequences, Output Directory

library(tidyverse)
library(seqinr)
library(Biostrings)
library(foreach)
library(doParallel)

#setwd("~/Dropbox (ASU)/poneglyph/Script")

# inFile = "../Data/geneId_update.txt"
# inDir  = "../Data/cds/"
# ouDir  = "../Data/Align_prep/"

main = function(inFile, inDir, ouDir) {
  numCores <- detectCores()
  nList = read_delim(inFile, "\t", col_names = FALSE)
  registerDoParallel(numCores)
  foreach (i = 1:length(nList[[1]])) %dopar% {
    geneA = paste0(inDir, nList[[1]][i], ".fa")
    geneB = paste0(inDir, nList[[2]][i], ".fa")
    dnaA  = readDNAStringSet(geneA, format = "fasta")
    dnaB  = readDNAStringSet(geneB, format = "fasta")
    namA  = names(dnaA)
    namB  = names(dnaB) 
    for (j in 1:length(dnaA)) {
      for (k in 1:length(dnaB)) {
        dnaAB = list(dnaA[j], dnaB[k])
        namAB = c(namA[j], namB[k])
        write.fasta(sequences = dnaAB, 
                    names = namAB, 
                    nbchar = 80,
                    open = "w", 
                    as.string = TRUE, 
                    file.out = paste0(ouDir, nList[[1]][i], ".", j, "_", k, ".fa"))
      }
    }
  }
}

args = commandArgs(trailingOnly = TRUE)
main(args[1], args[2], args[3])