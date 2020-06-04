# Poneglyph
Piecing all exons together to form an **imaginary isoform** is a "fairytale".  
This fairytale will come true in our **Poneglyph** project, just like piecing all stories together to form a "One piece". 

## Requirements
* R 3.4.4
* Bioconductor 3.6
* foeach 1.5.0
* doParallel 1.0.15

# Workflow
### 1. Find homologs in pairwise species (human/chimp) 
make homoCall
* 16471 orthologs

(Currently geneId.txt shrinks to 2000 samples to in geneID_cut.txt)

### 2. Extract all human/chimp alternative coding sequences
make cds_exon_call
* 3998 genes

### 3. Reorganize the human/chimp isoforms for alignment
make align_prep

### 4. Check for multiple (alternative isoforms) alignments using exon locations
make Align_mafft
* 29435 alignments

### 5. Exon order location determination
make exon_loc

### 6. Try different models to find the best homo isoform
* Parsimony
  * make fhi_parsi
* JC69
  * make fhi_jc69_hmm
* K2P
  * make fhi_k80_hmm

