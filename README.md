# Poneglyph

## Requirements
* R 3.4.4
* Bioconductor 3.6
* foeach 1.5.0
* doParallel 1.0.15

# Workflow
### 1. Find homologs in pairwise species (human/chimp) 
make homoCall

(Currently geneId.txt is cut to 200 samples to obtain geneID_cut.txt)

### 2. Extract all human/chimp alternative coding sequences
make cds_exon_call
### 3. Reorganize the human/chimp isoforms for alignment
make align_prep
### 4. Check for multiple (alternative isoforms) alignments using exon locations
make Align_mafft
### 5. Exon order location determination
make exon_loc
### 6. Try different models to find the best homo isoform
* Parsimony
  * make fhi_parsi
* JC69
  * make fhi_jc69

