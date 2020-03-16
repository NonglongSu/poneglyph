#!/bin/bash


##Pretest
RSCRIPT = Rscript --vanilla
MAFFT = ginsi --preservecase

default:all
all: 

.PHONY:default all clean
.PHONY:homoCall cds_exon_call align_prep align_mafft 		exon_loc fhi fhi_parsi fhi_jc69_hmm

namelist:=$(shell ls Align_prep/) 


Align_mafft:  align_mafft
align_mafft: $(patsubst %,Align_mafft/%,$(namelist))


################################################################

#Find homologs in pairwise species(human/chimp) 
homoCall:../Script/homo_Id_call_2.R 
	$(RSCRIPT) $< geneId.txt

#Extract all human/chimp alternative cds and move to cds/
cds_exon_call:../Script/cds_exon_call.R  geneId_cut.txt
	grep '$$*' $(word 2,$^) | parallel --colsep '\t' $(RSCRIPT) $< {1} {2}


#Reorganize the human/chimp isoforms for alignment
align_prep:../Script/align_prep.R geneId_cut.txt cds/ Align_prep/
	$(RSCRIPT) $< $(word 2,$^)  $(word 3,$^) $(word 4,$^) 


#Multiple (alternative isoforms) alignment 
Align_mafft/%.fa: Align_prep/%.fa
	$(MAFFT) $< > $@



################################################################

#Exon-order-loc determination (Randomly pick one gene)
exon_loc:../Script/exon_loc_map.R Exon_table/ geneId_cut.txt Results/From/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) 








#Try different models to find the best homo isoform

fhi:../Script/find_homo_isoform.R ../Script/fhi_gap_sub.R geneId_cut.txt Align_mafft/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) Results/iso_align.txt


fhi_parsi:../Script/find_homo_isoform.R ../Script/fhi_parsimony.R geneId_cut.txt Align_mafft/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) Results/iso_parsi.txt

fhi_jc69_hmm:../Script/find_homo_isoform.R ../Script/fhi_jc69.R geneId_cut.txt Align_mafft/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) Results/iso_jc69_hmm.txt





