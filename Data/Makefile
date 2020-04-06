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
file=Results/From/geneId.txt
f1  =Results/From/human_exonId.txt
f2  =Results/From/human_exonLoc.txt
f3  =Results/From/chimp_exonId.txt
f4  =Results/From/chimp_exonLoc.txt
inD =Align_mafft/
ouF1=exon_map_ref.fa
ouF2=exon_map_tes.fa
ouF3=exon_map.txt

fhi_parsi:../Script/exon_mapping_raw.R ../Script/fhi_parsimony.R $(file) $(f1) $(f2) $(f3) $(f4) $(inD)  
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) $(word 6,$^) $(word 7,$^) $(word 8,$^) Results/To/Parsi/$(ouF1) Results/To/Parsi/$(ouF2) Results/To/Parsi/$(ouF3)

fhi_jc69_hmm:../Script/exon_mapping_raw.R ../Script/fhi_jc69.R geneId_cut.txt Align_mafft/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) Results/iso_jc69_hmm.txt





