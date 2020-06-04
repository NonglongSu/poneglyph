#!/bin/bash

##Pretest
RSCRIPT = Rscript --vanilla
MAFFT   = ginsi --preservecase
Scp     = ../Script/One_test

default:all
all: 

.PHONY:default all clean
.PHONY:homoCall cds_exon_call align_prep align_mafft 			
.PHONY:exon_loc fhi fhi_parsi fhi_jc69_hmm


Database := one.txt
namelist := $(shell ls Align_prep/) 


Align_mafft:  align_mafft
align_mafft: $(patsubst %,Align_mafft/%,$(namelist))


################################################################

#Find 1 homolog in pairwise species(human/chimp) 
homoCall: $(Scp)/find_homo_Id.R 
	$(RSCRIPT) $< one.txt

#Extract all human/chimp alternative cds
cds_exon_call: $(Scp)/cdsExon_call.R  $(Database)
	grep '$$*' $(word 2,$^) | parallel --colsep '\t' $(RSCRIPT) $< {1} {2}


#Reorganize the human/chimp isoforms for alignment
align_prep:../Script/align_prep.R $(Database) 
	$(RSCRIPT) $< $(word 2,$^)  cds/ Align_prep/ 


#Multiple (alternative isoforms) alignment 
Align_mafft/%.fa: Align_prep/%.fa
	$(MAFFT) $< > $@


################################################################

#Exon-order-loc determination (Randomly pick one gene from the database we generated)
exon_loc:../Script/exon_loc_map.R Exon_table/ geneId_cut.txt Results/From/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) 


#Try different models/methods to approximate the best homo-isoform
file=Results/From/geneId.txt
f1  =Results/From/human_exonId.txt
f2  =Results/From/human_exonLoc.txt
f3  =Results/From/chimp_exonId.txt
f4  =Results/From/chimp_exonLoc.txt
inD =Align_mafft/
ouF1=exon_map_ref.fa
ouF2=exon_map_tes.fa
ouF3=exon_map.txt

fhi_parsi: $(Scp)/exon_mapping_raw.R ../Script/fhi_parsimony.R $(file) $(f1) $(f2) $(f3) $(f4) $(inD)  
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) $(word 6,$^) $(word 7,$^) $(word 8,$^) Results/To/Parsi/$(ouF1) Results/To/Parsi/$(ouF2) Results/To/Parsi/$(ouF3)


fhi_jc69_hmm: $(Scp)/exon_mapping_raw.R ../Script/fhi_jc69_hmm.R $(file) $(f1) $(f2) $(f3) $(f4) $(inD)  
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) $(word 6,$^) $(word 7,$^) $(word 8,$^) Results/To/JC_69/HMM/$(ouF1) Results/To/JC_69/HMM/$(ouF2) Results/To/JC_69/HMM/$(ouF3)


fhi_k80_hmm: $(Scp)/exon_mapping_raw.R ../Script/fhi_k80_hmm.R $(file) $(f1) $(f2) $(f3) $(f4) $(inD)  
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) $(word 6,$^) $(word 7,$^) $(word 8,$^) Results/To/K80/HMM/$(ouF1) Results/To/K80/HMM/$(ouF2) Results/To/K80/HMM/$(ouF3)


