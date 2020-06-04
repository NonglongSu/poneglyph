#!/bin/bash

##Pretest
RSCRIPT = Rscript --vanilla
MAFFT = ginsi --preservecase

default:all
all: 

.PHONY:default all clean
.PHONY:homoCall cds_exon_call align_prep align_mafft 			


Database := $(shell cut -f 1 geneId_cut.txt)
namelist := $(shell ls Align_prep/) 


cds: cds_exon_call
cds_exon_call: $(patsubst %, cds/%.fa, $(Database))

Align_mafft:  align_mafft
align_mafft: $(patsubst %, Align_mafft/%, $(namelist))


################################################################

#Find homologs in pairwise species(human/chimp) 
homoCall:../Script/homo_Id_call_2.R 
	$(RSCRIPT) $< geneId.txt


#update the geneId_cut.txt to obtain new dataset.
geneId_update.txt: assembGene.sh
	bash $< $@

#Extract all human/chimp alternative cds and move to cds/
#grep '$$*' $(word 2,$^) | parallel --colsep '\t' $(RSCRIPT) $< {1} {2}
cds/%.fa: ../Script/cds_exon_call.R  geneId_cut.txt cds/ Exon_table/
	$(RSCRIPT) $< $@ $(word 2,$^) $(word 3,$^) $(word 4,$^)


#Reorganize the human/chimp isoforms for alignment
align_prep:../Script/align_prep.R geneId_update.txt cds/ Align_prep/
	$(RSCRIPT) $< $(word 2,$^)  $(word 3,$^) $(word 4,$^) 


#Multiple (alternative isoforms) alignment 
Align_mafft/%.fa: Align_prep/%.fa
	$(MAFFT) $< > $@



