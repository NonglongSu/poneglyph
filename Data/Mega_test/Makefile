#!/bin/bash

##Megatest
RSCRIPT = Rscript --vanilla

default:all
all: 

.PHONY:default all clean
.PHONY:exon_loc fhi_parsi 

################################################################

Data = ../geneId_update.txt
D1   = ../Exon_table/
D2   = ../Align_mafft/

namelist := $(shell cut -f1 ../geneId_update.txt) 

Results/From/Exon_loc: exon_loc
Results/To/Iso_id:     fhi_parsi

exon_loc:  $(patsubst %, Results/From/Exon_loc/%.txt, $(namelist))
fhi_parsi: $(patsubst %, Results/To/Iso_id/%.txt,     $(namelist))


#Exon-order-loc determination (genome-wide database 2000)
Results/From/Exon_loc/%.txt: ../../Script/exonLoc_net.R $(Data) $(D1) $(D2)
	$(RSCRIPT) $< $@ $(word 2,$^) $(word 3,$^) $(word 4,$^) Results/From/Exon_id/ Results/From/Exon_loc/


#Try different models/methods to approximate the best homo-isoform
Results/To/Iso_id/%.txt: ../../Script/exon_mapping_net.R ../../Script/fhi_parsimony.R  $(Data) $(D2) Results/From/Exon_id/ Results/From/Exon_loc/
	$(RSCRIPT) $< $(word 2,$^) $@ $(word 3,$^) $(word 4,$^) $(word 5,$^) $(word 6,$^) $(word 7,$^) Results/To/Iso_id/ Results/To/Iso_seq/ 

