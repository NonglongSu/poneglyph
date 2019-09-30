#!/bin/bash

RSCRIPT = Rscript --vanilla




default:all
all: 

.PHONY:default all clean
.PHONY:geneCall cdsCall cdsFCall refList

geneCall:../Script/homo_geneId.R 
	$(RSCRIPT) $< geneId_with_MANE.txt transId_wtih_MANE.txt geneId_no_MANE.txt


NAMELIST:=$(shell test -f transId.1.txt && cut -f3 transId.1.txt )
cdsCall:  $(patsubst %,ref_cds/%.fa,$(NAMELIST)) 


ref_cds/%.fa:../Script/ref_cds_call.R | transId_with_MANE.txt
	grep '$*' $(word 2,$^) | parallel --colsep '\t' $(RSCRIPT) $< {3}

cdsFCall:../Script/mul_cds_call.R geneId_with_MANE.txt 
	grep '$*' $(word 2,$^) | parallel --colsep '\t' $(RSCRIPT) $< {2}

###########################################

#Usage: bash nameCheck.sh  trans.ref.txt ref_cds/
#Usage: bash emptyCheck.sh ref_cds/
#Usage: bash emptyCheck.sh focal_cds/
refList: 
	cat transId.1.txt | cut -f3 | awk 'NF{print $$0 ".fa"}' > trans.ref.txt

#Usage: bash appris_get.sh (Extrat all human appris_1 Ids)	
