
#Extract all geneId/transcriptId with the tag" appris_principal_1"
#Note: some genes owns more than one appris_1 transcripts. 

DATABASE=gencode.v31.annotation.gtf


for nums in 1 2 3 4 5;	do\
	pattern=appris_principal_${nums}
	cat ${DATABASE} | grep -v '^#' | grep  "$pattern" | cut -f9 | awk -F\" '{print $2 "\t" $4}'| sort | uniq\
        	                       | awk -F. '{print $1 "\t" $2}' | cut -f1,3  > Appris_tag/Appris_${nums}.txt
done 
