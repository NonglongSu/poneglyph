

#Extract all geneId/transcriptId with the tag" appris_principal_1"
#Note: some genes owns more than one appris_1 transcripts. 
cat gencode.v31.annotation.gtf | grep -v '^#' | grep  'appris_principal_1' | cut -f9 | awk -F\" '{print $2 "\t" $4}'| sort | uniq | awk -F. '{print $1 "\t" $2}' | cut -f1,3 > appris_Id.txt 
