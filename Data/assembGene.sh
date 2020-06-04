#
data=geneId_cut.txt
dir=Exon_table
#
ls ${dir}/ | grep 'ENSG' | cut -f 1 -d '.' > file1
#
cat ${data} | cut -f1 | sort > file2
#
grep -F -x -v -f file1 file2 > Diff
#
grep -v -f Diff ${data} > $1

#
rm file1 file2 Diff






