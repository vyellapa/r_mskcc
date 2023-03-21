rm -f *.ins.txt *.del.txt
for i in `ls *.indel.*`
do
INS=`echo ${i} |sed 's/.txt/.ins.txt/g'`
DEL=`echo ${i} |sed 's/.txt/.del.txt/g'`

less ${i} |awk -F'\t' '{OFS="\t";if($2=="I") {print "chr"$3,$6,$6+1000000,"-1"}}' > ${INS}
less ${i} |awk -F'\t' '{OFS="\t";if($2=="D") {print "chr"$3,$6,$6+1000000,"1"}}' > ${DEL}
less ${i} |awk -F'\t' '{OFS="\t";if($2=="DI") {print "chr"$3,$6,$6+1000000,"0"}}' >> ${DEL}



done
