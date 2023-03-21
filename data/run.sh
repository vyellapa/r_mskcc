SAMPLE=$1
NOS=`ls ${SAMPLE}*.vcf |wc -l|awk '{print $1}'`
cat ${SAMPLE}*.vcf|grep -v "^#" |awk -F'\t' '{OFS="\t";print $1,$2,$4,$5}'|sort |uniq -c |awk -v nos=${NOS} '{OFS="\t"; if($1==nos) {print $2,$3,$4,$5}}' > ${SAMPLE}_shared.temp 
cat ${SAMPLE}*.vcf | grep -v "^#" | awk -F'\t' '{OFS="\t";print $1,$2,$4,$5}'|sort |uniq -c |awk -v nos=${NOS} '{OFS="\t"; if($1!=nos) {print $2,$3,$4,$5}}' > ${SAMPLE}_uniqued.temp

cat ${SAMPLE}_shared.temp | awk '{if(($3=="C" && $4=="T") || ($3=="G" && $4=="A")) {print $0}}' > ${SAMPLE}_CT.txt
cat ${SAMPLE}_shared.temp | awk '{if(($3=="C" && $4=="G") || ($3=="G" && $4=="C")) {print $0}}' > ${SAMPLE}_CG.txt
cat ${SAMPLE}_shared.temp | awk '{if(($3=="C" && $4=="A") || ($3=="G" && $4=="T")) {print $0}}' > ${SAMPLE}_CA.txt
cat ${SAMPLE}_shared.temp | awk '{if(($3=="T" && $4=="A") || ($3=="A" && $4=="T")) {print $0}}' > ${SAMPLE}_TA.txt
cat ${SAMPLE}_shared.temp | awk '{if(($3=="T" && $4=="C") || ($3=="A" && $4=="G")) {print $0}}' > ${SAMPLE}_TC.txt
cat ${SAMPLE}_shared.temp | awk '{if(($3=="T" && $4=="G") || ($3=="A" && $4=="C")) {print $0}}' > ${SAMPLE}_TG.txt

mkdir -p unique
cat ${SAMPLE}_uniqued.temp | awk '{if(($3=="C" && $4=="T") || ($3=="G" && $4=="A")) {print $0}}' > unique/${SAMPLE}_CT.txt
cat ${SAMPLE}_uniqued.temp | awk '{if(($3=="C" && $4=="G") || ($3=="G" && $4=="C")) {print $0}}' > unique/${SAMPLE}_CG.txt
cat ${SAMPLE}_uniqued.temp | awk '{if(($3=="C" && $4=="A") || ($3=="G" && $4=="T")) {print $0}}' > unique/${SAMPLE}_CA.txt
cat ${SAMPLE}_uniqued.temp | awk '{if(($3=="T" && $4=="A") || ($3=="A" && $4=="T")) {print $0}}' > unique/${SAMPLE}_TA.txt
cat ${SAMPLE}_uniqued.temp | awk '{if(($3=="T" && $4=="C") || ($3=="A" && $4=="G")) {print $0}}' > unique/${SAMPLE}_TC.txt
cat ${SAMPLE}_uniqued.temp | awk '{if(($3=="T" && $4=="G") || ($3=="A" && $4=="C")) {print $0}}' > unique/${SAMPLE}_TG.txt

#less I-H-106917-T2-1-D1-2_vs_I-H-106917-N1-1-D1-2.flagged.annot.indeltype.output.txt |awk -F'\t' '{OFS="\t";if($2=="I") {a=1;} if($2=="D") {a=-1;} if($2=="DI") {a=0} if(a==1) print "chr"$3,$6,$6+3000000,a}'|grep -v "ACCE" > I-H-106917_shared.ins.txt
#less I-H-106917-T2-1-D1-2_vs_I-H-106917-N1-1-D1-2.flagged.annot.indeltype.output.txt |awk -F'\t' '{OFS="\t";if($2=="I") {a=1;} if($2=="D") {a=-1;} if($2=="DI") {a=0} if(a==-1 || a==0) print "chr"$3,$6,$6+3000000,a}'|grep -v "ACCE" > I-H-106917_shared.dels.txt

