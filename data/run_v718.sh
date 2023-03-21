#Usage
# ./run_v2.sh I-H-106917-T2
#(base) -bash-3.2$ ./run_v2.sh I-H-130719-T1
#(base) -bash-3.2$ ./run_v2.sh I-H-130720-T1
#(base) -bash-3.2$ ./run_v2.sh I-H-130718-T1

SAMPLE="I-H-130718-T1"
NOS=`ls ${SAMPLE}*.vcf |wc -l|awk '{print $1}'`
DP_SAMP=`echo ${SAMPLE} |awk '{a=substr($1,1,10); print a}'`
cat ${SAMPLE}*.vcf|grep -v "^#" |awk -F'\t' '{OFS="\t";print $1,$2,$4,$5}'|sort |uniq -c |awk -v nos=${NOS} '{OFS="\t"; {print $2,$3-1,$3,$4,$5}}' > ${SAMPLE}_all.bed 

cat ${SAMPLE}*.vcf|grep -v "^#" |awk -F'\t' '{OFS="\t";print $1,$2,$4,$5}'|sort |uniq -c |awk -v nos=${NOS} '{OFS="\t"; {print $1,$2,$3,$4,$5}}' > ${SAMPLE}_all.bedc 
less dp/${SAMPLE}*_10000iters_1000burnin_bestConsensusAssignments.bed |grep -v "^chr"|cut -f1,2,3,4 |sort -u > temp.bed
bedtools intersect -a ${SAMPLE}_all.bed -b temp.bed -wo |sort -u > ${SAMPLE}_all.txt

TRUNK=`cat ${SAMPLE}_all.txt |cut -f9 |sort |uniq -c |sort -nr|awk '{print $2}' |head -n1`


less ${SAMPLE}_all.txt |awk -v cl=${TRUNK} -F'\t' '{OFS="\t"; if($9==cl) {print $1,$3,$4,$5}}' > ${SAMPLE}.trunk
less ${SAMPLE}_all.bedc | awk '{if($2=="X" || $2=="Y") print $0}' | awk -v nos=${NOS} -F'\t' '{OFS="\t"; if($1==nos) {print $2,$3,$4,$5}}' >> ${SAMPLE}.trunk
less ${SAMPLE}_all.txt |awk -v cl=${TRUNK} -F'\t' '{OFS="\t"; if($9!=cl && $9!=2) {print $1,$3,$4,$5}}' > ${SAMPLE}.subclones
less ${SAMPLE}_all.bedc | awk '{if($2=="X" || $2=="Y") print $0}' | awk -v nos=${NOS} -F'\t' '{OFS="\t"; if($1!=nos) {print $2,$3,$4,$5}}' >> ${SAMPLE}.subclones

cat ${SAMPLE}.trunk | awk '{if(($3=="C" && $4=="T") || ($3=="G" && $4=="A")) {print $0}}' > ${SAMPLE}_CT.txt
cat ${SAMPLE}.trunk | awk '{if(($3=="C" && $4=="G") || ($3=="G" && $4=="C")) {print $0}}' > ${SAMPLE}_CG.txt
cat ${SAMPLE}.trunk | awk '{if(($3=="C" && $4=="A") || ($3=="G" && $4=="T")) {print $0}}' > ${SAMPLE}_CA.txt
cat ${SAMPLE}.trunk | awk '{if(($3=="T" && $4=="A") || ($3=="A" && $4=="T")) {print $0}}' > ${SAMPLE}_TA.txt
cat ${SAMPLE}.trunk | awk '{if(($3=="T" && $4=="C") || ($3=="A" && $4=="G")) {print $0}}' > ${SAMPLE}_TC.txt
cat ${SAMPLE}.trunk | awk '{if(($3=="T" && $4=="G") || ($3=="A" && $4=="C")) {print $0}}' > ${SAMPLE}_TG.txt

mkdir -p unique
#less I-H-106917-T2_all.txt |awk -v cl=${TRUNK} -F'\t' '{OFS="\t"; if($9!=cl) {print $0}}' > ${SAMPLE}.clones
#rm -f t
#for i in `cat dp/${DP_SAMP}.tsv |cut -f1|sort -u`
#do

#less ${SAMPLE}.clones |awk -v cl=${i} -F'\t' '{OFS="\t"; if($9==cl) {print $0}}' >> t

#done

#cat t | awk -F'\t' '{OFS="\t"; {print $1,$3,$4,$5}}' > ${SAMPLE}.clones


cat ${SAMPLE}.subclones | awk '{if(($3=="C" && $4=="T") || ($3=="G" && $4=="A")) {print $0}}' > unique/${SAMPLE}_CT.txt
cat ${SAMPLE}.subclones | awk '{if(($3=="C" && $4=="G") || ($3=="G" && $4=="C")) {print $0}}' > unique/${SAMPLE}_CG.txt
cat ${SAMPLE}.subclones | awk '{if(($3=="C" && $4=="A") || ($3=="G" && $4=="T")) {print $0}}' > unique/${SAMPLE}_CA.txt
cat ${SAMPLE}.subclones | awk '{if(($3=="T" && $4=="A") || ($3=="A" && $4=="T")) {print $0}}' > unique/${SAMPLE}_TA.txt
cat ${SAMPLE}.subclones | awk '{if(($3=="T" && $4=="C") || ($3=="A" && $4=="G")) {print $0}}' > unique/${SAMPLE}_TC.txt
cat ${SAMPLE}.subclones | awk '{if(($3=="T" && $4=="G") || ($3=="A" && $4=="C")) {print $0}}' > unique/${SAMPLE}_TG.txt

#less I-H-106917-T2-1-D1-2_vs_I-H-106917-N1-1-D1-2.flagged.annot.indeltype.output.txt |awk -F'\t' '{OFS="\t";if($2=="I") {a=1;} if($2=="D") {a=-1;} if($2=="DI") {a=0} if(a==1) print "chr"$3,$6,$6+3000000,a}'|grep -v "ACCE" > I-H-106917_shared.ins.txt
#less I-H-106917-T2-1-D1-2_vs_I-H-106917-N1-1-D1-2.flagged.annot.indeltype.output.txt |awk -F'\t' '{OFS="\t";if($2=="I") {a=1;} if($2=="D") {a=-1;} if($2=="DI") {a=0} if(a==-1 || a==0) print "chr"$3,$6,$6+3000000,a}'|grep -v "ACCE" > I-H-106917_shared.dels.txt

