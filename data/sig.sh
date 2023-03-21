mkdir -p vcf
for i in `ls *_shared.temp`
do

NAME=`echo ${i} | awk -F'_shared' '{print $1}'`
cat header > vcf/${NAME}_s.vcf
cat header > vcf/${NAME}_u.vcf
cat ${NAME}_shared.temp | awk '{OFS="\t";print $1,$2,".",$3,$4,".",".","."}'  >> vcf/${NAME}_s.vcf 
cat ${NAME}_uniqued.temp | awk '{OFS="\t";print $1,$2,".",$3,$4,".",".","."}'  >> vcf/${NAME}_u.vcf 


done
