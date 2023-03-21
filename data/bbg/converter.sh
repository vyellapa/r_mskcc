for i in `ls *ed.txt`
do

less ${i} | awk -F'\t' '{OFS="\t";if(NR==1) {$1=""} print $0}' > t

mv t ${i} 

done

