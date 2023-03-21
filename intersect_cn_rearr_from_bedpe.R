# Usage: Rscript /ifs/work/leukgen/home/gg10/soft/intersect_cn_rearr_from_bedpe.R /ifs/res/papaemme/users/gg10/renal/chromophobe/data/JHCHR6/brass/JHCHR6_vs_JHCHN6.r4 /ifs/res/papaemme/users/gg10/renal/chromophobe/data/JHCHR6/brass/JHCHR6.ngscn.segments.abs_cn.bg 10000 /ifs/res/papaemme/users/gg10/renal/chromophobe/data/JHCHR6/brass/extras
args = commandArgs(trailingOnly=TRUE)

rearrs = tools:::file_path_as_absolute(args[1])
cnsegs = tools:::file_path_as_absolute(args[2])
distance = args[3]
output_dir = args[4]
tumour_bam = tools:::file_path_as_absolute(args[5])
normal_bam = tools:::file_path_as_absolute(args[6])
#rearrs = "JHCHR3/brass/JHCHR3_vs_JHCHN3.r4"
#cnsegs = "JHCHR3/brass/JHCHR3.ngscn.segments.abs_cn.bg"
if(!dir.exists(output_dir)) {
	system(paste('mkdir ', output_dir, sep=''))
}
setwd(output_dir)
REARRS_LOWER_FILE = "tmp.rearrs_lower.txt"
REARRS_HIGHER_FILE = "tmp.rearrs_higher.txt"
SEGS_LOWER_FILE = "tmp.segs_lower.txt"
SEGS_HIGHER_FILE = "tmp.segs_higher.txt"
INTERSECT_FILE = "tmp.intersect"
REARRS_WITH_SEGS_FILE = "rearrs_with_cnsegs.txt"
ASSEMBLY_COMMANDS = "assembly_commands.sh"

cmd = paste("zless ", rearrs, " | awk 'OFS=\"\\t\"{print $1,$2-", distance, ",$3+", distance, ",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > ", REARRS_LOWER_FILE, sep='')
print(cmd)
system(cmd)
cmd = paste("zless ", rearrs, " | awk 'OFS=\"\\t\"{print $4,$5-", distance, ",$6+", distance, ",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > ", REARRS_HIGHER_FILE, sep='')
print(cmd)
system(cmd)
cmd = paste("zless ", cnsegs, " | awk 'OFS=\"\\t\"{print $1,$2-", distance, ",$2+", distance, ",$2,$3,$4,$5}' | awk 'OFS=\"\\t\"{if($2<0) print $1,1,$3,$4,$5,$6,$7; else {print $0}}' > ", SEGS_LOWER_FILE, sep='')
print(cmd)
system(cmd)
cmd = paste("zless ", cnsegs, " | awk 'OFS=\"\\t\"{print $1,$3-", distance, ",$3+", distance, ",$2,$3,$4,$5}' | awk 'OFS=\"\\t\"{if($2<0) print $1,1,$3,$4,$5,$6,$7; else {print $0}}' > ", SEGS_HIGHER_FILE)
print(cmd)
system(cmd)

cmd = paste("bedtools intersect -a ", SEGS_LOWER_FILE, " -b ", REARRS_LOWER_FILE, " -wa -wb | awk '{print $1,$4,$5,$6,$7,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' > ", INTERSECT_FILE, sep='')
print(cmd)
system(cmd)
cmd = paste("bedtools intersect -a ", SEGS_LOWER_FILE, " -b ", REARRS_HIGHER_FILE, " -wa -wb | awk '{print $1,$4,$5,$6,$7,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' >> ", INTERSECT_FILE, sep='')
print(cmd)
system(cmd)
cmd = paste("bedtools intersect -a ", SEGS_HIGHER_FILE, " -b ", REARRS_LOWER_FILE, " -wa -wb | awk '{print $1,$4,$5,$6,$7,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' >> ", INTERSECT_FILE, sep='')
print(cmd)
system(cmd)
cmd = paste("bedtools intersect -a ", SEGS_HIGHER_FILE, " -b ", REARRS_HIGHER_FILE, " -wa -wb | awk '{print $1,$4,$5,$6,$7,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' >> ", INTERSECT_FILE, sep='')
print(cmd)
system(cmd)
system(paste("sort -u", INTERSECT_FILE, ">", REARRS_WITH_SEGS_FILE, sep=" "))
system(paste('rm',REARRS_LOWER_FILE, REARRS_HIGHER_FILE, SEGS_LOWER_FILE, SEGS_HIGHER_FILE, INTERSECT_FILE, sep=" "))
selrearrs = read.table(file=REARRS_WITH_SEGS_FILE, header=F, sep=" ")
rearrs = read.table(file=rearrs, header=F, sep="\t")
selrearrs2 = rearrs[rearrs$V7 %in% selrearrs$V12,c(1:10)]
#selrearrs2[str_higher=='-','V10'] = '+'
#selrearrs2[str_higher=='+','V10'] = '-'
#rearrs2 = as.data.frame(cbind(selrearrs2$V1, selrearrs2$V11-1, selrearrs2$V11, selrearrs2$V4, selrearrs2$V12-1, selrearrs2[,4:8]))
rearrs2 = selrearrs2
s1 = seq(from=1,to=nrow(rearrs2),by=40)
if(s1[length(s1)]!=nrow(rearrs2)) {s1 = c(s1,nrow(rearrs2))}
idx = cbind(s1,c(s1[2:length(s1)],0))
for(i in 1:(nrow(idx)-2)) {
	write.table(rearrs2[idx[i,1]:(idx[i,2]-1),], file=paste('assembly_input_',i,'.txt',sep=''), row.names=F, col.names=F, quote=F, sep="\t")
}
i = nrow(idx)-1
write.table(rearrs2[idx[i,1]:(idx[i,2]),], file=paste('assembly_input_',i,'.txt',sep=''), row.names=F, col.names=F, quote=F, sep="\t")

for(infile in list.files('.', pattern='assembly_input')) {
	cmd = paste(
		'brass-assemble -X -m mem -O bedpe -r /ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta -T ', output_dir,
		' -o ', paste(output_dir, '/', sub('input','output',infile), sep=''), ' ', paste(output_dir, '/', infile, sep=''), ' ', tumour_bam,':',tumour_bam,'.bai ', normal_bam,':',normal_bam,'.bai', sep='')
	if(!file.exists(ASSEMBLY_COMMANDS)) {
		write.table(cmd, file=ASSEMBLY_COMMANDS, row.names=F, col.names=F, quote=F)
	} else {
		write.table(cmd, file=ASSEMBLY_COMMANDS, row.names=F, col.names=F, append=T, quote=F)
	}
}
