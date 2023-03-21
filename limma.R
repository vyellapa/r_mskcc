
packages <- c("readr", "tidyr", "splitstackshape", "plyr", "dplyr", "ggplot2", "ggpubr", "reshape2", "magrittr", "biomaRt", "tidyverse", "RColorBrewer", 
              "stringr", "BSgenome.Hsapiens.UCSC.hg19", "stringi", "tibble", "pander", "IRanges", "GenomicRanges", "plotly", "processx", "VennDiagram", 
              "mosaic", "survminer", "survival", "openxlsx", "vsn", "DESeq2", "AnnotationDbi", "limma", "lattice", "org.Hs.eg.db", "MASS", 
              "GenomicFeatures", "rtracklayer", "glmnet", "Hmisc", "ConsensusClusterPlus", "pheatmap", "edgeR", "affy", "rgl", "gridExtra", "valr")

invisible(suppressWarnings(suppressMessages(lapply(packages, library, character.only = TRUE))))

options(scipen=999)
options(repr.matrix.max.cols=70, repr.matrix.max.rows=50)
options(repr.plot.width=4, repr.plot.height=4)

set.seed(42)

source("./functions/utils.R")


### PATHS

# Hotspots
hotspots_matrix_path = "/Users/yellapav/Documents/MSKCC/Projects/SV_project/final_data/final_sv_hotspots.csv"

hotspots_path = "/Users/yellapav/Documents/MSKCC/Projects/SV_project/data/hotspots/190912_manual.hotspots.xlsx"

# SVs
SVs_path = "/Users/yellapav/Documents/MSKCC/Projects/SV_project/data/191120_delly_mapq_60_filtered_baseline.txt"

# CNV data
cnv_path = "../../CoMMpass/IA13/clean/commpass_cnv_new_nov_2019.txt"
minimal_path = '../ref_files/reference_cnv_prediction.txt'

# Canonical translocations and mutational signatures
canonical_tra_path = '/Users/yellapav/Documents/MSKCC/Projects/CoMMpass/IA13/clean/commpass_canonical_translocations_IA13.csv'
myc_nonIG_path = './output/MYC_non_immunoglobulin.txt'

# RNAseq
count_path = "/Users/yellapav/Documents/MSKCC/Projects/CoMMpass/IA13/MMRF_CoMMpass_IA13a_E74GTF_HtSeq_Gene_Counts.txt.gz"
gene_metadata_path = "../data/ENSG_to_ensembl74.RData"

# References
genes_path <- '../ref_files/UCSC_ref_genes.dms'

# load hotspots

hotspot_matrix = read.csv(hotspots_matrix_path, stringsAsFactors = F)

hot_manual = read.xlsx(hotspots_path)
hot_manual = hot_manual %>%
    filter(!grepl('artifact', hotspot_type),
           include_as_hotspot == 1, # remove artifacts
           !annotation %in% c("IGH", "IGL", "IGK")) %>% # remove IG
    mutate(hotspot_class = ifelse(CNV_TYPE == 'none' & limma_driver == 0 | hotspot.id %in% c("nonClustered_chr6_66.3Mb"), "unknown", 
                                     ifelse(CNV_TYPE == "AMP" | !is.na(limma_up) | hotspot.id %in% c("nonClustered_chr3_140.6Mb", "downsampling_chr11_78.1Mb"), "gain",
                                           ifelse(CNV_TYPE == "DEL" | !is.na(limma_down) | hotspot.id %in% c("nonClustered_chr17_39Mb"), "loss", NA)))) # annotate functional category


hot_limma = hot_manual %>%
    mutate(limma_name = paste0(gsub("chr", "", chr), cytoband, "_", hotspot_class),
           limma_name = ifelse(is.na(annotation), limma_name, paste(limma_name, annotation, sep = "_"))) %>%
    dplyr::select(hotspot.id, limma_name)

hot_limma_names = hot_limma$limma_name
names(hot_limma_names) <- hot_limma$hotspot.id

# load svs
SVs = read.delim(SVs_path, stringsAsFactors = F)
SVbp = splitSV(SVs)
SVbp = SVbp %>%
    dplyr::rename(chr = "chrom") %>%
    mutate(chr = paste0("chr", chr))

sv_samples = unique(SVbp$sample)

############################################################
#####           GENE EXPRESSION ANALYSIS
############################################################

##############################################################################
##
## Create desgin with clinical and biological info
##
##############################################################################

# IG translocations and non-IG-MYC
canonical_tra = read.csv(canonical_tra_path, stringsAsFactors = F)
canonical_tra = canonical_tra %>%
            dplyr::rename(t_MMSET = 't_WHSC1') %>%
            dplyr::select(sample, t_CCND1, t_MMSET, t_MAF, t_MAFB, t_MYC)

myc_nonig_tra = read.csv(myc_nonIG_path, stringsAsFactors = F)

canonical_tra = canonical_tra %>%
    left_join(myc_nonig_tra) %>%
    mutate(t_MYCnonIG = ifelse(is.na(t_MYCnonIG), 0, t_MYCnonIG)) %>%
    rename(t_MYC_IG = 't_MYC')

# CNV
cnv = read.delim(cnv_path, stringsAsFactors = F)

cnv_calls = cnvProcess(cnv)

# minimally deleted regions
minimal = read.delim(minimal_path, header=T, stringsAsFactors = F)

minimal = minimal %>%
    mutate(gene = ifelse(gene == "JAG1", band, gene),
           gene = ifelse(gene == "BIRCs", "BIRC", gene),
           code_gene = ifelse(gene != band, paste(code, gene, sep = "_"), code),
           arm = gsub("[0-9]|\\." ,"" , band),
           chrom = as.character(chrom))

# call recurrent CNVs

cnv_matrix = callRecurrent(cnv_calls, regions = minimal)

# Clinical matrix setup

clinical_matrix = cnv_matrix %>%
                dplyr::rename(gain1q21 = gain1q21.3_CKS1B) %>%
                dplyr::select(sample, HRD, gain1q21) %>%
                left_join(canonical_tra,
                             by = 'sample') %>%
                dplyr::select(sample, everything()) %>%
                left_join(hotspot_matrix, by = "sample")

# Remove cases with incomplete data and set all variables > 1 to 1
clinical_matrix = clinical_matrix[complete.cases(clinical_matrix),]

clinical_matrix[-1] = data.frame(apply(clinical_matrix[-1], 2, function(x) ifelse(x > 1, 1, x))) # set values above 1 to 1
              
                                        
##############################################################################
##
## Prepare hotspots
##
##############################################################################

hot_all = hot_manual %>%
                arrange(desc(analysis), chr, start.bp) %>%                        
                dplyr::select(hotspot.id, chr, start.bp, end.bp)

##############################################################################
##
## Upload data Gene expression based - raw counts
##
##############################################################################

RNAcounts = read.delim(count_path, stringsAsFactors = F)

RNAcounts = RNAcounts %>%
    dplyr::select(GENE_ID, which(names(.) %in% clinical_matrix$sample)) %>%
    column_to_rownames(var = 'GENE_ID')

# Gene ID/name conversion -- 57997 transcripts with ENSG ID  >>>> 33260 genes with resolved names
#load(gene_metadata_path) # Load G_list object
#names(G_list) = c('GENE_ID', 'gene')
#RNAcounts = left_join(RNAcounts, G_list)

##############################################################################
##
## definitive files intersecting RNAseq and design
##
##############################################################################
    
count_matrix = as.matrix(RNAcounts)

# Order clinical matrix rows by count materix columns
clinical_matrix = clinical_matrix[order(match(clinical_matrix$sample, colnames(count_matrix))),]
row.names(clinical_matrix) = NULL

# prepare design matrix with factor variables for normalization
design_matrix = clinical_matrix %>%
                        mutate(HRD = as.factor(HRD),
                               t_CCND1 = as.factor(t_CCND1),
                               t_MMSET = as.factor(t_MMSET))

design_matrix = design_matrix %>%
    column_to_rownames(var = 'sample')



##############################################################################
##
## DESeq normalization
##
##############################################################################

# Convert to DESeq matrix
# ddsMat = DESeqDataSetFromMatrix(countData = count_matrix,
#                                  colData = design_matrix,
#                                  design = ~  HRD + t_CCND1 + t_MMSET)

# Filter out genes without expression
# ddsMat = ddsMat[ rowSums(counts(ddsMat)) > 1, ]

# Perform normalization procedure using vst (varianceStabilizingTransformation) from the DESeq2 package
# RNAnorm = vst(ddsMat, blind = FALSE)

# saveRDS(RNAnorm, "../data/190626_vst_commpass.RDS")
RNAnorm = readRDS("../data/190626_vst_commpass.RDS")

exprMatrix = assay(RNAnorm)
exprMatrix = exprMatrix[,colnames(exprMatrix) %in% rownames(design_matrix)]

############################################################################################################################
############################################################################################################################
####                                                                                                                     ###
####                                                                                                                     ###
####                                         MORITZ LIMMA ANALYSIS                                                       ###
####                                                                                                                     ###
####                                                                                                                     ###
############################################################################################################################
############################################################################################################################

# Setup background
CNVs = cnv_calls %>%
    dplyr::rename(status = CNV) %>%
    mutate(status = gsub("DUP", "AMP", status)) %>%
    filter(status != 'normal') %>%
    mutate(chr = paste0('chr', chrom)) %>%
    dplyr::select(sample, chr, start, end, status)


# Setup genes
genes_path = '../ref_files/UCSC_ref_genes.dms'
gene= read.delim(genes_path, stringsAsFactors = F)
gene= dplyr::select(gene, hg19.knownGene.chrom, hg19.knownGene.txStart, hg19.knownGene.txEnd, hg19.kgXref.geneSymbol)
names(gene)=c("chr", "start_gene", "end_gene", "gene")
gene = gene %>%
        distinct() %>%
        group_by(gene) %>%
        filter((end_gene-start_gene) == max(end_gene-start_gene)) %>%
        ungroup()

### Setup gene reference
ensembl_75 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
hgnc_results = getBM(attributes=c('ensembl_gene_id','hgnc_symbol','hgnc_id'),filters = 'ensembl_gene_id', 
                      values = rownames(exprMatrix), mart = ensembl_75)

# Setup variables for analysis
clinical_matrix = clinical_matrix[clinical_matrix$sample %in% colnames(exprMatrix),]

core_clinical = clinical_matrix[c(1:9)]
hot_clinical = clinical_matrix[-c(2:9)]

# gain1q and SLAMF7
table(clinical_matrix$downsampling_chr1_160.2Mb,
      clinical_matrix$gain1q21)

window <- 10000000 # 10Mb
hotspots_with_genes_list <- list()
pdf("../data/hotspots/191125_hotspots_limma_clean.pdf", width = 10, height = 10)
for(n_hot in 1:nrow(hot_all[1:5,])){
    
    tryCatch(
    {
        # Setup
        h = hot_all$hotspot.id[n_hot]
        c = hot_all$chr[n_hot]
        st = hot_all$start.bp[n_hot] - window
        en = hot_all$end.bp[n_hot] + window

        sv_temp = SVbp[SVbp$chr == c & SVbp$pos > st & SVbp$pos < en,]
        cnv_temp = CNVs[CNVs$chr == c & CNVs$end > st & CNVs$start < en,]

        samples_exclude = unique(c(sv_temp$sample, cnv_temp$sample))

        design_sub = core_clinical %>%
                            left_join(hot_clinical %>%
                                         dplyr::select(sample, h), by = 'sample') %>%
                            filter(get(h) == 1 | !sample %in% samples_exclude)  %>% # Keep samples with involvement of the respective hotspot OR no SV or CNV within 10 Mb
                            mutate(offset = 1) %>%
                            column_to_rownames(var = 'sample') %>%
                            dplyr::select(offset, everything())

        # Limma prediction model
        lmOutput = lmFit(exprMatrix[,rownames(design_sub)], design = design_sub)
        lmOutput = eBayes(lmOutput)
        # Generating F-statistics for all but offset
        # Here we want to determine all genes which are associated with any covariate.
        # This will be based on an F-statistic. lmFit also tests wether the offset is different from zero (trivially true).
        F.stat = classifyTestsF(lmOutput[,-1], fstat.only=TRUE)
        lmOutput$F <- as.vector(F.stat)
        df1 = attr(F.stat, "df1") # degrees of freedom
        df2 = attr(F.stat, "df2")
        if(df2[1] > 1e6){ 
          lmOutput$F.p.value <- pchisq(df1*lmOutput$F,df1,lower.tail=FALSE)
        }else
          lmOutput$F.p.value <- pf(lmOutput$F,df1,df2,lower.tail=FALSE)

        ## Number of genes differentially expressed for each event
        testResults <- decideTests(lmOutput, method="hierarchical", adjust.method="BH", p.value=0.01)[,-1]
        significantGenes <- sapply(1:ncol(testResults), function(j){
          c <- lmOutput$coefficients[testResults[,j]!=0,j+1] # May be able to censor cases here -- 
          table(cut(c, breaks=c(-5,seq(-1.5,1.5,l=7),5)))
        })
        colnames(significantGenes) <- colnames(testResults)

        ## Genetic interactions                                        
        genomicData = design_sub[-1]
        interactions <- interactionsGenes <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(genomicData[,i], genomicData[,j]), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
        oddsRatio <- oddsGenes <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(genomicData[,i] + .5, genomicData[,j] +.5), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
        w <- p.adjust(lmOutput$F.p.value,"BH")<0.05
        oddsExpression <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(abs(testResults[w,i]), abs(testResults[w,j])), silent=TRUE); if(class(f)=="try-error") f=NA else f$estimate} ))
        interactionsExpression <- sapply(1:ncol(genomicData), function(i) sapply(1:ncol(genomicData), function(j) {f<- try(fisher.test(abs(testResults[w,i]), abs(testResults[w,j])), silent=TRUE); if(class(f)=="try-error") 0 else ifelse(f$estimate>1, -log10(f$p.val),log10(f$p.val))} ))
        oddsRatio[lower.tri(oddsRatio)] <- oddsExpression[lower.tri(oddsExpression)]
        interactions[lower.tri(interactions)] <- interactionsExpression[lower.tri(interactions)]

        diag(interactions) = NA
        diag(oddsRatio) = NA
        colnames(oddsRatio) <- rownames(oddsRatio) <- colnames(interactions) <- rownames(interactions) <- colnames(genomicData)
        oddsRatio[10^-abs(interactions) > 0.05] = 1
        oddsRatio[oddsRatio<1e-3] = 1e-4
        oddsRatio[oddsRatio>1e3] = 1e4
        logOdds=log10(oddsRatio)

        ############################################################
        #####       DIFFERENTIALLY EXPRESSED GENES BY HOTSPOT
        ############################################################

        i = which(colnames(testResults) == h)
        c = lmOutput$coefficients[testResults[,i]!=0,i+1]
        geneCoeffs = data.frame(hotspot.id = h,
                           ensembl_gene_id = names(c),
                           coeff = c,
                           stringsAsFactors = F)

        geneCoeffs = geneCoeffs %>%
            left_join(hgnc_results, by = "ensembl_gene_id") %>%
            dplyr::rename(gene = 'hgnc_symbol') %>%
            arrange(desc(abs(coeff))) %>%
            as.data.frame() %>%
            left_join(gene, by = 'gene') %>%
            left_join(hot_all %>%
                          dplyr::rename(chr_hot = 'chr'), by = 'hotspot.id') %>%
            mutate(chr_status = factor(ifelse(is.na(chr), "Gene not assigned", 
                                              ifelse(chr_hot != chr, "Different chromosome",
                                                     ifelse(start_gene < end.bp+5000000 & end_gene > start.bp-5000000, "<5 Mb from hotspot", "Same chromosome"))),
                                      levels = c("<5 Mb from hotspot", "Same chromosome", "Different chromosome", "Gene not assigned"))) %>%
            dplyr::select(hotspot.id, coeff, ensembl_gene_id, gene, chr, chr_status)


        if(nrow(geneCoeffs) == 0){
            next
        }

        ma <- max(geneCoeffs$coeff)+1

        if(min(geneCoeffs$coeff) < 0){
            vec <- c(min(geneCoeffs$coeff), min(geneCoeffs$coeff)/2, 0, max(geneCoeffs$coeff)/2, max(geneCoeffs$coeff))
            mi <- min(geneCoeffs$coeff)-1
        } else {
            vec <- c(0, max(geneCoeffs$coeff)/2, max(geneCoeffs$coeff))
            mi <- 0
        }

        mi <- round(mi, 2)
        vec <- round(vec, 2)

        temp_plot <- ggplot()+
            geom_point(data = filter(geneCoeffs, abs(coeff)<=0.5), aes(chr_status, coeff))+
            geom_text(data = filter(geneCoeffs, abs(coeff)>0.5), aes(chr_status, coeff, label = gene), size = 6)+
            scale_y_continuous(limits = c(mi, ma), breaks = vec)+
            labs(title = h,
                 y = 'Coefficient',
                 x = 'Distance between hotspot and gene')+
            theme(text = element_text(size = 16))

        print(temp_plot)

        hotspots_with_genes_list[[n_hot]] <- geneCoeffs  

        # Genetic interactions figure
        par(bty="n", mgp = c(2,.5,0), mar=c(20,20,4,4)+.1, las=2, tcl=-.33)
        m <- nrow(oddsRatio)
        n <- ncol(oddsRatio)
        r <- log10(oddsRatio)
        r[lower.tri(r)] <- NA
        image(x=1:n, y=1:m, r, col=brewer.pal(9,"PiYG"), breaks = c(-4:0-.Machine$double.eps,0:4), xaxt="n", yaxt="n", xlab="",ylab="", xlim=c(0, n+4), ylim=c(0, n+4))
        r <- log10(oddsRatio)
        r[upper.tri(r)] <- NA
        image(x=1:n, y=1:m, r, col=brewer.pal(9,"RdBu"), breaks = c(-4:0-.Machine$double.eps,0:4), add=TRUE)
        mtext(side=2, at=1:n, colnames(oddsRatio), font=ifelse(grepl('[[:lower:]]',colnames(oddsRatio)),1,3))
        axis(1, at=1:n,  labels = colnames(oddsRatio))

        abline(v=0:n+.5, col="white", lwd=.5)
        text(x=n/2, y=m+.5, "Genetic interactions", pos=3)
        text(x=n+1, y=m/2, "Overlap of expression targets", pos=3, srt=270)
        q <- p.adjust(10^-abs(interactions), method="BH")
        p <- p.adjust(10^-abs(interactions), method="holm")
        w = arrayInd(which(q < .1), rep(m,2))
        points(w, pch=".", col="white", cex=1.5)
        w = arrayInd(which(p < .05), rep(m,2))
        points(w, pch="*", col="white")
        image(y = 1:8 +6, x=rep(n,2)+c(2,2.5)+1, z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"PiYG"), add=TRUE)
        image(y = 1:8 +6, x=rep(n,2)+c(2.5,3)+1, z=matrix(c(1:8), nrow=1), col=brewer.pal(8,"RdBu"), add=TRUE)
        axis(side = 4, at = seq(1,7) + 6.5,  tcl=-.15, label=10^seq(-3,3), las=1, lwd=.5)
        mtext(side=4, at=10, "Odds ratio", las=3, line=3)
        par(xpd=NA)
        text(x=n+2.2, y=14, "Correlated", pos=4)
        text(x=n+2.2, y=6-.2, "Exclusive", pos=4)
        points(x=rep(n,2)+3.5, y=1:2, pch=c("*","."))
        image(x=rep(n,2)+c(2,3)+1, y=(3:4) -0.5, z=matrix(1), col=brewer.pal(3,"BrBG"), add=TRUE)
        mtext(side=4, at=1:3, c("P < 0.05", "Q < 0.1", "Not sig."), line=0.2)
        
    }, error=function(e){cat("ERROR for hotspot ", h, ":", conditionMessage(e), "\n")}
    )
}
dev.off()
                                                                             
hotspots_with_genes <- bind_rows(hotspots_with_genes_list) 

#write.csv(hotspots_with_genes, "../data/hotspots/191124_hotspots_limma_clean.csv", row.names = F)

#### UNADJUSTED ANALYSIS FOR THE DRIVER HOTSPOTS

window <- 10000000 # 10Mb

# subset driver hotspots
driver_ids <- c("nonClustered_chr4_1.8Mb", "nonClustered_chr11_69Mb", "nonClustered_chr16_78Mb", "nonClustered_chr8_126.2Mb", "nonClustered_chr8_128.5Mb")
hot_driver <- hot_all %>% filter(hotspot.id %in% driver_ids)
hotspots_drivers_with_genes_list <- list()

for(n_hot in 1:nrow(hot_driver)){
    
    tryCatch(
    {
        # Setup
        h <- hot_driver$hotspot.id[n_hot]
        c <- hot_driver$chr[n_hot]
        st <- hot_driver$start.bp[n_hot] - window
        en <- hot_driver$end.bp[n_hot] + window

        sv_temp <- SVs[SVs$chr == c & SVs$pos > st & SVs$pos < en,]
        cnv_temp <- CNVs[CNVs$chr == c & CNVs$end > st & CNVs$start < en,]

        samples_exclude <- unique(c(sv_temp$sample, cnv_temp$sample))

        design_sub <- core_clinical %>%
                            dplyr::select(sample, HRD, gain1q21) %>%
                            left_join(hot_clinical %>%
                                         dplyr::select(sample, h), by = 'sample') %>%
                            filter(get(h) == 1 | !sample %in% samples_exclude)  %>% # Keep samples with involvement of the respective hotspot OR no SV or CNV within 10 Mb
                            mutate(offset = 1) %>%
                            column_to_rownames(var = 'sample') %>%
                            dplyr::select(offset, everything())

        # Limma prediction model
        lmOutput = lmFit(exprMatrix[,rownames(design_sub)], design = design_sub)
        lmOutput = eBayes(lmOutput)
        # Generating F-statistics for all but offset
        # Here we want to determine all genes which are associated with any covariate.
        # This will be based on an F-statistic. lmFit also tests wether the offset is different from zero (trivially true).
        F.stat <- classifyTestsF(lmOutput[,-1], fstat.only=TRUE)
        lmOutput$F <- as.vector(F.stat)
        df1 <- attr(F.stat, "df1") # degrees of freedom
        df2 <- attr(F.stat, "df2")
        if(df2[1] > 1e6){ 
          lmOutput$F.p.value <- pchisq(df1*lmOutput$F,df1,lower.tail=FALSE)
        }else
          lmOutput$F.p.value = pf(lmOutput$F,df1,df2,lower.tail=FALSE)

        ## Number of genes differentially expressed for each event
        testResults <- decideTests(lmOutput, method="hierarchical", adjust.method="BH", p.value=0.01)[,-1]
        significantGenes <- sapply(1:ncol(testResults), function(j){
          c = lmOutput$coefficients[testResults[,j]!=0,j+1] # May be able to censor cases here -- 
          table(cut(c, breaks=c(-5,seq(-1.5,1.5,l=7),5)))
        })
        colnames(significantGenes) <- colnames(testResults)

        ############################################################
        #####       DIFFERENTIALLY EXPRESSED GENES BY HOTSPOT
        ############################################################

        i = which(colnames(testResults) == h)
        c = lmOutput$coefficients[testResults[,i]!=0,i+1]
        geneCoeffs = data.frame(hotspot.id = h,
                           ensembl_gene_id = names(c),
                           coeff = c,
                           stringsAsFactors = F)

        geneCoeffs = geneCoeffs %>%
            left_join(hgnc_results, by = "ensembl_gene_id") %>%
            dplyr::rename(gene = 'hgnc_symbol') %>%
            arrange(desc(abs(coeff))) %>%
            as.data.frame() %>%
            left_join(gene, by = 'gene') %>%
            left_join(hot_all %>%
                          dplyr::rename(chr_hot = 'chr'), by = 'hotspot.id') %>%
            mutate(chr_status = factor(ifelse(is.na(chr), "Gene not assigned", 
                                              ifelse(chr_hot != chr, "Different chromosome",
                                                     ifelse(start_gene < end.bp+5000000 & end_gene > start.bp-5000000, "<5 Mb from hotspot", "Same chromosome"))),
                                      levels = c("<5 Mb from hotspot", "Same chromosome", "Different chromosome", "Gene not assigned"))) %>%
            dplyr::select(hotspot.id, coeff, ensembl_gene_id, gene, chr, chr_status)


        if(nrow(geneCoeffs) == 0){
            next
        }

        ma = max(geneCoeffs$coeff)+1

        if(min(geneCoeffs$coeff) < 0){
            vec = c(min(geneCoeffs$coeff), min(geneCoeffs$coeff)/2, 0, max(geneCoeffs$coeff)/2, max(geneCoeffs$coeff))
            mi = min(geneCoeffs$coeff)-1
        } else {
            vec = c(0, max(geneCoeffs$coeff)/2, max(geneCoeffs$coeff))
            mi = 0
        }

        mi = round(mi, 2)
        vec = round(vec, 2)

        hotspots_drivers_with_genes_list[[n_hot]] <- geneCoeffs  
        
    }, error=function(e){cat("ERROR for hotspot ", h, ":", conditionMessage(e), "\n")}
    )
}
hotspots_drivers_with_genes <- bind_rows(hotspots_drivers_with_genes_list) 

hotspots_drivers_with_genes %>%
    filter(chr_status == "<5 Mb from hotspot") %>%
    filter(gene %in% c("MAF", "CCND1", "MMSET", "WHSC1", "MYC"))












