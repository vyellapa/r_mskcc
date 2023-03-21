setwd("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts")

muts_example_5cols <- read.delim("all.txt", stringsAsFactors = F, header=F) #%>% dplyr::select(-V3)
muts_example_5cols <- read.delim("all_subclones_oct2019.txt", stringsAsFactors = F, header=F) #%>% dplyr::select(-V3)
muts_example_5cols <- read.delim("all_subclones_oct2019.txt", stringsAsFactors = F, header=F) #%>% dplyr::select(-V3)




names(muts_example_5cols) <- c("sample", "chr", "pos", "ref", "alt")
muts_example_5cols$chr <- paste0("chr", muts_example_5cols$chr)

# Import signature reference
sig_ref <- read.delim("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts/../mmsig/data/mm_platinum_unix.txt", stringsAsFactors = F, header=T)
#sig_ref <- read.delim("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts/../mmsig/data/sbs.txt", stringsAsFactors = F, header=T)
sig_ref1 <- read.delim("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/mmsig/data/sbsLike_sig_extract_run.txt", stringsAsFactors = F, header=T)
sig_ref2 <- read.delim("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/mmsig/data/sig_exptract_run.txt", stringsAsFactors = F, header=T)
sig_ref3 <- read.delim("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/mmsig/data/sig_exptract_mm1_run.txt", stringsAsFactors = F, header=T)


                      
                      # Source functions
                      source("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts/../mmsig/R/util/main.R")
                      source("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts/../mmsig/R/util/fitting.R")
                      source("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts/../mmsig/R/util/helpers.R")
                      source("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts/../mmsig/R/util/plotting.R")
                      source("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts/../mmsig/R/util/strandbias.R")
                      source("/Users/yellapav/Desktop/p220_2019/signatures/mmsig/scripts/muts/../mmsig/R/util/bootstrap.R")
                      
                      sig_out <- mm_fit_signatures(muts.input=muts_example_5cols,
                      sig.input=sig_ref,
                      input.format = "vcf",
                      sample.sigt.profs = NULL,
                      strandbias = TRUE,
                      bootstrap = TRUE,
                      iterations = 20,
                      refcheck=TRUE,
                      dbg=FALSE)
                      
                      sig_out1 <- mm_fit_signatures(muts.input=muts_example_5cols,
                                                   sig.input=sig_ref1,
                                                   input.format = "vcf",
                                                   sample.sigt.profs = NULL,
                                                   strandbias = TRUE,
                                                   bootstrap = TRUE,
                                                   iterations = 20,
                                                   refcheck=TRUE,
                                                   dbg=FALSE)
                      
                     
                      
                      sig_out2 <- mm_fit_signatures(muts.input=muts_example_5cols,
                                                    sig.input=sig_ref2,
                                                    input.format = "vcf",
                                                    sample.sigt.profs = NULL,
                                                    strandbias = TRUE,
                                                    bootstrap = TRUE,
                                                    iterations = 20,
                                                    refcheck=TRUE,
                                                    dbg=FALSE)
                      
                      
                      sig_out3 <- mm_fit_signatures(muts.input=muts_example_5cols,
                                                    sig.input=sig_ref3,
                                                    input.format = "vcf",
                                                    sample.sigt.profs = NULL,
                                                    strandbias = TRUE,
                                                    bootstrap = TRUE,
                                                    iterations = 20,
                                                    refcheck=TRUE,
                                                    dbg=FALSE)
                      
                      #mm_col <- c(rep("lightskyblue",16), rep("black",16), rep("firebrick2",16),rep("gray88",16),rep("darkolivegreen3",16),
                     #             rep("lightpink1",16))
                     # barplot(as.numeric(sig_out$mutmatrix[,1]), col=mm_col,
                      #        names.arg = rownames(sig_out$mutmatrix), las=2, cex.names=0.7, border = F, space=rep(0,96), cex.axis=2.5)
                      
                      
                      
                      plot_signatures(sig_out$estimate)
                      plot_signatures(sig_out1$estimate)
                      plot_signatures(sig_out2$estimate)
                      
                      bootSigsPlot(filter(sig_out$bootstrap, sample %in% c("I-H-106917_2")))
                      
                      head(sig_out$strand_bias_all_3nt)
                      
                      head(sig_out$strand_bias_mm1)
                      
                      
                      library(gridExtra)
                      ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
                      library(ref_genome, character.only = TRUE)
                      library(MutationalPatterns)
                      library(viridis)
                      sample="I917_2"
                      
                      
                      
col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(8, "Paired"))   
mcol=c("#74a9cf","#969696","#de2d26","#fe9929","#2ca25f","#fbb4b9")
controls=read.table("/Users/yellapav/Desktop/p220_2019/signatures/melph/signatures_controls/extract/SBS96/All_solutions/SBS96_1_Signatures/SBS96_S1_Signatures.txt",sep="\t",header=T)
mel=read.table("/Users/yellapav/Desktop/p220_2019/signatures/melph/sig.serena/extract/SBS96/All_solutions/SBS96_1_Signatures/SBS96_S1_Signatures.txt",sep="\t",header=T)
mel=dplyr::left_join(controls,mel,by="MutationType") %>%
  dplyr::rename(control=Signature.A.x,mel=Signature.A.y) %>% dplyr::mutate(diff=ifelse((mel-control)<0,0,mel-control))
mel$change=substr(mel$MutationType,3,5)
m=melt(mel)
ggplot(m,aes(MutationType,value,fill=change))+geom_bar(stat="identity")+facet_grid(variable~change,scales = "free")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  scale_fill_manual(values = mcol)+coord_cartesian()

