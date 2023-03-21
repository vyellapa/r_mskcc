# An R file to be sourced.
#
# Rearrangement plot (aka chromothripsis paper type plot) for plotting 
# rearrangement breakpoints and copy number.
#
# After sourcing workspace will have a variable 'chr_lens', and three functions:
# arc(), window_means() and plot_rearrangements()
# 
#
# Author: Yilong Li
# Adapted by gg10 23.10.2013
# Adaptation allows you to plot rearrangements, CN from bed graph and segments from Battenberg algorithm
# $8 of the bedpe file should state which colour should be used while plotting the rearrangements)
#
###############################################################################

library(quantsmooth)  # For ideogram

chr_lens = read.table("/ifs/work/leukgen/ref/homo_sapiens/GRCh37d5/gr37.fasta.fai", header=F, sep="\t", row.names=1, colClasses = c("character", "numeric"))
temp = rownames(chr_lens)
chr_lens = chr_lens[,1]
names(chr_lens) = temp


arc = function(x0, x1, y, xr, yr, col, lwd) {
  x = (x0 + x1)/2  # Center of arc
  xr = x - x0    # x-radius of arc
  
  apply(
    cbind(x, y, xr, yr, col),
    1,
    function(z) {
      x   = as.numeric(z[1])
      y   = as.numeric(z[2])
      xr  = as.numeric(z[3])
      yr  = as.numeric(z[4])
      col = z[5]
      x_points = seq(x - xr, x + xr, length.out = 200)
      y_points = y + yr * sqrt( 1  -  ( (x_points-x) / xr )^2 )
      
      lines(
        x_points,
        y_points,
        col = col,
        lwd = lwd
      )
    }
  )
  
  return()
}

window_means = function (coords, cns, min_pos, max_pos, win_size) {
  cut_levels = cut(
    coords,
    seq(min_pos, max_pos + win_size, win_size),
    labels = F,
    include_lowest = T,
    right = F
  )
  
  cut_values = by(
    cns,
    cut_levels,
    mean
  )
  
  return(
    cut_values[ as.character(1:ceiling((max_pos-min_pos+1)/win_size)) ]
  )
}

plot_rearrangements = function(
  bedpe, chrs, cn_bedgraph = NULL, segments = NULL, yrange = NULL,
  ideogram=T, cn_cex=0.5, lwd = 0.75, cn_win_size = 1e5, BFB_ids = c(),
  arrow_ln = 0.15, xlim = NULL, chr_lim = NULL, legend_pos=NULL, rearr_col='type', cn_log_scale = F
  ) {
  
  chrs = as.character(sort(chrs))
  if (!is.null(chr_lim)) {
    chr_lim = as.character(chr_lim)
  }
  chr_cum_lns = c(0, cumsum(chr_lens[chrs])[-length(chrs)])
  names(chr_cum_lns) = chrs
  xrange_size = cumsum(chr_lens[chrs])
  
  if (!is.null(cn_bedgraph)) {
    cn = cn_bedgraph[cn_bedgraph[,1] %in% chrs, ]
    if(cn_log_scale){cn$V4 = log2(cn$V4)}
    
    if (is.null(yrange)) {
      yrange = quantile(cn[,4], p = c(0.01, 0.995))
      yrange = c(floor(yrange[1]), ceiling(yrange[2]))
    }
  }
#  else {
#    stop("cn_bedgraph must be provided")
#  }
  if (!is.null(segments)) {
    if (is.null(yrange)) {
      yrange = quantile(segments[segments[,'chr'] %in% chrs, 'LogR'], p=c(0.01, 0.99))
      yrange = c(floor(yrange[1]), ceiling(yrange[2]))
    }
  }
  
  yrange_size = yrange[2] - yrange[1]
  
  td_col         = rgb(238/255, 118/255, 0) # darkorange2
  del_col        = rgb(0, 178/255, 238/255) # deepskyblue2
  inter_chrs_col = rgb(191/255, 62/255, 255/255) #darkorchid1
  tail_tail_col  = rgb(205/255, 205/255, 0/255) # yellow3
  head_head_col  = rgb(85/255, 107/255, 47/255) # darkolivegreen
  if(rearr_col=='type' & nrow(bedpe)!=0)  {
   # Colour rearrangements by rearrangement type
   bedpe[,8] = as.character(bedpe[,8])
   bedpe[bedpe[,9]=='+' & bedpe[,10]=='-',8] = del_col
   bedpe[bedpe[,9]=='+' & bedpe[,10]=='+',8] = tail_tail_col
   bedpe[bedpe[,9]=='-' & bedpe[,10]=='-',8] = head_head_col
   bedpe[bedpe[,9]=='-' & bedpe[,10]=='+',8] = td_col
   bedpe[as.character(bedpe[,1])!=as.character(bedpe[,4]),8] = inter_chrs_col
  }
  col.vec = as.vector(bedpe[,8])
  # Create the plot
#  par(mar = c(5, 4, 2, 2) + .1)
  par(mai=c(1,1,0.82,0.42))
  
  if (is.null(xlim)) {
    xlim = c(1, sum(chr_lens[chrs]))
  }
  else if (length(chrs) > 1) {
    if (is.null(chr_lim)) {
      stop()
    }
    
    xlim = chr_cum_lns[chr_lim] + xlim
  }
  
  plot(
    c(),
    ylim = c(yrange[1] - .2*yrange_size, yrange[2] + 1.4*yrange_size),
    xlim = xlim,
    bty  = "n",
    yaxt = "n",
    xaxt = "n",
    xlab = "",
    ylab = "",
    yaxs = "i",
    xaxs = "i",
  )
   if(!is.null(legend_pos)) {
      legend.data = unique(bedpe[,c('color','presence.summary')])
      legend.data = legend.data[order(table(bedpe$presence.summary), decreasing=T),]
      legend(legend_pos, legend=as.character(legend.data[,2]), col=as.character(legend.data[,1]), lty=rep(1,nrow(legend.data)), lwd=2, bg='white', text.font=2)
   }
#  legend("topright", inset=c(-0.2,0), legend=unique(bedpe[,11]), col=unique(col.vec), lty=rep(1,length(unique(col.vec))))
  # X axis names and ticks
  par(mgp = par("mgp") + c(0,1,0))
  if (all(xlim == c(1, sum(chr_lens[chrs])))) {
    axis(
      1,
      at = (cumsum(chr_lens[chrs]) + chr_cum_lns)/2,
#      labels = paste("chr", chrs, " position (Mb)", sep=""),
#      labels = paste("chr", chrs, " (Mb)", sep=""),
      labels = paste("chr", chrs, sep=""),
      tick = F,
#      cex.lab = 1.5
      cex.lab=2
    )
  }
  else if (!is.null(chr_lim)) {
    axis(
      1,
      at = mean(xlim),
#      labels = paste("chr", chr_lim, " position (Mb)", sep=""),
      labels = paste("chr", chrs, " (Mb)", sep=""),
      tick = F,
#      cex.lab = 1.5
      cex.lab=2
    )
  }
  else {
    axis(
      1,
      at = mean(xlim),
#      labels = paste("chr", chrs, " position (Mb)", sep=""),
      labels = paste("chr", chrs, " (Mb)", sep=""),
      tick = F,
#      cex.lab = 1.5
      cex.lab=2
    )
  }
  par(mgp = par("mgp") - c(0,1,0))
  
  if (length(chrs) > 1) {
    if (all(xlim == c(1, sum(chr_lens[chrs])))) {
      for (c in chrs) {
        pretty_ticks = pretty(c(1, chr_lens[c]))
        pretty_ticks = pretty_ticks[which(pretty_ticks < chr_lens[c])]
        axis(1, at = pretty_ticks + chr_cum_lns[c], labels = pretty_ticks/1e6)
      }
    }
    else {
      if (is.null(chr_lim) || !(chr_lim %in% names(chr_lens))) {
        stop()
      }
      
      pretty_ticks = pretty(xlim - chr_cum_lns[chr_lim])
      
      axis(1, at = pretty_ticks + chr_cum_lns[chr_lim], labels = pretty_ticks/1e6)
    }
  }
  else {
    if (all(xlim == c(1, sum(chr_lens[chrs])))) {
      axis(1, at = axisTicks(usr=c(1, chr_lens[chrs]), log=F), labels = axisTicks(usr=c(1, chr_lens[chrs]), log=F)/1e6)
    }
    else {
      pretty_ticks = pretty(xlim)
      axis(1, at = pretty_ticks, labels = pretty_ticks / 1e6)
    }
  }
  
  
  # Shaded grid
  for (i in yrange[1]:yrange[2]) {
    polygon(
      c(1, sum(chr_lens[chrs]), sum(chr_lens[chrs]), 1),
      i - 0.5 + c(0, 0, 1, 1),
      col=rgb(.1, .1, .1, ifelse(i %% 2 == 0, 0.1, 0.05)),
      lty=0
    )
  }
  
  
  # Line to separate chromosomes
  if (length(chrs) > 1) {
    segments(
      x0 = cumsum(chr_lens[chrs])[-length(chrs)],
      y0 = yrange[1] - 0.5,
      y1 = yrange[2] + 0.5
    )
  }
  
  
  # Plot CN
  win_size = cn_win_size
  for (c in chrs) {
    if (!is.null(cn_bedgraph)) {
      # sel = cn[,1] == c
      # When plotting only one chromome with zoomed-in image, only plot
      # the data that will be visible in the graph.
      if (length(chrs) == 1 && !all(xlim == c(1, chr_lens[chrs]))) {
        sel = cn[,1] == c & cn[,2] > xlim[1] - 1e5 & cn[,3] < xlim[2] + 1e5
        x = win_size/2 + seq(xlim[1]-1e5, xlim[2]+1e5, win_size)
        y = window_means(rowMeans(cn[sel, 2:3]), cn[sel, 4], xlim[1]-1e5, xlim[2]+1e5, win_size)
        points(
          x = x,
          y = y,
          pch = 16,
          cex = cn_cex,
          col = "black"
        )
      }
      else {
        sel = cn[,1] == c
        points(
          x = chr_cum_lns[c] + win_size/2 + seq(1, chr_lens[c], win_size),
          y = window_means(rowMeans(cn[sel, 2:3]), cn[sel, 4], 1, chr_lens[c], win_size),
          pch = 16,
          cex = cn_cex,
          col = "black"
        )
      }
    }
    
    if (!is.null(segments)) {
      sel = segments[,'chr'] == c
      
      segments(
        x0 = chr_cum_lns[c] + segments[sel, 'startpos'] + 1,
        x1 = chr_cum_lns[c] + segments[sel, 'endpos'],
        y0 = segments[sel, 'LogR'],
        lwd = 2,
        col = segments[sel, 'colour']
      )
#      text(
#        x = chr_cum_lns[c] + segments[sel, 'endpos'],
#        y = segments[sel, 'LogR'] - 0.1,
#        labels = round(segments[sel, 'ntot'], 1),
#        col = "darkgray",
#        cex = 0.8
#      )
    }
  }
  aticks = axisTicks(usr=yrange, log=F)
  if(cn_log_scale) {
    axis(2, at=aticks, labels=2^aticks, las=2)
    title(ylab="Copy number log2-transformed", cex.lab=2)
  } else {
    axis(2, at=aticks, las=2)
    title(ylab="Copy number", cex.lab=2)
  }
  
  # Plot rearrangements: First dotted lines
  segments(
    x0 = c(1, 1),
    x1 = c(sum(chr_lens[chrs]), sum(chr_lens[chrs])),
    y0 = c(yrange[2] + 0.3*yrange_size, yrange[2] + 0.75*yrange_size),
    lty = 3
  )
  abline(h = yrange[2] + 0.5)
  
  # Then 'intra-chromosomal' translocations
  sel = bedpe[,1] %in% chrs & bedpe[,4] %in% chrs
  col = col.vec[sel]
  if (sum(sel) > 0) {
    arc(
      x0  = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
      x1  = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
      y   = ifelse(bedpe[sel,9] == bedpe[sel, 10], yrange[2]+.3*yrange_size, yrange[2]+.75*yrange_size),
      yr  = ifelse(bedpe[sel, 10] == "-", 1, -1) * 0.2 * yrange_size,
      col = col,
      lwd = lwd
    )
    segments(
      x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
      y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75*yrange_size, 0.3*yrange_size),
      y1 = yrange[1],
      col = rgb(t(col2rgb(col)), alpha = 127, max=255),
#      col = col,
      lwd = lwd
    )
    segments(
      x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
      y0 = yrange[2] + ifelse(col %in% c(del_col, td_col), 0.75*yrange_size, 0.3*yrange_size),
      y1 = yrange[1],
      col = rgb(t(col2rgb(col)), alpha = 127, max=255),
#      col = col,
      lwd = lwd
    )
  }
  
 
  
  
  # Then rearrangments where low end %in% chrs and !(high end %in% chrs) 
  sel = bedpe[,1] %in% chrs & !(bedpe[,4] %in% chrs)
  col = col.vec[sel]
  if (sum(sel) > 0) {
    arrows(
      x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
      y0 = yrange[1],
      y1 = yrange[2] + yrange_size,
      col = col,
      lwd = lwd,
      length = arrow_ln
    )
    text(
      as.character(bedpe[sel, 4]),
      x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
      y = yrange[2] + 1.2 * yrange_size
    )
  }
  
  # Then rearrangments where high end %in% chrs and !(low end %in% chrs) 
  sel = !(bedpe[,1] %in% chrs) & bedpe[,4] %in% chrs
  if (sum(sel) > 0) {
    arrows(
      x0 = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
      y0 = yrange[1],
      y1 = yrange[2] + yrange_size,
       col = col.vec[sel],
      lwd = lwd,
      length = arrow_ln
    )
    text(
      as.character(bedpe[sel, 1]),
      x = chr_cum_lns[as.character(bedpe[sel, 4])] + rowMeans(bedpe[sel, 5:6]),
      y = yrange[2] + 1.2 * yrange_size
    )
  }
  
  
  # Then BFBs, currently only supporting 2 BFB events max
  # Also, in case of BFBs, only plotting 1st end
  if (length(BFB_ids) > 2) {
    stop()
  }
  if (length(BFB_ids) == 2) {
    sel = bedpe[,7] == BFB_ids[2]
    col = col.vec[sel]
    segments(
      x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
      y0 = yrange[2] + 1.2*yrange_size,
      y1 = yrange[1],
      col = col
    )
    segments(
      x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
      y0 = yrange[2] + 1.2*yrange_size,
      y1 = yrange[2] + (1.2 + 0.1)*yrange_size,
      col = col
    )
    
    # The curved arrow
    lines(
      x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(seq(-pi/2, pi/2, length.out=20)) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50,
      y = yrange[2] + (1 + 0.3)*yrange_size +  0.05*yrange_size*sin(seq(-pi/2, pi/2, length.out=20)),
      col = col,
      lwd = 2 * lwd
    )
    x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(-pi/2) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50
    arrows(
      x0 = x0,
      x1 = x0 + ifelse(bedpe[sel,9] == "+", -1, 1) * xrange_size / 50,
      y0 = yrange[2] + (1 + 0.3)*yrange_size +  0.05*yrange_size*sin(-pi/2),
      col = col,
      lwd = 2 * lwd,
      length = arrow_ln
    )
  }
  if (length(BFB_ids) >= 1) {
    abline(h = yrange[2] + c(1, 1.2)*yrange_size)
    if (length(BFB_ids) > 1) abline(h = yrange[2] + 1.4*yrange_size)
    
    sel = bedpe[,7] == BFB_ids[1]
    col = 
      ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "+", head_head_col,
              ifelse( bedpe[sel, 9] == "+" & bedpe[sel, 10] == "-", del_col,
                      ifelse( bedpe[sel, 9] == "-" & bedpe[sel, 10] == "+", td_col,  tail_tail_col)))
    segments(
      x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
      y0 = yrange[2] + 1*yrange_size,
      y1 = yrange[1],
      col = rgb(t(col2rgb(col)), alpha = 127, max=255)
#      col = col
    )
    segments(
      x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3]),
      y0 = yrange[2] + 1*yrange_size,
      y1 = yrange[2] + (1 + 0.1)*yrange_size,
      col = col
    )
    
    # The curved arrow
    lines(
      x = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(seq(-pi/2, pi/2, length.out=20)) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50,
      y = yrange[2] + (1 + 0.1)*yrange_size +  0.05*yrange_size*sin(seq(-pi/2, pi/2, length.out=20)),
      col = col,
      lwd = 2 * lwd
    )
    x0 = chr_cum_lns[as.character(bedpe[sel, 1])] + rowMeans(bedpe[sel, 2:3])  +  xrange_size/50*cos(-pi/2) * ifelse(bedpe[sel, 9] == "+", 1, -1) + ifelse(bedpe[sel, 9] == "+", -1, 1)*xrange_size/50
    arrows(
      x0 = x0,
      x1 = x0 + ifelse(bedpe[sel,9] == "+", -1, 1) * xrange_size / 50,
      y0 = yrange[2] + (1 + 0.1)*yrange_size +  0.05*yrange_size*sin(-pi/2),
      col = col,
      lwd = 2 * lwd,
      length = arrow_ln
    )
  }
  
  # Finally ideogram
  if (xlim[2] - xlim[1] < 10e6) {
    print("Ideogram plotting disabled because xlim[2] - xlim[1] < 10e6")
    
    ideogram = F
  }
  if (ideogram) {
    for (c in chrs) {
      paintCytobands(
        c,
        pos = c(1 + chr_cum_lns[c], yrange[1]-.1*yrange_size),
        units="bases",
        width = 0.1*yrange_size,
        length.out = chr_lens[c],
        legend = F
      )
    }
  }
}

#Plot graphs

#setwd("/Volumes//team78pc15/ly2/MBC/REARRANGEMENTS/PLOT_REGIONS/PDFs")
#for (x in c("PD9694a", "PD9694c","PD9694d")){
# #pdf(paste(x,"_v2.pdf", sep=''))
#  bedpe<-read.delim(paste("/Volumes/team78pc15/ly2/MBC/REARRANGEMENTS/PLOT_REGIONS/BEDpe/IGV_REVIEWED/",x,".bedpe", sep=''), header=FALSE)
#  bg<-read.delim(paste("/Volumes/team78pc15/ly2/IGV/RAW_BIN_COUNTS/",x,".ngscn_raw_bin_counts.lrrs.gc_norm.runmed.bg", sep=''), header=FALSE)
#  for(i in c(1:22,"X")){
#  plot_rearrangements(bedpe=bedpe, chrs=i, cn_bedgraph=bg)
#  }
# dev.off()
#}
