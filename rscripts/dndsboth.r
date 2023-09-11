rm(list = ls())
library(ggplot2)
library(reshape2)
setwd('/home/bharat/Documents/yolanda/')

colorpal2 = c("#CC004D","#004D80")

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16), strip.background = element_rect(colour="white", fill="white"), 
        legend.position = "bottom",  axis.title = element_text(size=16), 
        legend.text = element_text(size=16), legend.title = element_blank(),
        axis.text.y = element_text(size=16), axis.text.x = element_text(size=16), 
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), 
        panel.grid.minor = element_blank())

nsdata = read.table("dnds_combined.tsv", stringsAsFactors = F)

colnames(nsdata) = c("Sample","Copy","dnds")

nsdata$Round = sub("[12][ABCDEF][XYZ]","",nsdata$Sample)
nsdata$Gate = substr(nsdata$Sample,2,2)
nsdata$Replicate = substr(nsdata$Sample,3,3)
nsdata$ltype = strtoi(substr(nsdata$Sample,1,1))

nsdata$Sample = NULL

gates = unique(c("A","B","C","E","F"))

xdata = nsdata[nsdata$Round %in% c(1:5) & nsdata$Gate!="D",]
xdata$Round = strtoi(xdata$Round)

plotter = function(dfX, clrv){
  mXX = do.call(rbind,lapply(c(1:5), function(r){
    do.call(rbind,lapply(gates, function(z){
      do.call(rbind, lapply(unique(dfX[,clrv]), function(p){
        mm = median(dfX[dfX$Round==r & dfX$Gate==z & dfX[,clrv]==p ,"dnds"])
        df = data.frame(dnds = mm, Gate = z, Round = r)
        df[,clrv] = p
        return(df)
      }))
    }))
  }))
  replines <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_line(data = rdata, aes_string(x = "Round", colour = clrv, y = "dnds"), linetype = 'dotted', size = 0.5)
  })
  
  repoints <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_point(data = rdata, aes_string(x = "Round", colour = clrv, y = "dnds"), size = 0.5)
  })
  fx = ggplot(mXX, aes_string(x= "Round", y = "dnds", color = clrv)) + ggopt + geom_line() + geom_line(size=1) + geom_point() +
    replines + repoints + labs(y = "dN/dS") + facet_wrap(~Gate, ncol = 3, scales = "free")
  return(fx)
}

xdata$Copy = sub("A","Active",xdata$Copy)
xdata$Copy = sub("I","Inactive",xdata$Copy)
xdata$ltype = paste0(xdata$ltype,"-copy")

xc = xdata[xdata$Copy!="Inactive",]
xc$Replicate = paste0(xc$Replicate,xc$Copy)

f1 = plotter(xdata[xdata$ltype=="1-copy",], "Copy") + scale_color_manual(values = colorpal2)
f2 = plotter(xdata[xdata$ltype=="2-copy",], "Copy") + scale_color_manual(values = colorpal2)
fB = plotter(xc, "ltype") + scale_color_manual(values = colorpal2)

ggsave(f1, filename = "Figures/dnds-1libs.pdf", device = "pdf", width  = 20, height = 15, units = "cm")
ggsave(f2, filename = "Figures/dnds-2libs.pdf", device = "pdf", width  = 20, height = 15, units = "cm")
ggsave(fB, filename = "Figures/dnds-both.pdf", device = "pdf", width  = 20, height = 15, units = "cm")

dmax =max(xdata$dnds)*1.1

wcx = do.call(rbind, lapply(gates, function(s){
  dd1 = xdata[xdata$Gate==s,]
  do.call(rbind, lapply(unique(dd1$Round), function(r){
    xdata = dd1[dd1$Round==r & dd1$ltype==1,"dnds"]
    ydata = dd1[dd1$Round==r & dd1$ltype==2,"dnds"]
    a = t.test(xdata,ydata)$p.value
    z = if(a<0.05) dmax else NA
    wxx = if(median(xdata)>median(ydata))  1 else 2
    data.frame(Gate = s, Round = strtoi(r), ltype = wxx, dnds = z)
  }))
}))


