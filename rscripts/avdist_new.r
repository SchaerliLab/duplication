rm(list = ls())
library(ggplot2)
library(reshape2)

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16), strip.background = element_rect(colour="white", fill="white"), 
        legend.position = "bottom",  axis.title = element_text(size=16), 
        legend.text = element_text(size=16), legend.title = element_blank(),
        axis.text.y = element_text(size=16), axis.text.x = element_text(size=16), 
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), 
        panel.grid.minor = element_blank())


setwd("~/Documents/yolanda/")

colorpal2 = c("#CC004D","#004D80")
colorpal3 = c("black","darkorange1","navyblue")

toplot = "avmuts"
# toplot = "avdist"

d1 = read.table(sprintf("%s_all_1.txt",toplot), header = T)

d2 = read.table(sprintf("%s_all_2.txt",toplot), header = T)

X1 = melt(d1[grep("AF",d1$Sample, invert = T),], variable.name = "Copy", value.name = "dist")
X1$Round = sub("[12][ABCDEF][XYZ]","",X1$Sample)
X1$Gate = substr(X1$Sample,2,2)
X1$Replicate = substr(X1$Sample,3,3)
X1$Sample = NULL

X1 = X1[X1$Round %in% c(1:5),]
X1$Round = strtoi(X1$Round)
ss = unique(X1[,c(1,4,5)])
ss$Round = 0
ss$dist = 0
X1 = rbind(X1,ss)

X2 = melt(d2[grep("AF",d2$Sample, invert = T),], variable.name = "Copy", value.name = "dist")
X2$Round = sub("[12][ABCDEF][XYZ]","",X2$Sample)
X2$Gate = substr(X2$Sample,2,2)
X2$Replicate = substr(X2$Sample,3,3)
X2$Sample = NULL

X2 = X2[X2$Round %in% c(1:5),]
X2$Round = strtoi(X2$Round)
ss = unique(X2[,c(1,4,5)])
ss$Round = 0
ss$dist = 0
X2 = rbind(X2,ss)


gates = c("A","B","C","E","F")

yl = if(toplot == "avmuts") "Average number of mutations" else "Average pairwise distance"

plotter = function(dfX){
  mXX = do.call(rbind,lapply(c(0:5), function(r){
    do.call(rbind,lapply(unique(dfX$Gate), function(z){
      do.call(rbind, lapply(unique(dfX$Copy), function(p){
        mm = median(dfX[dfX$Round==r & dfX$Gate==z & dfX$Copy==p ,"dist"])
        data.frame(dist = mm, Gate = z, Round = r, Copy = p)
      }))
    }))
  }))
  replines <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_line(data = rdata, aes(x = Round, colour = Copy, y = dist), linetype = 'dotted', size = 0.5)
  })
  
  repoints <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_point(data = rdata, aes(x = Round, colour = Copy, y = dist), size = 0.5)
  })
  fx = ggplot(mXX, aes(x= Round, y = dist, color = Copy)) + ggopt + geom_line() + geom_line(size=1) + geom_point() +
    replines + repoints + labs(y = yl) # + facet_wrap(~Gate, ncol = 3)
  return(fx)
}



X1$Copy = sub("Copy.A","Active",X1$Copy)
X1$Copy = sub("Copy.I","Inactive",X1$Copy)

XC = rbind(X1[X1$Copy=="Active",], X2)
XC$Replicate = paste0(XC$Replicate,".",XC$Copy)
XC$Copy = sub("Active","1-copy",XC$Copy)
XC$Copy = sub("Copy.[12]","2-copy",XC$Copy)

fx1 = plotter(X1) + scale_color_manual(values = colorpal2)
fx2 = plotter(X2) + scale_color_manual(values = colorpal2)
fxc = plotter(XC[XC$Gate=="A",]) + scale_color_manual(values = rev(colorpal2))

ggsave(fx1, filename = sprintf("Figures/%s-1libs.pdf",toplot), device = "pdf", width  = 20, height = 15, units = "cm")
ggsave(fx2, filename = sprintf("Figures/%s-2libs.pdf",toplot), device = "pdf", width  = 20, height = 15, units = "cm")
ggsave(fxc, filename = sprintf("Figures/%s-both.pdf",toplot), device = "pdf", width  = 11, height = 10, units = "cm")
