rm(list=ls())
library(ggplot2)
library(ggfortify)
library(stringr)

aas = data.frame(numcode = c(c(1:20),-10))
rownames(aas) = c("W","F","Y","I","V","L","M","C","D","E","G","A","P","H","K","R","S","T","N","Q","O")  # Ordered based on amino acid similarity https://doi.org/10.1186/1471-2105-10-394
rprot = "MGHHHHHHSIPENSGLTEEMPAQMNLEGVVNGHAFSMEGIGGGNILTGIQKLDIRVIEGDPLPFSFDILSVAFQYGNRTYTSYPAKIPDYFVQSFPEGFTFERTLSFEDGAIVKVESDISIEDGKFVGKIKYNGEGFPEDGPVMKKEVTKLEPGSESMYVSDGTLVGEVVLSYKTQSTHYTCHMKTIYRSKKPVENLPKFHYVHHRLEKKIVEEGYYYEQHETAIAKPOO"
gfprotseq = unlist(str_split(rprot,""))
gfpnumseq=sapply(gfprotseq, function(x) return(aas[x,]))

numseq <- function(haptype){
  numhap = gfpnumseq
  if(haptype!=""){
    pvmat = matrix(unlist(strsplit(haptype, "(?=[A-Za-z])(?<=[0-9])|(?=[0-9])(?<=[A-Za-z])", perl=TRUE)), nrow=2)
    for(i in seq(1,dim(pvmat)[2]))
      numhap[strtoi(pvmat[1,i])] = aas[pvmat[2,i],]
  }
  return(numhap)
}

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16, face = "bold"), strip.background = element_rect(colour="white", fill="white"), legend.position = 'none',  axis.title = element_text(size=16), axis.text.y = element_text(size=13), axis.text.x = element_text(size=13)) + 
  theme(plot.title = element_text(size=16, hjust = 0.5, face = "bold"), panel.grid.minor = element_blank())

colorpal2 = c("#CC004D","#004D80")

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16), strip.background = element_rect(colour="white", fill="white"), 
        legend.position = "bottom",  axis.title = element_text(size=16), 
        legend.text = element_text(size=16), legend.title = element_blank(),
        axis.text.y = element_text(size=16), axis.text.x = element_text(size=16), 
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), 
        panel.grid.minor = element_blank())


basepath="/home/bharat/Documents/yolanda/haptypesAll/"
setwd(basepath)

gates = c("A","B","C","E","F")

samplesize = 500

colorpal5 = c("#5BAC00","#0B8CCF","gray61","#060045","#586500","red")

for(ltype in c(1,2)){
  for(r in c(2,5)){
    hapfiles = dir(pattern = sprintf("%d[ABCEF][XYZ]%d_haptype",ltype,r))
    print(hapfiles)
    M <- do.call('cbind', lapply(hapfiles, function(V){
      if(ltype == 2){
        ht = read.table(V, fill = T, sep = "\t", stringsAsFactors = F)[,1:2]
        pvs = ht[sample(c(1:dim(ht)[1]), samplesize, replace = F),]
        mV = rbind(sapply(pvs$V1,numseq),sapply(pvs$V2,numseq))
      }
      if(ltype == 1){
        ht = read.table(V, fill = T, sep = "\t", stringsAsFactors = F)
        ht = c(ht$V1[ht$V3=="A"],ht$V2[ht$V4=="A"])
        pvs = sample(ht, samplesize, replace = F)
        mV = sapply(pvs, numseq)
      }
      return(mV)
    }))
    M2 = as.data.frame(t(M))
    fnames = rep(hapfiles, each = samplesize)
    M2$Replicate = factor(substr(fnames,3,3),levels = c("X","Y","Z"))
    M2$Gate = factor(substr(fnames,2,2),level= c(gates,"anc-coGFP"))
    sbst = sample(seq(1,dim(M2)[1]),dim(M2)[1], replace = F)
    M2$subset = sbst
    
    rdim = dim(M2)[2]-3
    WT = if(ltype == 2) data.frame(t(c(gfpnumseq,gfpnumseq))) else data.frame(t(gfpnumseq))
    colnames(WT) = colnames(M2)[1:rdim]
    WT$Replicate = "Q"
    WT$Gate = "anc-coGFP"
    WT$subset = 0
    zdata = rbind(M2,WT)
    zdata2 = zdata[,which(apply(zdata[,1:rdim],2,var)!=0)]
    zdata2$subset = zdata$subset
    zdata2$Replicate = zdata$Replicate
    zdata2$Gate = zdata$Gate
    zdata2[order(zdata2$subset),]
    lr = paste0(ltype,r)
    if(lr=="12"){
      xl = c(-0.1,0.1)
      yl = xl
    }
    if(lr=="15"){
      xl = c(-0.05,0.05)
      yl = xl
    }
    if(lr=="22"){
      xl = c(-0.1,0.1)
      yl = xl
    }
    if(lr=="25"){
      xl = c(-0.075,0.05)
      yl = c(-0.08,0.04)
    }
      
    f1 = autoplot(prcomp(zdata2[,1:(dim(zdata2)[2]-3)], scale. = T), data = zdata2, colour = "Gate", shape= "Replicate", labels=T) +
      ggopt + scale_color_manual(values = colorpal5) + ylim(yl) + xlim(xl)
    f1
    ggsave(f1, filename = sprintf("../Figures/PCA_%d-libs_%d-round.pdf",ltype,r), device = "pdf", units = "cm", height = 15, width = 15)
  }
}
