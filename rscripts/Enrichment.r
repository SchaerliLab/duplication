rm(list=ls())

library(reshape2)
library(ggplot2)
library(stringr)

basepath="/home/bharat/Documents/yolanda/mutFreqAll/"
setwd(basepath)

colorpal2 = c("#CC004D","#004D80")

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16), strip.background = element_rect(colour="white", fill="white"), 
        legend.position = "bottom",  axis.title = element_text(size=16), 
        legend.text = element_text(size=16), legend.title = element_blank(),
        axis.text.y = element_text(size=16), axis.text.x = element_text(size=16), 
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), 
        panel.grid.minor = element_blank())


glmresults <- function(dfx,mfac){
  dfx = dfx[dfx$Mutation %in% unique(dfx[dfx$Frequency>=5,"Mutation"]),]
  glf = as.formula(sprintf("cbind(Success,Failure)~%s",mfac))
  glmodels = lapply(unique(dfx$Mutation), function(Z){
    glm(formula = glf, family = binomial, data = dfx[dfx$Mutation==Z,])
  })
  
  names(glmodels) = unique(dfx$Mutation)
  
  results = do.call('rbind', lapply(names(glmodels), function(Z){
    amodel=glmodels[[Z]]
    d1 = data.frame(Mutation = Z, pval = anova(amodel,test='LRT')[2,"Pr(>Chi)"], row.names=NULL)
    d1[,names(coef(amodel))] = coef(amodel)
    return(d1)
  }))
  
  results$pval = p.adjust(results$pval,method='bonf')
  return(results[results$pval<0.01,])
}

gates = c("A","B","C","E","F")
reps = c("X","Y","Z")

cutoff = 5;
alldata = read.table(sprintf("mutations_geq%d.tsv",cutoff), header = F, stringsAsFactors = F)
colnames(alldata) = c("Mutation","ltype","Gate","Replicate","Round","ctype","Frequency","Success")
ncounts = read.table("ncounts.txt", stringsAsFactors = F, row.names = 1)
alldata$total = ncounts[apply(alldata[,2:5],1, paste0, collapse =""),]
alldata[alldata$ltype==2,"total"] = 2*alldata[alldata$ltype==2,"total"]
alldata$total = alldata$total+1
alldata$Success = alldata$Success+1
alldata$Frequency = 100*alldata$Success/alldata$total
alldata$Failure = alldata$total - alldata$Success
alldata$Gate = factor(alldata$Gate, levels = c("C","A","B","E","F"))

data5 = alldata[alldata$Round==5,]

glmbyltype = do.call(rbind,lapply(gates, function(G){
  m1 = glmresults(data5[data5$Gate==G,],"ltype")
  m1$Gate = G
  return(m1)
}))

glmbygate = do.call(rbind,lapply(unique(data5$ltype), function(L){
  m1 = glmresults(data5[data5$ltype==L & data5$Gate %in% c("E","F"),],"Gate")
  m1$ltype = L
  return(m1)
}))

data5.L = data5[paste(data5$Mutation,data5$Gate) %in% paste(glmbyltype$Mutation,glmbyltype$Gate),]

data5.L$Replicate = paste0(data5.L$Replicate, data5.L$ctype)
data5.L$ltype = paste0(data5.L$ltype,"-copy")
data5.L$ctype = NULL

for(G in gates){
  dfX = data5.L[data5.L$Gate==G,]
  ss = glmbyltype[glmbyltype$Gate==G,]
  kk = ss[order(abs(ss$beta), decreasing = T),"Mutation"]
  dfX = dfX[dfX$Mutation %in% kk[1:10],]
  mXX = do.call(rbind,lapply(c(1:5), function(r){
    do.call(rbind,lapply(unique(dfX$Mutation), function(m){
      do.call(rbind, lapply(unique(dfX$ltype), function(p){
        mm = median(dfX[dfX$Mutation==m & dfX$Round==r & dfX$ltype==p ,"Frequency"])
        data.frame(Frequency = mm, Round = r, Mutation = m, ltype = p)
      }))
    }))
  }))
  
  mXX$Frequency[is.na(mXX$Frequency) ] = 0.5
  replines <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_line(data = rdata, aes_string(x = "Round", colour = "ltype", y = "Frequency"), linetype = 'dotted', size = 0.5)
  })
  
  repoints <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_point(data = rdata, aes_string(x = "Round", colour = "ltype", y = "Frequency"), size = 0.5)
  })
  
  fx = ggplot(mXX, aes_string(x= "Round", y = "Frequency", color = "ltype")) + ggopt + geom_line() + geom_line(size=1) + geom_point() +
    replines + repoints + labs(y = "Frequency") + facet_wrap(~Mutation, scales = "free") + scale_color_manual(values = colorpal2) 
  
  ggsave(fx, filename = sprintf("../Figures/DE-%s-ltype.pdf",G), device = "pdf", width  = 30, height = 30, units = "cm")
  
}

data5.G = data5[paste(data5$Mutation,data5$ltype) %in% paste(glmbygate$Mutation,glmbygate$ltype),]

data5.G$Replicate = paste0(data5.G$Replicate, data5.L$ctype)
data5.G$ltype = paste0(data5.G$ltype,"-copy")
data5.G$ctype = NULL

for(L in unique(data5.G$ltype)){
  dfX = data5.L[data5.G$ltype==L,]
  ss = glmbygate[glmbygate$ltype==L,]
  kk = ss[order(abs(ss$beta), decreasing = T),"Mutation"]
  dfX = dfX[dfX$Mutation %in% kk[1:10],]
  mXX = do.call(rbind,lapply(c(1:5), function(r){
    do.call(rbind,lapply(unique(dfX$Mutation), function(m){
      do.call(rbind, lapply(unique(dfX$ltype), function(p){
        mm = median(dfX[dfX$Mutation==m & dfX$Round==r & dfX$ltype==p ,"Frequency"])
        data.frame(Frequency = mm, Round = r, Mutation = m, ltype = p)
      }))
    }))
  }))
  
  mXX$Frequency[is.na(mXX$Frequency) ] = 0.5
  replines <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_line(data = rdata, aes_string(x = "Round", colour = "ltype", y = "Frequency"), linetype = 'dotted', size = 0.5)
  })
  
  repoints <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_point(data = rdata, aes_string(x = "Round", colour = "ltype", y = "Frequency"), size = 0.5)
  })
  
  fx = ggplot(mXX, aes_string(x= "Round", y = "Frequency", color = "ltype")) + ggopt + geom_line() + geom_line(size=1) + geom_point() +
    replines + repoints + labs(y = "Frequency") + facet_wrap(~Mutation, scales = "free") + scale_color_manual(values = colorpal2) 
  
  ggsave(fx, filename = sprintf("../Figures/DE-%s-ltype.pdf",G), device = "pdf", width  = 30, height = 30, units = "cm")
  
}


# # is.wholenumber = function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
# 
# obgenr <- glmresults(onlyBG,ncount)
# 
# CBenr <- glmresults(BothB,ncount)
# 
# CGenr <- glmresults(BothG,ncount)
# 
# Cenr = merge(CBenr,CGenr,by = "Mutation", suffixes = c(".Blue",".Green"))
# 
# obgfreq = read.table(sprintf("freq_1i_%s_3.0+.txt",rnd), header = T, stringsAsFactors = F)[obgenr$Mutation,grep("[EF]", colnames(alldata))]
# outdata = cbind(Cenr,obgfreq)
# 
# obgfreq$Mutation = rownames(obgfreq)
# # obgfreq$Mutation=factor(obgfreq$Mutation, levels = obgfreq$Mutation[str_order(unique(substr(obgfreq$Mutation,2,5)),numeric = T)])
# dmat = dist(obgfreq, method = "euclidean")
# hcl = hclust(dmat, method = "average")
# 
# obgfreq$Mutation=factor(obgfreq$Mutation, levels = obgfreq$Mutation[hcl$order])
normfreq = do.call(rbind,lapply(obgfreq$Mutation, function(x){
  d = obgfreq[obgfreq$Mutation==x,2:10]
  d/max(d)
}))
# 
# 
# 
 normfreq$Mutation = obgfreq$Mutation
# 
freqdata = melt(normfreq, variable.name = 'Sample', value.name = 'Frequency')
freqdata$ltype = substr(freqdata$Sample,1,6)
freqdata$Replicate = substr(freqdata$Sample,8,10)
# 
# legendlabel = expression("Frequency (log"["2"] * ")")
# 
ggplot(freqdata, aes(x = Sample, y = Mutation, fill = log2(Frequency))) + geom_tile() +
  theme(legend.position = 'bottom', axis.title.x = element_blank()) + guides(fill = guide_colorbar(title = "legendlabel", direction = "horizontal", barheight = 0.5, title.hjust = 0.5, title.vjust = 1, label.position = 'bottom', title.position = "top"))
# 
# ggsave(fig1, filename = sprintf("../Figures/DE_oBG-abs-%s.pdf",rnd), device = "pdf", width = 7, height = 21, units = "cm")
# 
# 
# 
# normdata = melt(normfreq, variable.name = 'Sample', value.name = 'Frequency')
# normdata$Selection = substr(freqdata$Sample,1,1)
# normdata$Replicate = substr(freqdata$Sample,2,2)
# 
# legendlabel = "Relative frequency"
# 
# fig2 <-  ggplot(normdata, aes(x = Sample, y = Mutation, fill = Frequency)) + geom_tile() + 
#   theme(legend.position = 'bottom', axis.title.x = element_blank()) + guides(fill = guide_colorbar(title = legendlabel, direction = "horizontal", barheight = 0.5, title.hjust = 0.5, title.vjust = 1, label.position = 'bottom', title.position = "top"))
# 
# ggsave(fig2, filename = sprintf("../Figures/DE_Both-norm-%s.pdf",rnd), device = "pdf", width = 7, height = 21, units = "cm")
# 
# 
# 
# 
# write.table(outdata, file =  sprintf("../Enriched/BGonly_%s.tsv",rnd), row.names = F, quote = F, sep = "\t")
