rm(list=ls())
library(flowCore)
#library(outliers)
library(ggplot2)
library(reshape2)
library(lme4)

basewdir = '/home/bharat/Documents/yolanda/FCS/AllGates/GateA_remeasured/FCS/'
colorPal = c("#004D80","#CC004D") 
ereps = c("X","Y","Z") 


scaleFUN <- function(x) sprintf("%.1f", x)

scientific_10 <- function(x) {
  parse(text=gsub("1e[+]*", "10^", scales::scientific_format()(x)))
}

gchan = "FL7-H"
bchan = "FL6-H"
gcdes = "KO525-H"
bcdes = "PB450-H"

channels = c(gchan,bchan)

figdir = "/home/bharat/Documents/yolanda/Figures/"

cv1 = function(data){
  return(sd(data)/mean(data))
}

cv2 = function(data){
  return(mad(data)/median(data))
}

rnds = c(0:5)

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16, face = "bold"), strip.background = element_rect(colour="white", fill="white"), legend.position = 'none',  axis.title = element_text(size=16), axis.text.y = element_text(size=16), axis.text.x = element_text(size=16)) + 
  theme(plot.title = element_text(size =16, face = "bold", hjust = 0.5), panel.grid.minor = element_blank())


modx <- function(x,y){
  z <- sapply(x, function(xx) a = if(xx%%y==0) y else xx%%y)
  return(z)
}

setwd(file.path(basewdir))

# remoutlier <- function(data){
#   outliers <-  NULL
#   test <- data
#   gt <- grubbs.test(test)
#   pv <- gt$p.value
#   outliers <- strsplit(gt$alternative," ")[[1]][3]
#   while(pv<0.05){
#     outliers <- c(outliers,as.numeric(strsplit(gt$alternative," ")[[1]][3]))
#     test <- data[!data %in% outliers]
#     gt <- grubbs.test(test)
#     pv <- gt$p.value
#   }
#   return(data[!data %in% outliers])
# }

files = apply(expand.grid(c(1:2),ereps,rnds,c(1:3)), 1, function(x) paste0(x[1],"A",x[2],x[3],"_P",x[4],".fcs"))
wtfiles = c(dir(pattern ="WT[1-2]*both1_P[1-3]"),dir(pattern ="2WTboth2_P[1-3]"))

ndata = 100000;
nfile = length(files)
nchaps = 3
nconds = 3
nreps = 4

plotraj = function(mf,mmf,ym,FL){
  
  replines <- lapply(c(1:3), function(x){
    rdata = mf[mf$erep==x,]
    geom_line(data = rdata, aes_string(x = "Round", colour = "ltype", y = ym), linetype = 'dotted', size = 0.5)
  })
  
  repoints <- lapply(c(1:3), function(x){
    rdata = mf[mf$erep==x,]
    geom_point(data = rdata, aes_string(x = "Round", colour = "ltype", shape = "erep", y = ym), size = 0.5)
  })
  
  plt <-  ggplot(data = mmf, aes_string(x = "Round", y = ym, color = "ltype")) + geom_line(size=1) + geom_point() + 
    replines + repoints+ ggopt +
    scale_color_manual(values = colorPal) + 
    labs(x = "Rounds", y = FL)
  return(plt)
}

dmin = 20

gdata <- do.call(rbind,lapply(c(1:2), function(l){
  do.call(rbind,lapply(ereps, function(erep){
    do.call(rbind,lapply(c(1:3), function(trep){
      do.call(rbind,lapply(c(1:5), function(r){
        fname = paste0(l,"A",erep,r,"_P",trep,".fcs")
        paste(fname)
        fobj = read.FCS(fname)
        fkey = keyword(fobj)
        gcx = sub("S","N",names(which(fkey==gcdes)))
        bcx = sub("S","N",names(which(fkey==bcdes)))
        channels = unlist(c(fkey[gcx],fkey[bcx]))
        flr = exprs(fobj)[,channels]
        flr2 = flr[flr[,channels[1]]>dmin,]
        data.frame(Green = flr2[,channels[1]], Blue = flr2[,channels[2]], ltype = l, erep = erep, trep = trep, Round = r)
      }))
    }))
  }))
}))

wtdata = do.call(rbind,lapply(wtfiles, function(fname){
  paste(fname)
  fobj = read.FCS(fname)
  fkey = keyword(fobj)
  gcx = sub("S","N",names(which(fkey==gcdes)))
  bcx = sub("S","N",names(which(fkey==bcdes)))
  channels = unlist(c(fkey[gcx],fkey[bcx]))
  flr = exprs(fobj)[,channels]
  flr2 = flr[flr[,channels[1]]>dmin,]
  l = strtoi(substr(fname,1,1))
  trep = strtoi(gsub("[^0-9]","",substr(fname,10,16)))
  ndatax = dim(flr2)[1]
  do.call(rbind,lapply(ereps, function(erep){
    a = sample(ndatax,5000)
    data.frame(Green = flr2[a,channels[1]], Blue = flr2[a,channels[2]], ltype = l, erep = erep, trep = trep, Round = 0)
  }))
}))

allGdata = rbind(gdata,wtdata)

allGdata$ltype = factor(allGdata$ltype)

mFluor = do.call(rbind,lapply(c(1:2), function(l){
  do.call(rbind,lapply(ereps, function(e){
      do.call(rbind,lapply(rnds, function(r){
        do.call(rbind, lapply(c(1:3), function(t){
          data1 =  allGdata[allGdata$Round==r & allGdata$ltype==l & allGdata$erep==e & allGdata$trep==t,]
          data.frame(Green = median(data1$Green), Blue = median(data1$Blue), ltype = l, erep = e, trep = t, Round = r)
        }))
    }))
  }))
}))
mFluor$ltype = factor(mFluor$ltype, levels = c(1:2))


mmFluor = do.call(rbind,lapply(c(1:2), function(l){
  do.call(rbind,lapply(ereps, function(e){
    do.call(rbind,lapply(rnds, function(r){
      data1 =  allGdata[allGdata$Round==r & allGdata$ltype==l & allGdata$erep==e,]
      data.frame(Green = median(data1$Green), Blue = median(data1$Blue), ltype = l, erep = e, Round = r)
    }))
  }))
}))
mmFluor$ltype = factor(mmFluor$ltype, levels = c(1:2))

mmmFluor = do.call(rbind,lapply(c(1:2), function(l){
    do.call(rbind,lapply(rnds, function(r){
    data1 =  allGdata[allGdata$Round==r & allGdata$ltype==l,]
    data.frame(Green = median(data1$Green), Blue = median(data1$Blue), ltype = l, Round = r)
  }))
}))
mmmFluor$ltype = factor(mmmFluor$ltype, levels = c(1:2))

# Normalization #

for(l in c(1:2)){
  nfac = mmmFluor[mmmFluor$ltype==l & mmmFluor$Round==0,"Green"]
  mFluor[mFluor$ltype==l,"Green"] = mFluor[mFluor$ltype==l,"Green"]/nfac
  mmFluor[mmFluor$ltype==l,"Green"] = mmFluor[mmFluor$ltype==l,"Green"]/nfac
  mmmFluor[mmmFluor$ltype==l,"Green"] = mmmFluor[mmmFluor$ltype==l,"Green"]/nfac
}


replines <- lapply(ereps, function(x){
  rdata = mmFluor[mmFluor$erep==x,]
  geom_line(data = rdata, aes(x = Round, colour = ltype, y = log10(Green)), linetype = 'dotted', size = 1)
})

repoints <- lapply(ereps, function(x){
  rdata = mmFluor[mmFluor$erep==x,]
  geom_point(data = rdata, aes(x = Round, colour = ltype, y = log10(Green), shape = erep), size = 2)
})

medplt = ggplot(data = mmmFluor, aes(x = Round, y = log10(Green), color = ltype)) + geom_line(size = 0.75) + replines + ggopt + 
  scale_color_manual(values = colorPal) + repoints + 
  labs(y = "log10(Green fluorescence)")

ggsave(medplt, filename = file.path(figdir,"normedFluorTraj.pdf"), device = 'pdf', units = 'cm', width = 12, height = 10)

lmx = lmer(log10(Green)~ Round*ltype + (1|erep), data = mFluor)
nmx = lmer(log10(Green)~ Round + (1|erep), data = mFluor)
anova(nmx,lmx)

wcx = do.call(rbind,lapply(rnds, function(d){
  wdata = mFluor
  clr = "Green"
  xdata = wdata[wdata$Round == d & wdata$ltype == 1, clr]
  ydata = wdata[wdata$Round == d & wdata$ltype == 2, clr]
  p = wilcox.test(xdata,ydata,alternative = "g")$p.value
  data.frame(Round = d, p.val = p)
}))

# Variation #


vFluor = do.call(rbind,lapply(c(1:2), function(l){
  do.call(rbind,lapply(ereps, function(e){
    do.call(rbind,lapply(rnds, function(r){
      do.call(rbind, lapply(c(1:3), function(t){
        data1 =  allGdata[allGdata$Round==r & allGdata$ltype==l & allGdata$erep==e & allGdata$trep==t,]
        data.frame(Green = cv2(data1$Green), Blue = cv2(data1$Blue), ltype = l, erep = e, trep = t, Round = r)
      }))
    }))
  }))
}))
vFluor$ltype = factor(mFluor$ltype, levels = c(1:2))


vvFluor = do.call(rbind,lapply(c(1:2), function(l){
  do.call(rbind,lapply(ereps, function(e){
    do.call(rbind,lapply(rnds, function(r){
      data1 =  allGdata[allGdata$Round==r & allGdata$ltype==l & allGdata$erep==e,]
      data.frame(Green = cv2(data1$Green), Blue = cv2(data1$Blue), ltype = l, erep = e, Round = r)
    }))
  }))
}))
vvFluor$ltype = factor(vvFluor$ltype, levels = c(1:2))

vvvFluor = do.call(rbind,lapply(c(1:2), function(l){
  do.call(rbind,lapply(rnds, function(r){
    data1 =  allGdata[allGdata$Round==r & allGdata$ltype==l,]
    data.frame(Green = cv2(data1$Green), Blue = cv2(data1$Blue), ltype = l, Round = r)
  }))
}))
vvvFluor$ltype = factor(vvvFluor$ltype, levels = c(1:2))

lmx = lmer(log10(Green)~ Round*ltype + (1|erep), data = mFluor)
nmx = lmer(log10(Green)~ Round + (1|erep), data = mFluor)
anova(nmx,lmx)

wcx = do.call(rbind,lapply(rnds, function(d){
  wdata = mFluor
  clr = "Green"
  xdata = wdata[wdata$Round == d & wdata$ltype == 1, clr]
  ydata = wdata[wdata$Round == d & wdata$ltype == 2, clr]
  if(median(xdata)<median(ydata))
    alt = "l"
  else
    alt = "g"
  p = wilcox.test(xdata,ydata,alternative = alt)$p.value
  data.frame(Round = d, p.val = p, fc = median(xdata)/median(ydata), tail = alt)
}))

wcx[wcx$p.val<0.05,]

replines <- lapply(ereps, function(x){
  rdata = vvFluor[vvFluor$erep==x,]
  geom_line(data = rdata, aes(x = Round, colour = ltype, y = Green), linetype = 'dotted', size = 1)
})

repoints <- lapply(ereps, function(x){
  rdata = vvFluor[vvFluor$erep==x,]
  geom_point(data = rdata, aes(x = Round, colour = ltype, y = Green, shape = erep), size = 2)
})

varplt = ggplot(data = vvvFluor, aes(x = Round, y = Green, color = ltype)) + geom_line(size = 0.75) + replines + ggopt + 
  scale_color_manual(values = colorPal) + repoints + 
  labs(y = "Median absolute deviation \n (Green fluorescence)")

ggsave(varplt, filename = file.path(figdir,"mad_FluorTraj.pdf"), device = 'pdf', units = 'cm', width = 12.5, height = 10)
