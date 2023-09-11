rm(list=ls())

library(reshape2)
library(ggplot2)
library(stringr)

basepath="/home/bharat/Documents/yolanda/mutFreqAll"
setwd(basepath)

colorpal2 = c("#CC004D","#004D80")

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16), strip.background = element_rect(colour="white", fill="white"), 
        legend.position = "bottom",  axis.title = element_text(size=16), 
        legend.text = element_text(size=16), legend.title = element_blank(),
        axis.text.y = element_text(size=16), axis.text.x = element_text(size=16), 
        plot.title = element_text(size=16, face = "bold", hjust = 0.5), 
        panel.grid.minor = element_blank())

eallm = read.table("mutations_geq5_inAllGates.tsv", stringsAsFactors = F)
colnames(eallm) = c("Mutation", "ltype", "Gate", "Replicate","Round", "ctype", "Frequency", "Count")

gates = c("A","B","C","E","F")
reps = c("X","Y","Z")

eallm$Replicate = paste0(eallm$Replicate, eallm$ctype)
eallm$ltype = paste0(eallm$ltype,"-copy")
eallm$ctype = NULL

plotter = function(dfX, clrv){
  mXX = do.call(rbind,lapply(c(1:5), function(r){
    do.call(rbind,lapply(gates, function(z){
      do.call(rbind,lapply(unique(dfX$Mutation), function(m){
        do.call(rbind, lapply(unique(dfX[,clrv]), function(p){
          mm = median(dfX[dfX$Mutation==m & dfX$Round==r & dfX$Gate==z & dfX[,clrv]==p ,"Frequency"])
          df = data.frame(Frequency = mm, Gate = z, Round = r, Mutation = m)
          df[,clrv] = p
          return(df)
        }))
      }))
    }))
  }))
  mXX$Frequency[is.na(mXX$Frequency) ] = 0.005
  replines <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_line(data = rdata, aes_string(x = "Round", colour = clrv, y = "Frequency"), linetype = 'dotted', size = 0.5)
  })
  
  repoints <- lapply(unique(dfX$Replicate), function(x){
    rdata = dfX[dfX$Replicate==x,]
    geom_point(data = rdata, aes_string(x = "Round", colour = clrv, y = "Frequency"), size = 0.5)
  })
  fx = ggplot(mXX, aes_string(x= "Round", y = "Frequency", color = clrv)) + ggopt + geom_line() + geom_line(size=1) + geom_point() +
    replines + repoints + labs(y = "Frequency") + facet_grid(vars(Mutation), vars(Gate))
  return(fx)
}

f1 = plotter(eallm,"ltype") + scale_color_manual(values = colorpal2) + 
  scale_y_continuous(trans = 'log2', limits = c(0.005,100), breaks = c(0.06,0.25,1,4,16,64))
ggsave(f1, filename = "../Figures/allEnrichedMutations.pdf", device = "pdf", width  = 30, height = 30, units = "cm")


highfreq = unique(which(freqdata4>=10, arr.ind = T)[,1])
hfdata = freqdata4[highfreq,]
hfdata$Mutation = rownames(hfdata)
hfdata = hfdata[grep(x = hfdata$Mutation, pattern = "H[1-9][A-Z]",invert = T),]
freqsmdata$Mutation = rownames(freqsmdata)
mdata = merge(hfdata,freqsmdata, by = "Mutation", all.x = T, suffixes = c(4,1))
mdata[is.na(mdata)]=0
smdata = melt(mdata, value.name = "Frequency", variable.name = "Sample", id.vars = "Mutation")
smdata$Frequency[is.na(smdata$Frequency)]=0
smdata$Selection = substr(smdata$Sample,1,1)
smdata$Replicate = substr(smdata$Sample,2,2)
smdata$Round = factor(substr(smdata$Sample,3,3), levels = c(1,4))
# 
# mds = do.call(rbind, lapply(unique(smdata$Selection), function(s){
#   do.call(rbind,lapply(unique(smdata$Mutation), function(m){
#     dd1 = smdata[smdata$Mutation==m & smdata$Selection==s,]
#     do.call(rbind, lapply(c(1,4), function(r){
#       a = median(dd1[dd1$Round==r,"Frequency"])
#       data.frame(Selection = s, Round = factor(r), Mutation = m, Frequency = a)
#     }))
#   }))
# }))

# 
# replines <- lapply(unique(smdata$Replicate), function(x){
#   rdata = smdata[smdata$Rep==x,]
#   geom_line(data = rdata, aes(x = Round, color = Selection, y = Frequency), linetype = 'dotted', size = 0.5)
# })
# 
# repoints <- lapply(unique(smdata$Replicate), function(x){
#   rdata = smdata[smdata$Rep==x,]
#   geom_point(data = rdata, aes(x = Round, color = Selection, y = Frequency), size = 0.5)
# })

smdata$Mutation=factor(smdata$Mutation, levels = smdata$Mutation[str_order(unique(substr(smdata$Mutation,2,5)),numeric = T)])

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16, face = "bold"), strip.background = element_rect(colour="white", fill="white"), legend.position = "bottom", legend.text =  element_text(size = 16), legend.title =  element_text(size = 16), axis.title = element_text(size=16), axis.text.y = element_text(size=16), axis.text.x = element_text(size=16)) + 
  theme(plot.title = element_text(size=16, face = "bold", hjust = 0.5), panel.grid.minor = element_blank())

fig2 <- ggplot(smdata, aes(x = Selection, y = Frequency, color = Round)) + geom_point(position = position_dodge2(width=0.3)) + 
  facet_wrap(~Mutation) + ggopt +  scale_y_continuous(trans='log10')
fig2


ggsave(fig2, filename = "../Figures/Mutations.pdf", device = "pdf", width = 25, height = 25, units = "cm")

write.table(hfdata, file = "../Enriched/10+.txt", quote = F, row.names = F, sep = "\t")
