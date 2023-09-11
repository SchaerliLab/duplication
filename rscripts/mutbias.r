rm(list = ls())
library(ggplot2)
setwd('/home/bharat/Documents/yolanda/mutFreq/')
bias = read.table("mutbias.txt",stringsAsFactors = F, header = F)

colorpal = c("#56B4E9","black","#E69F00","#064E3B","#D55E00","#833292")

colnames(bias) = c("Type","Selection","Replicate","ncopy","Round","Frequency")
bias$Round = factor(bias$Round)
bias$Type = factor(bias$Type, levels = c("A:G","A:C","A:T","G:A","G:C","G:T"))
bias$Selection = factor(bias$Selection, levels = c("AF","A","B","C","E","F"))

ggopt <-  theme_bw() +
  theme(strip.text = element_text(size=16, face = "bold"), strip.background = element_rect(colour="white", fill="white"), legend.position = "bottom", legend.text =  element_text(size = 16), legend.title = element_blank(), axis.title = element_text(size=16), axis.text.y = element_text(size=16), axis.text.x = element_text(size=16)) + 
  theme(plot.title = element_text(size=16, face = "bold", hjust = 0.5), panel.grid.minor = element_blank())

p1 <- ggplot(bias[bias$Round==1 & bias$ncopy==1,], aes(x = Replicate, y = Frequency, fill = Type)) + geom_bar(stat = "identity", position = "stack") + facet_wrap(~Selection) +
  ggopt + scale_fill_manual(values = colorpal) + guides(fill = guide_legend(nrow = 1))

ggsave(p1, filename = "../Figures/mutbias_1.pdf", device = "pdf", width = 30, height = 20, units = "cm")

p2 <- ggplot(bias[bias$Round==1 & bias$ncopy==2,], aes(x = Replicate, y = Frequency, fill = Type)) + geom_bar(stat = "identity", position = "stack") + facet_wrap(~Selection) +
  ggopt + scale_fill_manual(values = colorpal) + guides(fill = guide_legend(nrow = 1))

ggsave(p2, filename = "../Figures/mutbias_2.pdf", device = "pdf", width = 30, height = 20, units = "cm")

