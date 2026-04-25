#!/usr/bin/env Rscript

args=commandArgs(TRUE)
if (length(args) != 3) {
  print ("usage: <plink.lmiss>  <plink.imiss>  <out.prefix> ")
  q()
}

lfile <- args[1]
ifile <- args[2]
prefix <- args[3]

#lfile <- "./202403.GWAS/GWAS_genek/data/plink.lmiss"
#ifile <- "./202403.GWAS/GWAS_genek/data/plink.imiss"
lmiss <- read.table(lfile, header = T)
imiss <- read.table(ifile, header = T)


library(ggplot2)

slable <- paste("Total:",  nrow(imiss), 
                "  Mean:", round(mean(imiss$F_MISS),3),
                "  SD:", round(sd(imiss$F_MISS),3), sep = "")
psample <- ggplot(imiss, aes(x = F_MISS )) +
  geom_histogram(color = "black", fill = "orange",boundary=0 ) + 
  #xlim(c(0,1))+
  xlab( "Missing rate") +
  ylab("Sample number") +
  annotate("text",  x=Inf, y = Inf, label = slable, vjust=1.5, hjust=1.1) + 
  theme_bw() 


vlable <- paste("Total:",  nrow(lmiss), 
                "  Mean:", round(mean(lmiss$F_MISS),3),
                "  SD:", round(sd(lmiss$F_MISS),3), sep = "")
pvariant <- ggplot(lmiss, aes(x = F_MISS )) +
  geom_histogram(color = "black", fill = "orange",boundary=0 ) + 
  xlab( "Missing rate") +
  ylab(" Variant number") +
  annotate("text",  x=Inf, y = Inf, label = vlable, vjust=1.5, hjust=1.1) + 
  theme_bw() 

ggsave(file=paste(prefix,"sample_missing.pdf",sep="."), psample, width = 10, height = 7)
ggsave(file=paste(prefix,"Variant_missing.pdf",sep="."), pvariant, width = 10, height = 7 )
