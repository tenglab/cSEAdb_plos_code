setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(ggplotify)
library(pheatmap)
source("~/Projects/super_enhancer/gene_se/scripts/se_reduce_functions.R")


enhancer <- read.table("all_peak_filtered.bed",sep="\t",header=F)
e_gr <- GRanges(
  seqnames=enhancer$V1,
  ranges=IRanges(as.integer(enhancer$V2),as.integer(enhancer$V3))
)

se_ce <- read.table("results/all_ce.bed",sep="\t",header = F)
# ce dataframe

se_ce$width <- as.integer(se_ce$V3)-as.integer(se_ce$V2)+1

ce_10k <- se_ce

ce_10k$group <- NA
ce_10k$group[which(ce_10k$width>=1 & ce_10k$width <1000)] <- "1-1k"
ce_10k$group[which(ce_10k$width>=1000 & ce_10k$width <2000)] <- "1k-2k"
ce_10k$group[which(ce_10k$width>=2000 & ce_10k$width <3000)] <- "2k-3k"
ce_10k$group[which(ce_10k$width>=3000 & ce_10k$width <4000)] <- "3k-4k"
ce_10k$group[which(ce_10k$width>=4000 & ce_10k$width <5000)] <- "4k-5k"
ce_10k$group[which(ce_10k$width>=5000 & ce_10k$width <6000)] <- "5k-6k"
ce_10k$group[which(ce_10k$width>=6000 & ce_10k$width <7000)] <- "6k-7k"
ce_10k$group[which(ce_10k$width>=7000 & ce_10k$width <8000)] <- "7k-8k"
ce_10k$group[which(ce_10k$width>=8000 & ce_10k$width <9000)] <- "8k-9k"
ce_10k$group[which(ce_10k$width>=9000 & ce_10k$width <100000)] <- "9k-10k"
ce_10k$group[which(ce_10k$width>=10000)] <- "10k+"



plot_df <- data.frame(table(ce_10k$group))
plot_df$Var1 <- factor(plot_df$Var1,levels=c("1-1k","1k-2k","2k-3k","3k-4k","4k-5k",
                                             "5k-6k","6k-7k","7k-8k","8k-9k","9k-10k",
                                             "10k+"))
write.table(plot_df,"results/all_peak_filtered_summary.txt",sep="\t",quote=F,row.names = F)


ggplot(plot_df, aes(x=Var1,y=Freq,fill=Var1))+
  geom_bar(stat="identity")+
  geom_text(aes(label=paste0(Freq,"\n(",round(Freq/147307*100,3),"%)")),
            vjust=0.5,size=3)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

#------------------------------------------------------------
# find overlap percent
enhancer <- read.table("all_peak_filtered.bed",sep="\t",header=F)

#reduce CE with at least 10% overlap of 120 samples
#make granges
e_gr <- GRanges(
  seqnames=enhancer$V1,
  ranges=IRanges(as.integer(enhancer$V2),as.integer(enhancer$V3))
)

mem.maxVSize(vsize = Inf)

overlap_percent <- seq(0.1,0.9,0.1)

group_count <- data.frame()
for (p in 1:length(overlap_percent)) {
  print(p)
  e_reduce_tmp <- as.data.frame(group_se_by_overlap(e_gr,overlap_percent[p]))

  e_reduce_gr_tmp <- GRanges(
    seqnames=e_reduce_tmp$seqnames,
    ranges=IRanges(as.integer(e_reduce_tmp$start),as.integer(e_reduce_tmp$end))
  )

  overlap_1 <- findOverlaps(e_reduce_gr_tmp,se_gr)
  ce_df_tmp <- data.frame(e_reduce_gr_tmp[queryHits(overlap_1)])

  # group
  ce_df_tmp$group <- NA
  ce_df_tmp$group[which(ce_df_tmp$width>=1 & ce_df_tmp$width <1000)] <- "1-1k"
  ce_df_tmp$group[which(ce_df_tmp$width>=1000 & ce_df_tmp$width <2000)] <- "1k-2k"
  ce_df_tmp$group[which(ce_df_tmp$width>=2000 & ce_df_tmp$width <3000)] <- "2k-3k"
  ce_df_tmp$group[which(ce_df_tmp$width>=3000 & ce_df_tmp$width <4000)] <- "3k-4k"
  ce_df_tmp$group[which(ce_df_tmp$width>=4000 & ce_df_tmp$width <5000)] <- "4k-5k"
  ce_df_tmp$group[which(ce_df_tmp$width>=5000 & ce_df_tmp$width <6000)] <- "5k-6k"
  ce_df_tmp$group[which(ce_df_tmp$width>=6000 & ce_df_tmp$width <7000)] <- "6k-7k"
  ce_df_tmp$group[which(ce_df_tmp$width>=7000 & ce_df_tmp$width <8000)] <- "7k-8k"
  ce_df_tmp$group[which(ce_df_tmp$width>=8000 & ce_df_tmp$width <9000)] <- "8k-9k"
  ce_df_tmp$group[which(ce_df_tmp$width>=9000 & ce_df_tmp$width <100000)] <- "9k-10k"
  ce_df_tmp$group[which(ce_df_tmp$width>=10000)] <- "10k+"

  group_count_tmp <- data.frame(table(ce_df_tmp$group))
  group_count_tmp$percentage <- overlap_percent[p]

  group_count <- rbind(group_count,group_count_tmp)
}

write.table(group_count,"results/ce_width_overlap_percentage.txt",sep="\t",quote=F,row.names = F)

#-----------------------
# percentage plot
#-----------------------
plot_df <- read.table("results/all_peak_filtered_summary.txt",sep="\t",header=T)
group_count <- read.table("results/ce_width_overlap_percentage.txt",header=T)

plot_df$percentage <- "0.0"
group_count <- rbind(group_count,plot_df)

group_count$Var1 <- factor(group_count$Var1,levels=c("1-1k","1k-2k","2k-3k","3k-4k","4k-5k",
                                                            "5k-6k","6k-7k","7k-8k","8k-9k","9k-10k",
                                                            "10k+"))

group_count$Freq <- as.integer(group_count$Freq)
group_count$percentage <- as.numeric(group_count$percentage)

ggplot(group_count,
       aes(x=as.numeric(percentage),y=Freq,color=Var1))+
  geom_line()+
  theme_classic()+
  labs(color='CE width')+
  xlab("Overlap percentage")
  #facet_wrap(~percentage)+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



ggplot(group_count[which(group_count$Var1=="10k-20k"),],
       aes(x=percentage,y=Freq,group=percentage))+
  geom_bar(stat='identity')+
  #facet_wrap(~percentage)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

total_ce <- aggregate(Freq~percentage,data=group_count,sum)

ce_10k_20k <- aggregate(Freq~percentage,data=group_count[which(group_count$Var1=="10k-20k"),],sum)


ggplot(total_ce,aes(x=percentage,y=Freq))+
  geom_bar(stat='identity')+
  facet_wrap(~percentage)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# overlap with se_reduce, extract CE
ce_df <- data.frame()
for (se in 1:nrow(se_reduce)) {
  if(se %% 500 ==0 ){
    print(se)
  }
  se_gr_tmp <- GRanges(
    seqnames=se_reduce$V1[se],
    ranges=IRanges(se_reduce$V2[se],se_reduce$V3[se],
                   names=se_reduce$V4[se]))

  hit_tmp <- findOverlaps(e_gr,se_gr_tmp)

  #ce_name <- names(e_gr[queryHits(hit_tmp)])
  ce_tmp <- data.frame(e_gr[queryHits(hit_tmp)])

  ce_df <- rbind(ce_df,ce_tmp)
}
















