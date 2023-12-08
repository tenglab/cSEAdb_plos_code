setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(Rsubread)

source("~/Projects/super_enhancer/se_data_portal/new_analysis/scripts/se_reduce_functions.R")
mem.maxVSize(vsize = Inf)

all_enhancer <- read.table("all_peak_filtered.bed",sep="\t",header=F)
all_se <- read.table("all_se_no_pro.bed",sep="\t",header=F)

#------------------------------------------------------
# reduce SE with at least 20% overlap of 120 samples
# make granges
se_gr <- GRanges(
  seqnames=all_se$V1,
  ranges=IRanges(all_se$V2,all_se$V3,
                 names=paste(all_se$V1,all_se$V2,all_se$V3,sep="_"))
)

# test merge percentage

percentage <- seq(0,1,0.1)

se_test <- data.frame()
se_test_summary <-data.frame()
for (p in 1:length(percentage)) {
  print(p)
  se_reduce <- as.data.frame(group_se_by_overlap(se_gr,percentage[p]))
  se_reduce_bed <- se_reduce[,c(1:3)]
  se_reduce_bed$name <- paste(se_reduce_bed$seqnames,se_reduce_bed$start,se_reduce_bed$end,sep="_")
  se_reduce_bed$signal <- 1
  se_reduce_bed$strand <- "."
  se_reduce_bed$group <- as.character(percentage[p])
  se_reduce_bed$width <- se_reduce_bed$end-se_reduce_bed$start+1
  se_test <- rbind(se_test,se_reduce_bed)

  se_tmp <- data.frame(min=min(se_reduce_bed$width),
                       max=max(se_reduce_bed$width),
                       median=median(se_reduce_bed$width),
                       mean=mean(se_reduce_bed$width),
                       n=nrow(se_reduce_bed),
                       group=percentage[p])

  se_test_summary <- rbind(se_test_summary,se_tmp)
}

write.table(se_test_summary,"results/se_test_summary.txt",sep="\t",row.names = F,quote=F)

se_test_summary <- read.table("results/se_test_summary.txt",sep="\t",header=T)
# test
plot(se_test_summary$group,log2(se_test_summary$median),type="l",
     ylab="log2 scale of Median width of SE",xlab="SE overlap percentage")


# violin plot
ggplot(se_test,aes(x=group,y=log2(width),color=group))+
  geom_violin()

# set overlap to 20%
write.table(se_reduce_bed,"se_reduce_r1.bed",
            sep="\t",row.names = F,quote=F,col.names = F)

#-----------------------------------------------------------
# reduce CE with at least 10% overlap of 120 samples
# make granges
e_gr <- GRanges(
  seqnames=all_enhancer$V1,
  ranges=IRanges(as.integer(all_enhancer$V2),as.integer(all_enhancer$V3))
)


# test overlap percentage
percentage <- seq(0.6,1,0.1)
ce_test_summary <-data.frame()
for (p in 1:length(percentage)) {
  print(p)
  e_reduce_tmp <- as.data.frame(group_se_by_overlap(e_gr,0.1))
  e_bed <- e_reduce_tmp[,c(1:3)]
  e_bed$group <- as.character(percentage[p])
  e_bed$width <- e_bed$end-e_bed$start+1

  ce_tmp <- data.frame(min=min(e_bed$width),
                       max=max(e_bed$width),
                       median=median(e_bed$width),
                       mean=mean(e_bed$width),
                       n=nrow(e_bed),
                       group=as.character(0))

  ce_test_summary <- rbind(ce_test_summary,ce_tmp)

}

e_bed$name <- paste(e_bed$seqnames,e_bed$start,e_bed$end,sep="_")


ce_test_summary$group <- as.numeric(ce_test_summary$group )
write.table(ce_test_summary,"ce_test_summary.txt",
            sep="\t",row.names = F,quote=F,col.names = F)


ce_test_summary <- ce_test_summary[c(17,1:16),]
# test
par(mfrow=c(2,2))
plot(ce_test_summary$group,ce_test_summary$median,type="l",ylab="median",xlab="percentage")
plot(ce_test_summary$group,ce_test_summary$mean,type="l",ylab="mean",xlab="percentage")
plot(ce_test_summary$group,ce_test_summary$n,type="l",ylab="n",xlab="percentage")



write.table(e_bed,"all_enhancer_reduce_filtered.bed",
            sep="\t",row.names = F,quote=F,col.names = F)



