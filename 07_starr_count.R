setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(Rsubread)
library(grid)
library(gridExtra)
library(rtracklayer)

encode_clist <- c("A549","HCT-116","K-562","MCF7")
cSEAdb <- readRDS("results/cSEAdb_r1.rds")
ce_spec_table <- cSEAdb$se_specificity
# expand ce_spec_table
ce_gain_only <- ce_spec_table[which(ce_spec_table$object_type=="cell"),]

se_ce <- cSEAdb$ce_signal[,c(1,62)]

all_ce <- cSEAdb$ce_bed

ce_gr <- GRanges(
  seqnames=all_ce$chr,
  ranges=IRanges(as.integer(all_ce$start),as.integer(all_ce$end))
)
all_se <- cSEAdb$se_bed

#---------------------------------------
# make saf file for featurecount
#---------------------------------------
encode_clist <- c("A549","HCT-116","K-562","MCF7")
for (c in 1:length(encode_clist)) {
  # gain
  cell_spec_gain<- ce_gain_only[which(ce_gain_only$spec_object %like% encode_clist[c] &
                                        ce_gain_only$specifity=="gain"),]
  cell_spec_gain_ce <- separate(unique(cell_spec_gain[,c("spec_ce"),drop=F]),spec_ce,c("V1","V2","V3"),remove = F)

  cell_spec_gain_ce$Strand <- "."
  colnames(cell_spec_gain_ce) <- c("GeneID","Chr","Start","End","Strand")
  write.table(cell_spec_gain_ce,paste0("enhancer_reporter_assay/saf/",encode_clist[c],"_gain.saf"),quote=F,row.names = F,sep="\t")

  # loss
  cell_spec_loss<- ce_gain_only[which(ce_gain_only$spec_object %like% encode_clist[c] &
                                        ce_gain_only$specifity=="loss"),]

  cell_spec_loss_ce <- separate(unique(cell_spec_loss[,c("spec_ce"),drop=F]),spec_ce,c("V1","V2","V3"),remove = F)

  cell_spec_loss_ce$Strand <- "."
  colnames(cell_spec_loss_ce) <- c("GeneID","Chr","Start","End","Strand")
  write.table(cell_spec_loss_ce,paste0("enhancer_reporter_assay/saf/",encode_clist[c],"_loss.saf"),quote=F,row.names = F,sep="\t")


  # non
  cell_not_spec_ce <- all_ce[!(all_ce$ce_name %in% ce_gain_only$spec_ce),]

  cell_not_spec_ce$Strand <- "."
  cell_not_spec_ce <- cell_not_spec_ce[,c(4,1,2,3,5)]
  colnames(cell_not_spec_ce) <- c("GeneID","Chr","Start","End","Strand")
  write.table(cell_not_spec_ce,paste0("enhancer_reporter_assay/saf/",encode_clist[c],"_non_spec.saf"),quote=F,row.names = F,sep="\t")
}

#-------------------------------
# make bed file
for (c in 1:4) {
  saf_1 <- read.table(paste0("enhancer_reporter_assay/saf/",encode_clist[c],"_gain.saf"),sep="\t",header=T)
  saf_2 <- read.table(paste0("enhancer_reporter_assay/saf/",encode_clist[c],"_loss.saf"),sep="\t",header=T)
  saf_3 <- read.table(paste0("enhancer_reporter_assay/saf/",encode_clist[c],"_non_spec.saf"),sep="\t",header=T)
  saf_out_1 <- saf_1[,c(2,3,4,1)]
  saf_out_2 <- saf_2[,c(2,3,4,1)]
  saf_out_3 <- saf_3[,c(2,3,4,1)]

  write.table(saf_out_1,paste0("enhancer_reporter_assay/",encode_clist[c],"_gain.bed"),
              quote=F,row.names = F,sep="\t",col.names = F)
  write.table(saf_out_2,paste0("enhancer_reporter_assay/",encode_clist[c],"_loss.bed"),
              quote=F,row.names = F,sep="\t",col.names = F)
  write.table(saf_out_3,paste0("enhancer_reporter_assay/",encode_clist[c],"_non_spec.bed"),
              quote=F,row.names = F,sep="\t",col.names = F)

}

# feature count command line

#while read bed bigwig;do multiBigwigSummary BED-file --BED $bed
#-b /Users/4472271/Projects/super_enhancer/se_data_portal/new_analysis/enhancer_reporter_assay/bigwig/$bigwig
#-bs 10 --outRawCounts assay_count/"$bed"_"$bigwig"_count -o assay_count/"$bed"_"$bigwig"_count.txt; done < assay_list

# make signal boxplot
plot_df <- data.frame()
for (c in 1:4) {
  gain <- read.table(paste0("enhancer_reporter_assay/assay_count/",encode_clist[c],"_gain.bed_count"),sep="\t",header=F)
  loss <- read.table(paste0("enhancer_reporter_assay/assay_count/",encode_clist[c],"_loss.bed_count"),sep="\t",header=F)
  non_spec <- read.table(paste0("enhancer_reporter_assay/assay_count/",encode_clist[c],"_non_spec.bed_count"),sep="\t",header=F)

  plot_df_tmp <- data.frame(ce_name=c(paste(gain$V1,gain$V2,gain$V3,sep="_"),
                                      paste(loss$V1,loss$V2,loss$V3,sep="_"),
                                      paste(non_spec$V1,non_spec$V2,non_spec$V3,sep="_")),
                            width=c(gain$V3-gain$V2+1,
                                    loss$V3-loss$V2+1,
                                    non_spec$V3-non_spec$V2+1),
                            signal=c(gain$V4,
                                     loss$V4,
                                     non_spec$V4),
                            group=c(rep("gain",nrow(gain)),
                                    rep("loss",nrow(loss)),
                                    rep("non_spec",nrow(non_spec))),
                            cell=encode_clist[c])

  #print(encode_clist[c])
  #print(summary(plot_df_tmp$width))
  #print(summary(plot_df_tmp$signal))
  plot_df <- rbind(plot_df,plot_df_tmp)
}

plot_df$signal[is.na(plot_df$signal)] <- 0
write.table(plot_df,"figs_df/fig4b_reporter_assay_r1.txt",sep="\t",quote=F,row.names = F)


#---------------
# count ce spec file with multiBigwigSummary
#---------------

# command line
# while read bed bigwig;
# do multiBigwigSummary BED-file --BED ../enhancer_reporter_assay/saf/$bed
# -b bigwig/$bigwig -bs 10 --outRawCounts bw_count/"$bed"_"$bigwig"_count
# -o bw_count/"$bed"_"$bigwig"_count.txt;
# done < count_bw_list_r1


#-----------------------
# eRNA
#-----------------------
# bigwig count
sample_list <- read.table("enhancer_reporter_assay/eRNA_count_list")
plot_df <- data.frame()
for (c in 1:length(sample_list$V1)) {
  gain_p <- read.table(paste0("enhancer_reporter_assay/eRNA_counts_r1/",sample_list$V1[c],"_gain_p.bw_count"),sep="\t",header=F)
  loss_p <- read.table(paste0("enhancer_reporter_assay/eRNA_counts_r1/",sample_list$V1[c],"_loss_p.bw_count"),sep="\t",header=F)
  non_spec_p <- read.table(paste0("enhancer_reporter_assay/eRNA_counts_r1/",sample_list$V1[c],"_non_spec_p.bw_count"),sep="\t",header=F)

  gain_m <- read.table(paste0("enhancer_reporter_assay/eRNA_counts_r1/",sample_list$V1[c],"_gain_m.bw_count"),sep="\t",header=F)
  loss_m <- read.table(paste0("enhancer_reporter_assay/eRNA_counts_r1/",sample_list$V1[c],"_loss_m.bw_count"),sep="\t",header=F)
  non_spec_m <- read.table(paste0("enhancer_reporter_assay/eRNA_counts_r1/",sample_list$V1[c],"_non_spec_m.bw_count"),sep="\t",header=F)


  gain_p$V4[is.na(gain_p$V4)] <- 0
  loss_p$V4[is.na(loss_p$V4)] <- 0
  non_spec_p$V4[is.na(non_spec_p$V4)] <- 0
  gain_m$V4[is.na(gain_m$V4)] <- 0
  loss_m$V4[is.na(loss_m$V4)] <- 0
  non_spec_m$V4[is.na(non_spec_m$V4)] <- 0

  gain_p$V5 <- gain_p$V4+abs(gain_m$V4)
  loss_p$V5 <- loss_p$V4+abs(loss_m$V4)
  non_spec_p$V5 <- non_spec_p$V4+abs(non_spec_m$V4)

  plot_df_tmp <- data.frame(ce_name=c(paste(gain_p$V1,gain_p$V2,gain_p$V3,sep="_"),
                                      paste(loss_p$V1,loss_p$V2,loss_p$V3,sep="_"),
                                      paste(non_spec_p$V1,non_spec_p$V2,non_spec_p$V3,sep="_")),
                            width=c(gain_p$V3-gain_p$V2,
                                    loss_p$V3-loss_p$V2,
                                    non_spec_p$V3-non_spec_p$V2),
                            signal=c(gain_p$V5,loss_p$V5,non_spec_p$V5),
                            group=c(rep("gain",nrow(gain_p)),
                                    rep("loss",nrow(loss_p)),
                                    rep("non_spec",nrow(non_spec_p))),
                            cell=sample_list$V1[c])

  plot_df <- rbind(plot_df,plot_df_tmp)
}

write.table(plot_df,"figs_df/fig4c_eRNA_r1.txt",sep="\t",quote=F,row.names = F)







