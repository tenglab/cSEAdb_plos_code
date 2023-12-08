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
library(plyr)
#---------------------------
# fig2d examples
#---------------------------
cSEAdb <- readRDS("results/cSEAdb_r1.rds")
ce_bed <- cSEAdb$ce_bed
mix_mod <- cSEAdb$ce_mixtrue_model_cell
gene_promoter <- cSEAdb$gene_promotor

diff_long_signal <- read.table("results/binary_model_diff_r1.txt",sep="\t",header=T)

gene_gr <- GRanges(seqnames = gene_promoter$V1,
                   IRanges(gene_promoter$V2,gene_promoter$V3))

#------
cell_line <- c("A549","HCT.116","K.562","MCF7")
chiapet <- c("A549","HCT116","K562","MCF7")

all_cell_example <- data.frame()
for (c in 1:4) {
  print(c)
  cell_chia <- read.table(paste0("chiapet/",chiapet[c],".bedpe"),sep="\t",header=F)
  # link grange
  left_link <- GRanges(seqnames = cell_chia$V1,
                       IRanges(cell_chia$V2,cell_chia$V3))
  right_link <- GRanges(seqnames = cell_chia$V4,
                        IRanges(cell_chia$V5,cell_chia$V6))

  #----------------------------------------
  # extract all ce~gene links
  #----------------------------------------
  # cell ce gr

  all_ce <- diff_long_signal[which(diff_long_signal$cell==cell_line[c]),]

  all_ce_gr <- GRanges(seqnames=sapply(strsplit(all_ce$ce_name,"_"),"[[",1),
                       IRanges(start=as.integer(sapply(strsplit(all_ce$ce_name,"_"),"[[",2)),
                               end=as.integer(sapply(strsplit(all_ce$ce_name,"_"),"[[",3))))

  # left ce right gene
  l_c_all <- findOverlaps(all_ce_gr,left_link)
  r_g_all <- findOverlaps(gene_gr,right_link)

  # left gene right ce
  r_c_all <- findOverlaps(all_ce_gr,right_link)
  l_g_all <- findOverlaps(gene_gr,left_link)

  # intersect
  l_inter <- unique(intersect(subjectHits(l_c_all),subjectHits(r_g_all)))
  r_inter <- unique(intersect(subjectHits(r_c_all),subjectHits(l_g_all)))

  l_pair <- l_c_all[which(subjectHits(l_c_all) %in% l_inter)]
  r_pair <- r_c_all[which(subjectHits(r_c_all) %in% r_inter)]

  # all ce include model only ce out dataframe
  l_out <- cbind(all_ce[queryHits(l_pair),],cell_chia[subjectHits(l_pair),])
  r_out <- cbind(all_ce[queryHits(r_pair),],cell_chia[subjectHits(r_pair),])

  all_out_tmp <- rbind(l_out,r_out)

  write.table(all_out_tmp[,c(8:14)],paste0("figs_df/fig3_igv_",cell_line[c],"_link_r1.bedpe"),sep="\t",row.names = F,quote = F)


  #----------------------------------------
  # find peak only ce list don't have links
  #----------------------------------------
  no_link_ce <- all_ce[-unique(c(queryHits(l_c_all),queryHits(r_c_all))),]

  no_link_ce <- no_link_ce[which(no_link_ce$group=="peak only"),]

  # make final out put
  final_out_tmp <- rbind.fill(all_out_tmp,no_link_ce)
  all_cell_example <- rbind(all_cell_example,final_out_tmp)
}

# add SE names
all_cell_example <- read.table("figs_df/fig3_examples_link_r1.txt",sep="\t",header=T)

write.table(all_cell_example_2,"figs_df/fig3_examples_link_r1.txt",sep="\t",row.names = F,quote = F)

# find peak only and no link examples

se_peak_no_link <- unique(all_cell_example$se_name[is.na(all_cell_example$V1)])
peak_tmp_1 <- all_cell_example[which(all_cell_example$se_name %in% se_peak_no_link),]

cell_line <- c("A549","HCT.116","K.562","MCF7")

peak_no_link <- data.frame()
for (c in 1:4) {
  print(c)
  tmp_1 <- peak_tmp_1[which(peak_tmp_1$cell==cell_line[c]),]
  for (se in 1:length(se_peak_no_link)) {
    tmp_2 <- tmp_1[which(tmp_1$se_name==se_peak_no_link[se]),]
    if (length(unique(tmp_2$group))>1 & length(which(is.na(tmp_2$V1)))!=0) {
      peak_no_link <- rbind(peak_no_link,tmp_2)
    } else {
      next
    }
  }
}

# remove peak and model both are 0
peak_no_link_2 <- unique(peak_no_link[!(peak_no_link$binary==0&peak_no_link$mod==0),])

write.table(peak_no_link_2,"figs_df/fig3_examples_peak_no_link_r1.txt",sep="\t",row.names = F,quote = F)

# examples
# peak only
# chr1_244435122_244438603 A549
# chr13_110554793_110558664 A549
# chr14_102529921_102533138 A549
# chr5_150049468_150055442 HCT116
# chr1_41944464_41949258 MCF7

# model only
# chr14_61188638_61189975 A549
# chr1_109617105_109618081 HCT116
# chr17_74988638_74988956 HCT116
# chr14_60655908_60656970 A549



