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

cSEAdb <- readRDS("results/cSEAdb_r1.rds")
all_ce <- cSEAdb$ce_bed

# extract link and genes
source("scripts/extract_gene_link_functions.R")
# loop cell
link_clist <- c("A549","HCT-116","K-562","MCF7")
link_file_name <- c("A549","HCT116","K562","MCF7")
encode_clist <- c("A549","HCT-116","K-562","MCF7")
comp_c <- c("A549","HCT.116","K.562","MCF7")
# model peak diff
compare_df <- read.table("results/binary_model_diff_r1.txt",sep="\t",header=T)

compare_df <- compare_df[which(compare_df$cell %in% c("A549","HCT.116","K.562","MCF7")),]


# CE
link_summary <- data.frame()
for (c in 1:length(link_clist)) {
  print(c)
  link_bed <- read.table(paste0("chiapet/",link_file_name[c],".bedpe"),sep='\t', header =F)

  diff_tmp <- compare_df[which(compare_df$cell==comp_c[c]),]
  ce_tmp <- all_ce[which(all_ce$ce_name %in% unique(diff_tmp$ce_name)),]
  ce_gr <- GRanges(
    seqnames=ce_tmp$chr,
    ranges=IRanges(as.integer(ce_tmp$start),as.integer(ce_tmp$end))
  )

  # count links
  link_gr_1 <- GRanges(seqnames = link_bed$V1,IRanges(link_bed$V2,link_bed$V3))
  link_gr_2 <- GRanges(seqnames = link_bed$V4,IRanges(link_bed$V5,link_bed$V6))

  # overlap
  overlap_1 <- as.data.frame(findOverlaps(ce_gr,link_gr_1))
  overlap_2 <- as.data.frame(findOverlaps(ce_gr,link_gr_2))

  link_out_tmp_1 <- cbind(diff_tmp[overlap_1$queryHits,],link_bed[overlap_1$subjectHits,])
  link_out_tmp_1$link_group <- "left_overlap"
  link_out_tmp_2 <- cbind(diff_tmp[overlap_2$queryHits,],link_bed[overlap_2$subjectHits,])
  link_out_tmp_2$link_group <- "right_overlap"

  link_out_tmp_all <- rbind(link_out_tmp_1,link_out_tmp_2)

  link_out_tmp_all <- link_out_tmp_all[,c(1,2,6,7,14)]

  ce_name <- unique(link_out_tmp_all$ce_name)
  link_summary_tmp <- data.frame()
  for (ce in 1:length(ce_name)) {
    if (ce %% 5000 ==0){
      print(ce)
    }
    link_summary_tmp_2 <- link_out_tmp_all[which(link_out_tmp_all$ce_name==ce_name[ce]),]
    link_summary_tmp_3 <- data.frame(ce_name=ce_name[ce],
                                     cell=comp_c[c],
                                     signal=unique(link_summary_tmp_2$signal),
                                     group=unique(link_summary_tmp_2$group),
                                     mean_link_score=mean(link_summary_tmp_2$V7),
                                     n_link=nrow(link_summary_tmp_2))
    link_summary_tmp <- rbind(link_summary_tmp,link_summary_tmp_3)
  }


  link_summary <- rbind(link_summary,link_summary_tmp)

}

write.table(link_summary,"results/fig3b_ce_link_count_r1.txt",sep="\t",row.names = F,quote = F)

# number of links boxplot

link_summary <- read.table("results/fig3b_ce_link_count_r1.txt",sep="\t",header=T)

# calculate p-value
link_summary_2 <- link_summary[which(link_summary$signal!=0),]

link_summary_2$mean_link_score <- log2(link_summary_2$mean_link_score+1)
link_summary_2$n_link <- log2(link_summary_2$n_link+1)

cell_line <- unique(link_summary_2$cell)
pvalue <- data.frame()
for (c in 1:4) {
  tmp_link <- link_summary_2[which(link_summary_2$cell==cell_line[c]),]
  ttest_1 <- wilcox.test(tmp_link$n_link[which(tmp_link$group=="model only")],
                         tmp_link$n_link[which(tmp_link$group=="peak only")],
                         alternative="greater")

  ttest_2 <- wilcox.test(tmp_link$n_link[which(tmp_link$group=="model only")],
                         tmp_link$n_link[which(tmp_link$group=="common")],
                         alternative="greater")

  ttest_3 <- wilcox.test(tmp_link$n_link[which(tmp_link$group=="peak only")],
                         tmp_link$n_link[which(tmp_link$group=="common")],
                         alternative="less")
  pvalue_tmp <- data.frame(cell=cell_line[c],
                           model_peak=ttest_1$p.value,
                           model_common=ttest_2$p.value,
                           peak_common=ttest_3$p.value)
  pvalue <- rbind(pvalue,pvalue_tmp)
}
pvalue

table(link_summary_2$cell,link_summary_2$group)

link_summary_2$group <- factor(link_summary_2$group,levels = c("model only","peak only","common"))
ggplot(link_summary_2,aes(x=cell,y=n_link,fill=group))+
  #geom_boxplot(outlier.shape = NA,width=0.5,position=position_dodge(0.6))+
  geom_boxplot(outlier.shape = NA,width=0.5,position=position_dodge(0.6))+
  #geom_violin()+
  #geom_point(position=position_jitterdodge(jitter.width=0.1,dodge.width=0.6),
  #           cex=0.3)+
  theme_classic()+
  ylim(c(0,10))+
  ylab("log2(n_link+1)")+
  scale_fill_manual(values = c(brewer.pal(8,"Paired")[8],
                               brewer.pal(8,"Blues")[8],
                               "grey50"))+
  scale_color_manual(values = c(brewer.pal(8,"Paired")[8],
                                brewer.pal(8,"Blues")[8],
                                "grey50"))

#------------------------------------------------------------------
#----------------------------------------------
# create gene reference grange
#----
gtf <- rtracklayer::import("~/Projects/gencode_ref/gencode.v43.annotation.gtf")
gtf_df <- as.data.frame(gtf)

# protein coding only
# gtf_p_region <- unique(gtf_df[which(gtf_df$gene_type=="protein_coding" &
#                                       gtf_df$type=="gene"),c(1,2,3,5,7,10,11,12)])

gtf_p_region <- unique(gtf_df[,c(1,2,3,5,7,10,11,12)])


# remove chrY and chrM
gtf_p_region <- gtf_p_region[!(gtf_p_region$seqnames %in% c("chrY","chrM")),]

gtf_p_region <- gtf_p_region[which(gtf_p_region$gene_type=="protein_coding" & gtf_p_region$type=="gene"),]

# add promotor to body
gtf_p_region_w_promoter <- gtf_p_region

gtf_p_region_w_promoter$start[which(gtf_p_region_w_promoter$strand=="+")] <-
  gtf_p_region_w_promoter$start[which(gtf_p_region_w_promoter$strand=="+")]-3000

gtf_p_region_w_promoter$end[which(gtf_p_region_w_promoter$strand=="-")] <-
  gtf_p_region_w_promoter$end[which(gtf_p_region_w_promoter$strand=="-")]+3000

# gene premoter gr
gene_gr_tmp <- GRanges(
  seqnames=gtf_p_region_w_promoter$seqnames,
  ranges=IRanges(gtf_p_region_w_promoter$start,gtf_p_region_w_promoter$end,names=gtf_p_region_w_promoter$gene_name),
  strand = gtf_p_region_w_promoter$strand
)

gene_gr <- gene_gr_tmp

#gene_gr <- promoters(gene_gr_tmp, upstream=3000, use.names=TRUE)



link_gene <- data.frame()
link_gene_summary <- data.frame()

for (c in 1:length(link_clist)) {
  print(c)
  link_bed <- read.table(paste0("chiapet/",link_file_name[c],".bedpe"),sep='\t', header =F)

  diff_tmp <- compare_df[which(compare_df$cell==comp_c[c]),]

  ce_tmp <- all_ce[which(all_ce$ce_name %in% unique(diff_tmp$ce_name)),]
  ce_gr <- GRanges(
    seqnames=ce_tmp$chr,
    ranges=IRanges(as.integer(ce_tmp$start),as.integer(ce_tmp$end))
  )

  # count links
  link_gr_1 <- GRanges(seqnames = link_bed$V1,IRanges(link_bed$V2,link_bed$V3))
  link_gr_2 <- GRanges(seqnames = link_bed$V4,IRanges(link_bed$V5,link_bed$V6))

  # overlap
  overlap_1 <- as.data.frame(findOverlaps(ce_gr,link_gr_1))
  overlap_2 <- as.data.frame(findOverlaps(ce_gr,link_gr_2))

  # overlap gene
  gene_1 <- as.data.frame(findOverlaps(gene_gr,link_gr_1))
  gene_2 <- as.data.frame(findOverlaps(gene_gr,link_gr_2))

  # find common link with 1-2,and 2-1
  link_gene_1 <- merge(overlap_1,gene_2,by="subjectHits")
  link_gene_2 <- merge(overlap_2,gene_1,by="subjectHits")


  link_out_tmp_1 <- cbind(diff_tmp[link_gene_1$queryHits.x,],link_bed[link_gene_1$subjectHits,],
                          gtf_p_region[link_gene_1$queryHits.y,8,drop=F])
  link_out_tmp_1$link_group <- "left_overlap"
  link_out_tmp_2 <- cbind(diff_tmp[link_gene_2$queryHits.x,],link_bed[link_gene_2$subjectHits,],
                          gtf_p_region[link_gene_2$queryHits.y,8,drop=F])
  link_out_tmp_2$link_group <- "right_overlap"

  link_out_tmp_all <- rbind(link_out_tmp_1,link_out_tmp_2)

  link_out_tmp_all <- link_out_tmp_all[,c(1,2,6,7,15)]

  link_gene <- rbind(link_gene,link_out_tmp_all)

  ce_name <- unique(link_out_tmp_all$ce_name)
  link_gene_tmp <- data.frame()
  for (ce in 1:length(ce_name)) {
    if (ce %% 5000 ==0){
      print(ce)
    }
    link_summary_tmp_2 <- link_out_tmp_all[which(link_out_tmp_all$ce_name==ce_name[ce]),]
    link_summary_tmp_3 <- data.frame(ce_name=ce_name[ce],
                                     cell=comp_c[c],
                                     signal=unique(link_summary_tmp_2$signal),
                                     group=unique(link_summary_tmp_2$group),
                                     n_gene=length(unique(link_summary_tmp_2$gene_name)))
    link_gene_tmp <- rbind(link_gene_tmp,link_summary_tmp_3)
  }

  link_gene_summary <- rbind(link_gene_summary,link_gene_tmp)

}

write.table(link_gene,"figs_df/fig3c_ce_link_gene_w_promoter_r1.txt",sep="\t",row.names = F,quote = F)
write.table(link_gene_summary,"figs_df/fig3c_ce_link_gene_w_promoter_summary_r1.txt",sep="\t",row.names = F,quote = F)


# calculate p-value
link_summary_2 <- link_gene_summary[which(link_gene_summary$signal!=0),]

cell_line <- unique(link_summary_2$cell)
pvalue <- data.frame()
for (c in 1:4) {
  tmp_link <- link_summary_2[which(link_summary_2$cell==cell_line[c]),]
  ttest_1 <- wilcox.test(tmp_link$n_gene[which(tmp_link$group=="model only")],
                         tmp_link$n_gene[which(tmp_link$group=="peak only")],
                         alternative="greater")

  ttest_2 <- wilcox.test(tmp_link$n_gene[which(tmp_link$group=="model only")],
                         tmp_link$n_gene[which(tmp_link$group=="common")],
                         alternative="greater")

  ttest_3 <- wilcox.test(tmp_link$n_gene[which(tmp_link$group=="peak only")],
                         tmp_link$n_gene[which(tmp_link$group=="common")],
                         alternative="less")
  pvalue_tmp <- data.frame(cell=cell_line[c],
                           model_peak=ttest_1$p.value,
                           model_common=ttest_2$p.value,
                           peak_common=ttest_3$p.value)
  pvalue <- rbind(pvalue,pvalue_tmp)
}

pvalue

table(link_summary_2$cell,link_summary_2$group)

link_summary_2$group <- factor(link_summary_2$group,levels = c("model only","peak only","common"))


ggplot(link_summary_2,aes(x=cell,y=n_gene,fill=group))+
  #geom_boxplot(outlier.shape = NA,width=0.5,position=position_dodge(0.6))+
  geom_boxplot(outlier.shape = NA,width=0.5,position=position_dodge(0.6))+
  #geom_violin()+
  #geom_point(position=position_jitterdodge(jitter.width=0.1,dodge.width=0.6),
  #           cex=0.3)+
  theme_classic()+
  ylim(c(0,25))+
  ylab("n_gene")+
  scale_fill_manual(values = c(brewer.pal(8,"Paired")[8],
                               brewer.pal(8,"Blues")[8],
                               "grey50"))+
  scale_color_manual(values = c(brewer.pal(8,"Paired")[8],
                                brewer.pal(8,"Blues")[8],
                                "grey50"))

compare_df[which(compare_df$ce_name=="chr1_66443354_66447012"),]
