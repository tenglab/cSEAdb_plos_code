setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(data.table)
library(RColorBrewer)


pair_list <- read.table("../pair_list",sep=" ",header=F)

se_binary_diff <- data.frame()
se_diff_count <- data.frame()
for (p in 1:length(pair_list$V1)) {
  if (p %% 100==0) {
    print(p)
  }

  pair_name <- paste0(pair_list$V1[p],"_",pair_list$V2[p])
  file_name <- paste0("results/DASE_pair/",pair_name,"/se_category.txt")
  se_dase <- read.table(file_name,sep="\t",header=T)

  se_dase_gr <- GRanges(
    seqnames=sapply(strsplit(se_dase$se_merge_name,"_"),"[[",1),
    ranges=IRanges(as.integer(sapply(strsplit(se_dase$se_merge_name,"_"),"[[",2)),
                   as.integer(sapply(strsplit(se_dase$se_merge_name,"_"),"[[",3)))
  )


  c1 <- pair_list$V1[p]
  c2 <- pair_list$V2[p]

  # read SE file
  se_c1_r1 <- read.table(paste0("rose_no_promotor/",c1,"_rep1_no_promotor_peak_Gateway_SuperEnhancers.bed"),sep="\t",header=F)
  se_c1_r2 <- read.table(paste0("rose_no_promotor/",c1,"_rep2_no_promotor_peak_Gateway_SuperEnhancers.bed"),sep="\t",header=F)

  se_c2_r1 <- read.table(paste0("rose_no_promotor/",c2,"_rep1_no_promotor_peak_Gateway_SuperEnhancers.bed"),sep="\t",header=F)
  se_c2_r2 <- read.table(paste0("rose_no_promotor/",c2,"_rep2_no_promotor_peak_Gateway_SuperEnhancers.bed"),sep="\t",header=F)

  # merge replicates
  se_c1_gr <- GRanges(
    seqnames=c(se_c1_r1$V1,se_c1_r2$V1),
    ranges=IRanges(c(se_c1_r1$V2,se_c1_r2$V2),
                   c(se_c1_r1$V3,se_c1_r2$V3))
  )
  se_c1_gr_reduce <- reduce(se_c1_gr)

  se_c2_gr <- GRanges(
    seqnames=c(se_c2_r1$V1,se_c2_r2$V1),
    ranges=IRanges(c(se_c2_r1$V2,se_c2_r2$V2),
                   c(se_c2_r1$V3,se_c2_r2$V3))
  )
  se_c2_gr_reduce <- reduce(se_c2_gr)

  # overlap with se_dase
  c1_overlap <- findOverlaps(se_c1_gr_reduce,se_dase_gr)
  c2_overlap <- findOverlaps(se_c2_gr_reduce,se_dase_gr)

  # exclude commom se in c1 and c2
  similar_se <- unique(intersect(unique(subjectHits(c1_overlap)),
                                     unique(subjectHits(c2_overlap))))
  binary_similar <- se_dase$se_merge_name[similar_se]
  binary_diff <- se_dase$se_merge_name[-similar_se]

  dase_similar <- unique(se_dase$se_merge_name[which(se_dase$category == "Similar")])
  dase_diff <- se_dase$se_merge_name[which(se_dase$category != "Similar")]

  bi_similar_in_dase_similar <- dase_similar[which(dase_similar %in% binary_similar)]
  bi_diff_in_dase_similar <- dase_similar[which(dase_similar %in% binary_diff)]

  bi_similar_in_dase_diff <- dase_diff[which(dase_diff %in% binary_similar)]
  bi_diff_in_dase_diff <- dase_diff[which(dase_diff %in% binary_diff)]


  #binary_diff$data_pair <- pair_name

  #se_binary_diff <- rbind(se_binary_diff,se_diff)

  se_diff_count_tmp <- data.frame(data_pair=pair_name,
                                  count=c(length(unique(bi_similar_in_dase_similar)),
                                          length(unique(bi_diff_in_dase_similar)),
                                          length(unique(bi_diff_in_dase_diff)),
                                          length(unique(bi_similar_in_dase_diff))),
                                  approch=c("Common non-differential",
                                          "Binary only",
                                          "Common differential",
                                          "DASE only"),
                                  group=c("Non-differential",
                                          "Non-differential",
                                          "Differential",
                                          "Differential"))

  se_diff_count <- rbind(se_diff_count,se_diff_count_tmp)

}

write.table(se_diff_count,"figs_df/fig1a_binary_dase_pairwise_diff.txt",sep="\t",quote = F,row.names = F)

#----------------------------------------------------------
# boxplot of diff SE count
#----------------------------------------------------------
se_binary_diff_count <- data.frame(table(se_binary_diff$data_pair))

dase_diff <- se_all_pair[which(se_all_pair$category != "Similary"),]
dase_diff_count <- data.frame(table(dase_diff$data_pair))

dase_diff_count$group <- "DASE_diff"
se_binary_diff_count$group <- "Binary_diff"

plot_df <- rbind(dase_diff_count,se_binary_diff_count)

plot_df_2 <- data.frame(group=c("DASE_diff","Binary_diff"),
                        mean=c(mean(dase_diff_count$Freq),
                               mean(se_binary_diff_count$Freq)),
                        sd=c(sd(dase_diff_count$Freq),
                             sd(se_binary_diff_count$Freq)))
ggplot(plot_df,aes(x=group,y=Freq,fill=group))+
  geom_boxplot(width=0.5)+
  #geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd))+
  theme_classic()

#-------------------------
# reduce all diff SE

binary_diff_gr <- reduce(GRanges(
  seqnames=se_binary_diff$seqnames,
  ranges=IRanges(se_binary_diff$start,
                 se_binary_diff$end))
  )

dase_diff_gr <- reduce(GRanges(
  seqnames=dase_diff$chr,
  ranges=IRanges(dase_diff$start,
                 dase_diff$end))
)

total_se_diff <- data.frame(count=c(length(binary_diff_gr),
                                    length(dase_diff_gr)),
                            group=c("Binary_diff","DASE_diff"))
ggplot(total_se_diff,aes(x=group,y=count,fill=group))+
  geom_bar(stat="identity")+
  #geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd))+
  theme_classic()


