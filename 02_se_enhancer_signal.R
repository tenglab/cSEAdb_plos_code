setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(Rsubread)
library(grid)
library(gridExtra)
library(ggplotify)
library(pheatmap)
source("~/Projects/super_enhancer/gene_se/scripts/se_reduce_functions.R")


e_reduce <- read.table("all_enhancer_reduce_filtered.bed",header=F,sep="\t")
e_reduce$merge_e_name <- paste(e_reduce$V1,e_reduce$V2,e_reduce$V3,sep="_")
cell_meta <- read.table("../cell_meta.txt",header=T,sep="\t")
sample_list <- read.table("../sample_list",sep="\t",header=F)

e_reduce_gr <- GRanges(
  seqnames=e_reduce$V1,
  ranges=IRanges(e_reduce$V2,e_reduce$V3,
                 names=e_reduce$merge_e_name)
)



merge_count_list <- list()
for (s in 1:nrow(sample_list)) {
  print(s)
  count_tmp <- read.table(paste0("../../fc_count/",sample_list$V1[s]),sep="\t",header=T)
  colnames(count_tmp)[7] <- sample_list$V1[s]
  count_tmp$merge_e_name <- NA

  e_gr_tmp <- GRanges(
    seqnames=count_tmp$Chr,
    ranges=IRanges(count_tmp$Start,count_tmp$End,
                   names=paste(count_tmp$Chr,count_tmp$Start,count_tmp$End,sep="_"))
  )

  # find overlap
  hit <- findOverlaps(e_gr_tmp,e_reduce_gr)

  count_tmp$merge_e_name[queryHits(hit)] <- e_reduce$merge_e_name[subjectHits(hit)]

  # aggregate count
  merge_count_list[[s]] <- aggregate(.~merge_e_name,data=count_tmp[,c(7,8)],sum)

}

# merge all merge_e_count to matrix
count_matrix_df <- e_reduce
for (a in 1:length(merge_count_list)) {
  print(a)
  count_matrix_df <- merge(count_matrix_df,merge_count_list[[a]],by="merge_e_name",all=T)
  count_matrix_df[is.na(count_matrix_df)] <- 0
}

# save matrix
count_matrix_df <- count_matrix_df[,-c(2,3,4)]
write.table(count_matrix_df,"results/all_enhancer_10p_overlap_signal_raw_120cell.txt",
            sep="\t",row.names = F,quote=F)

#-----------------------------------------------------
# DESeq2 normalize
# make count matrix
count_matrix <- count_matrix_df

duplicate_ce_idx <- which(duplicated(count_matrix_df$merge_e_name))
row.names(count_matrix) <- count_matrix$merge_e_name
count_matrix <- count_matrix[,-1]

# remove enhancers with 0 value for all sample
not_zero_in_se <- apply(count_matrix, 1, function(row) sum(row) != 0)
count_matrix <- count_matrix[not_zero_in_se,]

# make sample data for DESeq2
sample_data <- data.frame(row.names = colnames(count_matrix),
                          condition = gsub("_rep1|_rep2","",colnames(count_matrix)))

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_data,
                              design = ~ condition)

# save raw counts and normalized counts
dds <- estimateSizeFactors(dds)
normalized_count <- as.data.frame(round(counts(dds,normalized=TRUE),digits=3))
normalized_count$merge_e_name <- row.names(normalized_count)
write.table(normalized_count,"results/all_enhancer_signal_normalized_120cell.txt",
            sep="\t",row.names = F,quote=F)

# 120 sample to 60 sample
colnames(normalized_count) <- gsub("_rep1|_rep2","",colnames(normalized_count))
row.names(normalized_count) <- normalized_count$merge_e_name
normalized_count <- as.data.frame(t(rowsum(t(normalized_count[,-121]),
                                           group = colnames(normalized_count)[-121], na.rm = T)))

normalized_count$merge_e_name <- rownames(normalized_count)
write.table(normalized_count,"results/all_enhancer_signal_normalized_60cell.txt",sep="\t",row.names = F,quote=F)


#-----------------------------------------------------
# read SEs
#-----------------------------------------------------
# extract CE count
se_reduce <- read.table("se_reduce_r1.bed",sep="\t",header=F)
se_gr <- GRanges(
  seqnames=se_reduce$V1,
  ranges=IRanges(se_reduce$V2,se_reduce$V3,
                 names=se_reduce$V4)
)

ce_hit <- findOverlaps(e_reduce_gr,se_gr)

hit_df <- data.frame(merge_e_name=names(e_reduce_gr[queryHits(ce_hit)]),
                     se_name=names(se_gr[subjectHits(ce_hit)]))

normalized_count <- read.table("results/all_enhancer_signal_normalized_60cell.txt",sep="\t",header=T)

# merge signal
ce_signal <- merge(normalized_count,hit_df,by="merge_e_name")

write.table(ce_signal,"results/ce_signal_normalized_60cell.txt",sep="\t",quote = F,row.names = F)

#-----------------------------------------------------
# create binary
#-----------------------------------------------------
ce_signal <- read.table("results/ce_signal_normalized_60cell.txt",sep="\t",header=T)

# all selected ce region
ce_bed <- unique(data.frame(chr=sapply(strsplit(ce_signal$merge_e_name,"_"),"[[",1),
                            start=sapply(strsplit(ce_signal$merge_e_name,"_"),"[[",2),
                            end=sapply(strsplit(ce_signal$merge_e_name,"_"),"[[",3),
                            ce_name=ce_signal$merge_e_name,
                            se_name=ce_signal$se_name))

ce_gr <- GRanges(
  seqnames=ce_bed$chr,
  ranges=IRanges(as.integer(ce_bed$start),
                 as.integer(ce_bed$end),
                 names=ce_bed$ce_name)
)

ce_binary_matrix <- data.frame(merge_e_name=names(ce_gr))
for (c in 1:60) {
  print(c)
  cell_name <- colnames(ce_signal)[1+c]
  if (cell_name == "X786.0") {
    cell_name <- "786-0"
  } else {
    cell_name <- gsub("\\.","-",cell_name)
  }
  cell_peak_1 <- read.table(paste0("peak_filter/",cell_name,"_rep1_filtered.narrowPeak"),sep="\t",header=F)
  cell_peak_2 <- read.table(paste0("peak_filter/",cell_name,"_rep2_filtered.narrowPeak"),sep="\t",header=F)

  gr_tmp <- GRanges(
    seqnames=c(cell_peak_1$V1,
               cell_peak_2$V1),
    ranges=IRanges(c(cell_peak_1$V2,
                     cell_peak_2$V2),
                   c(cell_peak_1$V3,
                     cell_peak_2$V3))
  )

  overlap_tmp <- findOverlaps(ce_gr,gr_tmp)

  tmp_df <- data.frame(V1=rep(0,nrow(ce_signal)))
  tmp_df$V1[queryHits(overlap_tmp)] <- 1

  colnames(tmp_df)<- colnames(ce_signal)[1+c]
  ce_binary_matrix <- cbind(ce_binary_matrix,tmp_df)
}

ce_binary <- unique(merge(ce_binary_matrix,ce_signal[,c(1,62)],by="merge_e_name"))

write.table(ce_binary,"results/ce_binary_mod.txt",sep="\t",quote = F,row.names = F)


