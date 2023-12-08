setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/")
library(GenomicRanges)
library(rtracklayer)


#--------------------------------
# reference gene
#--------------------------------
gtf <- rtracklayer::import("~/Projects/gencode_ref/gencode.v43.annotation.gtf")
gtf_df <- as.data.frame(gtf)

# protein coding only
gtf_p_region <- unique(gtf_df[which(gtf_df$gene_type=="protein_coding" &
                               gtf_df$type=="gene"),c(1,2,3,5,7,10,11,12)])

# remove chrY and chrM
gtf_p_region <- gtf_p_region[!(gtf_p_region$seqnames %in% c("chrY","chrM")),]

# gene premoter gr
gene_gr_tmp <- GRanges(
  seqnames=gtf_p_region$seqnames,
  ranges=IRanges(gtf_p_region$start,gtf_p_region$end,names=gtf_p_region$gene_name),
  strand = gtf_p_region$strand
)

gene_gr <- promoters(gene_gr_tmp, upstream=3000, downstream=1000, use.names=TRUE)

gene_pro_bed <- data.frame(gene_gr)
gene_pro_bed$gene_name <- names(gene_gr)
write.table(gene_pro_bed[,c(1,2,3,6)],"cSEAdb_revision/gene_promotor.bed",sep="\t",row.names = F,col.names = F,quote = F)

#--------------------------------
# remove enhancer and SE overlap with gene premotor region
#--------------------------------
blacklist <- read.table("cSEAdb_revision/ENCFF356LFX_blacklist.bed",sep="\t",header=F)
gene_pro_bed <- read.table("cSEAdb_revision/gene_promotor.bed",sep="\t",header=F)

filter_gr <- GRanges(
  seqnames=c(blacklist$V1,gene_pro_bed$V1),
  ranges=IRanges(c(blacklist$V2,gene_pro_bed$V2),c(blacklist$V3,gene_pro_bed$V3)),
)

sample_list <- read.table("../cell_list",sep="\t",header=F)

for (c in 1:nrow(sample_list)) {
  print(c)
  #-----------------
  # enhancer
  f_name_1 <- paste0("../enhancer_bed/",sample_list$V1[c],"_rep1_peaks.narrowPeak")
  f_name_2 <- paste0("../enhancer_bed/",sample_list$V1[c],"_rep2_peaks.narrowPeak")
  old_peak_1 <- read.table(f_name_1,sep="\t",header=F)
  old_peak_2 <- read.table(f_name_2,sep="\t",header=F)

  e_gr_1 <- GRanges(
    seqnames=old_peak_1$V1,
    ranges=IRanges(old_peak_1$V2,old_peak_1$V3),
    names=old_peak_1$V4
  )

  e_gr_2 <- GRanges(
    seqnames=old_peak_2$V1,
    ranges=IRanges(old_peak_2$V2,old_peak_2$V3),
    names=old_peak_2$V4
  )

  # overlap
  overlap_1 <- findOverlaps(e_gr_1,filter_gr)
  overlap_2 <- findOverlaps(e_gr_2,filter_gr)

  # extract non-overlap
  keep_name_1 <- mcols(e_gr_1[-queryHits(overlap_1),])
  keep_name_2 <- mcols(e_gr_2[-queryHits(overlap_2),])

  new_peak_1 <- old_peak_1[which(old_peak_1$V4 %in% keep_name_1$names),]
  new_peak_2 <- old_peak_2[which(old_peak_2$V4 %in% keep_name_2$names),]

  write.table(new_peak_1,paste0("peak_no_promotor/",sample_list$V1[c],"_rep1_no_promotor.narrowPeak"),
              sep="\t",quote=F,row.names = F,col.names = F)
  write.table(new_peak_2,paste0("peak_no_promotor/",sample_list$V1[c],"_rep2_no_promotor.narrowPeak"),
              sep="\t",quote=F,row.names = F,col.names = F)

}








