setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(RColorBrewer)


se_ce <- read.table("results/ce_signal_normalized_60cell.txt",sep="\t",header=T)

#------------------------------------------------------------------
# calculate ce weight of max_coverage for each se
matrix_cat_final <- se_ce

#ce_signl_all <- as.vector(se_ce[,-c(61,62)])[as.vector(se_ce[,-c(61,62)])!=0]

se_name <- unique(matrix_cat_final$se_name)
ce_maxcover <- data.frame()
for (se in 1:length(se_name)) {
  if(se %% 500 ==0 ){
    print(se)
  }

  se_example <- matrix_cat_final[which(matrix_cat_final$se_name==se_name[se]),]

  # add rowname
  rownames(se_example) <- se_example$merge_e_name

  # add CE max coverage and percentage
  se_example$max <- rowMax(as.matrix(se_example[,-c(1,62)]))
  total_max_cover <- sum(se_example$max)

  ce_weight_percent <- se_example[,c(1,62,63)]
  ce_weight_percent$max_percent <- round(ce_weight_percent$max/total_max_cover*100,2)
  colnames(ce_weight_percent)[1] <- "ce_name"

  ce_maxcover <- rbind(ce_maxcover,ce_weight_percent)
}

write.table(ce_maxcover,"results/ce_maxcover_60cell.txt",quote = F,row.names = F,sep="\t")

# keep >3% CE
ce_keep <- ce_maxcover[which(ce_maxcover$max_percent > 3),]

se_ce_keep <- se_ce[which(se_ce$merge_e_name %in% ce_keep$ce_name),]

write.table(se_ce_keep,"results/ce_signal_normalized_60cell_filtered.txt",sep="\t",row.names = F,quote=F)

