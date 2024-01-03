setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/")
library(GenomicRanges)
library(tidyr)
library(data.table)

# make final database table

meta_data <- read.table("results/cell_meta.txt",sep="\t",header=T)
ce_matrix <- read.table("results/ce_signal_normalized_60cell_filtered.txt",sep="\t",header=T)
ce_mix_model_cell <- read.table("results/filtered_ce_mixture_model_r1.txt",sep="\t",header=T)


final_cancer_specific <- read.table("results/final_se_cancer_spec_r1.txt",sep="\t",header=T)
final_cell_specific <- read.table("results/final_se_cell_spec_r1.txt",sep="\t",header=T)


# cancer and cell together
all_se_specific <- rbind(final_cell_specific,final_cancer_specific)

write.table(all_se_specific,"results/final_all_se_spec_r1.txt",quote = F,row.names = F,sep="\t")

# make ce region
ce_region <- unique(ce_matrix[,1,drop=F]) %>%
  separate(merge_e_name,c("V1","V2","V3"))

ce_region$V4 <- paste(ce_region$V1,ce_region$V2,ce_region$V3,sep="_")

# rename meta cancer
meta_data$Cancer <- gsub(" |,|-",".",meta_data$Cancer)

# all se region
se_bed <- unique(data.frame(chr=sapply(strsplit(ce_matrix$se_name,"_"),"[[",1),
                        start=sapply(strsplit(ce_matrix$se_name,"_"),"[[",2),
                        end=sapply(strsplit(ce_matrix$se_name,"_"),"[[",3),
                        se_name=ce_matrix$se_name))

# all selected ce region
ce_bed <- unique(data.frame(chr=sapply(strsplit(ce_matrix$merge_e_name,"_"),"[[",1),
                            start=sapply(strsplit(ce_matrix$merge_e_name,"_"),"[[",2),
                            end=sapply(strsplit(ce_matrix$merge_e_name,"_"),"[[",3),
                            ce_name=ce_matrix$merge_e_name))

# gene file
gene_promotor <- read.table("gencode.v43_gene_promotor.bed",sep="\t",header=F)

# cSEAdb list rds
out_rds <- list(se_bed=se_bed,
                ce_bed=ce_bed,
                cell_cancer_meta=meta_data,
                ce_signal=ce_matrix,
                ce_mixtrue_model_cell=ce_mix_model_cell,
                se_specificity=all_se_specific,
                gene_promotor=gene_promotor)

saveRDS(out_rds,"results/cSEAdb_r1.rds")






