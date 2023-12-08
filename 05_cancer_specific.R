setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(Rsubread)
library(igraph)

#-----------------------------------------------------------------
cell_spec <- read.table("results/final_se_cell_spec_r1.txt",sep="\t",header=T)

meta_data <- read.table("results/cell_meta.txt",sep="\t",header=T)
meta_data$Cell_line[which(meta_data=="786-0")] <- "X786-0"
meta_data$Cancer <- gsub("-|,| ",".",meta_data$Cancer)

cancer_freq <- data.frame(table(meta_data$Cancer))

# loop spec ce
spec_ce <- unique(cell_spec[which(cell_spec$specifity !="none"),c("se_name","spec_ce","spec_object","specifity")])

se_name <- unique(spec_ce$se_name)
se_spec_cancer_final <- data.frame()
for (c in 1:length(se_name)) {
  if (c %% 1000 ==0) {
    print(c)
  }
  se_tmp_1 <- spec_ce[which(spec_ce$se_name==se_name[c]),]

  se_spec_cancer <- data.frame()
  for (i in 1:nrow(se_tmp_1)){
    se_tmp <- se_tmp_1[i,]
    spec_cell <- data.frame(Cell_line=unlist(strsplit(se_tmp$spec_object,",")))
    spec_cell_cancer <- merge(spec_cell,meta_data,by="Cell_line")

    # more than 2 cancer is specific
    cancer_freq_tmp <- data.frame(table(spec_cell_cancer$Cancer))
    cancer_freq_tmp_2 <- merge(cancer_freq_tmp,cancer_freq,by="Var1")
    cancer_freq_tmp_3 <- cancer_freq_tmp_2[which((cancer_freq_tmp_2$Freq.x>=2 & cancer_freq_tmp_2$Freq.y>=2) |
                                                   (cancer_freq_tmp_2$Freq.x==1 & cancer_freq_tmp_2$Freq.y==1)),]

    if (nrow(cancer_freq_tmp_3)>0) {
      spec_cancer_tmp <- data.frame(se_name=se_tmp$se_name,
                                    spec_ce=se_tmp$spec_ce,
                                    number_spec_object=nrow(cancer_freq_tmp_3),
                                    spec_object=paste(cancer_freq_tmp_3$Var1,collapse=","),
                                    specifity=se_tmp$specifity,
                                    object_type="cancer")
    } else {
      spec_cancer_tmp <- data.frame(se_name=se_tmp$se_name,
                                    spec_ce="none",
                                    number_spec_object=0,
                                    spec_object="none",
                                    specifity="none",
                                    object_type="cancer")
    }

    se_spec_cancer <- rbind(se_spec_cancer,spec_cancer_tmp)

  }

  se_spec_cancer$number_spec_ce <- nrow(se_spec_cancer[which(se_spec_cancer$number_spec_object>0),])

  se_spec_cancer <- se_spec_cancer[,c(1,7,2:6)]
  se_spec_cancer_final <- rbind(se_spec_cancer_final,se_spec_cancer)
}

se_spec_cancer_final_2 <- se_spec_cancer_final[-which(se_spec_cancer_final$spec_object=="none"),]
write.table(se_spec_cancer_final_2,"results/final_se_cancer_spec_r1.txt",sep="\t",quote=F,row.names = F)







