setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/")
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
library(mixtools)
library(igraph)

#-----------------------------------------------------------------
mixture_keep <- read.table("results/filtered_ce_mixture_model_r1.txt",sep="\t",header=T)

# remove all 0
mixture_keep <- mixture_keep[which(mixture_keep$rowmean!=0),]

#------------------------------------------------------------------
# compare cell specific SE cluster with cancer type using VI metric
#------------------------------------------------------------------

# ground cancer type cluster
sample_meta <- read.table("../results/cell_meta.txt",sep="\t",header=T)

# ture cancer_cluster
cancer_cluster <- sample_meta[,c(1,2)]
cancer_cluster$Cell_line <- gsub("-",".",cancer_cluster$Cell_line)
cancer_cluster$Cell_line[which(cancer_cluster$Cell_line=="786.0")] <- "X786.0"

cancer_cluster$Cancer <- as.numeric(factor(cancer_cluster$Cancer))

ture_cluster <- cancer_cluster$Cancer
names(ture_cluster) <- cancer_cluster$Cell_line

#------------------------------------------------------------------
# set cell specific cutoff 0.05,0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5
#------------------------------------------------------------------

cell_spc_cut <- seq(0.05,0.5,0.01)
se_name <- unique(mixture_keep$se_name)
se_cluster_score <- data.frame()
for (cut in 1:length(cell_spc_cut)){

  print(cut)
  cutoff <- cell_spc_cut[cut]
  mixture_keep_filter <- mixture_keep[which(mixture_keep$rowmean <= cutoff |
                                              (mixture_keep$rowmean >= 1-cutoff &
                                              mixture_keep$rowmean < 1)),]

  # calculate distance
  mixture_keep_filter <- unique(mixture_keep_filter[,-c(61,62,63)])
  mixture_matrix <- as.data.frame(t(as.matrix(mixture_keep_filter)))
  dist_tmp <- dist(mixture_matrix,method="euclidean")

  #----------------------------------
  # hclust cluster
  #----------------------------------
  hc_tmp <- hclust(dist_tmp,method="complete")

  # cancer cluster compare
  cluster_tmp <- cutree(hc_tmp, k = 28)
  out_tmp <- data.frame(cut=as.character(cell_spc_cut[cut]),
                        vi_score=round(compare(cluster_tmp,ture_cluster,method="vi"),3),
                        n_ce=nrow(mixture_keep_filter),
                        ce_percent=round(nrow(mixture_keep_filter)/nrow(unique(mixture_keep[,-61]))*100,2))

  se_cluster_score <- rbind(se_cluster_score,out_tmp)

}

write.table(se_cluster_score,"figs_df/fig4a_se_cell_cutoff_vi_score_r1.txt",quote = F,row.names = F,sep="\t")


# find cutoff ggplot
p <- ggplot(se_cluster_score,aes(x=cut,y=vi_score,group=1),alpha = 0.5)+
  geom_point()+
  geom_smooth(method="loess",se=F,span = 0.8,color="red")+
  theme_classic()+
  #ylim(c(2.15,2.3))+
  theme(axis.text.x = element_text(size=10,angle=90),
        axis.title = element_text(size=16))


lowes_fit <- ggplot_build(p)$data[[2]]
lowes_fit$x <- as.numeric(se_cluster_score$cut)

# find slope==1
slop_xy <- abs((lowes_fit$y[2:46]-lowes_fit$y[1:45])/(max(lowes_fit$y)-min(lowes_fit$y)))/
  (abs(lowes_fit$x[2:46]-lowes_fit$x[1:45])/(max(lowes_fit$x[2:46])-min(lowes_fit$x[2:46])))

cutoff <- lowes_fit$x[which.min(abs(slop_xy-1))]


#------------------------------------------------------------------
# set cutoff = 0.21
#------------------------------------------------------------------
se_name <- unique(mixture_keep$se_name)

se_cell_spec_gain <- data.frame()
se_cell_spec_loss <- data.frame()

cutoff <- 0.21
# loop se
for (se in 1:length(se_name)) {
  if(se %% 500 ==0 ){
    print(se)
  }

  se_example <- mixture_keep[which(mixture_keep$se_name==se_name[se]),]

  # make cell specific ce
  lower_ce <- se_example[which(se_example$rowmean <= cutoff),]
  upper_ce <- se_example[which(se_example$rowmean >= 1-cutoff &
                               se_example$rowmean < 1),]

  lower_cut <- nrow(lower_ce)
  upper_cut <- nrow(upper_ce)

  # cell specific gain
  if (lower_cut == 0) {
    se_gain_tmp <- data.frame(se_name=se_name[se],
                              specific_gain="no",
                              number_ce_gain=0,
                              number_cell_gain=0,
                              gain_ce="none",
                              gain_cell="none")

  } else if (lower_cut > 0){

    gain_ce <- c()
    gain_cell <- c()
    gain_cell_number <- c()
    for (i in 1:lower_cut) {
      ce_tmp <- lower_ce[i,]
      # add CE name
      gain_ce_tmp <- ce_tmp$ce_name

      # specific cell
      cell_gain <- colnames(ce_tmp)[which(ce_tmp==1)]
      gain_cell_number_tmp <- length(cell_gain)
      gain_cell_tmp <- paste0(gsub("\\.","-",cell_gain),
                          collapse = ",")

      gain_ce <- c(gain_ce,gain_ce_tmp)
      gain_cell <- c(gain_cell,gain_cell_tmp)
      gain_cell_number <- c(gain_cell_number,gain_cell_number_tmp)
    }

    se_gain_tmp <- data.frame(se_name=se_name[se],
                              specific_gain="yes",
                              number_ce_gain=lower_cut,
                              number_cell_gain=gain_cell_number,
                              gain_ce=gain_ce,
                              gain_cell=gain_cell)
  }


  # cell specific loss
  if (upper_cut == 0 ) {
    se_loss_tmp <- data.frame(se_name=se_name[se],
                              specific_loss="no",
                              number_ce_loss=0,
                              number_cell_loss=0,
                              loss_ce="none",
                              loss_cell="none")

  } else if (upper_cut > 0) {

    loss_ce <- c()
    loss_cell <- c()
    loss_cell_number <- c()
    for (i in 1:upper_cut) {
      ce_tmp <- upper_ce[i,]
      # add CE name
      loss_ce_tmp <- ce_tmp$ce_name

      # specific cell
      cell_loss <- colnames(ce_tmp)[which(ce_tmp==0)]
      loss_cell_number_tmp <- length(cell_loss)
      loss_cell_tmp <- paste0(gsub("\\.","-",cell_loss),
                              collapse = ",")

      loss_ce <- c(loss_ce,loss_ce_tmp)
      loss_cell <- c(loss_cell,loss_cell_tmp)
      loss_cell_number <- c(loss_cell_number,loss_cell_number_tmp)
    }

    se_loss_tmp <- data.frame(se_name=se_name[se],
                              specific_loss="yes",
                              number_ce_loss=upper_cut,
                              number_cell_loss=loss_cell_number,
                              loss_ce=loss_ce,
                              loss_cell=loss_cell)
  }

  se_cell_spec_gain <- rbind(se_cell_spec_gain,se_gain_tmp)
  se_cell_spec_loss <- rbind(se_cell_spec_loss,se_loss_tmp)

}
write.table(se_cell_spec_gain,"results/se_cell_spec_gain_020_r1.txt",quote = F,row.names = F,sep="\t")
write.table(se_cell_spec_loss,"results/se_cell_spec_loss_020_r1.txt",quote = F,row.names = F,sep="\t")

# cell specific table
cell_gain <- read.table("results/se_cell_spec_gain_020_r1.txt",sep="\t",header=T)
cell_loss <- read.table("results/se_cell_spec_loss_020_r1.txt",sep="\t",header=T)

# reformat
cell_gain$specifity <- NA
cell_gain$specifity <- ifelse(cell_gain$specific_gain=="yes","gain","no_gain")
cell_gain <- cell_gain[,-c(2)]
colnames(cell_gain) <- c("se_name","number_spec_ce","number_spec_object","spec_ce","spec_object","specifity")

cell_loss$specifity <- NA
cell_loss$specifity <- ifelse(cell_loss$specific_loss=="yes","loss","no_loss")
cell_loss <- cell_loss[,-c(2)]
colnames(cell_loss) <- c("se_name","number_spec_ce","number_spec_object","spec_ce","spec_object","specifity")

# add specific SE
final_cell_specific_1 <- rbind(cell_gain[which(cell_gain$specifity=="gain"),],
                               cell_loss[which(cell_loss$specifity=="loss"),])

# add non specific SE
final_cell_specific_2 <- rbind(cell_gain,cell_loss)

final_cell_specific_3 <- final_cell_specific_2[!(final_cell_specific_2$se_name %in% final_cell_specific_1$se_name),]
final_cell_specific_3$specifity <- "none"
final_cell_specific_3 <- unique(final_cell_specific_3)

final_cell_specific <- rbind(final_cell_specific_1,final_cell_specific_3)
final_cell_specific$object_type <- "cell"

write.table(final_cell_specific,"results/final_se_cell_spec_r1.txt",quote = F,row.names = F,sep="\t")


