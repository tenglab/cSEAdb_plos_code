setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(data.table)
library(RColorBrewer)
library(chromoMap)
library(htmltools)
cSEAdb <- readRDS("results/cSEAdb_v4.rds")

# chrom background
cancer_cell_meta <- cSEAdb$cell_cancer_meta
#------------------------
# all specific SE
#------------------------

# all cancer
chrom_size <- read.table("hg38_chrom_loc_w_centromere.txt",sep="\t",header=F)

cancer_sepc_tmp<- cSEAdb$se_specificity

all_cancer_sepc_1 <- cancer_sepc_tmp[which(cancer_sepc_tmp$object_type=="cancer"),]
all_cancer_sepc <- all_cancer_sepc_1[which(all_cancer_sepc_1$spec_object %like% "Lung.adenocarcinoma"),]


all_cancer_sepc <- unique(all_cancer_sepc[,c("se_name","specifity")])

chr_plot_df <-  data.frame(se_name=all_cancer_sepc$se_name,
                           V1=sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",1),
                           V2=as.integer(sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",2)),
                           V3=as.integer(sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",3)),
                           group=all_cancer_sepc$specifity)

chr_plot_df <-  data.frame(se_name=all_cancer_sepc$se_name,
                           V1=sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",1),
                           V2=(as.integer(sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",2))+
                             as.integer(sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",3)))/2-10,
                           V3=(as.integer(sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",2))+
                                 as.integer(sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",3)))/2+10,
                           group=all_cancer_sepc$specifity)

chr_plot_df <- chr_plot_df[which(chr_plot_df$group!="none"),]

se_both <- intersect(chr_plot_df$se_name[which(chr_plot_df$group=="gain")],
                     chr_plot_df$se_name[which(chr_plot_df$group=="loss")])

chr_plot_df$group[which(chr_plot_df$se_name %in% se_both)] <- "both"


chrom_size_2 <-  read.table("hg38.chrom.sizes",sep="\t",header=F)

chrom_size_2$V3 <- 1
chrom_size_2 <- chrom_size_2[,c(1,3,2)]
chromoMap(list(chrom_size_2),list(chr_plot_df),
          n_win.factor = 3,
          #fixed.window=T,
          #window.size=20000,
          win.summary.display = T,
          # chr_width = 20,
          chr_length = 10,
          chr_color = "grey",
          data_based_color_map = T,
          data_type = "categorical",
          legend = T,
          lg_x = 200,
          lg_y = 500,
          #title=title,
          text_font_size = c(16),
          data_colors = list(c("lightskyblue","firebrick3","gold3")))


#-------------------------
# loop cancer name
cancer <- unique(cancer_cell_meta$Cancer)

cancer_sepc <- cSEAdb$se_specificity

output_list <- list()
for (c in 1: length(cancer)) {
  print(c)
  cancer_tmp <- cancer[c]

  cell_sepc_gain <- cancer_sepc[which(cancer_sepc$spec_object %like% cancer_tmp &
                                        cancer_sepc$specifity=="gain"),]
  cell_sepc_loss <- cancer_sepc[which(cancer_sepc$spec_object %like% cancer_tmp &
                                        cancer_sepc$specifity=="loss"),]

  both_se <- intersect(unique(cell_sepc_gain$se_name),unique(cell_sepc_loss$se_name))
  gain_se <- cell_sepc_gain$se_name[!(cell_sepc_gain$se_name %in%both_se)]
  loss_se <- cell_sepc_loss$se_name[!(cell_sepc_loss$se_name %in%both_se)]
  #none_se <- cancer_sepc$se_name[!(cancer_sepc$se_name %in% c(both_se,gain_se,loss_se))]

  if (length(both_se) >0) {
    cell_sepc <- data.frame(se_name=c(gain_se,loss_se,both_se),
                            group=c(rep("gain",length(gain_se)),
                                    rep("loss",length(loss_se)),
                                    rep("both",length(both_se))))
    se_spec_df <- separate(data = cell_sepc,se_name, into=c("V1","V2","V3"), sep='_',remove = F)
    se_spec_df$V2 <- as.integer(se_spec_df$V2)
    se_spec_df$V3 <- as.integer(se_spec_df$V3)
  } else {
    cell_sepc <- data.frame(se_name=c(gain_se,loss_se),
                            group=c(rep("gain",length(gain_se)),
                                    rep("loss",length(loss_se))))
    se_spec_df <- separate(data = cell_sepc,se_name, into=c("V1","V2","V3"), sep='_',remove = F)
    se_spec_df$V2 <- as.integer(se_spec_df$V2)
    se_spec_df$V3 <- as.integer(se_spec_df$V3)
  }


  output_list[[c]] <- se_spec_df
  names(output_list)[[c]] <- cancer[c]
}
saveRDS(output_list,"results/cancer_spec_for_chr_plot.rds")

# plot chromosome map

cancer_spec_se <- readRDS("results/cancer_spec_for_chr_plot.rds")
chrom_size <- read.table("hg38_chrom_loc_w_centromere.txt",sep="\t",header=F)


for (i in 1:length(cancer_spec_se)) {

  se_spec_df <- cancer_spec_se[[i]]
  title <- names(cancer_spec_se)[[i]]

  if (length(unique(cancer_spec_se[[i]]$group))==2) {
    plot_tmp <- chromoMap(list(chrom_size),list(se_spec_df),
                          n_win.factor = 2,
                          # chr_width = 20,
                          chr_length = 6,
                          chr_color = "grey",
                          data_based_color_map = T,
                          data_type = "categorical",
                          legend = T,
                          lg_x = 200,
                          lg_y = 500,
                          title=title,
                          text_font_size = c(16),
                          data_colors = list(c("firebrick3","lightskyblue")))
  } else {
    plot_tmp <- chromoMap(list(chrom_size),list(se_spec_df),
                          n_win.factor = 2,
                          # chr_width = 20,
                          chr_length = 6,
                          chr_color = "grey",
                          data_based_color_map = T,
                          data_type = "categorical",
                          legend = T,
                          lg_x = 200,
                          lg_y = 500,
                          title=title,
                          text_font_size = c(16),
                          data_colors = list(c("firebrick3","lightskyblue","gold3")))
  }


  out_name <- paste0("results/cancer_spec_chromap/",title,".html")
  htmltools::save_html(plot_tmp, file = out_name)

}

cancer_spec_se <- readRDS("results/cell_spec_for_chr_plot.rds")
for (i in 1:length(cancer_spec_se)) {

  se_spec_df <- cancer_spec_se[[i]]
  title <- names(cancer_spec_se)[[i]]

  if (length(unique(cancer_spec_se[[i]]$group))==2) {
    plot_tmp <- chromoMap(list(chrom_size),list(se_spec_df),
                          n_win.factor = 2,
                          # chr_width = 20,
                          chr_length = 6,
                          chr_color = "grey",
                          data_based_color_map = T,
                          data_type = "categorical",
                          legend = T,
                          lg_x = 200,
                          lg_y = 500,
                          title=title,
                          text_font_size = c(16),
                          data_colors = list(c("firebrick3","lightskyblue")))
  } else {
    plot_tmp <- chromoMap(list(chrom_size),list(se_spec_df),
                          n_win.factor = 2,
                          # chr_width = 20,
                          chr_length = 6,
                          chr_color = "grey",
                          data_based_color_map = T,
                          data_type = "categorical",
                          legend = T,
                          lg_x = 200,
                          lg_y = 500,
                          title=title,
                          text_font_size = c(16),
                          data_colors = list(c("firebrick3","lightskyblue","gold3")))
  }


  out_name <- paste0("results/cell_spec_chromap/",title,".html")
  htmltools::save_html(plot_tmp, file = out_name)

}




cell_sepc <- cSEAdb$se_cancer_specificity[,c(1,2,5)]



cell_sepc$group <- NA
cell_sepc$group[which(cell_sepc$cell_specific_gain=="yes" &
                        cell_sepc$cell_specific_loss=="no")] <- "gain"
cell_sepc$group[which(cell_sepc$cell_specific_gain=="no" &
                        cell_sepc$cell_specific_loss=="yes")] <- "loss"
cell_sepc$group[which(cell_sepc$cell_specific_gain=="yes" &
                        cell_sepc$cell_specific_loss=="yes")] <- "both"
cell_sepc$group[which(cell_sepc$cell_specific_gain=="no" &
                        cell_sepc$cell_specific_loss=="no")] <- "none"

# make se_df for plot

se_spec_df <- separate(data = cell_sepc[,-c(2,3)],se_name, into=c("V1","V2","V3"), sep='_',remove = F)

se_reduce <- read.table("../interface_table/se_reduce.bed",sep="\t",header=F)
se_reduce <- se_reduce[,c(4,1,2,3)]
colnames(se_reduce)[1] <- "se_name"


test_1 <- merge(se_reduce,se_spec_df,by="se_name")
test_1 <- test_1[,c(1,2,3,4,8)]



chromoMap(list(chrom_size),list(test_1))

chromoMap(list(chrom_size),list(test_1),
          #n_win.factor = 2,
          # chr_width = 20,
          chr_length = 6,
          #chr_color = "white",
          data_based_color_map = T,
          data_type = "categorical",
          legend = T,
          lg_x = 200,
          lg_y = 500,
          text_font_size = c(16),
          data_colors = list(c("firebrick3","green3","gold2","lightskyblue")))







#----------------------
# make SE annotation file
#----------------------
# cell specific example:MCF7

cell_sepc <- cSEAdb$se_cell_specificity

# extract all MCF7 spec SE
se_spec_gain <- cell_sepc$se_name[which(cell_sepc$gain_ce %like% "MCF7")]
se_spec_loss <- cell_sepc$se_name[which(cell_sepc$loss_ce %like% "MCF7")]

# se have gain and loss
se_gain_loss <- intersect(se_spec_gain,se_spec_loss)
se_gain_only <- se_spec_gain[!(se_spec_gain %in% se_gain_loss)]
se_loss_only <- se_spec_loss[!(se_spec_loss %in% se_gain_loss)]

# make se_df for plot
se_spec_df <- data.frame(se_name=c(se_gain_only,se_loss_only,se_gain_loss),
                         group=c(rep("gain",length(se_gain_only)),
                                 rep("loss",length(se_loss_only)),
                                 rep("both",length(se_gain_loss))))

se_spec_df_2 <- separate(data = se_spec_df,se_name, into=c("V1","V2","V3"), sep='_',remove = F)

chromoMap(list(chrom_size),list(se_spec_df_2),
          n_win.factor = 2,
          # chr_width = 20,
          # chr_length = 5,
          data_based_color_map = T,
          data_type = "categorical",
          data_colors = list(c("firebrick3","lightskyblue","green3")))




# reformat
se_df <- se_df[,c(4,1,2,3)]

se_test <- head(se_df)

chromoMap(list(chrom_size),list(se_df),
          n_win.factor = 3)


chromoMap(list(chrom_size),list(se_test))








