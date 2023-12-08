setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(data.table)
library(RColorBrewer)
library(ggchicklet)

cSEAdb <- readRDS("results/cSEAdb_r1.rds")

# chrom background
cancer_cell_meta <- cSEAdb$cell_cancer_meta
#------------------------
# all specific SE
#------------------------

# all cancer
chrom_size_3 <- read.table("hg38.chrom.sizes",sep="\t",header=F)

chrom_size <- read.table("hg38_chrom_loc_w_centromere.txt",sep="\t",header=F)
chrom_size$V5<- chrom_size$V3-chrom_size$V4
chrom_size_2 <- chrom_size[,c(1,4)]

colnames(chrom_size_2) <- c("V1","V2")
chrom_size_2$group <- "p1"

chrom_size <- chrom_size[,c(1,5)]
colnames(chrom_size) <- c("V1","V2")
chrom_size$group <- "p2"

chrom_size_final <- rbind(chrom_size_2,chrom_size)

chrom_size_final$V1 <- factor(chrom_size_final$V1,levels=rev(unique(chrom_size_final$V1)))

chrom_bg <- ggplot(chrom_size_final,aes(x=V1,y=V2,color=group))+
  geom_chicklet(radius = grid::unit(0.6, "mm"),
                width=0.5,
                fill="NA",
                color="black",size=0.2
                )+
  theme_classic()+
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank())+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, NA),
                     labels = c("0","50M","100M","150M","200M","NA"))

chrom_bg_2 <- ggplot(chrom_size_final,aes(x=V1,y=V2,fill=group))+
  geom_chicklet(radius = grid::unit(0.6, "mm"),
                width=0.5,
                fill="black",
                #color="black",size=0.2
  )+
  theme_classic()+
  theme(axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank())+
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, NA),
                     labels = c("0","50M","100M","150M","200M","NA"))



#-------------------------------
cancer_sepc_tmp<- cSEAdb$se_specificity
chrom_size <- read.table("hg38_chrom_loc_w_centromere.txt",sep="\t",header=F)
#----------------------------
# lung cancer
all_cancer_sepc_1 <- cancer_sepc_tmp[which(cancer_sepc_tmp$object_type=="cancer"),]
all_cancer_sepc <- all_cancer_sepc_1[which(all_cancer_sepc_1$spec_object %like% "Lung.adenocarcinoma"),]

# all cancer
all_cancer_sepc <- cancer_sepc_tmp[which(cancer_sepc_tmp$object_type=="cancer"),]

# A549 cancer
all_cancer_sepc_1 <- cancer_sepc_tmp[which(cancer_sepc_tmp$object_type=="cell"),]
all_cancer_sepc <- all_cancer_sepc_1[which(all_cancer_sepc_1$spec_object %like% "A549"),]



all_cancer_sepc <- unique(all_cancer_sepc[,c("se_name","specifity")])

chr_plot_df <-  data.frame(se_name=all_cancer_sepc$se_name,
                           V1=sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",1),
                           V2=as.integer(sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",2)),
                           V3=as.integer(sapply(strsplit(all_cancer_sepc$se_name,"_"),"[[",3)),
                           group=all_cancer_sepc$specifity)

chr_plot_df <- chr_plot_df[which(chr_plot_df$group!="none"),]

se_both <- intersect(chr_plot_df$se_name[which(chr_plot_df$group=="gain")],
                     chr_plot_df$se_name[which(chr_plot_df$group=="loss")])

chr_plot_df$group[which(chr_plot_df$se_name %in% se_both)] <- "both"

chr_plot_df_2 <- merge(chr_plot_df,chrom_size,by="V1")

chr_plot_df_2$V4 <- ifelse(chr_plot_df_2$V3.x+500000<=chr_plot_df_2$V3.y,chr_plot_df_2$V3.x+500000,chr_plot_df_2$V3.y)


#--------------
# diff location
chrom_bg+
  geom_segment(data=chr_plot_df_2[which(chr_plot_df_2$group=="gain"),],
               aes(x=V1,xend=V1,y=V2.x,yend=V4,color="gain"),
                      size = 4,position=position_nudge(x = 0.25))+
  geom_segment(data=chr_plot_df_2[which(chr_plot_df_2$group=="both"),],
               aes(x=V1,xend=V1,y=V2.x,yend=V4,color="both"),
               size = 4,position=position_nudge(x = 0))+
  geom_segment(data=chr_plot_df_2[which(chr_plot_df_2$group=="loss"),],
               aes(x=V1,xend=V1,y=V2.x,yend=V4,color="loss"),
               size = 4,position=position_nudge(x = -0.25))+
  scale_color_manual(values = c("gain"="#762a83",
                                "both"="#1b7837",
                                "loss"="#e08214"))+
  coord_flip()

#--------------
# V1
#
# chrom_bg+
#   geom_segment(data=chr_plot_df_2,aes(x=V1,xend=V1,y=V2.x,yend=V4,color=group),
#                size = 4)+
#   scale_color_manual(values = c(gain="#762a83",
#                                 both="#1b7837",
#                                 loss="#e08214"))+
#   coord_flip()
#
#
# # dark bg
# chrom_bg_2+geom_segment(data=chr_plot_df_2,aes(x=V1,xend=V1,y=V2.x,yend=V4,color=group),
#                       size = 4)+
#   scale_color_manual(values = c(gain="#762a83",
#                                 both="#4d9221",
#                                 loss="#e08214"))+
#   coord_flip()








