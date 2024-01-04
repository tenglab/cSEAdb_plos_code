setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/")
library(ggplot2)
library(tidyr)
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(grid)
library(ggrepel)
library(gridExtra)
library(factoextra)
library(pheatmap)
library(VennDiagram)
library(dendextend)
library(ggvenn)
library(matrixStats)

mem.maxVSize(vsize = Inf)
#-----------------------------------------------------
# fig1a pairwise
#------------------------

pairsiwe_diff <- read.table("figs_df/fig1a_binary_dase_pairwise_diff.txt",sep="\t",header=T)

pairsiwe_diff_2 <- aggregate(count~data_pair+group,pairsiwe_diff[,c(1,2,4)],FUN=sum)

# calculate paired p-value

ttest <- t.test(pairsiwe_diff_2$count[which(pairsiwe_diff_2$group=="Non-differential")],
                pairsiwe_diff_2$count[which(pairsiwe_diff_2$group=="Differential")],
                paired = T)

plot_df <- data.frame(group=c("Non-differential",
                              "Differential"),
                      mean=c(mean(pairsiwe_diff_2$count[which(pairsiwe_diff_2$group=="Non-differential")]),
                             mean(pairsiwe_diff_2$count[which(pairsiwe_diff_2$group=="Differential")])),
                      sd=c(sd(pairsiwe_diff_2$count[which(pairsiwe_diff_2$group=="Non-differential")]),
                           sd(pairsiwe_diff_2$count[which(pairsiwe_diff_2$group=="Differential")])))

plot_df$group <- factor(plot_df$group,levels=c("Non-differential","Differential"))

ggplot(plot_df,
       aes(x=group,y=mean,fill=group))+
  geom_bar(stat="identity",color="black",width=0.5, position=position_dodge(0.6))+
  geom_errorbar(aes(ymin=mean,ymax=mean+sd),
                width=0.1,
                position=position_dodge(.6))+
  theme_classic()+
  ylab("Number of SE")+
  scale_fill_manual(values = c("grey70","purple4"))

write.table(pairsiwe_diff_2,"figs_df/supplementary_data/fig1a_1.txt",sep = "\t",row.names = F,quote = F)
write.table(plot_df,"figs_df/supplementary_data/fig1a_2.txt",sep = "\t",row.names = F,quote = F)

#-----------------------------------------------------
# fig1c CE count
#-----------------------------------------------------
ce_coverage<- read.table("results/ce_maxcover_60cell.txt",sep="\t",header=T)
ce_binary <- read.table("results/ce_binary_mod.txt",sep="\t",header=T)

# change name
colnames(ce_binary)[1] <- "ce_name"


ce_peak_w_percent <- merge(ce_coverage[,-3],ce_binary,by=c("ce_name","se_name"))
#write.table(ce_peak_w_percent,"results/s1_ce_binary_w_percent.txt",sep="\t",row.names = F,quote=F)

ce_larger_3_percent <- unique(ce_peak_w_percent$ce_name[which(ce_peak_w_percent$max_percent > 3)])
ce_peak_plot_df <- data.frame(ce_name=ce_peak_w_percent$ce_name,
                              max_percent=ce_peak_w_percent$max_percent,
                              varance=round(matrixStats::rowVars(as.matrix(ce_peak_w_percent[,c(4:63)])),digits = 3),
                              n_1=rowSums(as.matrix(ce_peak_w_percent[,c(4:63)])))

lable_x <- c()
for (i in 1:60) {
  lable_x<- c(lable_x,unique(ce_peak_plot_df$varance[which(ce_peak_plot_df$n_1==i)]))
}

ce_peak_plot_df$n_1 <- factor(ce_peak_plot_df$n_1,levels = sort(unique(ce_peak_plot_df$n_1)))

ce_peak_plot_df_3 <- ce_peak_plot_df %>%
  group_by(n_1) %>%
  tally()


ce_peak_plot_df_keep <- ce_peak_plot_df[which(ce_peak_plot_df$ce_name %in% ce_larger_3_percent),]

ce_peak_plot_df_4 <- ce_peak_plot_df_keep %>%
  group_by(n_1) %>%
  tally()

ce_peak_plot_df_4$group <- "Coverage more than 3%"
ce_peak_plot_df_5 <- data.frame(n_1=ce_peak_plot_df_4$n_1,
                                n=ce_peak_plot_df_3$n-ce_peak_plot_df_4$n,
                                group="Coverage less than 3%")
ce_peak_plot_df_6 <- rbind(ce_peak_plot_df_4,ce_peak_plot_df_5)
ce_peak_plot_df_6$group <- as.factor(ce_peak_plot_df_6$group)


ggplot(ce_peak_plot_df_6, aes(x=as.factor(n_1), y=n,fill=group)) +
  geom_bar(position="stack", stat="identity")+
  theme_classic()+
  scale_fill_manual(values = brewer.pal(2,"Blues"))+
  # geom_text(aes(label=n),
  #           position=position_stack(vjust = 0),size=3)+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank())+
  xlab("Number of cell CE present")+
  ylab("Number of CE")

write.table(ce_peak_plot_df_6,"figs_df/supplementary_data/fig1c.txt",sep = "\t",row.names = F,quote = F)


#-------
# fig1d
binary_more_3 <- ce_peak_plot_df_6[which(ce_peak_plot_df_6$group=="Coverage more than 3%"),]
ggplot(binary_more_3, aes(x=as.factor(n_1), y=n)) +
  geom_bar(stat="identity",fill=brewer.pal(2,"Blues")[2])+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  # theme(axis.text.x = element_blank(),
  #       axis.text.y = element_blank(),
  #       axis.ticks.x = element_blank())+
  xlab("Number of cell CE present")+
  ylab("Number of CE")
#theme(axis.text.x = element_blank())

write.table(binary_more_3,"figs_df/supplementary_data/fig1d.txt",sep = "\t",row.names = F,quote = F)


#-----------------------------------------------------
# fig1e heatmap
#-----------------------------------------------------
# CE PCA
all_ce_signal <- read.table("results/ce_signal_normalized_60cell.txt",sep="\t",header=T)
selected_ce_signal <- read.table("results/ce_signal_normalized_60cell_filtered.txt",sep="\t",header=T)
all_enhancer_signal <- read.table("results/all_enhancer_signal_normalized_60cell.txt",sep="\t",header=T)
cell_meta <- read.table("results/cell_meta.txt",sep="\t",header=T)

non_ce <- all_enhancer_signal[!(all_enhancer_signal$merge_e_name %in%all_ce_signal$merge_e_name),]

cell_meta$Cell_line <- gsub("-","\\.",cell_meta$Cell_line)
cell_meta$Cell_line[1] <-"X786.0"

# chose cancer
cancer_cell_count <- data.frame(table(cell_meta$Cancer))

# keep cancer have at least 3 cells
cancer_keep <- cancer_cell_count$Var1[which(cancer_cell_count$Freq>=3)]
cell_keep <- cell_meta$Cell_line[which(cell_meta$Cancer %in% cancer_keep)]

ce_binary <- read.table("results/s1_ce_binary_w_percent.txt",sep="\t",header=T)


# subset kept cell
selected_ce_signal_2 <- unique(selected_ce_signal[,c("merge_e_name",cell_keep)])
all_enhancer_signal_2 <- unique(all_enhancer_signal[,c("merge_e_name",cell_keep)])
non_ce_2 <- unique(non_ce[,c("merge_e_name",cell_keep)])

# all enhancer
ce_keep_signal <- all_enhancer_signal_2
ce_keep_signal <- selected_ce_signal_2
ce_keep_signal <- non_ce_2

#ce_keep_signal[,-1] <- log2(ce_keep_signal[,-1]+1)
ce_keep_signal[,-1] <- ce_keep_signal[,-1]

# correlation heatmap
heatmap_df <- ce_keep_signal
# raw signal use correlation
cor_matrix <- cor(heatmap_df[,-1])

write.table(cor_matrix,"figs_df/supplementary_data/fig1e_non_ce.txt",sep = "\t",quote = F)

cancer_df <- cell_meta[which(cell_meta$Cancer %in% cancer_keep),c(1,2)]
rownames(cancer_df) <- cancer_df$Cell_line
cancer_df <- cancer_df[,-1,drop=F]


cancer_color <- c(brewer.pal(8,"Dark2"))

names(cancer_color) <- unique(cancer_df$Cancer)

#tissue_color <- c(brewer.pal(7,"Dark2"))
#names(tissue_color) <- unique(cancer_df$Tissue)

# annot_color <- list(Cancer=cancer_color,
#                      Tissue=tissue_color)

annot_color <- list(Cancer=cancer_color)

# for distance
ph <- pheatmap(cor_matrix,
               color=brewer.pal(9,"Blues"),
               #breaks = c(0.2,0.3,0.4,0.5,0.7,0.9,1),
               treeheight_row=0,
               treeheight_col=0,
               annotation_col = cancer_df,
               annotation_colors= annot_color
)
h_row_name <- ph$tree_row$labels[ph$tree_row$order]
#write.table(h_row_name,"figs_df/fig1e_select_CE_mixture_heat_rowname.txt",sep="\t",quote = F,row.names = F,col.names = F)

#-----------------------------------------------------
# fig2b density band
#-----------------------------------------------------
se_ce_tmp <- read.table("results/ce_signal_normalized_60cell_filtered.txt",sep="\t",header=T)
mod2 <- read.table("results/filtered_ce_mixture_model_r1.txt",sep="\t",header=T)

cat1_ce <- unique(mod2$ce_name[which(mod2$rowmean==0)])
cat2_ce <- unique(mod2$ce_name[which(mod2$rowmean>0 & mod2$rowmean<1)])
cat3_ce <- unique(mod2$ce_name[which(mod2$rowmean==1)])


cat1_signal <- se_ce_tmp[which(se_ce_tmp$merge_e_name %in% cat1_ce),]
cat2_signal <- se_ce_tmp[which(se_ce_tmp$merge_e_name %in% cat2_ce),]
cat3_signal <- se_ce_tmp[which(se_ce_tmp$merge_e_name %in% cat3_ce),]

ce_name <- c(cat2_ce)
plot_df <- c()
for (ce in 1:length(ce_name)){
  if (ce %% 2000 ==0){
    print(ce)
  }

  tmp_signal <- unique(se_ce_tmp[which(se_ce_tmp$merge_e_name ==ce_name[ce]),c(2:61)])
  den_tmp <-density(as.numeric(log2(tmp_signal+1)),from=0,n=250)
  tmp_1 <- c(den_tmp$x,den_tmp$y)
  plot_df <- rbind(plot_df,tmp_1)
}

plot_df_t1 <- as.data.frame(plot_df)


#write.table(plot_df_t1,"figs_df/fig2b_density_cat2.txt",sep="\t",row.names = F,quote = F)


# cat 2 only
plot_df <- read.table("figs_df/fig2b_density_cat2.txt",sep="\t",header=T)

densities.qtiles <- data.frame(seq=seq(1:250),
                               median_y=colMedians(as.matrix(plot_df[,c(251:500)])),
                               q5=colQuantiles(as.matrix(plot_df[,c(251:500)]),
                                               probs = 0.05),
                               q95=colQuantiles(as.matrix(plot_df[,c(251:500)]),
                                                probs=0.95),
                               x_median=colMedians(as.matrix(plot_df[,c(1:250)])))
write.table(densities.qtiles,"figs_df/supplementary_data/fig2b_cat2.txt",sep="\t",row.names = F,quote = F)

# cat 1 and 3
plot_df <- read.table("figs_df/fig2b_density_cat3.txt",sep="\t",header=T)
densities.qtiles <-
  plot_df %>%
  #rename(log_signal=x, dens = y) %>%
  ungroup() %>%
  group_by(seq) %>%
  summarise(median_y = median(y),
            q5 = quantile(y,0.05),
            q95=quantile(y,0.95),
            x_median=median(x))

write.table(densities.qtiles,"figs_df/supplementary_data/fig2b_cat3.txt",sep="\t",row.names = F,quote = F)


# plot
ggplot(densities.qtiles, aes(x=x_median, y=median_y)) +
  geom_line(color=brewer.pal(5,"Dark2")[2])+
  geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.5, fill = brewer.pal(5,"Dark2")[2])+
  theme_classic()+
  scale_color_manual(values = brewer.pal(5,"Dark2")[2])+
  xlab("log2 scale signal")


#----------------
# fig2c pie chart
#----------------
pie_df <- data.frame(count=c(length(cat1_ce),
                             length(cat2_ce),
                             length(cat3_ce)),
                     group=c("Cat1","Cat2","Cat3"))


pie_df$percent <- round(pie_df$count/sum(pie_df$count)*100,digits = 2)

ggplot(pie_df,aes(x="",y=count,fill=group))+
  geom_bar(stat="identity",width=0.5)+
  coord_polar("y")+
  theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_blank(),
    axis.text.x=element_blank()
  )+
  geom_label_repel(aes(label = paste0(count,"(",percent,"%)")),
            position=position_stack(vjust = 0.5),size=4)+
  scale_fill_manual(values = brewer.pal(5,"Dark2"))

#--------------------
# fig2d venn diagram
#

diff_long_signal <- read.table("results/binary_model_diff_r1.txt",sep="\t",header=T)

model_activate <-  unique(diff_long_signal[which(diff_long_signal$mod==1),c(1,2)])
peak_activate <- unique(diff_long_signal[which(diff_long_signal$binary==1),c(1,2)])

model_activate$pair <- paste(model_activate$ce_name,model_activate$cell,sep=":")
peak_activate$pair <- paste(peak_activate$ce_name,peak_activate$cell,sep=":")

v <- venn.diagram(x=list(model_activate$pair,peak_activate$pair),
             category.names = c("model",paste("peak", "calling",sep="\n")),
             fill = c(brewer.pal(8,"Paired")[7],brewer.pal(8,"Blues")[5]),
             alpha=c(1,1),
             cex=1,
             cat.cex=0,
             filename=NULL,
             #filename="figs_df/fig2d_venn_no_text.png",
             disable.logging=T)

grid.newpage()
grid.draw(v)

53056/(53056+727962+39980)
727962/(53056+727962+39980)
39980/(53056+727962+39980)

#------------------------
# fig3 compare of binary
#------------------------
diff_long_signal <- read.table("results/binary_model_diff_r1.txt",sep="\t",header=T)

# signal boxplot
# remove signal= 0
signal_diff_plot <- diff_long_signal[which(diff_long_signal$signal!=0),]
signal_diff_plot$pair <- paste(signal_diff_plot$ce_name,signal_diff_plot$cell,sep="_")
signal_diff_plot$group <- factor(signal_diff_plot$group,levels=c("model only","peak only","common"))
ggplot(signal_diff_plot,aes(x=group,y=log2(signal+1),fill=group))+
  #geom_boxplot(width=0.5)+
  geom_violin()+
  theme_classic()+
  scale_fill_manual(values=c(brewer.pal(8,"Paired")[8],
                             brewer.pal(8,"Blues")[8],
                             "grey50"))

write.table(signal_diff_plot,"figs_df/supplementary_data/fig3a.txt",sep="\t",row.names = F,quote = F)


nrow(signal_diff_plot[which(signal_diff_plot$group=="model only"),])/nrow(signal_diff_plot)*100
nrow(signal_diff_plot[which(signal_diff_plot$group=="peak only"),])/nrow(signal_diff_plot)*100
nrow(signal_diff_plot[which(signal_diff_plot$group=="common"),])/nrow(signal_diff_plot)*100

#--------------
# fig3b: link count
#--------------
link_summary <- read.table("results/fig3b_ce_link_count_r1.txt",sep="\t",header=T)
# calculate p-value
link_summary_2 <- link_summary[!(link_summary$signal==0&link_summary$group=="common"),]


link_summary_2$mean_link_score <- log2(link_summary_2$mean_link_score+1)
link_summary_2$n_link <- log2(link_summary_2$n_link+1)

(mean(link_summary_2$n_link[which(link_summary_2$cell=="A549"&link_summary_2$group=="model only")])+
    mean(link_summary_2$n_link[which(link_summary_2$cell=="HCT.116"&link_summary_2$group=="model only")])+
    mean(link_summary_2$n_link[which(link_summary_2$cell=="K.562"&link_summary_2$group=="model only")])+
    mean(link_summary_2$n_link[which(link_summary_2$cell=="MCF7"&link_summary_2$group=="model only")]))/4

(mean(link_summary_2$n_link[which(link_summary_2$cell=="A549"&link_summary_2$group=="peak only")])+
    mean(link_summary_2$n_link[which(link_summary_2$cell=="HCT.116"&link_summary_2$group=="peak only")])+
    mean(link_summary_2$n_link[which(link_summary_2$cell=="K.562"&link_summary_2$group=="peak only")])+
    mean(link_summary_2$n_link[which(link_summary_2$cell=="MCF7"&link_summary_2$group=="peak only")]))/4

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
                           gain_loss=ttest_1$p.value,
                           gain_not=ttest_2$p.value,
                           loss_not=ttest_3$p.value)
  pvalue <- rbind(pvalue,pvalue_tmp)
}
pvalue
table(link_summary_2$cell,link_summary_2$group)

link_summary_2$group <- factor(link_summary_2$group,levels = c("model only","peak only","common"))
ggplot(link_summary_2,aes(x=cell,y=n_link,fill=group))+
  geom_boxplot(outlier.shape = NA,width=0.5,position=position_dodge(0.6))+
  #geom_violin()+
  #geom_point(position=position_jitterdodge(jitter.width=0.1,dodge.width=0.6),
  #           cex=0.3)+
  theme_classic()+
  ylim(c(0,10))+
  #ylab("log2(mean_link_score+1)")+
  scale_fill_manual(values = c(brewer.pal(8,"Paired")[8],
                               brewer.pal(8,"Blues")[8],
                               "grey50"))+
  scale_color_manual(values = c(brewer.pal(8,"Paired")[8],
                                brewer.pal(8,"Blues")[8],
                                "grey50"))

write.table(link_summary_2,"figs_df/supplementary_data/fig3b.txt",sep="\t",row.names = F,quote = F)

#--------------
# fig3c: target gene count
#--------------
link_gene_summary <- read.table("figs_df/fig3c_ce_link_gene_w_promoter_summary_r1.txt",sep="\t",header=T)
# calculate p-value
link_summary_2 <- link_gene_summary[!(link_gene_summary$signal==0 &link_gene_summary$group=="common"),]

(mean(link_summary_2$n_gene[which(link_summary_2$cell=="A549"&link_summary_2$group=="model only")])+
  mean(link_summary_2$n_gene[which(link_summary_2$cell=="HCT.116"&link_summary_2$group=="model only")])+
  mean(link_summary_2$n_gene[which(link_summary_2$cell=="K.562"&link_summary_2$group=="model only")])+
  mean(link_summary_2$n_gene[which(link_summary_2$cell=="MCF7"&link_summary_2$group=="model only")]))/4

(mean(link_summary_2$n_gene[which(link_summary_2$cell=="A549"&link_summary_2$group=="peak only")])+
    mean(link_summary_2$n_gene[which(link_summary_2$cell=="HCT.116"&link_summary_2$group=="peak only")])+
    mean(link_summary_2$n_gene[which(link_summary_2$cell=="K.562"&link_summary_2$group=="peak only")])+
    mean(link_summary_2$n_gene[which(link_summary_2$cell=="MCF7"&link_summary_2$group=="peak only")]))/4

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
                           gain_loss=ttest_1$p.value,
                           gain_not=ttest_2$p.value,
                           loss_not=ttest_3$p.value)
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
  ylab("log2(mean_link_score+1)")+
  scale_fill_manual(values = c(brewer.pal(8,"Paired")[8],
                               brewer.pal(8,"Blues")[8],
                               "grey50"))+
  scale_color_manual(values = c(brewer.pal(8,"Paired")[8],
                                brewer.pal(8,"Blues")[8],
                                "grey50"))

write.table(link_summary_2,"figs_df/supplementary_data/fig3c.txt",sep="\t",row.names = F,quote = F)

#------------
# fig_3d
#------------

# make plot data frame
depmap_g_ranking_ratio <- read.table("results/depmap_gene_ranking_ratio.txt",sep="\t",header=T)
depmap_g_ranking_ratio_2 <- depmap_g_ranking_ratio[which(depmap_g_ranking_ratio$cell %in%
                                                           c("ACH-000681","ACH-000971","ACH-000551","ACH-000019")),]

depmap_rank_2 <- as.data.frame(t(depmap_g_ranking_ratio_2))
colnames(depmap_rank_2) <- c("HCT.116","MCF7","A549","K.562")
depmap_rank_2 <- depmap_rank_2[-1,]
depmap_rank_2$gene_name <- row.names(depmap_rank_2)

link_gene <- read.table("../figs_df/fig3c_ce_link_gene_w_promoter.txt",sep="\t",header=T)
depmap_2_long <- gather(depmap_rank_2,cell,ratio,HCT.116:K.562)
link_gene_summary <- merge(link_gene,depmap_2_long,by=c("cell","gene_name"))

# calculate p-value
link_summary_2 <- unique(link_gene_summary[!(link_gene_summary$signal==0 &link_gene_summary$group=="common"),])

link_summary_2$ratio <- as.numeric(link_summary_2$ratio)
link_summary_2$group <- factor(link_summary_2$group,levels=c("b0_m1","b1_m0","common"))
ggplot(link_summary_2,aes(x=cell,y=ratio,fill=group))+
  #geom_boxplot(outlier.shape = NA,width=0.5,position=position_dodge(0.6))+
  geom_boxplot(outlier.shape = NA,width=0.5,position=position_dodge(0.6))+
  #geom_violin(outlier.shape = NA,width=0.5,position=position_dodge(0.6))+
  #geom_point(position=position_jitterdodge(jitter.width=0.1,dodge.width=0.6),
  #           cex=0.3)+
  theme_classic()+
  #ylim(c(0,25))+
  ylab("depmap ranking ratio")+
  scale_fill_manual(values = c(brewer.pal(8,"Paired")[8],
                               brewer.pal(8,"Blues")[8],
                               "grey50"))+
  scale_color_manual(values = c(brewer.pal(8,"Paired")[8],
                                brewer.pal(8,"Blues")[8],
                                "grey50"))

write.table(link_summary_2,"figs_df/supplementary_data/fig3d.txt",sep="\t",row.names = F,quote = F)


(mean(link_summary_2$ratio[which(link_summary_2$cell=="A549"&link_summary_2$group=="b0_m1")])+
    mean(link_summary_2$ratio[which(link_summary_2$cell=="HCT.116"&link_summary_2$group=="b0_m1")])+
    mean(link_summary_2$ratio[which(link_summary_2$cell=="K.562"&link_summary_2$group=="b0_m1")])+
    mean(link_summary_2$ratio[which(link_summary_2$cell=="MCF7"&link_summary_2$group=="b0_m1")]))/4

(mean(link_summary_2$ratio[which(link_summary_2$cell=="A549"&link_summary_2$group=="b1_m0")])+
    mean(link_summary_2$ratio[which(link_summary_2$cell=="HCT.116"&link_summary_2$group=="b1_m0")])+
    mean(link_summary_2$ratio[which(link_summary_2$cell=="K.562"&link_summary_2$group=="b1_m0")])+
    mean(link_summary_2$ratio[which(link_summary_2$cell=="MCF7"&link_summary_2$group=="b1_m0")]))/4



cell_line <- unique(link_summary_2$cell)
pvalue <- data.frame()
for (c in 1:4) {
  tmp_link <- link_summary_2[which(link_summary_2$cell==cell_line[c]),]
  ttest_1 <- wilcox.test(tmp_link$ratio[which(tmp_link$group=="b0_m1")],
                         tmp_link$ratio[which(tmp_link$group=="b1_m0")],
                         alternative="greater")
  
  ttest_2 <- wilcox.test(tmp_link$ratio[which(tmp_link$group=="b0_m1")],
                         tmp_link$ratio[which(tmp_link$group=="common")],
                         alternative="greater")
  
  ttest_3 <- wilcox.test(tmp_link$ratio[which(tmp_link$group=="b1_m0")],
                         tmp_link$ratio[which(tmp_link$group=="common")],
                         alternative="less")
  pvalue_tmp <- data.frame(cell=cell_line[c],
                           gain_loss=ttest_1$p.value,
                           gain_not=ttest_2$p.value,
                           loss_not=ttest_3$p.value)
  pvalue <- rbind(pvalue,pvalue_tmp)
}
pvalue
table(link_summary_2$cell,link_summary_2$group)

#------------------------
# fig4a
#------------------------
se_cluster_score <- read.table("figs_df/fig4a_se_cell_cutoff_vi_score_r1.txt",sep="\t",header=T)
ggplot(se_cluster_score,aes(x=cut,y=vi_score,group=1),alpha = 0.6)+
  geom_point()+
  geom_smooth(method="loess",se=F,span = 0.8,color="red")+
  theme_classic()+
  #ylim(c(2.15,2.3))+
  theme(axis.text.x = element_text(size=10,angle=90),
        axis.title = element_text(size=16))

write.table(se_cluster_score,"figs_df/supplementary_data/fig4a.txt",sep="\t",row.names = F,quote = F)


# barplot
se_spec_cancer_final <- read.table("results/final_se_cancer_spec_r1.txt",sep="\t",header=T)

# gain
gain_df <- se_spec_cancer_final[which(se_spec_cancer_final$specifity=="gain"),c("se_name","number_spec_object","spec_object")]
gain_se <- unique(gain_df$se_name)
gain_plot <- data.frame()
for (se in 1:length(gain_se)) {
  g_tmp <- gain_df[which(gain_df$se_name==gain_se[se]),]
  uniq_cancer <- unique(unlist(strsplit(g_tmp$spec_object,",")))
  gain_plot_tmp <- data.frame(se_name=gain_se[se],
                              n_unique_cancer=length(uniq_cancer))
  gain_plot <- rbind(gain_plot,gain_plot_tmp)
}

# loss
loss_df <- se_spec_cancer_final[which(se_spec_cancer_final$specifity=="loss"),c("se_name","number_spec_object","spec_object")]
loss_se <- unique(loss_df$se_name)
loss_plot <- data.frame()
for (se in 1:length(loss_se)) {
  l_tmp <- loss_df[which(loss_df$se_name==loss_se[se]),]
  uniq_cancer <- unique(unlist(strsplit(l_tmp$spec_object,",")))
  loss_plot_tmp <- data.frame(se_name=loss_se[se],
                              n_unique_cancer=length(uniq_cancer))
  loss_plot <- rbind(loss_plot,loss_plot_tmp)
}

# plot df
gp_df <- data.frame(table(gain_plot$n_unique_cancer))
lp_df <- data.frame(table(loss_plot$n_unique_cancer))

se_both <- intersect(gain_plot$se_name,loss_plot$se_name)

length(se_both)
length(gain_plot$se_name[!(gain_plot$se_name %in% se_both)])
length(loss_plot$se_name[!(loss_plot$se_name %in% se_both)])+length(se_both)

write.table(gp_df,"figs_df/fig4d_bar_gain.txt",sep="\t",quote = F,row.names = F)
write.table(lp_df,"figs_df/fig4e_bar_loss.txt",sep="\t",quote = F,row.names = F)

gp_df <- read.table("figs_df/fig4d_bar_gain.txt",sep="\t",header=T)
lp_df <- read.table("figs_df/fig4e_bar_loss.txt",sep="\t",header=T)

gp_df$Var1 <- factor(gp_df$Var1,levels=gp_df$Var1)
ggplot(gp_df,aes(x=Var1,y=Freq))+
  geom_bar(stat="identity",fill="firebrick")+
  theme_classic()

lp_df$Var1 <- factor(lp_df$Var1,levels=lp_df$Var1)
ggplot(lp_df,aes(x=Var1,y=Freq))+
  geom_bar(stat="identity",fill="skyblue")+
  theme_classic()


#------------------------
# fig4 assay and boxplot
#------------------------
# fig4b-c

# signal boxplot
plot_df <- read.table("figs_df/fig4b_reporter_assay_r1.txt",sep="\t",header=T)

plot_df_erna <- read.table("figs_df/fig4c_eRNA_r1.txt",sep="\t",header=T)

# pvalue
# calculate p-value
cell_line <- c("A549","MCF7","HCT-116","K-562")
plot_df_2 <- plot_df

cell_line <- c("A549_pro_seq","K562_PROseq_GEO")
plot_df_2 <- plot_df_erna

pvalue <- data.frame()
for (c in 1:length(cell_line)) {
  tmp_link <- plot_df_2[which(plot_df_2$cell==cell_line[c]),]
  ttest_1 <- wilcox.test(tmp_link$signal[which(tmp_link$group=="gain")],
                         tmp_link$signal[which(tmp_link$group=="loss")])

  ttest_2 <- wilcox.test(tmp_link$signal[which(tmp_link$group=="gain")],
                         tmp_link$signal[which(tmp_link$group=="non_spec")])

  ttest_3 <- wilcox.test(tmp_link$signal[which(tmp_link$group=="loss")],
                         tmp_link$signal[which(tmp_link$group=="non_spec")])
  pvalue_tmp <- data.frame(cell=cell_line[c],
                           gain_loss=ttest_1$p.value,
                           gain_not=ttest_2$p.value,
                           loss_not=ttest_3$p.value)
  pvalue <- rbind(pvalue,pvalue_tmp)
}
pvalue
table(plot_df_2$cell,plot_df_2$group)

ggplot(plot_df_2[which(plot_df_2$cell %in% cell_line),],aes(x=cell,y=log2(signal+1)/width*1000,fill=group))+
  geom_boxplot(width=0.5,position=position_dodge(0.6),outlier.shape = NA)+
  # geom_text(aes(label=paste0(rate_ce,"%")),
  #           vjust=-0.5,position = position_dodge(0.9),size=3)+
  theme_classic()+
  ylim(c(0,3))+
  scale_fill_manual(values=c("firebrick","lightskyblue","grey50"))

plot_df_3 <- plot_df_2[which(plot_df_2$cell %in% cell_line),]

plot_df_3$signal <- log2(plot_df_3$signal+1)/plot_df_3$width*1000

write.table(plot_df_3,"figs_df/supplementary_data/fig4b.txt",sep="\t",row.names = F,quote = F)


#------------------------
# fig4f SE summary
#------------------------
cSEAdb <- readRDS("results/cSEAdb_r1.rds")
se_sepc <- cSEAdb$se_specificity
se_sepc_cancer <- se_sepc[which(se_sepc$object_type=="cancer"),]
se_all <- cSEAdb$se_bed
# summary SE based on specificity
se_gain <- unique(se_sepc_cancer$se_name[which(se_sepc_cancer$specificity=="active")])
se_loss <- unique(se_sepc_cancer$se_name[which(se_sepc_cancer$specificity=="inactive")])
se_none <- unique(se_sepc_cancer$se_name[which(se_sepc_cancer$specificity=="none")])
se_both <- intersect(se_gain,se_loss)


plot_d <- data.frame(group=c("both","active","inactive","none"),
                     value=c(length(se_both),
                             length(se_gain)-length(se_both),
                             length(se_loss)-length(se_both),
                             length(unique(se_all$se_name))-length(unique(se_sepc_cancer$se_name))))
plot_d$group <- factor(plot_d$group,levels=c("both","active","inactive","none"))
ggplot(plot_d,aes(x=group,y=value,fill=group))+
  geom_bar(stat="identity",color="black")+
  theme_classic()+
  scale_fill_manual(values = c("gold3",
                               "firebrick",
                               "lightskyblue",
                               "grey50"))

write.table(plot_d,"figs_df/supplementary_data/fig4f.txt",sep="\t",row.names = F,quote = F)


#-----------------------------------
# supplementary
#-----------------------------------
#------------
# fig_s1

se_width <- read.table("results/se_test_summary_0.01.txt",sep="\t",header=T)
ce_width <- read.table("results/ce_width_summary_0.01.txt",sep="\t",header=T)

plot_df <- se_width[which(se_width$group %in% seq(0,0.26,0.01)),]

plot_df <- ce_width[which(ce_width$group %in% seq(0.09,0.31,0.01)),]

ggplot(se_width,aes(x=group,y=log2(median)))+
  geom_line()+
  #geom_smooth(span=0.9)+
  theme_classic()

ggplot(se_width,aes(x=group,y=log2(n)))+
  geom_line()+
  #geom_smooth(span=0.9)+
  theme_classic()

ggplot(ce_width,aes(x=group,y=log2(median)))+
  geom_line()+
  #geom_smooth(span=0.9)+
  theme_classic()

ggplot(ce_width,aes(x=group,y=log2(n)))+
  geom_line()+
  #geom_smooth(span=0.9)+
  theme_classic()


#------------
# fig_s2
mix_per_cell <- read.table("figs_df/fig_s1_mix_piror.txt",sep="\t",header=T)

plot_df <- gather(mix_per_cell,group,mean,m_1:m_2)
ggplot(plot_df,aes(x=group,y=mean,fill=group))+
  geom_violin()+
  geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_classic()

write.table(plot_df,"figs_df/supplementary_data/figS2.txt",sep="\t",row.names = F,quote = F)


#------------
# fig_s4

cSEAdb <- readRDS("results/cSEAdb_r1.rds")
se_df <- cSEAdb$se_bed

se_gr <- GRanges(seqnames = se_df$chr,IRanges(as.integer(se_df$start),as.integer(se_df$end)))

meta_data <- read.table("results/cell_meta.txt",sep="\t",header=T)
meta_data$Cell_line[which(meta_data=="786-0")] <- "X786-0"
meta_data$Cancer <- gsub("-|,| ",".",meta_data$Cancer)

cancer_n <- unique(meta_data$Cancer)

se_cancer <- data.frame()
for (i in 1:length(cancer_n)) {
  cell_tmp <- meta_data$Cell_line[which(meta_data$Cancer==cancer_n[i])]
  cell_tmp[which(cell_tmp=="X786-0")] <- "786-0"

  s_cancer <- data.frame()
  for (j in 1:length(cell_tmp)) {
    s1 <- read.table(paste0("rose_no_promotor/",cell_tmp[j],"_rep1_no_promotor_peak_Gateway_SuperEnhancers.bed"),sep="\t",header=F)
    s2 <- read.table(paste0("rose_no_promotor/",cell_tmp[j],"_rep2_no_promotor_peak_Gateway_SuperEnhancers.bed"),sep="\t",header=F)
    s_all <- rbind(s1,s2)

    print(paste0(cell_tmp[j],":",nrow(s_all)))
    s_cancer <- rbind(s_cancer,s_all)
  }

  s_tmp_gr <- GRanges(seqnames = s_cancer$V1,IRanges(s_cancer$V2,s_cancer$V3))

  ov_lap <- findOverlaps(se_gr,s_tmp_gr)

  se_cancer_tmp <- se_df[unique(queryHits(ov_lap)),]
  se_cancer_tmp$cancer <- cancer_n[i]
  se_cancer <- rbind(se_cancer,se_cancer_tmp)
}

write.table(se_cancer,"figs_df/fig_s2_cancer_se_count.txt",sep="\t",quote=F,row.names = F)

se_cancer_count <- read.table("figs_df/fig_s2_cancer_se_count.txt",sep="\t",header=T)

se_cout_freq <- as.data.frame(table(se_cancer_count$cancer))

se_cout_freq <- se_cout_freq[order(-se_cout_freq$Freq),]
se_cout_freq$Var1 <- factor(se_cout_freq$Var1,levels=se_cout_freq$Var1)
ggplot(se_cout_freq,aes(x=Var1,y=Freq))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

write.table(se_cout_freq,"figs_df/supplementary_data/figS4.txt",sep="\t",row.names = F,quote = F)


#------------
# fig_s6
mixture_model <- read.table("results/filtered_ce_mixture_model_r1.txt",sep="\t",header=T)
ce_peak_w_percent <- read.table("results/s1_ce_binary_w_percent.txt",sep="\t",header=T)

model_plot_df <- data.frame(ce_name=mixture_model$ce_name,
                            number_of_1 = apply(mixture_model[,c(1:60)], 1, function(x) sum(x==1)))

model_plot_df$group <- "Mixture model"
binary_plot_df <- data.frame(ce_name=ce_peak_w_percent$ce_name[which(ce_peak_w_percent$max_percent>3)],
                             number_of_1 = apply(ce_peak_w_percent[which(ce_peak_w_percent$max_percent>3),c(4:63)],
                                                 1, function(x) sum(x==1)))
binary_plot_df$group <- "Binary"

final_plot_df <- rbind(model_plot_df,binary_plot_df)
ggplot(final_plot_df, aes(x=number_of_1,fill=group,color=group)) +
  geom_density(alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c(brewer.pal(8,"Blues")[8],
                             brewer.pal(8,"Paired")[8]))+
  scale_color_manual(values=c(brewer.pal(8,"Blues")[8],
                              brewer.pal(8,"Paired")[8]))+
  xlab("Number of cell CE present")

write.table(final_plot_df,"figs_df/supplementary_data/figS6.txt",sep="\t",row.names = F,quote = F)


#------------
# fig_s9
plot_df_erna <- read.table("figs_df/fig4c_eRNA_r1.txt",sep="\t",header=T)

# pvalue
# calculate p-value
cell_line <- c("K562_gro_cap_encode","K562_GROcap_GEO","K562_GROseq_GEO","K562_pro_cap_encode")
plot_df_2 <- plot_df_erna[which(plot_df_erna$cell %in% cell_line),]

plot_df_2$cell <- factor(plot_df_2$cell,levels=cell_line)
ggplot(plot_df_2,aes(x=cell,y=log2(signal+1)/width*1000,fill=group))+
  geom_boxplot(width=0.5,position=position_dodge(0.6),outlier.shape = NA)+
  # geom_text(aes(label=paste0(rate_ce,"%")),
  #           vjust=-0.5,position = position_dodge(0.9),size=3)+
  theme_classic()+
  ylim(c(0,6))+
  scale_fill_manual(values=c("firebrick","lightskyblue","grey50"))


plot_df_2$signal <- log2(plot_df_2$signal+1)/plot_df_2$width*1000

write.table(plot_df_2,"figs_df/supplementary_data/figS9.txt",sep="\t",row.names = F,quote = F)


#------------
# fig_s10
plot_df <- read.table("figs_df/fig4b_reporter_assay_r1.txt",sep="\t",header=T)
plot_df_erna <- read.table("figs_df/fig4c_eRNA_r1.txt",sep="\t",header=T)

# pvalue
# calculate p-value
cell_line <- c("A549_pro_seq","K562_PROseq_GEO")
plot_df_2 <- plot_df_erna[which(plot_df_erna$cell %in% cell_line),]

cell_line <- c("A549","MCF7","HCT-116","K-562")
plot_df_2 <- plot_df


#plot_df_2$cell <- factor(plot_df_2$cell,levels=cell_line)
ggplot(plot_df_2,aes(x=cell,y=log2(signal+1),fill=group))+
  geom_boxplot(width=0.5,position=position_dodge(0.6),outlier.shape = NA)+
  # geom_text(aes(label=paste0(rate_ce,"%")),
  #           vjust=-0.5,position = position_dodge(0.9),size=3)+
  theme_classic()+
  ylim(c(0,5))+
  scale_fill_manual(values=c("firebrick","lightskyblue","grey50"))

plot_df_2$signal <- log2(plot_df_2$signal+1)

write.table(plot_df_2,"figs_df/supplementary_data/figS10_assay.txt",sep="\t",row.names = F,quote = F)



#------------
# paper number summary
#------------
cSEAdb <- readRDS("results/cSEAdb_r1.rds")
se_bed <- cSEAdb$se_bed
ce_bed <- cSEAdb$ce_bed
mix_mod <- cSEAdb$ce_mixtrue_model_cell
se_spec <- cSEAdb$se_specificity

ce_binary <- read.table("results/s1_ce_binary_w_percent.txt",sep="\t",header=T)

ce_binary <- ce_binary[which(ce_binary$max_percent>3),]
ce_binary$rowmean <- apply(ce_binary[,c(4:63)],1,sum)
length(unique(ce_binary$ce_name[which(ce_binary$rowmean==60)]))/length(unique(ce_binary$ce_name))

length(unique(ce_binary$se_name[which(ce_binary$rowmean==60)]))/length(unique(ce_binary$se_name))

se_cancer_spec <- se_spec[which(se_spec$object_type=="cancer"),]

length(unique(se_cancer_spec$se_name[which(se_cancer_spec$number_spec_object>=2)]))

# se cancer cell line
se_cancer_count <- read.table("figs_df/fig_s2_cancer_se_count.txt",sep="\t",header=T)

# count how many cancer each se have
se_tmp <- table(se_cancer_count$se_name)


length(se_tmp[which(se_tmp==28)])




