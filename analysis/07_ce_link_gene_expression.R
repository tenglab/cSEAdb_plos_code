setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(Rsubread)
library(grid)
library(gridExtra)

link_gene <- read.table("../figs_df/fig3c_ce_link_gene_w_promoter.txt",sep="\t",header=T)
depmap <- read.table("OmicsExpressionProteinCodingGenesTPMLogp1.csv",sep=",",header=T)

depmap_2 <- depmap
depmap_2[depmap_2==0] <- NA

depmap_test <- depmap_2[1:50,1:5]


replace_with_rankings <- function(x) {
  rank(x, na.last = "keep",ties.method="average")
}

rank(depmap_2$X,na.last = "keep",ties.method="average")[1:5]

# get gene ranking
depmap_g_ranking <- as.data.frame(apply(depmap_2[], 2, replace_with_rankings))

colnames(depmap_g_ranking) <- c("cell",gsub("\\...*","",colnames(depmap_g_ranking)[-1]))
depmap_g_ranking[is.na(depmap_g_ranking)] <-0
depmap_g_ranking$cell <- depmap_2$X
#write.table(depmap_g_ranking,"results/depmap_gene_ranking.txt",quote=F,row.names = F,sep="\t")


depmap_g_ranking_ratio <- depmap_g_ranking
link_summary <- unique(link_gene[!(link_gene$signal==0 &link_gene$group=="common"),])

cell_map <- data.frame(dep_id=c("ACH-000681","ACH-000971","ACH-000551","ACH-000019"),
                       cell_name=c("A549","HCT.116","K.562","MCF7"))

depmap_g_ranking_ratio[,-1] <- apply(depmap_g_ranking_ratio[,-1], 2, function(x) round(x/max(x),digits = 2))

#write.table(depmap_g_ranking_ratio,"results/depmap_gene_ranking_ratio.txt",quote=F,row.names = F,sep="\t")


depmap_g_ranking_ratio <- read.table("results/depmap_gene_ranking_ratio.txt",sep="\t",header=T)
depmap_g_ranking_ratio_2 <- depmap_g_ranking_ratio[which(depmap_g_ranking_ratio$cell %in%
                                                           c("ACH-000681","ACH-000971","ACH-000551","ACH-000019")),]

depmap_rank_2 <- as.data.frame(t(depmap_g_ranking_ratio_2))
colnames(depmap_rank_2) <- c("HCT.116","MCF7","A549","K.562")
depmap_rank_2 <- depmap_rank_2[-1,]
depmap_rank_2$gene_name <- row.names(depmap_rank_2)






