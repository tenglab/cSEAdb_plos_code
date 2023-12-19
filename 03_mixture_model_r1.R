setwd("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision")
library(GenomicRanges)
library(ggplot2)
library(tidyr)
library(dplyr)
library(data.table)
library(DESeq2)
library(RColorBrewer)
library(Rsubread)
library(mixtools)
library(mclust)
mem.maxVSize(vsize = Inf)

#chr1_66443354_66447012


#--------------------------------------
# all enhancer mixture prior
#--------------------------------------
se_ce_tmp <- read.table("results/ce_signal_normalized_60cell_filtered.txt",sep="\t",header=T)
# get prior distribution
se_ce_long<- gather(se_ce_tmp,cell,signal,X786.0:UO.31)
se_ce_long <- se_ce_long[se_ce_long$signal!=0,]
se_ce_long$log_signal <- log2(se_ce_long$signal)

prior_2 <- normalmixEM(se_ce_long$log_signal,k=2)
prior_df <- data.frame(u1=prior_2$mu[1],
                       u2=prior_2$mu[2],
                       sd1=prior_2$sigma[1],
                       sd2=prior_2$sigma[2],
                       p1=prior_2$lambda[1],
                       p2=prior_2$lambda[2])

write.table(prior_df,"results/all_ce_prior_r1.txt",sep="\t",quote = F,row.names = F)

#------------------
# all CE, from 0 mixture modeling
#------------------
# prior
prior_df <- read.table("results/all_ce_prior_r1.txt",sep="\t",header=T)
#prior_df2 <- read.table("results/all_ce_prior_r2.txt",sep="\t",header=T)

u1 <- prior_df$u1
u2 <- prior_df$u2
sd1 <- prior_df$sd1
sd2 <- prior_df$sd2
p <- prior_df$p2

# mixture model
se_ce_tmp <- read.table("results/ce_signal_normalized_60cell_filtered.txt",sep="\t",header=T)

# 0 imputation with half of min value
se_ce <- se_ce_tmp %>% mutate_at(c(2:61),~replace(., . == 0, min(.[.>0], na.rm = TRUE)/2))
#se_ce <- se_ce_tmp %>% mutate_at(c(2:61),~replace(., . == 0, rnorm(47594,mean=u1,sd=sd1)))


#A549 example: chr1_66443354_66447012



se_ce[,c(2:61)] <- log2(se_ce[,c(2:61)])
se_name <- unique(se_ce$se_name)

source("~/Projects/super_enhancer/se_data_portal/new_analysis/cSEAdb_revision/scripts/normalmixEM2comp_r1.R")

#se=9
#ce=4
matrix_cat_final <- data.frame()
for (se in 1:length(se_name)) {
  if(se %% 500 ==0 ){
    print(se)
  }
  se_tmp <- se_ce[which(se_ce$se_name==se_name[se]),]
  se_raw <- se_ce_tmp[which(se_ce_tmp$se_name==se_name[se]),]
  rownames(se_tmp) <- se_tmp$merge_e_name
  matrix_tmp <- t(as.matrix(se_tmp[,-c(1,62)]))
  matrix_out <- matrix_tmp
  ce <- which(se_tmp$merge_e_name=="chr1_66443354_66447012")

  # check each CE within SE
  for (ce in 1:ncol(matrix_tmp)) {
    ce_tmp <- matrix_tmp[,ce]
    ce_raw <- se_raw[ce,2:61]
    ce_0_name <- colnames(ce_raw)[which(ce_raw==0)]

    #----------------------
    # check all small
    if (all(ce_tmp <= u1+sd1) ) {
      matrix_out[,ce] <- 0
    } else if (all(ce_tmp >= u2-sd2)) {
      matrix_out[,ce] <- 1
    } else {
      #----------------------
      # use mixture model
      # remove outlier >u2+sd2
      ce_1 <- ce_tmp[ce_tmp >= u2+sd2]
      ce_1_name <- names(ce_1)

      ce_mod <- ce_tmp[ce_tmp<u2+sd2]
      ce_mod_name <- names(ce_mod)

      mix_mod <- normalmixEM2comp_r1(ce_mod,
                                     lambda = c(prior_df$p1,prior_df$p2),
                                     mu=c(prior_df$u1,prior_df$u2),
                                     sigsqrd=c(prior_df$sd1,prior_df$sd2))

      #plot(ce_mod,mix_mod$posterior[,2],pch=16)
      #------------
      # find cutoff:which min_up0.5_signal greater than max_low0.5_signal,
      # less than cutoff is 0, greater is 1

      # check which mu is greater
      if (mix_mod$mu[1]<=mix_mod$mu[2]) {
        post_up05_name <- names(which(mix_mod$posterior[,2]>0.5))
        post_low05_name <- names(which(mix_mod$posterior[,2]<0.5))
      } else {
        post_up05_name <- names(which(mix_mod$posterior[,1]>0.5))
        post_low05_name <- names(which(mix_mod$posterior[,1]<0.5))
      }

      #-----------
      # min max signal of two groups
      signal_up05_min <- min(ce_mod[names(ce_mod) %in% post_up05_name])
      signal_up05_max <- max(ce_mod[names(ce_mod) %in% post_up05_name])
      signal_low05_min <- min(ce_mod[names(ce_mod) %in% post_low05_name])
      signal_low05_max <- max(ce_mod[names(ce_mod) %in% post_low05_name])

      #----------
      # check shape and assign groups
      # linear
      if (signal_low05_min<=signal_up05_min & signal_low05_max<=signal_up05_max) {
        # cutoff = signal_low05_max
        idx_0 <- which(names(matrix_out[,ce]) %in% names(which(ce_mod<=signal_low05_max)))
        idx_1 <- which(names(matrix_out[,ce]) %in% names(which(ce_mod>signal_low05_max)))
        matrix_out[idx_0,ce] <- 0
        matrix_out[idx_1,ce] <- 1
      } else if (signal_low05_min < signal_up05_min & signal_low05_max > signal_up05_max) {
        # bump above 0.5
        # cutoff= signal_up05_min
        idx_0 <- which(names(matrix_out[,ce]) %in% names(which(ce_mod<signal_up05_min)))
        idx_1 <- which(names(matrix_out[,ce]) %in% names(which(ce_mod>=signal_up05_min)))
        matrix_out[idx_0,ce] <- 0
        matrix_out[idx_1,ce] <- 1
      } else {
        # bump below 0.5,
        # cutoff=signal_low05_max
        idx_0 <- which(names(matrix_out[,ce]) %in% names(which(ce_mod<=signal_low05_max)))
        idx_1 <- which(names(matrix_out[,ce]) %in% names(which(ce_mod>signal_low05_max)))
        matrix_out[idx_0,ce] <- 0
        matrix_out[idx_1,ce] <- 1
      }

      # set all 0 signal with 0  and outlier with 1
      matrix_out[(names(matrix_out[,ce]) %in% ce_0_name),ce] <-0
      matrix_out[(names(matrix_out[,ce]) %in% ce_1_name),ce] <- 1
    }
  }

  #cbind(mix_mod$posterior[,2],ce_mod)
  #cbind(t(ce_raw),ce_tmp,matrix_out[,ce])
  matrix_out <- data.frame(t(matrix_out))

  # add se_name
  matrix_out$ce_name <- row.names(matrix_out)
  matrix_out$se_name <- se_name[se]

  # combine all
  matrix_cat_final <- rbind(matrix_cat_final,matrix_out)
}

write.table(matrix_cat_final,"results/ce_mixture_model_r1.txt",quote = F,row.names = F,sep="\t")

#------------------------------------------------------------------
# CE rowMean distribution after remove < 3% weight
#------------------------------------------------------------------
matrix_cat_final <- read.table("results/ce_mixture_model_r1.txt",sep="\t",header=T)
ce_width <- read.table("results/ce_maxcover_60cell.txt",sep="\t",header=T)
ce_keep <- unique(ce_width$ce_name[which(ce_width$max_percent > 3)])
mixture_keep <- matrix_cat_final[which(matrix_cat_final$ce_name %in% ce_keep),]
mixture_keep$rowmean <- rowMeans(mixture_keep[,-c(61,62)])

write.table(mixture_keep,"results/filtered_ce_mixture_model_r1.txt",quote = F,row.names = F,sep="\t")

mixture_keep <- read.table("results/filtered_ce_mixture_model_r1.txt",sep="\t",header=T)
se_ce_tmp <- read.table("results/ce_signal_normalized_60cell_filtered.txt",sep="\t",header=T)

length(unique(mixture_keep$ce_name[which(mixture_keep$rowmean==0)]))
length(unique(mixture_keep$ce_name[which(mixture_keep$rowmean==1)]))

mod_tmp <- unique(mixture_keep[,-c(62,63)])
mixture_long <- gather(mod_tmp, cell, mod, X786.0:UO.31, factor_key=TRUE)

# add signal
se_ce_tmp_2_long <- gather(unique(se_ce_tmp[,-62]), cell, signal, X786.0:UO.31, factor_key=TRUE)
colnames(se_ce_tmp_2_long)[1] <-"ce_name"
compare_df <- merge(mixture_long,se_ce_tmp_2_long,by=c("cell","ce_name"))

summary(compare_df$signal[which(compare_df$mod==0)])
summary(compare_df$signal[which(compare_df$mod==1)])

compare_df[which(compare_df$mod==0&compare_df$signal>1000),]

#------------------------------------------------------------------
# comapre model and binary df
#------------------------------------------------------------------

mod2 <- read.table("results/filtered_ce_mixture_model_r1.txt",sep="\t",header=T)
binary_ce <- read.table("results/s1_ce_binary_w_percent.txt",sep="\t",header=T)


# peak calling--------
binary_selected <- binary_ce[which(binary_ce$max_percent > 3),-2]

# remove duplicate CE
binary_selected <- unique(binary_selected[,-2])

binary_long <- gather(binary_selected, cell, binary, X786.0:UO.31, factor_key=TRUE)

# model--------
# remove duplicate CE
mod_tmp <- unique(mod2[,-c(62,63)])
mod_long <- gather(mod_tmp, cell, mod, X786.0:UO.31, factor_key=TRUE)

# diff data frame
diff_long <- merge(binary_long,mod_long,by=c("ce_name","cell"))

diff_long$bi_minor_mod <- diff_long$binary-diff_long$mod

# add signal
# signal difference
se_ce_tmp <- read.table("results/ce_signal_normalized_60cell_filtered.txt",sep="\t",header=T)
se_ce <- unique(se_ce_tmp[,-c(62)])

se_ce_long <- gather(se_ce, cell, signal, X786.0:UO.31, factor_key=TRUE)
colnames(se_ce_long)[1] <- "ce_name"

# merge diff and signal
diff_long_signal <- merge(diff_long,se_ce_long,by=c("ce_name","cell"))

# boxplot

diff_long_signal$group <- NA
diff_long_signal$group[which(diff_long_signal$bi_minor_mod==0)] <- "common"
diff_long_signal$group[which(diff_long_signal$bi_minor_mod==-1)] <- "model only"
diff_long_signal$group[which(diff_long_signal$bi_minor_mod==1)] <- "peak only"

write.table(diff_long_signal,"results/binary_model_diff_r1.txt",sep="\t",quote=F,row.names = F)

