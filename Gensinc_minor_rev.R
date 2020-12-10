#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se


#eLife minor revision - recalibrating P-values and re-plotting
#Figure 3
library(tidyverse)
library(ggplot2)
BA_AA_data <- process_chi("~/Projects/Herring/data/British_populations/gensinc_revision/BalticAutumn_AtlanticAutumn.chi.out", lambda = BAut_v_AAut_chi$reg$coefficients[1])
BS_AS_data <- process_chi("~/Projects/Herring/data/British_populations/gensinc_revision/BalticSpring_AtlanticSpring.chi.out", lambda = BSpr_v_ASpr_chi$reg$coefficients[1])

BvA_join <- BS_AS_data %>% full_join(BA_AA_data, by=c("CHR", "POS"))
rm(BS_AS_data, BA_AA_data)

BvA_join <- BvA_join[order(BvA_join$CHR, BvA_join$POS),]
rownames(BvA_join) <- paste(BvA_join$CHR, BvA_join$POS, sep = "_")
names(BvA_join) <- c("CHR", "POS",  "log10P_spr", "index_spr", "adj_P_spr", "log10P_aut", "index_aut", "adj_P_aut")
BvA_join$chr_col <- c("grey30", "grey70")[(as.integer(BvA_join$CHR) %% 2)+1]

#Tagging SNPs with -logP > 10 in both contrasts
p_cutoff <- 10
BvA_join$dual_signal <- F
BvA_join$dual_signal[which(BvA_join$adj_P_spr > p_cutoff  & BvA_join$adj_P_aut > p_cutoff)] <- T
BvA_join[,"cumulative_pos"] <- BvA_join$POS + Ch_v2.0.2_sizes[match(BvA_join$CHR, Ch_v2.0.2_sizes$name),"offset"]

png("~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/Baltic_Atlantic.manhattan.png", width=1500, height = 1000)
plot(x= 1, y = 1, type = "n", xlim = c(0, max(BvA_join$cumulative_pos, na.rm = T)), ylim = c(-max(BvA_join$adj_P_aut, na.rm = T), max(BvA_join$adj_P_spr, na.rm = T)), xlab = "", ylab = "", axes = F)
lines(x = c(-1e6, max(BvA_join$cumulative_pos, na.rm = T)+1e6), y = c(0,0), lwd = 2)
points(x = BvA_join$cumulative_pos[!BvA_join$dual_signal], y = 1 +BvA_join$adj_P_spr[!BvA_join$dual_signal], col = BvA_join$chr_col[!BvA_join$dual_signal], pch = 16, cex = 0.5)
points(x = BvA_join$cumulative_pos[!BvA_join$dual_signal], y = -1 -BvA_join$adj_P_aut[!BvA_join$dual_signal], col = BvA_join$chr_col[!BvA_join$dual_signal], pch = 16, cex = 0.5)
points(x = BvA_join$cumulative_pos[BvA_join$dual_signal], y = 1 +BvA_join$adj_P_spr[BvA_join$dual_signal], col = "red", pch = 16, cex = 0.9)
points(x = BvA_join$cumulative_pos[BvA_join$dual_signal], y = -1 -BvA_join$adj_P_aut[BvA_join$dual_signal], col = "red", pch = 16, cex = 0.9)
axis(2, pos =-1e6, at = seq(from = 1, to = max(BvA_join$adj_P_spr, na.rm = T), by = 10), labels = seq(from = 1, to = max(BvA_join$adj_P_spr, na.rm = T), by = 10)-1, lwd = 2, lwd.ticks = 2)
axis(2, pos =-1e6, at = -seq(from = 1, to = max(BvA_join$adj_P_aut, na.rm = T), by = 10), labels = -seq(from = 1, to = max(BvA_join$adj_P_aut, na.rm = T), by = 10)+1, lwd = 2, lwd.ticks = 2)

dev.off()

BA_BS_data <- process_chi("~/Projects/Herring/data/British_populations/gensinc_revision/SpringBaltic_AutumnBaltic.chi.out", lambda = BAut_v_BSpr_chi$reg$coefficients[1])
AA_AS_data <- process_chi("~/Projects/Herring/data/British_populations/gensinc_revision/SpringAtlantic_AutumnAtlantic.chi.out", lambda = AAut_v_ASpr_chi$reg$coefficients[1])

AvS_join <- BA_BS_data %>% full_join(AA_AS_data, by=c("CHR", "POS"))
rm(BA_BS_data, AA_AS_data)

AvS_join <- AvS_join[order(AvS_join$CHR, AvS_join$POS),]
rownames(AvS_join) <- paste(AvS_join$CHR, AvS_join$POS, sep = "_")
names(AvS_join) <- c("CHR", "POS",  "log10P_Balt", "index_Balt", "adj_P_Balt", "log10P_Atl", "index_Atl", "adj_P_Atl")
AvS_join$chr_col <- c("grey30", "grey70")[(as.integer(AvS_join$CHR) %% 2)+1]
p_cutoff <- 10
AvS_join$dual_signal <- F
AvS_join$dual_signal[which(AvS_join$adj_P_Balt > p_cutoff  & AvS_join$adj_P_Atl > p_cutoff)] <- T
AvS_join[,"cumulative_pos"] <- AvS_join$POS + Ch_v2.0.2_sizes[match(AvS_join$CHR, Ch_v2.0.2_sizes$name),"offset"]

save(AvS_join, BvA_join, file = "~/Projects/Herring/data/British_populations/gensinc_revision/joint_chi_dfs.RData")

NvS_data <- process_chi("~/Projects/Herring/data/British_populations/gensinc_revision/South_North.chi.out", lambda = S_v_N_chi$reg$coefficients[1])
NvS_data[,"cumulative_pos"] <- NvS_data$POS + Ch_v2.0.2_sizes[match(NvS_data$CHR, Ch_v2.0.2_sizes$name),"offset"]
NvS_data$chr_col <- c("grey30", "grey70")[(as.integer(NvS_data$CHR) %% 2)+1]

NSvNEatl_data <- process_chi("~/Projects/Herring/data/British_populations/gensinc_revision/NorthSea_otherAtlantic.chi.out", lambda = NS_v_NEatl_chi$reg$coefficients[1])
NSvNEatl_data[,"cumulative_pos"] <- NSvNEatl_data$POS + Ch_v2.0.2_sizes[match(NSvNEatl_data$CHR, Ch_v2.0.2_sizes$name),"offset"]
NSvNEatl_data$chr_col <- c("grey30", "grey70")[(as.integer(NSvNEatl_data$CHR) %% 2)+1]


#3a BvA
joint_manhattan(joint_df = BvA_join, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/Baltic_Atlantic_manhattan.png", up_col = which(names(BvA_join) == "adj_P_spr"), down_col = which(names(BvA_join) == "adj_P_aut"), bg_cex = 1, hl_cex = 1.2)
#3b chr12 BvA
joint_manhattan(joint_df = BvA_join, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/Baltic_Atlantic_manhattan.png", up_col = which(names(BvA_join) == "adj_P_spr"), down_col = which(names(BvA_join) == "adj_P_aut"), target_chr = "chr12", bg_cex = 1, hl_cex = 1.2)
#3c AvS
joint_manhattan(joint_df = AvS_join, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/Spring_Autumn_manhattan.png", up_col = which(names(AvS_join) == "adj_P_Balt"), down_col = which(names(AvS_join) == "adj_P_Atl"), bg_cex = 1, hl_cex = 1.2)
#3d chr15 AvS 7-11.5 Mb
joint_manhattan(joint_df = AvS_join, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/Spring_Autumn_manhattan.png", up_col = which(names(AvS_join) == "adj_P_Balt"), down_col = which(names(AvS_join) == "adj_P_Atl"), target_chr = "chr15", reg_lim = c(7e6, 11.5e6), bg_cex = 1, hl_cex = 1.2)


#4a NSvNEatl
single_manhattan(chi_df = NSvNEatl_data, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/NS_v_NEatl_manhattan.png", up_col = "adj_P", bg_cex = 1)
#4b chr6 NSvNEatl
single_manhattan(chi_df = NSvNEatl_data, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/NS_v_NEatl_manhattan.png", up_col = "adj_P", target_chr = "chr6", bg_cex = 1)
#4c chr 12 NSvNEatl
single_manhattan(chi_df = NSvNEatl_data, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/NS_v_NEatl_manhattan.png", up_col = "adj_P", target_chr = "chr12", bg_cex = 1)
#4d chr 17 NSvNEatl
single_manhattan(chi_df = NSvNEatl_data, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/NS_v_NEatl_manhattan.png", up_col = "adj_P", target_chr = "chr17", bg_cex = 1)
#4e chr23 NSvNEatl
single_manhattan(chi_df = NSvNEatl_data, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/NS_v_NEatl_manhattan.png", up_col = "adj_P", target_chr = "chr23", bg_cex = 1)


#5a NvS
single_manhattan(chi_df = NvS_data, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/N_v_S_manhattan.png", up_col = "adj_P", bg_cex = 1)
#5b chr2 NvS
single_manhattan(chi_df = NvS_data, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/N_v_S_manhattan.png", up_col = "adj_P", target_chr = "chr2", bg_cex = 1)
#5c chr2 AHR region NvS
anno <- read.table("~/Projects/Herring/data/British_populations/gensinc_revision/chr2.log20.missense.list", header=F, sep="\t")
colnames(anno) <- c("CHR","POS", "log10P", "AA")
single_manhattan(chi_df = NvS_data, png_file = "~/Projects/Herring/doc/Gensinc/Revision/minor_rev_eLife/N_v_S_manhattan_AHR.png", up_col = "adj_P", target_chr = "chr2", reg_lim = c(15.8e6, 15.82e6), bg_cex = 2, hl_df = anno)

#Correcting values in the text
#LRRC8C
139/(BSpr_v_ASpr_chi$reg$coefficients)
exp_q[1:floor(length(exp_q) * 0.9)] 
72.28328 
62/(BAut_v_AAut_chi$reg$coefficients)
exp_q[1:floor(length(exp_q) * 0.9)]
37.62795
#TSHR
106/(BAut_v_BSpr_chi$reg$coefficients)
exp_q[1:floor(length(exp_q) * 0.9)] 
68.24975 
82/(AAut_v_ASpr_chi$reg$coefficients)
exp_q[1:floor(length(exp_q) * 0.9)] 
53.68841
#SOX11B
101/(BAut_v_BSpr_chi$reg$coefficients)
exp_q[1:floor(length(exp_q) * 0.9)] 
65.03042 
88/(AAut_v_ASpr_chi$reg$coefficients)
exp_q[1:floor(length(exp_q) * 0.9)] 
57.61683 
#AHR2B2
87/(S_v_N_chi$reg$coefficients)
exp_q[1:floor(length(exp_q) * 0.9)] 
45.68292 
#SLC12A2
47/(S_v_N_chi$reg$coefficients)
exp_q[1:floor(length(exp_q) * 0.9)] 
24.67928 

#Support functions
process_chi <- function(chi_file, lambda){
  require(tidyverse)
  require(ggplot2)
  chi_data <- read.table(chi_file, header=T, sep="\t", skipNul = TRUE)
  chi_data <- chi_data %>% filter(!grepl("unplaced", CHR))
  chr_order <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23","chr24","chr25","chr26")
  chi_data$CHR <- factor(chi_data$CHR, levels = chr_order)
  chi_data <- chi_data[order(chi_data$CHR, chi_data$POS),]
  #chi_data$index <- 1:nrow(chi_data)
  chi_data$index <- paste(chi_data$CHR, chi_data$POS, sep = "_")
  chi_data$adj_P <-  chi_data$log10P/lambda
  return(chi_data)
}

joint_manhattan <- function(joint_df, png_file, up_col, down_col, target_chr = NULL, reg_lim = NULL, hl_cex = 0.9, bg_cex = 0.5){
  if(is.null(target_chr)){
    png(file = png_file, width=1500, height = 1000)
    plot(x= 1, y = 1, type = "n", xlim = c(0, max(joint_df$cumulative_pos, na.rm = T)), ylim = c(-max(joint_df[,down_col], na.rm = T), max(joint_df[,up_col], na.rm = T)), xlab = "", ylab = "", axes = F)
    lines(x = c(-1e6, max(joint_df$cumulative_pos, na.rm = T)+1e6), y = c(0,0), lwd = 2)
    points(x = joint_df$cumulative_pos[!joint_df$dual_signal], y = 1 +joint_df[!joint_df$dual_signal,up_col], col = joint_df$chr_col[!joint_df$dual_signal], pch = 16, cex = bg_cex)
    points(x = joint_df$cumulative_pos[!joint_df$dual_signal], y = -1 -joint_df[!joint_df$dual_signal,down_col], col = joint_df$chr_col[!joint_df$dual_signal], pch = 16, cex = bg_cex)
    points(x = joint_df$cumulative_pos[joint_df$dual_signal], y = 1 +joint_df[joint_df$dual_signal, up_col], col = "red", pch = 16, cex = hl_cex)
    points(x = joint_df$cumulative_pos[joint_df$dual_signal], y = -1 -joint_df[joint_df$dual_signal,down_col], col = "red", pch = 16, cex = hl_cex)
    axis(2, pos =-1e6, at = seq(from = 0, to = max(joint_df[,up_col], na.rm = T), by = 10)+1, labels = seq(from = 0, to = max(joint_df[,up_col], na.rm = T), by = 10), lwd = 3, lwd.ticks = 3)
    axis(2, pos =-1e6, at = -seq(from = 0, to = max(joint_df[,down_col], na.rm = T), by = 10)-1, labels = seq(from = 0, to = max(joint_df[,down_col], na.rm = T), by = 10), lwd = 3, lwd.ticks = 3)
    dev.off()
  }
  if(!is.null(target_chr)){
    y_max <- max(joint_df[,up_col], na.rm = T)
    y_min <- -max(joint_df[,down_col], na.rm = T)
    joint_df <- joint_df[joint_df$CHR == target_chr, ]
    if(is.null(reg_lim)) reg_lim <- c(0, max(joint_df$POS, na.rm = T))
    png(file = paste0(sub("[.]png", "", png_file), "_", target_chr, ".png"), width=1500, height = 1000)
    plot(x= 1, y = 1, type = "n", xlim = reg_lim , ylim = c(y_min, y_max), xlab = "", ylab = "", axes = F)
    #lines(x = c(-1e6, max(joint_df$POS, na.rm = T)+1e6), y = c(0,0), lwd = 2)
    axis(1, pos = 0)
    points(x = joint_df$POS[!joint_df$dual_signal], y = 3 +joint_df[!joint_df$dual_signal,up_col], col = "black", pch = 16, cex = bg_cex)
    points(x = joint_df$POS[!joint_df$dual_signal], y = -5 -joint_df[!joint_df$dual_signal,down_col], col = "black", pch = 16, cex = bg_cex)
    points(x = joint_df$POS[joint_df$dual_signal], y = 3 +joint_df[joint_df$dual_signal, up_col], col = "red", pch = 16, cex = hl_cex)
    points(x = joint_df$POS[joint_df$dual_signal], y = -5 -joint_df[joint_df$dual_signal,down_col], col = "red", pch = 16, cex = hl_cex)
    axis(2, pos = reg_lim[1] - 0.05* diff(reg_lim), at = seq(from = 0, to = y_max, by = 10)+3, labels = seq(from = 0, to = y_max, by = 10), lwd = 3, lwd.ticks = 3)
    axis(2, pos = reg_lim[1] - 0.05* diff(reg_lim), at = -seq(from = 0, to = -y_min, by = 10)-5, labels = seq(from = 0, to = -y_min, by = 10), lwd = 3, lwd.ticks = 3)
    dev.off()
  }
}

single_manhattan <- function(chi_df, png_file, up_col, target_chr = NULL, reg_lim = NULL, bg_cex = 0.5, hl_df = NULL){
  if(is.null(target_chr)){
    y_max <- max(chi_df[,up_col], na.rm = T)
    png(file = png_file, width=1500, height = 1000)
    plot(x= 1, y = 1, type = "n", xlim = c(0, max(chi_df$cumulative_pos, na.rm = T)), ylim = c(0, y_max), xlab = "", ylab = "", axes = F)
    #lines(x = c(-1e6, max(chi_df$cumulative_pos, na.rm = T)+1e6), y = c(0,0), lwd = 2)
    points(x = chi_df$cumulative_pos, y = chi_df[,up_col], col = chi_df$chr_col, pch = 16, cex = bg_cex)
    axis(2, pos =-1e6, at = seq(from = 0, to = max(chi_df[,up_col], na.rm = T), by = 10), labels = seq(from = 0, to = max(chi_df[,up_col], na.rm = T), by = 10), lwd = 2, lwd.ticks = 2)
    dev.off()
  }
  if(!is.null(target_chr)){
    y_max <- max(chi_df[,up_col], na.rm = T)
    chi_df <- chi_df[chi_df$CHR == target_chr, ]
    if(is.null(reg_lim)) reg_lim <- c(0, max(chi_df$POS, na.rm = T))
    png(file = paste0(sub("[.]png", "", png_file), "_", target_chr, ".png"), width=1500, height = 1000)
    plot(x= 1, y = 1, type = "n", xlim = reg_lim, ylim = c(0, y_max), xlab = "", ylab = "", axes = F)
    #lines(x = c(-1e6, max(chi_df$cumulative_pos, na.rm = T)+1e6), y = c(0,0), lwd = 2)
    points(x = chi_df$POS, y = chi_df[,up_col], col = "black", pch = 16, cex = bg_cex)
    if(!is.null(hl_df)){
      points(x = chi_df$POS[chi_df$POS %in% hl_df$POS], y = chi_df[chi_df$POS %in% hl_df$POS,up_col], col = "red", pch = 16, cex = bg_cex*1.2)
    }
    axis(2, pos = reg_lim[1] - 0.05* diff(reg_lim), at = seq(from = 0, to = y_max, by = 10), labels = seq(from = 0, to = y_max, by = 10), lwd = 3, lwd.ticks = 3)
    axis(1, lwd = 3, lwd.ticks = 3)
    dev.off()
  }
}


