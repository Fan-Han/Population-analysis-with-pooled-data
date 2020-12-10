#Author: Mats Pettersson
#mats.pettersson@imbim.uu.se

#British population paper
pool_freq_df <- read.table("~/Projects/Herring/data/British_populations/60.Neff.freq.gz", sep = "\t", header = T, stringsAsFactors = F)
#alternative, to subset before importing
#zgrep -m5 -E 'chr19|CHROM' 60.Neff.freq.gz > chr19_60_Neff.freq
pool_freq_df$SNP_id <- paste(pool_freq_df$CHROM, pool_freq_df$POS, sep = "_")

North_v_South_DAF <- read.table("~/Projects/Herring/data/diagnostic_SNP_set/South_North.deltaAF", stringsAsFactors = F, sep = "\t", header = T)
plot(x = North_v_South_DAF$POS[North_v_South_DAF$CHR == "chr19"], y = abs(North_v_South_DAF$deltaAF[North_v_South_DAF$CHR == "chr19"]), type = "p", pch = 16, cex = 0.4)
chr19_high_diff <- North_v_South_DAF[North_v_South_DAF$CHR == "chr19" & North_v_South_DAF$POS > 5e6 & abs(North_v_South_DAF$deltaAF > 0.4),]
rownames(chr19_high_diff) <- paste(chr19_high_diff$CHR, chr19_high_diff$POS, sep = "_")

chr19_high_diff_freq <- pool_freq_df[pool_freq_df$CHROM == "chr19" & pool_freq_df$POS %in% chr19_high_diff$POS, ]
rownames(chr19_high_diff_freq) <- paste(chr19_high_diff_freq$CHROM, chr19_high_diff_freq$POS, sep = "_")
pool_order_vec <- c(58, 5, 6, 7, 8, 9, 10, 61, 15, 53, 37, 54, 52, 3, 51, 14, 44, 60, 57, 56, 55, 30, 31, 20, 29, 62, 4, 11, 50, 34, 35, 32, 33, 45, 13, 49, 59, 36, 48, 28, 27, 26, 25, 46, 16, 23, 24, 19, 21, 47, 22, 17, 18, 12)

chr19_6.4_Mb_freq <- pool_freq_df[pool_freq_df$CHROM == "chr19" & pool_freq_df$POS > 6.35e6 & pool_freq_df$POS < 6.39e6,]
rownames(chr19_6.4_Mb_freq) <- paste(chr19_6.4_Mb_freq$CHROM, chr19_6.4_Mb_freq$POS, sep = "_")
pool_freq_hm(chr19_6.4_Mb_freq, pdf_file = "~/Projects/Herring/doc/British_pop/chr_19_6.4_Mb_hm_full.pdf", pool_order_vec = rev(pool_order_vec))

#random region for checking
tmp_freq <- pool_freq_df[pool_freq_df$CHROM == "chr7" & pool_freq_df$POS > 8.35e6 & pool_freq_df$POS < 8.38e6,]
pool_freq_hm(tmp_freq, pdf_file = "~/Projects/Herring/doc/British_pop/ctest_hm_full.pdf", pool_order_vec = rev(pool_order_vec))

#Using annotation data frame from Polygenic.R
load(file = "~/Projects/Herring/data/polygenic/Salinity_M_values.RData")
chr19_high_diff_rows <- Baltic_EastAtl_ann_df$SNP_id %in% rownames(chr19_high_diff)
chr19_high_diff_B_v_EA <- Baltic_EastAtl_ann_df[chr19_high_diff_rows,]
rownames(chr19_high_diff_B_v_EA) <- chr19_high_diff_B_v_EA$SNP_id

chr19_6.4_Mb_rows <- Baltic_EastAtl_ann_df$SNP_id %in% rownames(chr19_6.4_Mb_freq)
chr19_6.4_Mb_B_v_EA <- Baltic_EastAtl_ann_df[chr19_6.4_Mb_rows,]
rownames(chr19_6.4_Mb_B_v_EA) <- chr19_6.4_Mb_B_v_EA$SNP_id

hm_col <- c("#004B87", "#FFCD00")

chr19_high_diff_freq$NonSyn <- chr19_high_diff_B_v_EA$NonSyn
pool_freq_hm(chr19_high_diff_freq, pdf_file = "~/Projects/Herring/doc/British_pop/chr_19_6.4_Mb_hm_update.pdf",pool_order_vec = rev(pool_order_vec), hm_col = hm_col)

chr19_6.4_Mb_freq$NonSyn <- chr19_6.4_Mb_B_v_EA$NonSyn[match(rownames(chr19_6.4_Mb_freq), rownames(chr19_6.4_Mb_B_v_EA))]
pool_freq_hm(chr19_6.4_Mb_freq, pdf_file = "~/Projects/Herring/doc/British_pop/chr_19_6.4_Mb_hm_full.pdf", pool_order_vec = rev(pool_order_vec))

#Including HWS samples
hws_pool_order_vec  <- c(38:43, pool_order_vec)
pool_freq_hm(chr19_high_diff_freq, pdf_file = "~/Projects/Herring/doc/British_pop/chr_19_6.4_Mb_hm_with_HWS.pdf", pool_order_vec = rev(hws_pool_order_vec))

#Individual haplotypes
which(h79_w$scaffold == "chr19" & h79_w$start == 6.3e6+1)
plot_top_regions_2(5542, win_df = h79_w, geno = herring_79$geno, sample_list = herring_79$sample_list, plot_dir = "~/Projects/Herring/doc/British_pop/", snp_col="blue", snp_lwd = 0.25, snp_alpha = 120)

hm_col <- c("#004B87", "#FFCD00")

#Revision for eLife
#Inflation & QQ plots


#BA_BS_chi <- read.table("~/Projects/Herring/data/British_populations/gensinc_revision/SpringBaltic_AutumnBaltic.chi.out", header = T)
#obs_q <- BA_BS_chi$log10P[order(BA_BS_chi$log10P)]
#exp_q <- -log10((length(BA_BS_chi$log10P):1)/length(BA_BS_chi$log10P))
#lambda_reg <- lm(obs_q[1:floor(length(exp_q)*0.9)]~0+exp_q[1:floor(length(exp_q)*0.9)])
#png("~/Projects/Herring/doc/Gensinc/Revision/BAvBS_qq.png", height = 1000, width = 1000)
#plot(x = exp_q, y = obs_q, pch = 16, cex = 0.5, xlim = c(0,7), ylim = c(0,7))
#abline(a = 0, b = 1, col = "red")
#abline(lambda_reg)
#dev.off()

#png("~/Projects/Herring/doc/Gensinc/Revision/BAvBS_full_qq.png", height = 1000, width = 1000)
#plot(x = exp_q, y = obs_q, pch = 16, cex = 0.5, xlim = c(0,7))
#abline(a = 0, b = 1, col = "red")
#abline(lambda_reg)
#dev.off()

#Salinity
BAut_v_AAut_chi <- chi_qq(chi_file = "~/Projects/Herring/data/British_populations/gensinc_revision/BalticAutumn_AtlanticAutumn.chi.out", contrast_name = "BaltAutvAtlAut_recalc")
BAut_v_AAut_chi$raw$SNP_id <-  paste(BAut_v_AAut_chi$raw$CHR, BAut_v_AAut_chi$raw$POS, sep = "_")
summary(BAut_v_AAut_chi$reg)
BAut_v_AAut_DAF <- daf_process("~/Projects/Herring/data/British_populations/gensinc_revision/BalticAutumn_AtlanticAutumn.deltaAF", BAut_v_AAut_chi$raw)
daf_and_chi_plots(BAut_v_AAut_DAF, file_prefix = "~/Projects/Herring/doc/Gensinc/Revision/BaltAutvAtlAut")

BSpr_v_ASpr_chi <- chi_qq(chi_file = "~/Projects/Herring/data/British_populations/gensinc_revision/BalticSpring_AtlanticSpring.chi.out", contrast_name = "BaltSprvAtlSpr")
BSpr_v_ASpr_chi$raw$SNP_id <-  paste(BSpr_v_ASpr_chi$raw$CHR, BSpr_v_ASpr_chi$raw$POS, sep = "_")
summary(BSpr_v_ASpr_chi$reg)

BSpr_v_ASpr_DAF <- daf_process("~/Projects/Herring/data/British_populations/gensinc_revision/BalticSpring_AtlanticSpring.deltaAF", BSpr_v_ASpr_chi$raw)
daf_and_chi_plots(BSpr_v_ASpr_DAF, file_prefix = "~/Projects/Herring/doc/Gensinc/Revision/BaltSprvAtlSpr")


#Spawning time
AAut_v_ASpr_chi <- chi_qq(chi_file = "~/Projects/Herring/data/British_populations/gensinc_revision/SpringAtlantic_AutumnAtlantic.chi.out", contrast_name = "AtlAutvAtlSpr")
AAut_v_ASpr_chi$raw$SNP_id <-  paste(AAut_v_ASpr_chi$raw$CHR, AAut_v_ASpr_chi$raw$POS, sep = "_")
summary(AAut_v_ASpr_chi$reg)
AAut_v_ASpr_DAF <- daf_process("~/Projects/Herring/data/British_populations/gensinc_revision/SpringAtlantic_AutumnAtlantic.deltaAF", AAut_v_ASpr_chi$raw)
daf_and_chi_plots(AAut_v_ASpr_DAF, file_prefix = "~/Projects/Herring/doc/Gensinc/Revision/AtlAutvAtlSpr")

BAut_v_BSpr_chi <- chi_qq(chi_file = "~/Projects/Herring/data/British_populations/gensinc_revision/SpringBaltic_AutumnBaltic.chi.out", contrast_name = "BaltAutvBaltSpr")
BAut_v_BSpr_chi$raw$SNP_id <-  paste(BAut_v_BSpr_chi$raw$CHR, BAut_v_BSpr_chi$raw$POS, sep = "_")
summary(BAut_v_BSpr_chi$reg)
save(BAut_v_AAut_chi, BSpr_v_ASpr_chi, AAut_v_ASpr_chi, BAut_v_BSpr_chi, file = "~/Projects/Herring/data/British_populations/gensinc_revision/contrast_chi_collection.RData")

BAut_v_BSpr_DAF <- daf_process("~/Projects/Herring/data/British_populations/gensinc_revision/SpringBaltic_AutumnBaltic.deltaAF", BAut_v_BSpr_chi$raw)
daf_and_chi_plots(BAut_v_BSpr_DAF, file_prefix = "~/Projects/Herring/doc/Gensinc/Revision/BaltAutvBaltSpr")


#British vs NE_atl
NS_v_NEatl_chi <- chi_qq(chi_file = "~/Projects/Herring/data/British_populations/gensinc_revision/NorthSea_otherAtlantic.chi.out", contrast_name = "NS_v_NEatl")
NS_v_NEatl_chi$raw$SNP_id <-  paste(NS_v_NEatl_chi$raw$CHR, NS_v_NEatl_chi$raw$POS, sep = "_")
summary(NS_v_NEatl_chi$reg)
NS_v_NEatl_DAF <- daf_process("~/Projects/Herring/data/British_populations/gensinc_revision/NorthSea_otherAtlantic.deltaAF", NS_v_NEatl_chi$raw)
daf_and_chi_plots(NS_v_NEatl_DAF, file_prefix = "~/Projects/Herring/doc/Gensinc/Revision/NS_v_NEatl")
summary(NS_v_NEatl_chi$reg)



#South vs North
S_v_N_chi <- chi_qq(chi_file = "~/Projects/Herring/data/British_populations/gensinc_revision/South_North.chi.out", contrast_name = "SouthvNorth")
S_v_N_chi$raw$SNP_id <-  paste(S_v_N_chi$raw$CHR, S_v_N_chi$raw$POS, sep = "_")
summary(S_v_N_chi$reg)
S_v_N_DAF <- daf_process("~/Projects/Herring/data/British_populations/gensinc_revision/South_North.deltaAF", S_v_N_chi$raw)
daf_and_chi_plots(S_v_N_DAF, file_prefix = "~/Projects/Herring/doc/Gensinc/Revision/SouthvNorth")
summary(S_v_N_chi$reg)

save(BAut_v_AAut_chi, BSpr_v_ASpr_chi, AAut_v_ASpr_chi, BAut_v_BSpr_chi, S_v_N_chi, NS_v_NEatl_chi, file = "~/Projects/Herring/data/British_populations/gensinc_revision/contrast_chi_collection.RData")




#Haplotypes around AHR
ahr_missense <- read.table("~/Projects/Herring/data/British_populations/gensinc_revision/chr2.log20.missense.list", stringsAsFactors = F, sep = "\t")
ahr_missense$NonSyn <- T
names(ahr_missense)[1:4] <- c("CHR", "POS", "log10_p", "ANN")
region_haplotype_vis(target_chr = "chr2", upper_bound = 15.9e6, lower_bound = 15.8e6, plot_dir = "~/Projects/Herring/doc/Gensinc/Revision/HapDist/", snp_df = ahr_missense)


#BSpr_v_ASpr_DAF <- read.table("~/Projects/Herring/data/British_populations/gensinc_revision/BalticSpring_AtlanticSpring.deltaAF", header = T, stringsAsFactors = F)
#BSpr_v_ASpr_DAF$SNP_id <-  paste(BSpr_v_ASpr_DAF$CHR, BSpr_v_ASpr_DAF$POS, sep = "_")
#BSpr_v_ASpr_DAF$chi <- BSpr_v_ASpr_chi$raw[match(BSpr_v_ASpr_DAF$SNP_id, BSpr_v_ASpr_chi$raw$SNP_id ),"log10P"]
#BSpr_v_ASpr_DAF$mean <- (BSpr_v_ASpr_DAF$BalticSpring + BSpr_v_ASpr_DAF$AtlanticSpring)/2
#BSpr_v_ASpr_DAF$Ftot <- BSpr_v_ASpr_DAF$mean*(1-BSpr_v_ASpr_DAF$mean)
#BSpr_v_ASpr_DAF$FA <- 0.5*(BSpr_v_ASpr_DAF$AtlanticSpring*(1-BSpr_v_ASpr_DAF$AtlanticSpring))
#BSpr_v_ASpr_DAF$FB <- 0.5*(BSpr_v_ASpr_DAF$BalticSpring*(1-BSpr_v_ASpr_DAF$BalticSpring))
#BSpr_v_ASpr_DAF$Fst <- (BSpr_v_ASpr_DAF$Ftot - (BSpr_v_ASpr_DAF$FA + BSpr_v_ASpr_DAF$FB))/BSpr_v_ASpr_DAF$Ftot
#BSpr_v_ASpr_DAF$ZFst <- (BSpr_v_ASpr_DAF$Fst - mean(BSpr_v_ASpr_DAF$Fst))/sd(BSpr_v_ASpr_DAF$Fst)

Z_t <- qnorm(p = 1-(0.05/length(BSpr_v_ASpr_DAF$ZFst))) #Bonferroni corrected Z threshold
Fst_zt <- mean(BSpr_v_ASpr_DAF$Fst) + Z_t*sd(BSpr_v_ASpr_DAF$Fst) #Back-transformed threshold

#png("~/Projects/Herring/doc/Gensinc/Revision/BSvAS_DAFvCHI.png", height = 1000, width = 1000)
#plot(x = 0, y = 1, type = "n", xlim = c(0,1), ylim = c(0, 150), xlab = "Delta Allele Frequency", ylab = "-log10(P)")
#abline(v = seq(0.1, 0.9, 0.1), col = "grey50")
#abline(h = 10, col = "red")
#rect(xleft = c(0, 0.15), xright = c(0.15, 1.2), ybottom = c(10, 0), ytop = c(160, 10), col = c("salmon", "steelblue"), border = NULL)
#points(x = jitter(abs(BSpr_v_ASpr_DAF$deltaAF)), y = BSpr_v_ASpr_DAF$chi, pch = 16, cex = 0.5)
#text(x = 0.05, y = 50, labels = paste0("n = ", sum(BSpr_v_ASpr_DAF$chi > 10 & abs(BSpr_v_ASpr_DAF$deltaAF) < 0.15, na.rm = T)))
#text(x = 0.5, y = 5, labels = paste0("n = ", sum(BSpr_v_ASpr_DAF$chi < 10 & abs(BSpr_v_ASpr_DAF$deltaAF) > 0.15, na.rm = T)))
#dev.off()
summary(abs(BSpr_v_ASpr_DAF$deltaAF)[which(BSpr_v_ASpr_DAF$chi > 10)])
hist(abs(BSpr_v_ASpr_DAF$deltaAF)[which(BSpr_v_ASpr_DAF$chi > 10)])

#png("~/Projects/Herring/doc/Gensinc/Revision/BSvAS_FSTvDAF.png", height = 1000, width = 1000)
#plot(x = 0, y = 1, type = "n", xlim = c(0,1), ylim = c(0, 1), xlab = "Delta Allele Frequency", ylab = "Fst")
#abline(v = seq(0.1, 0.9, 0.1), col = "grey50")
#abline(h = 10, col = "red")
#rect(xleft = c(0, 0.15), xright = c(0.15, 1.2), ybottom = c(10, 0), ytop = c(160, 10), col = c("salmon", "steelblue"), border = NULL)
#points(x = jitter(abs(BSpr_v_ASpr_DAF$deltaAF)), y = BSpr_v_ASpr_DAF$Fst, pch = 16, cex = 0.2)
#text(x = 0.05, y = 50, labels = paste0("n = ", sum(BSpr_v_ASpr_DAF$chi > 10 & abs(BSpr_v_ASpr_DAF$deltaAF) < 0.15, na.rm = T)))
#text(x = 0.5, y = 5, labels = paste0("n = ", sum(BSpr_v_ASpr_DAF$chi < 10 & abs(BSpr_v_ASpr_DAF$deltaAF) > 0.15, na.rm = T)))
dev.off()

#Fst_co <- 0.0849 # from Lamichhaney et al 2012
#Fst_co <- quantile(BSpr_v_ASpr_DAF$Fst, probs = 0.99) # assuming 1% non-neutral loci
Fst_co <- Fst_zt #Very conservative

#png("~/Projects/Herring/doc/Gensinc/Revision/BSvAS_FSTvCHI.png", height = 1000, width = 1000)
#plot(x = 0, y = 1, type = "n", xlim = c(0,1), ylim = c(0, 150), xlab = "Fst", ylab = "-log10(P)")
#abline(v = seq(0.1, 0.9, 0.1), col = "grey50")
#abline(h = 10, col = "red")
#rect(xleft = c(-1, Fst_co), xright = c(Fst_co, 1.2), ybottom = c(10, -10), ytop = c(160, 10), col = c("salmon", "steelblue"), border = NULL)
#points(x = jitter(abs(BSpr_v_ASpr_DAF$Fst)), y = BSpr_v_ASpr_DAF$chi, pch = 16, cex = 0.2)
#text(x = 0.025, y = 50, labels = paste0("n = ", sum(BSpr_v_ASpr_DAF$chi > 10 & BSpr_v_ASpr_DAF$Fst < Fst_co, na.rm = T)))
#text(x = 0.5, y = 5, labels = paste0("n = ", sum(BSpr_v_ASpr_DAF$chi < 10 & BSpr_v_ASpr_DAF$Fst > Fst_co, na.rm = T)))
#dev.off()

#png("~/Projects/Herring/doc/Gensinc/Revision/BSvAS_DAF_hist.png", height = 1000, width = 1000)
#hist_darkorchid <- col2rgb("darkorchid")/255
#hist_salmon <- col2rgb("salmon")/255
#par(mar = c(5,5, 2, 2))
#sub_10 <- hist(abs(BSpr_v_ASpr_DAF$deltaAF)[BSpr_v_ASpr_DAF$chi < 10], breaks = seq(from = 0, to = 1, by = 0.02), plot = F) 
#super_10 <- hist(abs(BSpr_v_ASpr_DAF$deltaAF)[BSpr_v_ASpr_DAF$chi > 10], breaks = seq(from = 0, to = 1, by = 0.02), plot = F)
#sub_10$density <- (sub_10$density/sum(sub_10$density))*100
#super_10$density <- (super_10$density/sum(super_10$density))*100
#plot(sub_10, col = rgb(t(hist_salmon), alpha = 0.5), freq = F, xlim = c(0,1), main = "", xlab = "Delta Allele Frequency", ylab = "Density (percent of markers)",cex.axis = 2, cex.lab = 2) 
#plot(super_10, col = rgb(t(hist_darkorchid), alpha = 0.5), add = T, freq = F, xlim = c(0,1))
#legend(x = "topright", legend = c("SNPs below threshold", "SNPs above threshold"), fill = c(rgb(t(hist_salmon), alpha = 0.5), rgb(t(hist_darkorchid), alpha = 0.5)), cex = 2)
#dev.off()


p_tot <- 0.5
daf <- 0.2
f_tot <- p_tot*(1-p_tot)
p1 <- p_tot - daf/2
f1 <- 0.5*(p1*(1-p1))
p2 <- p_tot + daf/2
f2 <- 0.5*(p2*(1-p2))
f_st <- (f_tot - (f1 + f2))/f_tot

#Support functions
daf_process <- function(daf_file, chi_df){
  DAF_df <- read.table(daf_file, header = T, stringsAsFactors = F)
  DAF_df $SNP_id <-  paste(DAF_df$CHR, DAF_df$POS, sep = "_")
  DAF_df$chi <- chi_df[match(DAF_df$SNP_id, chi_df$SNP_id ),"log10P"]
  DAF_df$mean <- (DAF_df[,3] + DAF_df[,4])/2
  DAF_df$Ftot <-  DAF_df$mean*(1- DAF_df$mean)
  DAF_df$F1 <- 0.5*(DAF_df[,3]*(1-DAF_df[,3]))
  DAF_df$F2 <- 0.5*(DAF_df[,4]*(1-DAF_df[,4]))
  DAF_df$Fst <- (DAF_df$Ftot - (DAF_df$F1 + DAF_df$F2))/DAF_df$Ftot
  DAF_df$ZFst <- (DAF_df$Fst - mean(DAF_df$Fst))/sd(DAF_df$Fst)
  return(DAF_df)
}

chi_qq <- function(chi_file, chi_data = NULL, contrast_name){
  chi_raw <- chi_data
  if(is.null(chi_data[1])){
    chi_raw <- read.table(chi_file, header = T, stringsAsFactors = F)
  }
  obs_q <- chi_raw$log10P[order(chi_raw$log10P)]
  exp_q <- -log10((length(chi_raw$log10P):1)/length(chi_raw$log10P))
  lambda_reg <- lm(obs_q[1:floor(length(exp_q)*0.9)]~0+exp_q[1:floor(length(exp_q)*0.9)])
  png(paste0("~/Projects/Herring/doc/Gensinc/Revision/", contrast_name,"_qq.png"), height = 1000, width = 1000)
  plot(x = exp_q, y = obs_q, pch = 16, cex = 0.5, xlim = c(0,7), ylim = c(0,7), main = contrast_name, xlab = "Expected", ylab = "Observed")
  abline(a = 0, b = 1, col = "red")
  abline(lambda_reg)
  dev.off()
  
  png(paste0("~/Projects/Herring/doc/Gensinc/Revision/", contrast_name,"_full_qq.png"), height = 1000, width = 1000)
  plot(x = exp_q, y = obs_q, pch = 16, cex = 0.5, xlim = c(0,7), main = contrast_name, xlab = "Expected", ylab = "Observed")
  abline(a = 0, b = 1, col = "red")
  abline(lambda_reg)
  dev.off()
  return(list(raw =  chi_raw, reg =  lambda_reg))
}

daf_and_chi_plots <- function(daf_df, file_prefix = ""){
  png(paste0(file_prefix,"DAFvCHI.png"), height = 1000, width = 1000)
  plot(x = 0, y = 1, type = "n", xlim = c(0,1), ylim = c(0, 150), xlab = "Delta Allele Frequency", ylab = "-log10(P)")
  points(x = jitter(abs(daf_df$deltaAF)), y = daf_df$chi, pch = 16, cex = 0.5)
  dev.off()

  png(paste0(file_prefix,"FSTvDAF.png"), height = 1000, width = 1000)
  plot(x = 0, y = 1, type = "n", xlim = c(0,1), ylim = c(0, 1), xlab = "Delta Allele Frequency", ylab = "Fst")
  points(x = jitter(abs(daf_df$deltaAF)), y = daf_df$Fst, pch = 16, cex = 0.2)
  dev.off()
  
  png(paste0(file_prefix,"FSTvCHI.png"), height = 1000, width = 1000)
  plot(x = 0, y = 1, type = "n", xlim = c(0,1), ylim = c(0, 150), xlab = "Fst", ylab = "-log10(P)")
  points(x = daf_df$Fst, y = daf_df$chi, pch = 16, cex = 0.2)
  dev.off()
  
  png(paste0(file_prefix,"DAF_hist.png"), height = 1000, width = 1000)
  rgb_darkorchid <- t(col2rgb("darkorchid")/255)
  rgb_salmon <- t(col2rgb("salmon")/255)
  par(mar = c(5,5, 2, 2))
  sub_10 <- hist(abs(daf_df$deltaAF)[daf_df$chi < 10], breaks = seq(from = 0, to = 1, by = 0.02), plot = F) 
  super_10 <- hist(abs(daf_df$deltaAF)[daf_df$chi > 10], breaks = seq(from = 0, to = 1, by = 0.02), plot = F)
  sub_10$density <- (sub_10$density/sum(sub_10$density))*100
  super_10$density <- (super_10$density/sum(super_10$density))*100
  
  plot(sub_10, col = rgb(rgb_salmon, alpha = 0.5), freq = F, xlim = c(0,1), main = "", xlab = "Delta Allele Frequency", ylab = "Density (percent of markers)",cex.axis = 2, cex.lab = 2) 
  plot(super_10, col = rgb(rgb_darkorchid, alpha = 0.5), add = T, freq = F, xlim = c(0,1))
  legend(x = "topright", legend = c(paste("SNPs below threshold; n =", sum(sub_10$counts)), paste("SNPs above threshold; n =", sum(super_10$counts))), fill = c(rgb(rgb_salmon, alpha = 0.5), rgb(rgb_darkorchid, alpha = 0.5)), cex = 2)
  dev.off()
}


pool_freq_hm <- function(freq_df, pdf_file = "./pool_hm.pdf", pool_order_vec, hm_col = c("blue4", "gold1"), ...){
  pdf(file = pdf_file, height = 10, width = 10)
  par(xpd = NA)
  row_lab_vec <- sub("[A-Z0-9]+_", "", colnames(freq_df)[pool_order_vec])
  target_snp_id <- rownames(freq_df)
  #pool_order_vec <- c(17, 12, 20, 6, 7, 8, 18, 5, 4, 21, 19, 10, 9, 3, 2, 11, 14, 15, 16, 13, 1)
  cr <- colorRamp(hm_col)
  crp = colorRampPalette(colors = hm_col)(500)
  col_colors <- rep("black", dim(freq_df)[1])
  col_colors[freq_df$NonSyn] <- "darkorchid"
  heatmap(t(as.matrix(freq_df[,pool_order_vec])), scale = "none", Colv = NA, Rowv = NA, labCol = target_snp_id, labRow = row_lab_vec, col = crp, margins =c(12,12), ColSideColors = col_colors)
  scale_x_vec <- seq(from = par("usr")[1], to = par("usr")[2]*0.5, length.out = 505)
  rect(ybottom = par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.09, xleft = scale_x_vec[-(501:505)], ytop = par("usr")[3] - (par("usr")[4] -par("usr")[3])*0.06, xright =  scale_x_vec[-(1:5)], col = rgb(cr((1:500)/500), maxColorValue=255), border = NA)
  text(x=scale_x_vec[c(10,253,505)], y = par("usr")[3] - (par("usr")[4]-par("usr")[3])*0.10, labels = paste0(c(0, 50, 100), "%"))
  dev.off()
  #return(par("usr"))
}
