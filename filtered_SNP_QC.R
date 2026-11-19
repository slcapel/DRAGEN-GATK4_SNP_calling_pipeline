#!/usr/bin/env Rscript

# parse arguements
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Error: Missing required arguments. Usage: Rscript raw_SNP_QC_metric_dists.R <figures_dir_path> <filtered_SNPs.stats.tsv> <filtered_SNPs.hetVAF.tsv> <filtered_SNPs.FILT_DIST.tsv> <MINMAF>", call. = FALSE)
}
dir <- args[1]
stats <- args[2]
hvaf <- args[3]
filt <- args[4]
mm <- args[5]

# Load and/or download required packages
required.pkg <- c("ggplot2", "ggpubr", "tidyr","data.table", "cowplot", "svglite")
pkgs.not.installed <- required.pkg[!sapply(required.pkg, function(p) require(p, character.only=T))]
if (length(pkgs.not.installed) > 0) install.packages(pkgs.not.installed, dependencies=TRUE)
lapply(required.pkg, library, character.only = TRUE)

# read in stats and calculate percentages
data <- read.csv(stats, header = T, sep = "\t")
data$prop.Rhom <- data$nRefHom/(data$nRefHom+data$nNonRefHom+data$nHets)
data$prop.NRhom <- data$nNonRefHom/(data$nRefHom+data$nNonRefHom+data$nHets)
data$prop.het <- data$nHets/(data$nRefHom+data$nNonRefHom+data$nHets)
data$perc.miss <- (data$nMissing/(data$nMissing+data$nRefHom+data$nNonRefHom+data$nHets))*100

# depth vs. Ho
fit <- lm(prop.het ~ average.depth, data)
r2adj <- round(summary(fit)$adj.r.squared, 3)
b <- round(coef(fit)[[1]], 1)
m <- round(coef(fit)[[2]], 4)
stats <- bquote(paste(italic(y) == .(m), italic(x)," + ", .(b), "; ",
                      R[adj]^2 == .(r2adj)))
ho_range <- range(data$prop.het)[2]-range(data$prop.het)[1]
ho_mid <- range(data$prop.het)[1]+ho_range/2
plot_ho_dp <- ggplot(data, aes(x = average.depth, y = prop.het)) +
  geom_vline(xintercept = mean(data$average.depth), linetype = "dashed",
             color = "#b8b8b8") +
  geom_point(size = 1.5, alpha = 0.6, color = "#157882") +
  geom_smooth(method = "lm", se = F, fullrange = TRUE,
              linewidth = 0.75, color = "darkgrey") +

  scale_x_continuous(expand = c(0,0), limits = c(0,max(data$average.depth)+3)) +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(ho_mid-(ho_range*1.2/2),
                                ho_mid+(ho_range*1.2/2))) +
  labs(x = "Mean SNP Depth", y = bquote(italic(H)[o]), subtitle = stats) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 9.5))

# missingness vs. Ho
fit <- lm(prop.het ~ perc.miss, data)
r2adj <- round(summary(fit)$adj.r.squared, 3)
b <- round(coef(fit)[[1]], 1)
m <- round(coef(fit)[[2]], 4)
stats <- bquote(paste(italic(y) == .(m), italic(x)," + ", .(b), "; ",
                      R[adj]^2 == .(r2adj)))
plot_ho_mis <- ggplot(data, aes(x = perc.miss, y = prop.het)) +
  geom_point(size = 1.5, alpha = 0.6, color = "#157882") +
  geom_smooth(method = "lm", se = F, fullrange = TRUE,
              linewidth = 0.75, color = "darkgrey") +
  scale_y_continuous(expand = c(0,0), 
                     limits = c(ho_mid-(ho_range*1.2/2),
                                ho_mid+(ho_range*1.2/2))) +
  labs(x = "% Missingness", y = bquote(italic(H)[o]), subtitle = stats) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 9.5))

# relative genotype ratios
pal=c("#D59583", "#964242", "#501C1C")
gt_long <- data[,c(1,13:15)] %>% pivot_longer(cols = starts_with("prop"),
                                              names_to = "GT_type",
                                              values_to = "perc_GT")
gt_long$GT_type <- gsub("prop.het", "Het", gt_long$GT_type)
gt_long$GT_type <- gsub("prop.NRhom", "AltHom", gt_long$GT_type)
gt_long$GT_type <- gsub("prop.Rhom", "RefHom", gt_long$GT_type)
gt_long$GT_type <- factor(gt_long$GT_type, levels = c("RefHom", "AltHom", "Het"))
plot_gt <- ggplot(gt_long, aes(x = GT_type, y = perc_GT,
                               fill = GT_type, color = GT_type)) +
  geom_violin(linewidth = 0.3) +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  labs(x = "Genotype", y = "% Called") +
  theme_classic() +
  theme(legend.position = "none")

# het VAF S.D.
vaf <- fread(hvaf, sep = "\t", header = F)
names(vaf) <- c("sample", "VAF", "DP")
inds <- unique(vaf$sample)

sd_vaf_trim <- data.frame()
n <- 25
a_tot <- c()
VAF_tot <- c()
for (i in inds) {
  print(i)
  vaf_ind <- vaf[(vaf$sample == i),2]
  a <- vaf_ind$VAF
  b <- data.frame(matrix(a[1:(floor(length(a)/n)*n)], ncol=n,byrow=TRUE))

  b <- unlist(lapply(transpose(b), function(x) sd(x)))
  b <- unlist(unname(lapply(b, function(x) rep(x,n))))
  b <- runmed(b,n)

  Q <- unname(quantile(b))
  upper.outlier.bound<- Q[2]+1.5*(Q[2]-Q[1])
  lower.outlier.bound <- Q[2]-1.5*(Q[2]-Q[1])

  VAF <- a[1:length(b)]
  VAF <- VAF[b < upper.outlier.bound & b > lower.outlier.bound]

  sd_vaf_trim <- rbind(sd_vaf_trim, data.frame("sample" = i,
                                               "sd_vaf" = sd(a),
                                               "sd_vaf_trim" = sd(VAF)))

  a_tot <- append(a, a_tot)
  VAF_tot <- append(VAF, VAF_tot)
}
data <- merge(data, sd_vaf_trim, by = "sample")
q <- quantile(data$sd_vaf_trim)

fit <- lm(sd_vaf_trim ~ average.depth, data)
r2adj <- round(summary(fit)$adj.r.squared, 3)
b <- round(coef(fit)[[1]], 1)
m <- round(coef(fit)[[2]], 4)
stats <- bquote(paste(italic(y) == .(m), italic(x)," + ", .(b), "; ",
                      R[adj]^2 == .(r2adj)))

vaf_dp <- ggplot(data, aes(x = average.depth, y = sd_vaf_trim)) +
  geom_point(size = 1.5, alpha = 0.6, color = "#157882") +
  geom_smooth(method = "lm", se = T, level = 0.99, fullrange = TRUE,
              linewidth = 0.75, color = "darkgrey") +
  labs(x = "Mean SNP Depth", y = expression(sigma[het]^","),
       subtitle = stats) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 9.5))

ggsave(paste(dir,"/filtered_SNP_QC_plots.svg", sep = ""),
       plot = plot_grid(plot_ho_dp, plot_ho_mis, vaf_dp, plot_gt, nrow = 2, labels = LETTERS[1:4]),
       width = 7, height = 7)


# SNP filtering by scaffold
data <- read.csv(filt, header = T, sep = "\t")
data$filter <- factor(data$filter, levels = c("DRAGENHardQUAL", "MISSING0.2", "INVARIABLE",
	       			   	      paste("MINMAF",mm, sep = ""), "MINMQ166.67", "MINQD2",
					      "MAXSOR3", "MAXFS60", "MINMQRankSum-12.5", "MINReadPosRankSum-8",
					      "PASS"))
pal <- c("#628146","#7695d5","#752c57","#b3ba58","#45bcce",
	 "#a35442","#50498d","#bf81c0","#54b48f","#c28d5a","#c86d88")

p_filt <- ggplot(data, aes(x = scaff, y = no.snps, group = filter, fill = filter)) +
  geom_bar(position="fill", stat="identity", width = 1) +
  scale_y_continuous(expand = c(0,0)) +
  ylab("Prop. Variants") +
  xlab("Scaffolds") +
  labs(fill = "Filter") +
  scale_fill_manual(values = pal) +
  theme(panel.background = element_blank(), plot.margin = margin(l = 8, t = 3, unit = "mm"),
        axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1))

ggsave(paste(dir,"/per-scaffold_SNP_filtering.svg", sep = ""),
       plot = p_filt,
       width = 6.5, height = 4)
