#!/usr/bin/env Rscript

# parse arguements
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Error: Missing required arguments. Usage: Rscript raw_SNP_QC_metric_dists.R <figures_dir_path> <input_file>", call. = FALSE)
}
dir <- args[1]
file <- args[2]

# Load and/or download required packages
required.pkg <- c("ggplot2", "dplyr", "scales", "svglite")
pkgs.not.installed <- required.pkg[!sapply(required.pkg, function(p) require(p, character.only=T))]
if (length(pkgs.not.installed) > 0) install.packages(pkgs.not.installed, dependencies=TRUE)
lapply(required.pkg, library, character.only = TRUE)


# read in data
d <- read.csv(file, header = T, sep = "\t")
d$DP <- as.numeric(d$DP)
d <- na.omit(d)
d <- d[(d$GT != "./."),]
d$GT <- factor(d$GT, levels = c("0/0", "1/1", "0/1"))

# calculate per-sample mean
dp_mean <- d %>%
  group_by(SAMPLE) %>%
  summarise(DPmean = mean(DP, na.rm = TRUE))

# potential upper bounds
cat(paste("mean depth range:", 
          round(max(dp_mean$DPmean)-min(dp_mean$DPmean), 1), "\n"))
sd6 <- mean(d$DP) + sd(d$DP)*6
max3 <- max(dp_mean$DPmean)*3

# calculate % GTs filtered
sd6p <- round(nrow(d[(d$DP > sd6),])/nrow(d)*100,2)
max3p <- round(nrow(d[(d$DP > max3),])/nrow(d)*100,2)
m10p <- round(nrow(d[(d$DP < 10),])/nrow(d)*100,0)

# base plot
sub <- d %>% slice_sample(n = 1000000)  # randomly subsample 1M rows for quicker plotting
p <- ggplot(sub, aes(x = DP)) +
  geom_histogram(aes(y = after_stat(count / sum(count))), 
                 binwidth = 1, alpha = 0.6) +
  labs(x = "Genotype read depth", y = "% Genotypes")
max_y <- layer_scales(p)$y$get_limits()[2]

# if min DP 10 filters < 23%, generate plot, otherwise calculate and plot min DP 6-9
if (m10p < 23) {
  pout <- p + geom_vline(xintercept = c(sd6, max3, 10), linetype = "dashed") +
    scale_x_continuous(expand = c(0.12,0), limits = c(0,max(c(sd6,max3))*.25+max(c(sd6,max3)))) +
    scale_y_continuous(expand = c(0,0)) +
    annotate("text", x = sd6, y = (max_y-max_y*.1), hjust = 0,
             label = paste(" SD6=", sd6p,"%", sep = "")) +
    annotate("text", x = max3, y = (max_y-max_y*.1), hjust = 0,
             label = paste(" max*3=", max3p,"%", sep = "")) +
    annotate("text", x = 10, y = (max_y-max_y*.1), hjust = 1,
             label = paste("10x=", m10p,"% ", sep = "")) +
    theme_classic()
} else {
  for (i in 6:9) {
    perc <- round(nrow(d[(d$DP < i),])/nrow(d)*100,1)
    assign(paste("m",i,"p", sep = ""), perc)
  }
  pout <- p + geom_vline(xintercept = c(sd6, max3, 6:9), linetype = "dashed") +
    scale_x_continuous(expand = c(0.13,0), limits = c(0,max(c(sd6,max3))*.25+max(c(sd6,max3)))) +
    scale_y_continuous(expand = c(0,0)) +
    annotate("text", x = sd6, y = (max_y-max_y*.1), hjust = 0,
             label = paste(" SD6: ", sd6p,"%", sep = "")) +
    annotate("text", x = max3, y = (max_y-max_y*.1), hjust = 0,
             label = paste(" max*3: ", max3p,"%", sep = "")) +
    annotate("text", hjust = 1, x = c(6:9), 
             y = c((max_y-max_y*.4), (max_y-max_y*.3),
                   (max_y-max_y*.2), (max_y-max_y*.1)), 
             label = c(paste("6x=", m6p,"% ", sep = ""),
                       paste("7x=", m7p,"% ", sep = ""),
                       paste("8x=", m8p,"% ", sep = ""),
                       paste("9x=", m9p,"% ", sep = ""))) +
    theme_classic()
}

setwd(dir)
ggsave("raw_SNP_GT_read_depth_distribution.svg", plot = pout, width = 5.7, height = 3)
