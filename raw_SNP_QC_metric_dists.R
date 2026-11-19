#!/usr/bin/env Rscript

# parse arguements
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Error: Missing required arguments. Usage: Rscript raw_SNP_QC_metric_dists.R <figures_dir_path> <input_file>", call. = FALSE)
}
dir <- args[1]
file <- args[2]

# Load and/or download required packages
required.pkg <- c("ggplot2", "cowplot", "svglite")
pkgs.not.installed <- required.pkg[!sapply(required.pkg, function(p) require(p, character.only=T))]
if (length(pkgs.not.installed) > 0) install.packages(pkgs.not.installed, dependencies=TRUE)
lapply(required.pkg, library, character.only = TRUE)

# run script
snp_qual <- read.csv(file, header = T, sep = "\t")

qd <- 2        # filter values < 2
fs <- 60       # filter > 60
sor <- 3       # filter > 3
mq <- 166.67   # filter < 166.67; analagous to phred < 40
mqrs <- -12.5  # filter < -12.5
rprs <- -8     # filter < -8
co <- c(qd, fs, sor, mq, mqrs, rprs)

i <- 1
for (col in colnames(snp_qual)) {
  print(col)
  print(co[i])
  d <- data.frame("var" = na.omit(as.numeric(snp_qual[[col]])))
  p <- ggplot(d, aes(x = var)) +
    geom_density(color = "grey", fill = "grey", alpha = 0.6) +
    geom_vline(xintercept = co[i], color = "red", linetype = "dashed") +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(x = col) +
    theme_classic()
  if (i >= 5) {
    p <- p + scale_x_continuous(limits = c((co[i]-.5),max(d$var)))
  }
  if (i == 2 || i == 3) {
    p <- p + annotate("rect", xmin = co[i], xmax = Inf, ymin = 0, ymax = Inf,
                      fill = "red", alpha = 0.2)
  } else {
    p <- p + annotate("rect", xmin = -Inf, xmax = co[i], ymin = 0, ymax = Inf,
                      fill = "red", alpha = 0.2)
  }
    assign(paste("p", i, sep = ""), p)
    i <- i+1
}
plots <- mget(ls(pattern = "^p[0-9]"))
out <- plot_grid(plotlist = plots, labels = LETTERS[1:6])

setwd(dir)
ggsave("raw_SNP_QC_metric_distributions.svg", plot = out, width = 6.24, height = 4.32)
