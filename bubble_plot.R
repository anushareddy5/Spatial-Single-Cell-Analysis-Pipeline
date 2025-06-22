#!/usr/bin/env Rscript
# bubble_plot.R

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "CSV input file"),
  make_option(c("-o", "--output"), type = "character", help = "PDF/PNG output")
)
opt <- parse_args(OptionParser(option_list = option_list))

message("[", Sys.time(), "] Reading input: ", opt$input)
df <- read.csv(opt$input, row.names = 1)

message("[", Sys.time(), "] Melting data frame")
melted <- melt(df,
  id.vars = c("gene", "sample", "region"),
  variable.name = "metric", value.name = "value"
)

message("[", Sys.time(), "] Building bubble plot")
p <- ggplot(melted, aes(x = region, y = gene, size = value1, color = value2)) +
  geom_point() +
  facet_wrap(~sample) +
  theme_minimal() +
  labs(size = "Size", color = "Color")

message("[", Sys.time(), "] Saving plot to ", opt$output)
ggsave(opt$output, plot = p, width = 6, height = 4)
message("[", Sys.time(), "] Done.")
