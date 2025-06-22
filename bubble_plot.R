# bubble_plot.R

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(reshape2)
})

option_list <- list(
  make_option(c("-i","--input"), type="character", help="CSV input"),
  make_option(c("-o","--output"), type="character", help="PDF/PNG output")
)
opt <- parse_args(OptionParser(option_list=option_list))

df <- read.csv(opt$input, row.names=1)
melted <- melt(df, id.vars=c("gene","sample","region"), 
               measure.vars=c("value1","value2"),
               variable.name="metric", value.name="val")

ggplot(melted, aes(x=region, y=gene, size=value1, color=value2)) +
  geom_point() + theme_minimal() +
  scale_size_continuous(name="Size") +
  scale_color_gradient(name="Color") +
  facet_wrap(~sample) +
  theme(axis.text.y=element_text(size=6))

ggsave(opt$output, width=6, height=4)
