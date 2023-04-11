set.seed(12345)

options(stringsAsFactors = FALSE)
options(warn = -1)
dataDir = "../data/" # assumes being run from within R-scripts directory as e.g. Rscript 00_script.R
figureDir = "../output-figures/" # for figures
outputDir = "../output-data/" # for output data

# SAVE FIGURE FUNCTION
saveFigure = function(some.ggplot, filename, width, height, with.date=TRUE, dpi=600){
  if (with.date==TRUE){
    ggsave(some.ggplot, file=paste0(figureDir, filename, "-", Sys.Date(), ".pdf"), width=width, height=height, dpi=dpi)
  }
  if (with.date==FALSE){
    ggsave(some.ggplot, file=paste0(figureDir, filename, ".pdf"), width=width, height=height, dpi=dpi)
  }
}

# Pangenome component colours
component.colours = c("#1f78b4", "#a6cee3", "#b2df8a")
names(component.colours) = c("Core", "Non-core", "Plasmid")

# Colours for size of plasmids
size.colours = c("#fee5d9","#fcae91","#fb6a4a","#cb181d")

# LIBRARIES
print("Loading libraries...")
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(formatR))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(phytools))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(tidyr))
print("Done!")
