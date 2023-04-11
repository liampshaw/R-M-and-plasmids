# R-M target analysis: taxonomic hierarchy"

source('setup.R')

# READ IN DATASETS
K = 4
SUBSAMPLING = 50000
main.df <- read.csv(paste0(dataDir, K, '-merged-genome-results-levels-inclusive.csv'), header=F)
colnames(main.df) <- c("species", "kmer_category", "genome", "section", "subsampling", "n", "rank", "score")
main.df$section <- ordered(main.df$section, 
                           levels=c("core", "accessory_chrom", "accessory_plas"))
main.df$kmer_category <- ordered(main.df$kmer_category,
                                 levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"))
summary.df.all.4 <- main.df[which(main.df$subsampling==SUBSAMPLING & main.df$score!="NaN" & main.df$kmer_category!="Palindromic"),] %>% 
  group_by(species, kmer_category, section) %>% 
  summarise(score=mean(score),
            rank=mean(rank))

K = 5
SUBSAMPLING = 50000
main.df <- read.csv(paste0(dataDir, K, '-merged-genome-results-levels-inclusive.csv'), header=F)
colnames(main.df) <- c("species", "kmer_category", "genome", "section", "subsampling", "n", "rank", "score")
main.df$section <- ordered(main.df$section, 
                           levels=c("core", "accessory_chrom", "accessory_plas"))
main.df$kmer_category <- ordered(main.df$kmer_category,
                                 levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"))
summary.df.all.5 <- main.df[which(main.df$subsampling==SUBSAMPLING & main.df$score!="NaN" & main.df$kmer_category!="Palindromic"),] %>% 
  group_by(species, kmer_category, section) %>% 
  summarise(score=mean(score),
            rank=mean(rank),
            .groups = "drop")

K= 6
SUBSAMPLING = 50000
main.df <- read.csv(paste0(dataDir, K, '-merged-genome-results-levels-inclusive.csv'), header=F)
colnames(main.df) <- c("species", "kmer_category", "genome", "section", "subsampling", "n", "rank", "score")
main.df$section <- ordered(main.df$section, 
                           levels=c("core", "accessory_chrom", "accessory_plas"),
                           labels=c("Core", "Non-core", "Plasmid"))
main.df$kmer_category <- ordered(main.df$kmer_category,
                                 levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"))
summary.df.all.6 <- main.df[which(main.df$subsampling==SUBSAMPLING & main.df$score!="NaN" & main.df$kmer_category!="Palindromic"),] %>% 
  group_by(species, kmer_category, section) %>% 
  summarise(score=mean(score),
            rank=mean(rank))

# SCORE PLOT FUNCTION
scorePlot <- function(df, panel.one=FALSE, legend.position="none"){
  p.panel = ggplot(df[,], aes(kmer_category, score, group=interaction(section, species)))+
    geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
    geom_line(alpha=0.2, aes(colour=section), size=0.5)+
    facet_wrap(~section, nrow=3)+ 
    theme_bw()+
    scale_colour_manual(values=as.character(component.colours))+
    theme(panel.grid = element_blank(), axis.text=element_text(colour="black"))+
    xlab("")+
    ylab("")+
    guides(fill="none")+
    theme(axis.text.x=element_text(angle=45, hjust=1, size=8))+
    theme(legend.position = "none")
  p.all = ggplot(df[,], aes(kmer_category, score))+
    geom_hline(yintercept = 0, linetype='dashed', colour='black', size=1)+
    stat_smooth(aes(colour=section, group=section), size=2, alpha=0.2, se=FALSE)+
    theme_bw()+
    scale_colour_manual(values=as.character(component.colours))+
    theme(panel.grid = element_blank(), axis.text=element_text(colour="black"))+
    xlab("")+
    ylab("Exceptionality score of targets")+
    guides(fill="none")+
    theme(axis.text.x=element_text(angle=45, hjust=1))+
    theme(legend.position = "none")
  
  p.l = cowplot::plot_grid(p.all+ggtitle("(d)"), p.panel+ggtitle("(e)"), nrow=1, 
                           rel_widths = c(1, 0.4),
                           align='h', axis='b')
  if (panel.one==TRUE){
    return(p.all+theme(legend.position = legend.position))
  }
  else{
    return(p.l)
  }
}

# SCORE PLOTS, SMOOTHED COMPONENTS
p.scores.4 = scorePlot(summary.df.all.4)
p.scores.5 = scorePlot(summary.df.all.5)
p.scores.6 = scorePlot(summary.df.all.6)

# Save figures
saveFigure(p.scores.4, 'FigureS7_smoothed-component-K-4', width=6, height=6)
saveFigure(p.scores.5, 'FigureS8_smoothed-component-K-5', width=6, height=6)
saveFigure(p.scores.6,'Figure2de_smoothed-component-K-6', width=6, height=6)



# COMBINE SMOOTHED COMPONENT PLOTS FOR ALL K
p.4 = scorePlot(summary.df.all.4, panel.one = TRUE)
p.5 = scorePlot(summary.df.all.5, panel.one = TRUE)
p.6 = scorePlot(summary.df.all.6, panel.one = TRUE, legend.position = c(0.8, 0.2))
p.scores.combined = cowplot::plot_grid(p.4+ggtitle("(a) k=4"),
                                       p.5+ggtitle("(b) k=5"),
                                       p.6+ggtitle("(c) k=6")+
                                         theme(legend.position = "none"),
                                       nrow=1)
saveFigure(p.scores.combined,
           'FigureX_smoothed-components-all', width=7, height=2.5)

# FIGURE 2 PLOT CONSTRUCTION
schematic_file <- system.file("../manuscript-figures/schematic-approach.png", package = "cowplot")
p2 <- ggdraw() + draw_image("../manuscript-figures/schematic-approach.png", scale = 0.9)
p.combined.figure2 = cowplot::plot_grid(p2, p.scores.6,
                   nrow=2, 
                   rel_heights = c(1,1.3))
png(paste0(figureDir, "Figure2_k-6-subsampling-", SUBSAMPLING, "-", Sys.Date()), width=3500, height=3500, res=600)
p.combined.figure2
dev.off()
saveFigure(p.combined.figure2, paste0("Figure2_", "k-6-subsampling-", SUBSAMPLING), width=7, height=8)