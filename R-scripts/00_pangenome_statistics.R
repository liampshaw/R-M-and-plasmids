# PANGENOME STATISTICS
  
source('setup.R')

 
# READ STATISTICS AND REARRANGE
pangenome.stats = read.csv(paste0(dataDir, "pangenome-stats.csv"),
                header=T,
                stringsAsFactors=F)
pangenome.stats$mean.genes.per.genome = pangenome.stats$total.genes/pangenome.stats$genomes
pangenome.stats$mean.plasmid.genes.per.genome = pangenome.stats$total.plasmid.genes/pangenome.stats$genomes
pangenome.stats$mean.nonplasmid.noncoregenes.per.genome = pangenome.stats$mean.genes.per.genome-pangenome.stats$mean.plasmid.genes.per.genome-pangenome.stats$core.genes

pangenome.stats.melt = melt(pangenome.stats, id.vars = "species", measure.vars = c("core.genes",
                                                                                   "mean.nonplasmid.noncoregenes.per.genome",
                                                                                   "mean.plasmid.genes.per.genome"))
pangenome.stats.melt$variable = ordered(pangenome.stats.melt$variable,
                                        levels=rev(c("core.genes",
                                                                                   "mean.nonplasmid.noncoregenes.per.genome",
                                                                                   "mean.plasmid.genes.per.genome")),
                                        labels=rev(c("Core", "Non-core", "Plasmid")))

# PANGENOME PLOT
p.pangenome = ggplot(pangenome.stats.melt, aes(species, value, fill=variable))+
  geom_bar(stat="identity", position="stack")+
  theme_bw()+
  coord_flip()+
  scale_fill_manual(values=component.colours)+
  theme(panel.grid = element_blank())
saveFigure(p.pangenome, "FigureX_pangenome-plot.pdf", 
           width=5.5, height=10)



# CONTIGS DATABASE
contigs.db = read.csv(paste0(dataDir, "2022-11-18-contig-database.csv"), 
                      header=T,
                      stringsAsFactors = F,
                      row.names = 1)
#View(contigs.db %>% filter(contig.number>1, plasmid!=TRUE))
