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


# STATISTICS OF R-M SYSTEMS
rms.db = read.csv(paste0(dataDir, "2023-04-06-rms-db.csv"))
rms.db$genome = gsub("\\.[0-9][0-9][0-9][0-9]$", "", rms.db$contig)
  
# Total systems (non-unique)
print("Total putative R-M systems:") 
print(nrow(rms.db))
print("Number of unique genomes containing at least one R-M system:")
print(length(unique(rms.db$genome)))

# Number containing >1 RM system
table(table(rms.db$genome)>1)
# Table of numbers of RM systems
table(table(rms.db$genome))

# Investigate table further
genome.counts = table(rms.db$genome)
genome.counts[genome.counts>20]

# To avoid overlapping systems, we summarise by genome and unique targets
genomes.and.unique.targets = rms.db %>% group_by(genome) %>% 
  summarise(targets=length(unique(sequence)),
            total.putative.systems=length(sequence))

table(genomes.and.unique.targets$targets>1)
table(genomes.and.unique.targets$targets)
as.character(genomes.and.unique.targets[genomes.and.unique.targets$targets>9,"genome"])

# Associated paper text:
#We detected 8,616 putative R-M systems with a predicted target motif, 
# with 2,592 genomes carrying at least one R-M system (30.3%). 
# Some putative R-M systems ‘overlapped’ e.g. the same REase could be included
# in multiple putative systems if there were multiple MTases in close proximity. 
# To avoid overcounting these systems, we considered the unique targets recognised 
# per genome. Of the R-M-containing genomes, 1,875/2,592 (72.3%) had R-M systems 
# recognising one motif (range: 0-18 putative R-M systems; Helicobacter pylori 
# genomes accounted for all those with >9 R-M systems). 
