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
rms.db$species = gsub("\\..*", "", rms.db$genome)
rms.db$type = ifelse(gsub(".*\\.", "", rms.db$contig)=="0001", 
                     "chromosome", "plasmid")
print("R-M system carriage on chromosome vs. plasmid:")
table(rms.db$type)

species.names = read.csv(paste0(dataDir, "species.csv"), header=F, stringsAsFactors = F)
species.long = species.names$V1
names(species.long) = species.names$V2
rms.db$species = species.long[rms.db$species]

# Write species and sequence occurrences
write.csv(rms.db[,c("species", "sequence")], 
          file=paste0(dataDir, "species-motifs.csv"),
          quote=F,
          row.names = F)

  
write.csv(rms.db[,c("species", "sequence")], 
          file=paste0(dataDir, "species-motifs.csv"),
          quote=F,
          row.names = F)
# Can then run compute-motifs.py in data-generation-scripts to 
# get the unambiguous base numbers for Table 1

# Write summary
rms.db.summary = rms.db %>% group_by(sequence) %>% 
  summarise(n=length(species))
  
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
genomes.and.unique.targets = rms.db %>% group_by(genome, species) %>% 
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

# Number of species
print("Species without any R-M systems:")
cat(species.long[!species.long %in% unique(genomes.and.unique.targets$species)], 
    sep=", ")

# Unique motifs
print("Number of unique motifs:")
length(unique(rms.db$sequence))

# Summarise for Table 1
rms.db$length.of.target = nchar(rms.db$sequence)
print("(most of) Table 1 in manuscript:")
rms.db %>% group_by(length.of.target) %>%
  summarise(rebase.motifs=length(unique(sequence)),
            genomes=length(unique(genome)),
            species=length(unique(species)))


# Write to file 
cat(unique(rms.db$sequence), sep="\n", file="../data/rms-motifs.txt")

# Targeting by species level
motifs.by.species = rms.db %>% group_by(sequence) %>% 
  summarise(species=length(unique(species)),
            genomes=length(unique(genome)))
print("Motifs targeted by only a single species:")
table(motifs.by.species$species==1)
print("Median, range for number of species targeting motif:")
median(motifs.by.species$species)
range(motifs.by.species$species)
print("Median, range for number of genomes targeting motif:")
median(motifs.by.species$genomes)
range(motifs.by.species$genomes)

# PLASMID PTUs
print("Plasmid PTU statistics")


# PTU stats
df = read.csv(paste0(dataDir, 'redondo-salvo-2020-plasmid-DB.csv'), header=T, row.names = 1)
# Group into PTUs
ptu.df = df %>% group_by(PTU, PTU.hostrange) %>% 
  summarise(n=length(PTU.hostrange),
            median.size=median(Size))
# Plasmids that are in PTUs with host range >= II
sum(ptu.df$n[which(! ptu.df$PTU.hostrange %in% c("-", "I"))])/sum(ptu.df$n) * 100
# Proportion of PTUs with host range >= II
hostrange.ptu.table = table(ptu.df$PTU.hostrange[which(ptu.df$PTU!="-")])
sum(hostrange.ptu.table[-1])/sum(hostrange.ptu.table) * 100

# We can compare these to statistics from the Acman dataset. 
acman.df = read.csv(paste0(dataDir, 'acman-plasmid-DB.csv'), header=T, row.names = 1)
acman.clique.df = acman.df %>% group_by(Clique) %>% 
  summarise(n=length(Clique),
            median.size=median(Length),
            species=length(unique(Species)))
table(acman.clique.df[which(!is.na(acman.clique.df$Clique)),]$species>1)
table(acman.clique.df[which(acman.clique.df$n>3),"species"]==1)

