#  "Plasmid size in species and R-M presence"

source('setup.R')

SPECIES.OF.INTEREST = c("Neisseria gonorrhoeae", 
                        "Helicobacter pylori")

# GET K-MER DB FOR k=6
db.k.6 = read.csv(paste0(dataDir,"6-kmer-db.csv"), header=T, stringsAsFactors = F, row.names = 1)
species.sums = colSums(db.k.6)

species.names = sort(colnames(db.k.6))
species.one = sapply(gsub("_.*", "", species.names), function(x) substr(x, 1, 2))
species.two = sapply(gsub(".*_", "", species.names), function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2,2)))
species.short = paste0(species.one,species.two)
names(species.short) = species.names
names(species.sums) = species.short

species.long = gsub("_", " ", names(species.short))
names(species.long) = species.short

# READ IN CONTIG DB
contigs.db = read.csv(paste0(dataDir, "2022-11-18-contig-database.csv"), 
                      header=T,
                      stringsAsFactors = F,
                      row.names = 1)
contigs.db$species = gsub("\\..*", "", contigs.db$genome)

species.total.bases = contigs.db %>% group_by(species) %>%
  summarise(total_bases=sum(size),
            chrom_bases=sum(size[type=="chromosome"]),
            plas_bases=sum(size[type=="plasmid"]),
            n_genomes=length(size[type=="chromosome"])) %>%
  mutate(plas_proportion_overall=plas_bases/chrom_bases,
         mean_plas_bases=plas_bases/n_genomes)
rownames(species.total.bases) = species.total.bases$species
# Get plasmid stats
plasmid.db = contigs.db[which(contigs.db$type=="plasmid"),]
# Exclude the secondary chromosome species
EXCLUDE = c("Vibrio parahaemolyticus", 
            "Vibrio cholerae", 
            "Burkholderia pseudomallei")
plasmid.db$species.long = species.long[plasmid.db$species]
plasmid.db = plasmid.db[which(!plasmid.db$species.long %in% EXCLUDE),]

plasmid.by.species = plasmid.db %>% group_by(species, species.long) %>% 
  summarise(median=median(size),
            max=max(size),
            min=min(size),
            mean=mean(size),
            n.plasmids.total=length(size),
            total.plasmid.bases=sum(size))



# RMS PER GENOME
rms = read.csv(paste0(dataDir, "2023-04-06-rms-db.csv"),
               header=T, stringsAsFactors = F)
rms$genome =gsub("\\.[^\\.]+$", "", rms$contig)

rms$genome2 = gsub("\\..*\\.", ".", gsub(".prt", "", rms$genome))
#rms$contig.type = sapply(gsub(".*\\.", "", rms$contig), 
#                         function(x) ifelse(x=="0001", "chromosome", "plasmid"))

rms$species = gsub("\\..*", "", rms$genome)
rms$genome = rms$genome2
rms$genome2 = NULL


rms$contig2 = gsub("\\.[0-9][0-9][0-9][0-9]\\.", ".", rms$contig)
rms$contig = rms$contig2
rms$contig2 = NULL

# Reduce to only R-M systems with different targets
# to avoid overcounting issue
# Summarise
rms$genome_and_target = paste(rms$genome, rms$sequence)
rms = rms[order(rms$genome_and_target),]
rms.unique = rms[!duplicated(rms$genome_and_target),]
n.with.rms = rms.unique %>% group_by(species) %>% 
  summarise(n.rms.unique.target=length(species), 
            n.genomes.with.rms=length(unique(genome)),)
n.with.rms = data.frame(n.with.rms)
rownames(n.with.rms) = n.with.rms$species

# n.rms.unique.target is the number of R-M systems in total within a species,
# not counting any additional R-M systems within a genome which recognise 
# the same target as another one
# n.genomes.with.rms is the number of genomes within a species that 
# carry at least one R-M system

# Pangenome stats
pangenome.stats = read.csv(paste0(dataDir, "pangenome-stats.csv"),
                           header=T,
                           stringsAsFactors=F)
pangenome.stats$species.short = species.short[pangenome.stats$species]
# get rid of e. faecium/faecalis
pangenome.stats = pangenome.stats[!duplicated(pangenome.stats$species.short),]
pangenome.stats = na.omit(pangenome.stats)
rownames(pangenome.stats) = pangenome.stats$species.short

# Add to plasmid stats
plasmid.by.species$n.genomes = pangenome.stats[plasmid.by.species$species,"genomes"]
plasmid.by.species$n.rms.unique.target = n.with.rms[plasmid.by.species$species,"n.rms.unique.target"]
plasmid.by.species$n.genomes.with.rms = n.with.rms[plasmid.by.species$species,"n.genomes.with.rms"]
plasmid.by.species$mean.plasmid.bases = species.total.bases[plasmid.by.species$species,]$mean_plas_bases
plasmid.by.species$prop.with.rms = plasmid.by.species$n.genomes.with.rms/plasmid.by.species$n.genomes
# ANALAYS: Look at them by sorting the maximum:
plasmid.by.species$mean.rms.unique.target.per.genome = plasmid.by.species$n.rms.unique.target/plasmid.by.species$n.genomes

plasmid.by.species = plasmid.by.species[order(plasmid.by.species$mean.rms.unique.target.per.genome, decreasing = TRUE),]

p.mean.rms.per.genome = ggplot(plasmid.by.species, aes(mean.rms.unique.target.per.genome, max))+
  geom_point()+
  scale_y_log10()+
  theme_bw()+
  ggrepel::geom_text_repel(aes(label=species.long),
                           data=plasmid.by.species[which(plasmid.by.species$mean.rms.unique.target.per.genome>6),],
                           size=3)+
  theme(panel.grid = element_blank())+
  ylab("Largest plasmid observed in species (bp)")+
  xlab("Mean R-M systems per genome (unique targets)")

p.mean.bases.plasmid = ggplot(plasmid.by.species, aes(mean.rms.unique.target.per.genome, mean.plasmid.bases))+
  geom_point()+
  scale_y_log10()+
  theme_bw()+
  ggrepel::geom_text_repel(aes(label=species.long),
                           data=plasmid.by.species[which(plasmid.by.species$mean.rms.unique.target.per.genome>6),],
                           size=3)+
  theme(panel.grid = element_blank())+
  ylab("Mean plasmid bases per genome (bp)")+
  xlab("Mean R-M systems per genome (unique targets)")
p.combined = cowplot::plot_grid(p.mean.rms.per.genome+ggtitle("(a)"),
                   p.mean.bases.plasmid+ggtitle("(b)"))
saveFigure(p.combined, "FigureS9", width=10, height=5)
pggplot(plasmid.by.species, aes(mean.rms.per.genome, mean))+
  geom_point()+
  scale_y_log10()
ggplot(plasmid.by.species, aes(mean.rms.per.genome, median))+
  geom_point()+
  scale_y_log10()
