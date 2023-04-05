#  "Plasmid size in species and R-M presence"

source('setup.R')

# READ IN CONTIG DB
contigs.db = read.csv(paste0(dataDir, "2022-11-18-contig-database.csv"), 
                      header=T,
                      stringsAsFactors = F,
                      row.names = 1)
# Get plasmid stats
plasmid.db = contigs.db[which(contigs.db$type=="plasmid"),]
plasmid.db$species = gsub("\\..*", "", plasmid.db$genome)
plasmid.by.species = plasmid.db %>% group_by(species) %>% 
  summarise(median=median(size),
            max=max(size),
            min=min(size),
            mean=mean(size))

# GET K-MER DB FOR k=6
db.k.6 = read.csv(paste0(dataDir,"6-kmer-db.csv"), header=T, stringsAsFactors = F, row.names = 1)
species.sums = colSums(db.k.6)

# COMBINE TOGETHER
species.names = sort(colnames(db.k.6))
species.one = sapply(gsub("_.*", "", species.names), function(x) substr(x, 1, 2))
species.two = sapply(gsub(".*_", "", species.names), function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2,2)))
species.short = paste0(species.one,species.two)
names(species.short) = species.names
names(species.sums) = species.short

plasmid.by.species$targets.k6 = species.sums[plasmid.by.species$species]
ggplot(plasmid.by.species, aes(targets.k6, max))+
  geom_point()


# RMS PER GENOME
rms = read.csv(paste0(dataDir, "2022-12-06-all-RMS.csv"),
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
# Summarise
n.with.rms = rms %>% group_by(species) %>% 
  summarise(n=length(species))
n.with.rms = data.frame(n.with.rms)
rownames(n.with.rms) = n.with.rms$species

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
plasmid.by.species$n.rms = n.with.rms[plasmid.by.species$species,"n"]

# ANALAYS: Look at them by sorting the maximum:
plasmid.by.species$mean.rms.per.genome = plasmid.by.species$n.rms/plasmid.by.species$n.genomes
plasmid.by.species = plasmid.by.species[order(plasmid.by.species$mean.rms.per.genome, decreasing = TRUE),]
ggplot(plasmid.by.species, aes(mean.rms.per.genome, max))+
  geom_point()+
  scale_y_log10()
ggplot(plasmid.by.species, aes(mean.rms.per.genome, mean))+
  geom_point()+
  scale_y_log10()
ggplot(plasmid.by.species, aes(mean.rms.per.genome, median))+
  geom_point()+
  scale_y_log10()
