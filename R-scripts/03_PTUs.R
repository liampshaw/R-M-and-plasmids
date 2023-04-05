# Plasmid taxonomic units (PTUs)"

source('setup.R')


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
