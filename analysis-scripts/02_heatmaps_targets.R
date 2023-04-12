# "Trees and heatmaps"

source('setup.R')


plotHeatmap <- function(K){
  
  # Database of k-mer presence/absence
  db.6mers = read.csv(paste0(dataDir, K, '-kmer-db.csv'),
                      header=T, row.names = 1)
  db.6mers = db.6mers[,species.tree$tip.label]
  
  db.6mers.nonzero = t(db.6mers[which(rowSums(db.6mers)!=0),])
  db.hclust = hclust(dist(t(db.6mers.nonzero)))
  
  # Store if they are unique in the dataset?
  unique.names = colnames(db.6mers.nonzero)[colSums(db.6mers.nonzero)==1]
  
  db.6mers.nonzero.rel.counts = apply(db.6mers.nonzero, MARGIN=2, 
                                      function(x) x/sum(x)) # store relative proportion of counts in each species
  #db.6mers.nonzero.rel.counts = apply(db.6mers.nonzero.rel.counts, MARGIN=1:2,
  #                                    function(x) ifelse(x==1, "unique", ifelse(x!=0, "common", "absent")))
  
  species.tree$tip.label = gsub("_", " ", species.tree$tip.label)
  
  rownames(db.6mers.nonzero) = gsub("_", " ", rownames(db.6mers.nonzero))
  #db.6mers.nonzero.rel.counts
  p = ggtree(species.tree)+
    geom_tiplab(align = TRUE, size=0.9, linesize = 0.1,fontface=3 , linetype=NULL)
  g <- gheatmap(p,
                db.6mers.nonzero, 
                colnames_angle=90, 
                hjust=1, 
                font.size = 0.6, 
                offset = 0.3, 
                width=1.2*ncol(db.6mers.nonzero)/128, 
                family = "mono")+
    scale_fill_continuous(low="white", high= "black")+
    theme(legend.position = "none")
  g = g  + ggtree::vexpand(.1, -1)
  return(g)
}


species.tree <- read.tree(paste0(dataDir, '16S_species.nwk'))

db.4mers = read.csv(paste0(dataDir, '4-kmer-db.csv'),
                    header=T, row.names = 1)
db.4mers = db.4mers[,species.tree$tip.label]
db.4mers.nonzero = t(db.4mers[which(rowSums(db.4mers)!=0),])


db.5mers =  read.csv(paste0(dataDir, '5-kmer-db.csv'),
                     header=T, row.names = 1)
db.5mers = db.5mers[,species.tree$tip.label]
db.5mers.nonzero = t(db.5mers[which(rowSums(db.5mers)!=0),])

db.6mers =  read.csv(paste0(dataDir, '6-kmer-db.csv'),
                     header=T, row.names = 1)
db.6mers = db.6mers[,species.tree$tip.label]
db.6mers.nonzero = t(db.6mers[which(rowSums(db.6mers)!=0),])

#db.palindromes = rbind(db.4mers[palindromes.4,],db.6mers[palindromes,])

#db.palindromes = db.4mers[palindromes.4,]
#db.palindromes.t = t(db.palindromes)

# Order colnames by clustering?
heatmap.4 = plotHeatmap("4")
heatmap.5 = plotHeatmap("5")
heatmap.6 = plotHeatmap("6")

plot.widths = c(ncol(db.4mers.nonzero), 
                ncol(db.5mers.nonzero), 
                ncol(db.6mers.nonzero))
plot.widths/sum(plot.widths)
p.heatmap.4 = heatmap.4+ggtitle("k=4")
p.heatmap.5 = heatmap.5+ggtitle("k=5")
p.heatmap.6 = heatmap.6+ggtitle("k=6")


species.tree$tip.label = gsub("_", " ", species.tree$tip.label)
p.tree = ggtree(species.tree)+
    geom_tiplab(align = TRUE, size=0.9, linesize = 0.1,fontface=3, linetype=NULL)+
  xlim(c(0,0.5))+
  geom_treescale(width=0.1, x = 0, y=-2)

saveFigure(p.tree, 
       'FigureS1_species-tree', width=6, height=5)

saveFigure(p.heatmap.4, 
  'FigureS4_heatmap-k4-targets', width=8, height=4)
saveFigure(p.heatmap.5, 
       'FigureS5_heatmap-k5-targets', width=8, height=4)
saveFigure(p.heatmap.6, 
       'FigureS6_heatmap-k6-targets', width=8, height=4)

