
source('setup.R')

# LOAD DATA
print("Loading data")
counts_4mers = read.csv(paste0(dataDir, '4mer_occurrences_plasmids.csv'), sep=' ', header=T)
normalising = rowSums(counts_4mers)
counts_4mers_normalised = counts_4mers/normalising


counts_5mers = read.csv(paste0(dataDir, '5mer_occurrences_plasmids.csv'), sep=' ', header=T)
normalising = rowSums(counts_5mers)
counts_5mers_normalised = counts_5mers/normalising


counts_6mers = read.csv(paste0(dataDir, '6mer_occurrences_plasmids.csv'), sep=' ', header=T)
normalising = rowSums(counts_6mers)
counts_6mers_normalised = counts_6mers/normalising



genomes = read.csv('../../old-analysis/2022-11-15-all-filtered-genomes.csv', 
                   row.names = 1, 
                   header=T, stringsAsFactors = F)
kmer_dict_4 = read.csv(paste0(dataDir, '4-kmer-db.csv'),
                     header=T, stringsAsFactors = F, row.names = 1)
kmer_dict_5 = read.csv(paste0(dataDir, '5-kmer-db.csv'),
                     header=T, stringsAsFactors = F, row.names = 1)
kmer_dict_6 = read.csv(paste0(dataDir, '6-kmer-db.csv'),
                     header=T, stringsAsFactors = F, row.names = 1)

# PLASMID PLOT FUNCTION
plotForSpecies = function(SPECIES, K){
  species_targets = rownames(get(paste0("kmer_dict_", K)))[which(get(paste0("kmer_dict_", K))[,SPECIES]!=0)]
  species_short = paste0(substr(gsub("_.*", "", SPECIES), 1, 2),
                         toupper(substr(gsub(".*_", "", SPECIES), 1, 1)),
                         substr(gsub(".*_", "", SPECIES), 2, 2))
  species_kmers = get(paste0("counts_", K, "mers_normalised"))[grepl(species_short, rownames(get(paste0("counts_", K, "mers_normalised")))),]
  
  species_kmers$genome = rownames(species_kmers)
  species_kmers$size = rowSums(get(paste0("counts_", K, "mers"))[grepl(species_short, rownames(get(paste0("counts_", K, "mers")))),])
  
  species_kmers_long = species_kmers %>% pivot_longer(cols=colnames(species_kmers)[1:(ncol(species_kmers)-2)])
  # If investigating beyond-species level:
  #species_kmers_long$type = sapply(species_kmers_long$name,
  #                               function(x) 
  #                                 ifelse(x %in% species_kmers,
  #                               "Within-species\nR-M target", ifelse(x %in% all_kmers, "R-M target from beyond-species level", "Other k-mer")))
  species_kmers_long$type = ifelse(species_kmers_long$name %in% species_targets, "Within-species\nR-M targets", "Others")
  species_kmers_long$type = ordered(species_kmers_long$type, 
                                    levels=c("Within-species\nR-M targets", "Others"))
  species_kmers_long_medians = species_kmers_long %>% group_by(genome, type, size) %>%
    summarise(value=mean(value), .groups = "drop")
  species_kmers_long_medians = species_kmers_long_medians %>% 
    mutate(size_bin=cut(size, breaks=c(-1, 10000, 50000, 100000, 5e10),
                        labels=c("<10", "10-50", "50-100", ">100")))
  
  p = ggplot(species_kmers_long_medians[which(species_kmers_long_medians$type=="Within-species\nR-M targets"),], aes(size_bin, value, colour=size_bin))+
    #geom_hline(yintercept = 1/(4**K), linetype='dashed')+
    scale_color_manual(values=size.colours)+
    geom_quasirandom(alpha=0.4, width=0.25)+
        stat_summary(fun = median, size=2, shape="-", colour="black")+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"), 
          axis.text=element_text(colour="black"))+
    xlab("Plasmid size (kb)")+
    ylab("Mean target density (per bp)")+
    ggtitle(gsub("_", " ", SPECIES))+
    theme(axis.title=element_text(face="plain"))+
    theme(legend.position = "none")
  return(p)
}

# PLOTTING FOR SPECIES
#p.ecoli.4 = plotForSpecies("Escherichia_coli", 4)+ggtitle("(a) k=4 ()")
p.ecoli.5 = plotForSpecies("Escherichia_coli", 5)+ggtitle("(a) k=5")
p.ecoli.6 = plotForSpecies("Escherichia_coli", 6)+ggtitle("(b) k=6")+ylab("")
p.ecoli = plot_grid(p.ecoli.5, p.ecoli.6, nrow=1)
saveFigure(p.ecoli, "FigureX_plasmids-escherichia-coli",  width=6, height=4)

p.kpneu = plotForSpecies("Klebsiella_pneumoniae", 6)
p.sal = plotForSpecies("Salmonella_enterica", 6)
saveFigure(p.kpneu, "FigureX_plasmids-klebsiella-pneumoniae", width=6, height=4)
saveFigure(p.sal, "FigureX_plasmids-salmonella-enterica", width=6, height=4)

# SUMMARISE FOR SPECIES FUNCTION
summariseForSpecies = function(SPECIES, K){
  species_targets = rownames(get(paste0("kmer_dict_", K)))[which(get(paste0("kmer_dict_", K))[,SPECIES]!=0)]
  species_short = paste0(substr(gsub("_.*", "", SPECIES), 1, 2),
                         toupper(substr(gsub(".*_", "", SPECIES), 1, 1)),
                         substr(gsub(".*_", "", SPECIES), 2, 2))
  species_kmers = get(paste0("counts_", K, "mers_normalised"))[grepl(species_short, rownames(get(paste0("counts_", K, "mers_normalised")))),]
    if (nrow(species_kmers)==0){
    return(NA)
  }
  
  species_kmers$genome = rownames(species_kmers)
  species_kmers$size = rowSums(get(paste0("counts_", K, "mers"))[grepl(species_short, rownames(get(paste0("counts_", K, "mers")))),])
  
  species_kmers_long = species_kmers %>% pivot_longer(cols=colnames(species_kmers)[1:(ncol(species_kmers)-2)])
  # If investigating beyond-species level:
  #species_kmers_long$type = sapply(species_kmers_long$name,
  #                               function(x) 
  #                                 ifelse(x %in% species_kmers,
  #                               "Within-species\nR-M target", ifelse(x %in% all_kmers, "R-M target from beyond-species level", "Other k-mer")))
  species_kmers_long$type = ifelse(species_kmers_long$name %in% species_targets, "Within-species\nR-M targets", "Others")
  species_kmers_long$type = ordered(species_kmers_long$type, 
                                    levels=c("Within-species\nR-M targets", "Others"))

  species_kmers_long_medians = species_kmers_long %>% group_by(genome, type, size) %>%
    summarise(value=mean(value), .groups="drop")
  species_kmers_long_medians = species_kmers_long_medians %>% 
    mutate(size_bin=cut(size, breaks=c(-1, 10000, 50000, 100000, 5e10),
                        labels=c("<10", "10-50", "50-100", ">100")))



  return.df = species_kmers_long_medians[which(species_kmers_long_medians$type=="Within-species\nR-M targets"),]
  if (nrow(return.df)==0){
    return(NA)
  } 
  return.df = data.frame(return.df %>% group_by(size_bin) %>%
    summarise(mean=mean(value),
              se=sd(value)/sqrt(length(value)),
              n=length(value), .groups = "drop"))
  return.df$species = SPECIES
  return(return.df)
}

print("Combining results for species with plasmids...")
species.with.plasmids =  names(table(genomes[which(genomes$nb_conts>1),"species"]))
species.with.plasmids = species.with.plasmids[which(species.with.plasmids!="")]
#species.with.rms = colnames(kmer_dict_6)[colSums(kmer_dict_6)!=0]
#species.with.both = species.with.plasmids[species.with.plasmids %in% species.with.rms]
#species.with.both = species.with.both[-which(species.with.both=="Mycobacterium_intracellulare")]
combined.results.4 = matrix(nrow=0, ncol=5)
combined.results.5 = matrix(nrow=0, ncol=5)
combined.results.6 = matrix(nrow=0, ncol=5)

for (species in species.with.plasmids){
  print(species)
  results.4 = summariseForSpecies(species, 4)
  results.5 = summariseForSpecies(species, 5)
  results.6 = summariseForSpecies(species, 6)
  if (!(is.na(results.4))){
      combined.results.4 = rbind(combined.results.4, results.4)
  }
    if (!(is.na(results.5))){
      combined.results.5 = rbind(combined.results.5, results.5)
    }
     if (!(is.na(results.6))){
      combined.results.6 = rbind(combined.results.6, results.6)
  }

}


print("Plotting results for those with >N.MIN results...")
N.MIN = 5
p.combined.4 = ggplot(combined.results.4[which(combined.results.4$n>N.MIN),], aes(size_bin, mean, colour=size_bin))+
    geom_hline(yintercept = 1/(4**4), linetype='dashed')+
    scale_color_manual(values=size.colours)+
geom_quasirandom(width=0.25)+    theme_bw()+        stat_summary(fun = median, size=3, shape="-", colour="black")+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"), 
          axis.text=element_text(colour="black"))+
    xlab("Plasmid size (kb)")+
    ylab("Mean target density (per bp)")+
    theme(axis.title=element_text(face="plain"))+
    theme(legend.position = "none")+ggsignif::geom_signif(comparisons=list( c("<10", ">100")),step=0.05 ,textsize = 2, colour="black")
p.combined.5 = ggplot(combined.results.5[which(combined.results.5$n>N.MIN),], aes(size_bin, mean, colour=size_bin))+
    geom_hline(yintercept = 1/(4**5), linetype='dashed')+
    scale_color_manual(values=size.colours)+
geom_quasirandom(width=0.25)+    theme_bw()+        stat_summary(fun = median, size=3, shape="-", colour="black")+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"), 
          axis.text=element_text(colour="black"))+
    xlab("Plasmid size (kb)")+
    ylab("Mean target density (per bp)")+
    theme(axis.title=element_text(face="plain"))+
    theme(legend.position = "none")+ggsignif::geom_signif(comparisons=list(c("<10", ">100")),step=0.05,textsize = 2, colour="black")
p.combined.6 = ggplot(combined.results.6[which(combined.results.6$n>N.MIN),], aes(size_bin, mean, colour=size_bin))+
    geom_hline(yintercept = 1/(4**6), linetype='dashed')+
    scale_color_manual(values=size.colours)+
geom_quasirandom(width=0.25)+    theme_bw()+        stat_summary(fun = mean, size=2, shape="-", colour="black")+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"), 
          axis.text=element_text(colour="black"))+
    xlab("Plasmid size (kb)")+
    ylab("Mean target density (per bp)")+
    theme(axis.title=element_text(face="plain"))+
    theme(legend.position = "none")+ggsignif::geom_signif(comparisons=list( c("<10", ">100")),step=0.05,textsize = 2, colour="black")
p.combined.boxplot = plot_grid(p.combined.4+ggtitle("(c) k=4"),
                               p.combined.5+ggtitle("(d) k=5"),
                               p.combined.6+ggtitle("(e) k=6"),
                               nrow=3)
saveFigure(p.combined.boxplot, "Figure3cde_plasmid-boxplots", width=3.5, height=9)

# COMBINE THE COMBINED PLOTS FOR EACH VALUE OF K INTO SINGLE PLOT
print("Combining plots for each value of k...")
p.combined.boxplot = plot_grid(p.combined.4+ggtitle("(c) k=4")+ylab("")+theme(axis.text = element_text(size=8)),
                               p.combined.5+ggtitle("(d) k=5")+theme(axis.text = element_text(size=8)),
                               p.combined.6+ggtitle("(e) k=6")+ylab("")+theme(axis.text = element_text(size=8)),
                               nrow=3)
p.ecoli = plot_grid(p.ecoli.5+theme(axis.title=element_text(size=14), 
                                              axis.text=element_text(size=12)), 
                    p.ecoli.6+theme(axis.title=element_text(size=14), 
                                              axis.text=element_text(size=12)), nrow=1)

p.combined.combined = plot_grid(p.ecoli,p.combined.boxplot,
                                nrow=1,
                                align='h', axis='bt',
                                rel_widths = c(1,0.35))
saveFigure(p.combined.combined, "Figure3_plasmid-target-densities", width=10, height=5.2)


# LOOK AT PALINDROMES (see figure 4)
print("Looking at palindromes...")
palindromes.4 = read.csv("../data/4mer_palindromes.txt", header=F, stringsAsFactors = F)$V1
df.4 = data.frame(genome=rownames(counts_4mers_normalised), 
                                  median.score=apply(counts_4mers_normalised[,palindromes.4], 1, median))
df.4$size = rowSums(counts_4mers)
df.4 = df.4  %>% 
    mutate(size_bin=cut(size, breaks=c(-1, 10000, 50000, 100000, 5e10),
                        labels=c("<10", "10-50", "50-100", ">100")))
  

p.palindromes.4 = ggplot(df.4, aes(size_bin, median.score, colour=size_bin))+
  geom_hline(yintercept = 1/(4**4), linetype='dashed')+
    scale_color_manual(values=size.colours)+
geom_quasirandom(width=0.25)+    theme_bw()+       
  stat_summary(fun = mean, size=2, shape="-", colour="black", aes(group=size_bin))+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"), 
          axis.text=element_text(colour="black"))+
    xlab("Plasmid size (kb)")+
    ylab("Mean target density (per bp)")+
    theme(axis.title=element_text(face="plain"))+
    theme(legend.position = "none")+ggsignif::geom_signif(comparisons=list( c("<10", ">100")),step=0.05,textsize = 2, colour="black")

saveFigure(p.palindromes.4, "Figure3_palindrome-version-k-4",
       width=9, height=6.5)


# K=6
palindromes.6 = read.csv("../data/6mer_palindromes.txt", header=F, stringsAsFactors = F)$V1
df.6 = data.frame(genome=rownames(counts_6mers_normalised), 
                                  median.score=apply(counts_6mers_normalised[,palindromes.6], 1, median))
df.6$size = rowSums(counts_6mers)
df.6 = df.6  %>% 
    mutate(size_bin=cut(size, breaks=c(-1, 10000, 50000, 100000, 5e10),
                        labels=c("<10", "10-50", "50-100", ">100")))
  

p.palindromes.6 = ggplot(df.6, aes(size_bin, median.score, colour=size_bin))+
  geom_hline(yintercept = 1/(4**6), linetype='dashed')+
    scale_color_manual(values=size.colours)+
geom_quasirandom(width=0.25)+    theme_bw()+       
  stat_summary(fun = mean, size=2, shape="-", colour="black", aes(group=size_bin))+
    theme(panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_line(colour = "black"), 
          axis.text=element_text(colour="black"))+
    xlab("Plasmid size (kb)")+
    ylab("Mean target density (per bp)")+
    theme(axis.title=element_text(face="plain"))+
    theme(legend.position = "none")+ggsignif::geom_signif(comparisons=list( c("<10", ">100")),step=0.05,textsize = 2, colour="black")
saveFigure(p.palindromes.6, "Figure3_palindrome-version-k6.pdf",
       width=9, height=6.5)