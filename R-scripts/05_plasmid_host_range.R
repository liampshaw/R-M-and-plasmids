# Look at avoidance relative to core genes
# And plot by taxonomic distribution
source('setup.R')

plasmid_db = read.csv(paste0(dataDir, 'redondo-salvo-2020-plasmid-DB.csv'), header=T, row.names = 1)
# Group into PTUs


coreGeneAvoidance <- function(K, SUBSAMPLING){
  main.df <- read.csv(paste0(dataDir, K, '-merged-genome-results-levels-inclusive.csv'), header=F)
  
  colnames(main.df) <- c("species", "kmer_category", "genome", "section", "subsampling", "n", "rank", "score")
  
  main.df$section <- ordered(main.df$section, 
                             levels=c("core", "accessory_chrom", "accessory_plas"))
  
  main.df$kmer_category <- ordered(main.df$kmer_category,
                                   levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"))
  
  
  summary.df.all <- main.df[which(main.df$subsampling==SUBSAMPLING & main.df$score!="NaN" & 
                                    main.df$section=="core"),] %>% 
    group_by(species, kmer_category, section) %>% 
    summarise(score=median(score),
              rank=median(rank))
  summary.scores = summary.df.all$score
  names(summary.scores) = paste0(summary.df.all$species, summary.df.all$kmer_category)
  return(summary.scores)
}

plasmidAvoidanceDF <- function(K, SUBSAMPLING){
  scores.k = coreGeneAvoidance(K, SUBSAMPLING)
  plasmid.df = read.csv(paste0(dataDir, 'plasmids/redondo-salvo-summary-k-', K, '-n-', prettyNum(SUBSAMPLING),'-True.csv'),
                        stringsAsFactors = FALSE,
                        header=T)
  
  #plasmid.df = plasmid.df[,]
  plasmid.df$Size = plasmid_db[plasmid.df$plasmid, "Size"]
  plasmid.df$PTU = plasmid_db[plasmid.df$plasmid, "PTU"]
  plasmid.df$species.rs = gsub(" ", "_", plasmid_db[plasmid.df$plasmid, "TaxSpecies"])
  plasmid.df$hostrange = plasmid_db[plasmid.df$plasmid, "PTU.hostrange"]
  
  #plasmid.df.modelling = plasmid.df[which(plasmid.df$PTU!="-"),]
  plasmid.df.modelling = plasmid.df
  plasmid.df.modelling$hostrange = ordered(plasmid.df.modelling$hostrange,
                                           levels=c("-", "I", "II", "III", "IV", "V", "VI"))
  plasmid.df.modelling$hostrange.numeric = as.numeric(plasmid.df.modelling$hostrange)
  # Remove NAs
  plasmid.df.modelling = na.omit(plasmid.df.modelling)
  plasmid.df.modelling$relative.score = sapply(1:nrow(plasmid.df.modelling),
                                               function(i) plasmid.df.modelling$score[i]-scores.k[paste0(plasmid.df.modelling$species[i], plasmid.df.modelling$kmer_category[i])])
  return(plasmid.df.modelling)
}


# Make dataframe
plasmid.df.modelling = plasmidAvoidanceDF(6, "10000")

# Merge by PTU
ptu.df.plot = plasmid.df.modelling %>% group_by(PTU, species, kmer_category, hostrange) %>%
  summarise(score=mean(score))
ptu.df.plot$hostrange = ordered(ptu.df.plot$hostrange,
                                levels=c("-", "I", "II", "III", "IV", "V", "VI"))
ptu.df.plot$kmer_category = ordered(ptu.df.plot$kmer_category,
                                    levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"),
                                    labels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"))
ptu.df.plot$host.range.numeric = as.numeric(ptu.df.plot$hostrange)


# add fake thing to track unassigned plasmids
plasmid.df.modelling$PTU.fake = plasmid.df.modelling$PTU
plasmid.df.modelling$PTU.fake[which(plasmid.df.modelling$PTU=="-")] = plasmid.df.modelling$plasmid[which(plasmid.df.modelling$PTU=="-")]


ptu.df.plot.palindromic = plasmid.df.modelling %>% filter(kmer_category=="Palindromic") %>% group_by(PTU.fake, PTU, hostrange) %>%
  summarise(score=mean(score), size=median(Size))
ptu.df.plot.palindromic$hostrange = ordered(ptu.df.plot.palindromic$hostrange, 
                                            levels=c("-", "I", "II", "III", "IV", "V", "VI"),
                                            labels=c("Unassigned\nplasmids", "I", "II", "III", "IV", "V", "VI"))

#ptu.df.plot.palindromic$PTU.fake = ptu.df.plot.palindromic$PTU 
#ptu.df.plot.palindromic$PTU.fake[which(ptu.df.plot.palindromic$PTU=="-")] = ptu.df.plot.palindromic$[which(ptu.df.plot.palindromic$PTU=="-")]

p.all.PTUs.palindromic = ggplot(ptu.df.plot.palindromic, aes(hostrange, score, colour=hostrange))+
  geom_hline(yintercept = 0, linetype='dashed', colour='black')+
  geom_quasirandom( width=0.25, size=2)+
  stat_summary(fun = median, size=4, shape="-", colour="black")+
  theme_bw()+scale_x_discrete(drop=FALSE)+theme(panel.grid=element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("PTU host range")+
  ylab("Mean exceptionality score of palindromes")+
  geom_vline(xintercept = 1.5, linetype='dashed')+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"))+
  scale_colour_manual(values=c("grey", RColorBrewer::brewer.pal(6, "Set2")))+
  theme(legend.position = "none")
saveFigure(p.all.PTUs.palindromic, "Figure4a_ptus-palindromic-k6-hostrange", width=5.5, height=4)


# Compute a statistical test of correlation
# Select only those without missing PTU
df.for.correlation = ptu.df.plot.palindromic[which(ptu.df.plot.palindromic$PTU!="-"),]
print("Spearman correlation of mean avoidance with PTU host range (excluding unassigned plasmids)")
cor.test(as.numeric(df.for.correlation$hostrange), df.for.correlation$score, method="spearman")
# Rho=-0.26, p=0.003

ptu.df.plot.palindromic = ptu.df.plot.palindromic %>%   mutate(new_bin = cut(log10(size), 
                                                   breaks=c(2, 3, 4, 4.39794, 4.69897, 5, 10),
                                                   labels=c( "<1kb", "1-10kb", "10-25kb", "25-50kb", "50-100kb", ">100kb")))
p.all.PTUs.palindromic.size = ggplot(ptu.df.plot.palindromic, aes(new_bin, score, colour=hostrange))+
  geom_hline(yintercept = 0, linetype='dashed', colour='black')+
  geom_quasirandom( width=0.25, size=2)+
  stat_summary(fun = median, size=4, shape="-", colour="black")+
  theme_bw()+scale_x_discrete(drop=FALSE)+theme(panel.grid=element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("PTU host range")+
  ylab("Mean exceptionality score of palindromes")+
  geom_vline(xintercept = 1.5, linetype='dashed')+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"))+
  scale_colour_manual(values=c("grey", RColorBrewer::brewer.pal(6, "Set2")))+
  theme(legend.position = "none")
saveFigure(p.all.PTUs.palindromic.size, "Figure4b_ptus-palindromic-k6-hostrange-size", width=5.5, height=4)


# REPEAT FOR K=4
# Make dataframe
plasmid.df.modelling = plasmidAvoidanceDF(4, "2500")

# Merge by PTU
ptu.df.plot = plasmid.df.modelling %>% group_by(PTU, species, kmer_category, hostrange) %>%
  summarise(score=mean(score))
ptu.df.plot$hostrange = ordered(ptu.df.plot$hostrange,
                                levels=c("-", "I", "II", "III", "IV", "V", "VI"))
ptu.df.plot$kmer_category = ordered(ptu.df.plot$kmer_category,
                                    levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"),
                                    labels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"))
ptu.df.plot$host.range.numeric = as.numeric(ptu.df.plot$hostrange)


# add fake thing to track unassigned plasmids
plasmid.df.modelling$PTU.fake = plasmid.df.modelling$PTU
plasmid.df.modelling$PTU.fake[which(plasmid.df.modelling$PTU=="-")] = plasmid.df.modelling$plasmid[which(plasmid.df.modelling$PTU=="-")]


ptu.df.plot.palindromic = plasmid.df.modelling %>% filter(kmer_category=="Palindromic") %>% group_by(PTU.fake, PTU, hostrange) %>%
  summarise(score=mean(score), size=median(Size))
ptu.df.plot.palindromic$hostrange = ordered(ptu.df.plot.palindromic$hostrange, 
                                            levels=c("-", "I", "II", "III", "IV", "V", "VI"),
                                            labels=c("Unassigned\nplasmids", "I", "II", "III", "IV", "V", "VI"))

#ptu.df.plot.palindromic$PTU.fake = ptu.df.plot.palindromic$PTU 
#ptu.df.plot.palindromic$PTU.fake[which(ptu.df.plot.palindromic$PTU=="-")] = ptu.df.plot.palindromic$[which(ptu.df.plot.palindromic$PTU=="-")]

p.all.PTUs.palindromic = ggplot(ptu.df.plot.palindromic, aes(hostrange, score, colour=hostrange))+
  geom_hline(yintercept = 0, linetype='dashed', colour='black')+
  geom_quasirandom( width=0.25, size=2)+
  stat_summary(fun = median, size=4, shape="-", colour="black")+
  theme_bw()+scale_x_discrete(drop=FALSE)+theme(panel.grid=element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("PTU host range")+
  ylab("Mean exceptionality score of palindromes")+
  geom_vline(xintercept = 1.5, linetype='dashed')+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"))+
  scale_colour_manual(values=c("grey", RColorBrewer::brewer.pal(6, "Set2")))+
  theme(legend.position = "none")
saveFigure(p.all.PTUs.palindromic, "FigureX_ptus-palindromic-k4-hostrange", width=5.5, height=4)

ptu.df.plot.palindromic = ptu.df.plot.palindromic %>%   mutate(new_bin = cut(log10(size), 
                                                                             breaks=c(2, 3, 4, 4.39794, 4.69897, 5, 10),
                                                                             labels=c( "<1kb", "1-10kb", "10-25kb", "25-50kb", "50-100kb", ">100kb")))

p.all.PTUs.palindromic.size = ggplot(ptu.df.plot.palindromic, aes(new_bin, score, colour=hostrange))+
  geom_hline(yintercept = 0, linetype='dashed', colour='black')+
  geom_quasirandom( width=0.25, size=2)+
  stat_summary(fun = median, size=4, shape="-", colour="black")+
  theme_bw()+scale_x_discrete(drop=FALSE)+theme(panel.grid=element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("PTU host range")+
  ylab("Mean exceptionality score of palindromes")+
  geom_vline(xintercept = 1.5, linetype='dashed')+
  theme(panel.border = element_blank(),
        axis.line=element_line(colour="black"))+
  scale_colour_manual(values=c("grey", RColorBrewer::brewer.pal(6, "Set2")))+
  theme(legend.position = "none")+
  facet_wrap(~hostrange, ncol=2)
saveFigure(p.all.PTUs.palindromic.size, "FigureX_ptus-palindromic-k4-size", width=5.5, height=4)
