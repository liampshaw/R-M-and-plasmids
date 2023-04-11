# PALINDROME ANALYSIS: K=4 and K=6
# USING MCMCGLMM

source('setup.R')

source('01_modelling_functions.R')

mean.10000.k6 = runModelMakePlot(subsampling_value  = "10000", K=6)
mean.50000.k6 = runModelMakePlot(subsampling_value = "50000", K=6)
mean.100000.k6 = runModelMakePlot(subsampling_value = "100000", K=6)

mean.10000.k4 = runModelMakePlot(subsampling_value  = "10000", K=4)
mean.50000.k4  = runModelMakePlot(subsampling_value = "50000", K=4)
mean.100000.k4 = runModelMakePlot(subsampling_value = "100000", K=4)


corr.df = data.frame(species=mean.50000.k6$species,
                     k6.plasmid=mean.50000.k6$p, 
                    k6.core=mean.50000.k6$c,
                    k4.plasmid=mean.50000.k4$p, 
                    k4.core=mean.50000.k4$c)
p.corr.2 = ggplot(corr.df, aes(k4.plasmid, k6.plasmid))+
  geom_point()+
  xlab("Plasmid avoidance (k=4)")+
  ylab("Plasmid avoidance (k=6)")+
  ggtitle("(b) Plasmid genes")+
  theme_bw()+
  geom_hline(yintercept = 0, linetype='dashed')+
    geom_vline(xintercept = 0, linetype='dashed')+
  ggrepel::geom_text_repel(aes(label=species), size=2 )

p.corr.1 = ggplot(corr.df, aes(k4.core, k6.core))+
  geom_point()+
  xlab("Core avoidance (k=4)")+
  ylab("Core avoidance (k=6)")+
    ggtitle("(a) Core genes")+
  theme_bw()+
    geom_hline(yintercept = 0, linetype='dashed')+
    geom_vline(xintercept = 0, linetype='dashed')+
    ggrepel::geom_text_repel(aes(label=species), size=2)
p.corr = cowplot::plot_grid(p.corr.1, p.corr.2, nrow=1)

saveFigure(p.corr, file=paste0("FigureS2_correlation-plasmid-corr"), width=9, height=6.5)

# VARIATION IN PLASMIDS (subplot of Figure 1) 
summaryDF = function(subsampling_value, K){
  main.df = get(paste0("main.df.", K))
  
  df = summariseDFpalindromeCoreAccessory(subsampling_value, K)
  df.plasmid = summariseDFpalindromeCorePlasmid(subsampling_value, K)
  df.plasmid.accessory = summariseDFpalindromeCorePlasmidAccessory(subsampling_value, K)
  summary.df <- main.df[which(main.df$subsampling==subsampling_value & main.df$score!="NaN" &
                                main.df$kmer_category=="Palindromic" & 
                                main.df$section %in% c("accessory_plas", "core", "accessory_chrom")),] %>% 
    group_by(genome, species, kmer_category) %>% 
    summarise(n=length(section),
              score_diff=score[section=="accessory_plas"]-score[section=="core"],
              plasmid_score=score[section=="accessory_plas"],
              chrom_score=score[section=="core"],
              accessory_score=score[section=="accessory_chrom"]) %>%
    filter(n==3) 
  #summary.df$animal <- summary.df$species
  # make into a data frame, with a numeric variable for kmer_category
  summary.df$kmer_category_ordinal <- as.numeric(summary.df$kmer_category)
  summary.df <- as.data.frame(summary.df)
  
  # Abundant species
  species.counts = table(summary.df$species)
  include.species = names(species.counts)[species.counts>N.MIN]
  
  # Remove non-abundant species from what follows
  summary.df = summary.df[which(summary.df$species %in% include.species),]
  return(summary.df)
}

summary.df.4 = summaryDF("50000", "4")
summary.df.6 = summaryDF("50000", "6")

# Summarise plasmid variation
summary.df.6.sd = summary.df.6 %>% 
  group_by(species) %>% 
  summarise(Plasmid=sd(plasmid_score)**2, Core=sd(chrom_score)**2, `Non-core`=sd(accessory_score)**2) %>%
  pivot_longer(cols=c("Plasmid", "Core", "Non-core"), names_to="component")
summary.df.6.sd$species = gsub("_", " ", summary.df.6.sd$species)
p.6.variation = ggplot(summary.df.6.sd, aes(component, value, colour=component))+
  geom_quasirandom()+
  scale_colour_manual(values=component.colours)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Within-species variation")+
  theme(legend.position = "none")+
  xlab("")+
  theme(axis.text=element_text(colour="black"))+
      ggsignif::geom_signif(colour="black", test="wilcox.test", test.args=list(paired=TRUE),
                          comparisons=list(c("Core", "Plasmid"),
                                           c("Non-core", "Plasmid")),
                          step=0.1, 
                          textsize=3.5)+
  theme(axis.title.y=element_text(size=18),
        axis.text.y=element_text(size=14),
        axis.text.x=element_text(size=14),
        panel.border = element_rect(size=2))
saveFigure(p.6.variation, file=paste0("Figure1inset_variation-pangenome-components"), width=4, height=4)





