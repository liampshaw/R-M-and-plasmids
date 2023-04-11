suppressPackageStartupMessages(library(MCMCglmm))

# FUNCTIONS
readDF <- function(K){
  main.df <- read.csv(paste0(dataDir, K, '-merged-genome-results-levels-inclusive.csv'), header=F)
  colnames(main.df) <- c("species", "kmer_category", "genome", "section", "subsampling", "n", "rank", "score")
  
  main.df$section <- ordered(main.df$section, 
                             levels=c("core", "accessory_chrom", "accessory_plas"))
  
  main.df$kmer_category <- ordered(main.df$kmer_category,
                                   levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"))
  
  return(main.df) 
}

summariseDFpalindromeCoreAccessory <- function(subsampling, K){
  df = get(paste0("main.df.", K))
  summary.df <- df[which(df$subsampling==subsampling  & df$kmer_category=="Palindromic"),] %>% 
    group_by(genome, species, kmer_category, n) %>% 
    summarise(diff.score=score[section=="core"]-score[section=="accessory_chrom"],
              diff.rank=rank[section=="core"]-rank[section=="accessory_chrom"],
              k=K)
  return(summary.df)
}

summariseDFpalindromeCorePlasmid <- function(subsampling, K){
  df = get(paste0("main.df.", K))
  summary.df <- df[which(df$subsampling==subsampling  & df$kmer_category=="Palindromic"),] %>% 
    group_by(genome, species, kmer_category, n) %>% 
    summarise(diff.score=score[section=="core"]-score[section=="accessory_plas"],
              k=K)
  return(summary.df)
}

summariseDFpalindromeCorePlasmidAccessory <- function(subsampling, K){
  df = get(paste0("main.df.", K))
  summary.df <- df[which(df$subsampling==subsampling  & df$kmer_category=="Palindromic"),] %>% 
    group_by(genome, species, kmer_category, n) %>% 
    summarise(diff.plasmid=score[section=="core"]-score[section=="accessory_plas"],
              diff.accessory=score[section=="core"]-score[section=="accessory_chrom"],
              n=length(genome),
              k=K)
  return(summary.df)
}

# MCMC MODELLING
# MCMC PARAMETERS - equivalent to effective sample N=1000
N_ITERATIONS = 1100000#1100000 #1100000# 1100000
BURNIN = 100000#100000 #1000#100000 #
THIN = 1000 #1000#1000#1000
N.MIN = 2 # N-1 where N is genomes needed for inclusion in plots 

# Read in tree
species.tree <- read.tree(paste0(dataDir, '16S_species.nwk'))
phylo <- force.ultrametric(species.tree)
# Collapse dichotomies into multichotomy
phylo_2 <- di2multi(phylo)
# Add node labels
dataTreeNode <- makeNodeLabel(phylo_2, method = "number")
# Inverse matrix of phylogeny, for use in MCMC modelling
inv.phylo <- inverseA(dataTreeNode)
inv.phylogenetic.matrix <- inv.phylo$Ainv

# MCMC modelling:
# Prior for three random effects
prior3 <- list(G=list(G1=list(V=1,nu=0.002), 
                      G2=list(V=1,nu=0.002), 
                      G3=list(V=1,nu=0.002))) 
prior2 <- list(G=list(G1=list(V=1,nu=0.002), 
                      G2=list(V=1,nu=0.002))) 

# Read in k=4,6 palindrome data
main.df.4 = readDF(4)
main.df.6 = readDF(6)



runModelMakePlot <- function(subsampling_value, K, file_prefix=""){
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
  
  #Conditional R2 (total variance explained by the model)
  
  # P-value of the intercept?
  summary.df.longer = pivot_longer(summary.df, c("plasmid_score", "chrom_score", "accessory_score"))
  
  # MCMC GLMM
  df = data.frame(summary.df.longer)
  df$name.ordered = relevel(as.factor(df$name), ref="chrom_score")
  
  
  # Summarise at level of species
  summary.results.per.species = data.frame(df %>% group_by(species, name) %>%
                                             summarise(value=mean(value),
                                                       n=length(unique(genome))))
  # Order wrt chromosome (so it goes chrom > noncore > plasmid)
  summary.results.per.species$name.ordered = relevel(as.factor(summary.results.per.species$name), ref="chrom_score")
  
  summary.results.per.species$phylogeny = summary.results.per.species$species
  
  # mcmc_model_complete = MCMCglmm(value ~ name.ordered, 
  #                                  random=~phylogeny+species+n, 
  #                                  data=summary.results.per.species,
  #                                  pl=TRUE, slice=TRUE,
  #                                  nitt = N_ITERATIONS, burnin = BURNIN,
  #                                  thin = THIN, 
  #                                  prior=prior3, 
  #                                  ginverse = list(phylogeny = inv.phylogenetic.matrix),
  #                                  verbose = TRUE)
  mcmc_model_complete = MCMCglmm(value ~ name.ordered, 
                                 random=~species+n, 
                                 data=summary.results.per.species,
                                 pl=TRUE, slice=TRUE,
                                 nitt = N_ITERATIONS, burnin = BURNIN,
                                 thin = THIN, 
                                 prior=prior2, 
                                 ginverse = list(species = inv.phylogenetic.matrix),
                                 verbose = TRUE)
  
  # CALCULATIONS FOR MODEL
  model = mcmc_model_complete
  
  mFixed_accessory = mean(model$Sol[,2])*model$X[,2] # first entry
  mFixed_plasmid = mean(model$Sol[,3])*model$X[,3]
  mVarF_accessory = var(mFixed_accessory) #variance of fixed effects
  mVarF_plasmid = var(mFixed_plasmid) #variance of fixed effects
  mean_random_effects = apply(model$VCV,2,mean)
  mRandomF <- sum(mean_random_effects)
  
  # mVarF_accessory/(mVarF_accessory+mVarF_plasmid+mRandomF) * 100
  # mVarF_plasmid/(mVarF_accessory+mVarF_plasmid+mRandomF) * 100
  # mean_random_effects/(mVarF_accessory+mVarF_plasmid+mRandomF) * 100
  
  # Looking at mean differences
  df.mean.summary = df %>% group_by(species) %>%
    summarise(mean=mean(score_diff),
              se=sd(score_diff)/sqrt(length(score_diff)),
              n=length(score_diff))
  df.mean.summary$species = ordered(df.mean.summary$species, 
                                    levels=df.mean.summary$species[order(df.mean.summary$mean)])
  
  
  # Test for difference in variance
  variance.test.results = matrix(nrow=length(unique(summary.df$species)), ncol=4)
  rownames(variance.test.results) = unique(summary.df$species)
  for (s in unique(summary.df$species)){
    d = summary.df[which(summary.df$species==s),]
    if (nrow(d)>N.MIN){
      v = var.test(d$plasmid_score, d$chrom_score)
      variance.test.results[s,] = c(v$estimate, v$conf.int[1], v$conf.int[2], v$p.value)
      
    }
  }
  
  variance.test.results = data.frame(variance.test.results)
  variance.test.results$species = rownames(variance.test.results)
  
  variance.test.results$p.corrected = p.adjust(variance.test.results$X4)
  table(variance.test.results$p.corrected<0.05)
  
  # Correlation of mean avoidance in Plasmid vs core
  df.sum = summary.df %>% group_by(species) %>%
    summarise(p=mean(plasmid_score),
              c=mean(chrom_score),
              a=mean(accessory_score),
              n=length(species))
  df.sum = data.frame(df.sum)
  cor.test(df.sum$p,
           df.sum$c,
           method="spearman")
  cor.test(df.sum$p,
           df.sum$a,
           method="spearman")
  df.sum$species = gsub("_", " ", df.sum$species)
  df.sum.plot = df.sum[which(df.sum$n>N.MIN),]
  p.correlation = ggplot(df.sum, aes(c, p))+
    geom_hline(yintercept = 0, linetype='dashed')+
    geom_vline(xintercept = 0, linetype='dashed')+
    geom_point(size=2)+
    xlab("Core avoidance")+
    ylab("Plasmid avoidance")+
    theme_bw()+
    theme(panel.grid = element_blank())+
    theme(axis.text=element_text(colour="black"))+
    xlim(c(min(c(df.sum$c, df.sum$p)),max(0.1, max(c(df.sum$c, df.sum$p)))))+
    ylim(c(min(c(df.sum$c, df.sum$p)),max(0.1, max(c(df.sum$c, df.sum$p)))))
  
  df.sum.sum = df.sum %>% pivot_longer(cols=c("p", "c", "a"))
  df.sum.sum$name = ordered(df.sum.sum$name, 
                            levels=c("c", "a", "p"),
                            labels=c("Core", "Non-core", "Plasmid"))
  #df.sum.sum = df.sum.sum[which(df.sum.sum$n>N.MIN),] # remove those with fewer observations
  p.boxplot = ggplot(df.sum.sum, aes(name, value, colour=name))+
    geom_hline(yintercept = 0, linetype='dashed')+
    geom_boxplot(outlier.shape = NA)+
    geom_quasirandom(alpha=0.7, width=0.25, height=0)+
    ggsignif::geom_signif(colour="black", test="wilcox.test", test.args=list(paired=TRUE),
                          comparisons=list(c("Core", "Plasmid"),
                                           c("Non-core", "Plasmid")),
                          step=0.1, 
                          textsize=2)+
    theme_bw()+
    theme(panel.grid=element_blank())+
    theme(axis.text=element_text(colour="black"))+
    xlab("")+
    ylab(expression(" "%<-%"Stronger avoidance"))+
    scale_colour_manual(values=component.colours)+
    theme(legend.position = "none")+
    theme(axis.title.y=element_text(size=10))
  
  
  
  df.plasmid$species.name = gsub("_", " ", df.plasmid$species)
  df.plasmid.abundant.mean.se = df.plasmid %>% 
    group_by(species.name) %>%
    summarise(mean=mean(diff.score),
              se=sd(diff.score)/sqrt(length(diff.score)))
  df.plasmid.abundant.mean.se$species.name = ordered(df.plasmid.abundant.mean.se$species.name,
                                                     levels=df.plasmid.abundant.mean.se$species.name[order(df.plasmid.abundant.mean.se$mean, decreasing = TRUE)])
  
  # Add modelling
  modelling.output.plasmid = summary(model)$solutions["name.orderedplasmid_score",]
  modelling.estimate = data.frame(species.name="Modelled effect",
                                  mean=modelling.output.plasmid[1],
                                  lower=modelling.output.plasmid[2],
                                  upper=modelling.output.plasmid[3])
  p.bars = ggplot(df.plasmid.abundant.mean.se, 
                  aes(y=species.name, x=mean, xmin=mean-se, xmax=mean+se))+
    theme_bw()+
    theme(panel.grid=element_blank())+
    geom_vline(xintercept = 0, linetype='dashed')+
    geom_bar(stat="identity", aes(fill=mean<0))+
    geom_errorbar(width=0)+
    geom_point()+
    theme(axis.text=element_text(colour="black"),
          axis.text.y=element_text(face="italic", size=8))+
    xlab(expression(""%<-%" "%->%""))+
    ylab("")+
    scale_y_discrete(position = "right")+
    theme(legend.position = "none")+
    theme(panel.border = element_blank(),
          axis.line = element_line(colour = "black"))+
    scale_fill_manual(values=as.character(c(component.colours["Plasmid"],
                                            component.colours["Core"])))+
    scale_x_continuous(position="top")
  
  p.model = ggplot(data=modelling.estimate, 
                   aes(x=-mean, y=species.name, xmin=-lower, xmax=-upper))+ 
    geom_bar(stat="identity", aes(fill=mean<0))+
    geom_errorbar(width=0)+
    geom_point()+
    xlim(layer_scales(p.bars)$x$range$range)+
    theme_bw()+
    theme(panel.border = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x =element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_text(colour="black", face="bold"))+
    xlab("")+
    ylab("")+
    scale_y_discrete(position="right")+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    scale_fill_manual(values=as.character(c(component.colours["Plasmid"],
                                            component.colours["Core"])))+
    theme(legend.position = "none")
  
  p.bars.2 = cowplot::plot_grid(p.bars+ggtitle("(c)"), 
                                p.model, 
                                nrow=2,
                                rel_heights = c(1, 0.1),
                                align="v", axis="b")
  
  p.left = cowplot::plot_grid(p.boxplot+ggtitle("(a)"),
                              p.correlation+ggtitle("(b)"), 
                              nrow=2, align='v', axis='b',
                              rel_widths=c(1,1))
  
  
  
  p.combined = cowplot::plot_grid(p.left, 
                                  p.bars.2,
                                  nrow=1,
                                  rel_heights = c(1.1, 1),
                                  rel_widths = c(1, 1.5))
  
  # Save the plot
  if (K==6){
    saveFigure(p.combined, file=paste0(file_prefix, "Figure1_k", K, "-subsampling-", subsampling_value),
               width=9, height=6.5)
  }
  if (K==4){
    saveFigure(p.combined, file=paste0(file_prefix, "FigureS1_k", K, "-subsampling-", subsampling_value),
               width=9, height=6.5)
  }
  
  
  # Write results
  N.genomes = length(unique(summary.df$genome))
  N.species = length(unique(summary.df$species))
  median.variance.plasmid = median(na.omit(variance.test.results$X1))
  
  results.table = data.frame(summary(model)$solutions)
  rownames(results.table) = c("Intercept", "Category (relative to core): Non-core", "Category (relative to core): Plasmid")
  colnames(results.table) = c("Mean estimate","Lower (95% CI)", "Upper (95%)", "Effective sample size", "pMCMC" )
  variance.explained =    
    variance.explained = c(mVarF_accessory/(mVarF_accessory+mVarF_plasmid+mRandomF) * 100,
                           mVarF_plasmid/(mVarF_accessory+mVarF_plasmid+mRandomF) * 100,
                           mean_random_effects/(mVarF_accessory+mVarF_plasmid+mRandomF) * 100)
  names(variance.explained) = c("Non-core", "Plasmid",
                                "Phylogeny", "No. of genomes", "Residuals")
  
  header.string = c(K, subsampling_value, N.genomes, N.species)
  names(header.string) = c("K", "Subsampling (bp)", "No. of genomes", "No. of species (n>=3 genomes)")
  
  output_file = paste0(outputDir, file_prefix, "_", Sys.Date(), "-K-", K, "-subsampling-", subsampling_value, ".tsv")
  write.table(header.string, file=output_file, append=FALSE, 
              quote=F, col.names =F, sep="\t")
  write.table("\nModel output", file=output_file, append = TRUE, 
              quote=F, row.names = F, col.names = F, sep="\t")
  write.table(results.table, file=output_file, append = TRUE, quote=F, row.names = T, sep="\t")
  write.table("\nVariance explained", file=output_file, append = TRUE, 
              quote=F, row.names = F, col.names = F, sep="\t")
  write.table(t(variance.explained), file=output_file, 
              append = TRUE, quote=F, row.names = F,sep="\t")
  
  return(df.sum) 
}
