# PLASMID TAXONOMIC UNITS
  
source('setup.R')

# READ DATA
# Redondo Salvo database
plasmid_db = read.csv(paste0(dataDir, 'redondo-salvo-2020-plasmid-DB.csv'),
                      header=T, 
                      stringsAsFactors = F,
                      row.names = 1)


# ANALYSIS
K = 4
SUBSAMPLING = 2500

plasmid.df.4 = read.csv(paste0(dataDir, 'plasmids/redondo-salvo-summary-k-', K, '-n-', prettyNum(SUBSAMPLING),'-True.csv'),
                        stringsAsFactors = FALSE,
                        header=T)
plasmid.df.4 = plasmid.df.4[which(plasmid.df.4$kmer_category=="Palindromic"),]
plasmid.df.4$Size = plasmid_db[plasmid.df.4$plasmid, "Size"]
plasmid.df.4$PTU = plasmid_db[plasmid.df.4$plasmid, "PTU"]
plasmid.df.4$hostrange = plasmid_db[plasmid.df.4$plasmid, "PTU.hostrange"]

plasmid.df.4$species.rs = gsub(" ", "_", plasmid_db[plasmid.df.4$plasmid, "TaxSpecies"])

plasmid.df.4$size.category = sapply(plasmid.df.4$Size,
                                    function(x) ifelse(x<10000, "2.5-10kb", 
                                                       ifelse(x<50000, "10-50kb", 
                                                              ifelse(x<100000, "50-100kb", ">100kb" ))))
plasmid.df.4$size.category = ordered(plasmid.df.4$size.category, 
                                     levels=c("2.5-10kb", "10-50kb", "50-100kb", ">100kb"))
p.4 = ggplot(plasmid.df.4, aes(score))+
  geom_histogram()+
  facet_wrap(~size.category, nrow=4)+
  xlab("Avoidance score")+
  ylab("Number of plasmids")+
  theme_bw()+
  ggtitle("(a) k = 4")

K = 6
SUBSAMPLING = 10000
plasmid.df.6 = read.csv(paste0(dataDir, 'plasmids/redondo-salvo-summary-k-', 6, '-n-', prettyNum(SUBSAMPLING),'-True.csv'),
                        stringsAsFactors = FALSE,
                        header=T)
plasmid.df.6 = plasmid.df.6[which(plasmid.df.6$kmer_category=="Palindromic"),]
plasmid.df.6$Size = plasmid_db[plasmid.df.6$plasmid, "Size"]
plasmid.df.6$PTU = plasmid_db[plasmid.df.6$plasmid, "PTU"]
plasmid.df.6$species.rs = gsub(" ", "_", plasmid_db[plasmid.df.6$plasmid, "TaxSpecies"])
plasmid.df.6$hostrange = plasmid_db[plasmid.df.6$plasmid, "PTU.hostrange"]

plasmid.df.6$size.category = sapply(plasmid.df.6$Size,
                                    function(x) ifelse(x<20000, "10-20kb", 
                                                       ifelse(x<50000, "20-50kb", 
                                                              ifelse(x<100000, "50-100kb", ">100kb" ))))
plasmid.df.6$size.category = ordered(plasmid.df.6$size.category, 
                                     levels=c("10-20kb", "20-50kb", "50-100kb", ">100kb"))
p.6 = ggplot(plasmid.df.6, aes(score))+
  geom_histogram()+
  facet_wrap(~size.category, nrow=4)+
  xlab("Avoidance score")+
  theme_bw()+
  ggtitle("(b) k=6")

# PLOTTING FUNCTIONS
# Try to summarise this information into a plot
modelForPTUs <- function(SUBSAMPLING, K, KMER_CATEGORY, SPECIES="", whole.model=FALSE){
  
  plasmid.df = read.csv(paste0(dataDir, 'plasmids/redondo-salvo-summary-k-', K, '-n-', prettyNum(SUBSAMPLING),'-True.csv'),
                        stringsAsFactors = FALSE,
                        header=T)
  if (SPECIES ==""){
    SPECIES = unique(plasmid.df$species)
  }
  plasmid.df = plasmid.df[which(plasmid.df$kmer_category==KMER_CATEGORY & plasmid.df$species %in% SPECIES),]
  plasmid.df$Size = plasmid_db[plasmid.df$plasmid, "Size"]
  plasmid.df$PTU = plasmid_db[plasmid.df$plasmid, "PTU"]
  plasmid.df$species.rs = gsub(" ", "_", plasmid_db[plasmid.df$plasmid, "TaxSpecies"])
  plasmid.df$hostrange = plasmid_db[plasmid.df$plasmid, "PTU.hostrange"]
  
  plasmid.df.modelling = plasmid.df[which(plasmid.df$PTU!="-"),]
  plasmid.df.modelling$hostrange = ordered(plasmid.df.modelling$hostrange,
                                           levels=c("I", "II", "III", "IV", "V", "VI"))
  plasmid.df.modelling$hostrange.numeric = as.numeric(plasmid.df.modelling$hostrange)
  # Remove NAs
  plasmid.df.modelling = na.omit(plasmid.df.modelling)
  # Summarise by PTU
  ptu.df.modelling = plasmid.df.modelling %>% group_by(PTU) %>%
    summarise(hostrange=unique(hostrange.numeric),
              score=median(score),
              size=log10(median(Size)),
              n=length(hostrange.numeric))
  lm.model = lm(score ~ hostrange + size + n, data=ptu.df.modelling)
  if (whole.model==FALSE){
    return(list(adj.r.squared=summary(lm.model)$adj.r.squared,
                coef=summary(lm.model)$coefficients,
                n.plasmids=nrow(plasmid.df.modelling),
                n.species=length(unique(plasmid.df.modelling$species)),
                n.PTUs = length(unique(plasmid.df.modelling$PTU))))
  }
  
  else{
    return(lm.model)
  }
  
}

# Sum of squres
sumOfSquares = function(SUBSAMPLING, K, level){
  m = modelForPTUs(SUBSAMPLING, K, level, whole.model = TRUE)
  af = anova(m)
  afss =af$"Sum Sq"
  return(afss/sum(afss) * 100)
}


makePlotEffectsErrorBars <- function(SUBSAMPLING, K, YLIM=0, TITLE=""){
  results.df = data.frame(rbind(c(modelForPTUs(SUBSAMPLING, K, "Species")$coef[,"Estimate"], 
                                  modelForPTUs(SUBSAMPLING, K, "Species")$coef[,"Std. Error"]),
                                c(modelForPTUs(SUBSAMPLING, K, "Genus")$coef[,"Estimate"],
                                  modelForPTUs(SUBSAMPLING, K, "Genus")$coef[,"Std. Error"]),
                                c(modelForPTUs(SUBSAMPLING, K, "Family")$coef[,"Estimate"],
                                  modelForPTUs(SUBSAMPLING, K, "Family")$coef[,"Std. Error"]),
                                c(modelForPTUs(SUBSAMPLING, K, "Order")$coef[,"Estimate"],
                                  modelForPTUs(SUBSAMPLING, K, "Order")$coef[,"Std. Error"]),
                                c(modelForPTUs(SUBSAMPLING, K, "Class")$coef[,"Estimate"],
                                  modelForPTUs(SUBSAMPLING, K, "Class")$coef[,"Std. Error"]),
                                c(modelForPTUs(SUBSAMPLING, K, "Phylum")$coef[,"Estimate"],
                                  modelForPTUs(SUBSAMPLING, K, "Phylum")$coef[,"Std. Error"]),
                                c(modelForPTUs(SUBSAMPLING, K, "Kingdom")$coef[,"Estimate"],
                                  modelForPTUs(SUBSAMPLING, K, "Kingdom")$coef[,"Std. Error"])))
                                #c(modelForPTUs(SUBSAMPLING, K, "Palindromic")$coef[,"Estimate"],
                                #  modelForPTUs(SUBSAMPLING, K, "Palindromic")$coef[,"Std. Error"])))
  #results.df$level = ordered(c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"),
  #                           levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Palindromic"))
  results.df$level = ordered(c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom" ),
                             levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom"))
  p.hostrange = ggplot(results.df, aes(level, hostrange, ymin=hostrange-hostrange.1, ymax=hostrange+hostrange.1))+
    geom_hline(yintercept = 0, linetype='dashed', colour='black')+
    geom_point(size=3, colour="#e41a1c")+
    ylab("Coefficient")+
    theme_bw()+
    theme(panel.grid = element_blank())+
    xlab("R-M targets within same")+
    geom_errorbar(width=0, colour="#e41a1c")+
    ggtitle("(a) PTU host range")+
     theme(axis.text=element_text(colour="black"),
                     title =element_text(size=8),
           axis.title=element_text(size=8))
  
  p.size = ggplot(results.df, aes(level, size, ymin=size-size.1, ymax=size+size.1))+
    geom_hline(yintercept = 0, linetype='dashed', colour='black')+
    geom_point(size=3, colour="#377eb8")+
    ylab("Coefficient")+
    theme_bw()+
    theme(panel.grid = element_blank())+
    xlab("")+
    geom_errorbar(width=0, colour="#377eb8")+
    ggtitle("(b) PTU median length (log10)")+
    theme(axis.text=element_text(colour="black"),
          title =element_text(size=8))+
      theme(axis.text.y=element_blank())

  
  
  p.n = ggplot(results.df, aes(level, n, ymin=n-n.1, ymax=n+n.1))+
    geom_hline(yintercept = 0, linetype='dashed', colour='black')+
    geom_point(size=3, colour="#4daf4a")+
    ylab("Coefficient")+
    theme_bw()+
    theme(panel.grid = element_blank())+
    xlab("")+
    geom_errorbar(width=0, colour="#4daf4a")+
    ggtitle("(c) No. of plasmids in PTU")+
    theme(axis.text=element_text(colour="black"),
                   title =element_text(size=8))+
      theme(axis.text.y=element_blank())

  
  
  if (YLIM!=0){
    p.hostrange = p.hostrange+ylim(YLIM)
    p.size = p.size+ylim(YLIM)
    p.n = p.n+ylim(YLIM)
  }
  if (YLIM==0){
    p.hostrange = p.hostrange+ylim(c(-max(abs(layer_scales(p.hostrange)$y$get_limits())),
                                     max(abs(layer_scales(p.hostrange)$y$get_limits()))))#+ylim(c(0-max(abs(results.df$hostrange))-2*max(abs(results.df$hostrange.1)),
    #       0+max(abs(results.df$hostrange))-2*max(abs(results.df$hostrange.1))))
    p.size = p.size+ylim(c(-max(abs(layer_scales(p.size)$y$get_limits())),
                           max(abs(layer_scales(p.size)$y$get_limits()))))#+ylim(c(0-max(abs(results.df$size))-2*max(abs(results.df$size.1)),
    #            0+max(abs(results.df$size))-2*max(abs(results.df$size.1))))
    p.n = p.n+ylim(c(-max(abs(layer_scales(p.n)$y$get_limits())),
                     max(abs(layer_scales(p.n)$y$get_limits()))))#+ylim(c(0-max(abs(results.df$size))-2*max(abs(results.df$size.1)),
    #+ylim(c(0-max(abs(results.df$n))-2*max(abs(results.df$n.1)),
    #               0+max(abs(results.df$n))-2*max(abs(results.df$n.1))))
    
  }
  title_theme <- ggdraw() +
    draw_label(TITLE, 
               fontfamily = "sans",  x = 0.05, hjust = 0)
  ss.df = data.frame(rbind(sumOfSquares(SUBSAMPLING, K, "Species"),
                           sumOfSquares(SUBSAMPLING, K, "Genus"),
                           sumOfSquares(SUBSAMPLING, K, "Family"),
                           sumOfSquares(SUBSAMPLING, K, "Order"),
                           sumOfSquares(SUBSAMPLING, K, "Class"),
                           sumOfSquares(SUBSAMPLING, K, "Phylum"),
                           sumOfSquares(SUBSAMPLING, K, "Kingdom")))
                          #sumOfSquares("10000", 6, "Palindromic")))
  colnames(ss.df) = c("hostrange", "size", "n", "residuals")
  #ss.df$level = c("Species", "Genus", "Family", "Order", "Class", "Phylum", 
  #                "Kingdom", "Palindromic")
  ss.df$level = c("Species", "Genus", "Family", "Order", "Class", "Phylum", 
                  "Kingdom")
  ss.df.melt = melt(ss.df, id.vars = "level")
  #ss.df.melt$level = ordered(ss.df.melt$level,
  #                           levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", 
  #                                    "Kingdom", "Palindromic"))
  ss.df.melt$level = ordered(ss.df.melt$level,
                             levels=c("Species", "Genus", "Family", "Order", "Class", "Phylum", 
                                      "Kingdom"))
  ss.plot = ggplot(ss.df.melt[which(ss.df.melt$variable!="residuals"),], aes(level, fill=variable, value))+
    geom_bar(position="stack", stat="identity")+
    theme_bw()+
    theme(panel.grid = element_blank())+
    xlab("")+
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    theme(axis.line.x=element_line(colour="black"))+
    theme(panel.border = element_blank())+
    ggtitle("(d) Model summary")+
    theme(title=element_text(size=8))+
    theme(axis.text=element_text(colour="black"))+
    scale_fill_brewer(palette="Set1")+
    theme(legend.position = "none")+
    ylab("Variance explained (%)")
  p.combined = cowplot::plot_grid(p.hostrange+coord_flip(), p.size+coord_flip(), p.n+coord_flip(), 
                                  ss.plot+coord_flip(), nrow=1,
                                  rel_widths = c(1.2, 1,1, 0.8))
  return(plot_grid(title_theme, p.combined, ncol=1, rel_heights=c(0.1,1)))
}


# MAKING PLOTS
p.6 = makePlotEffectsErrorBars("10000", 6,
                               TITLE="")



saveFigure(p.6, "Figure5_PTU-modelling-size-effect-10k-k6", 
       width=10, height=4)
ggsave(plot=p.6, "Figure5_PTU-modelling-size-effect-10k-k6.png", 
       width=25, height=8, unit="cm", dpi=300)

p.4 = makePlotEffectsErrorBars("10000", 4,
                               TITLE="")
saveFigure(p.4, 'FigureS10_PTU-modelling-size-effect-10k-k4', 
       width=10, height=4)
ggsave(plot=p.4, file=paste0(figureDir, Sys.Date(), 'FigureS9_PTU-modelling-size-effect-10k-k4.png'), 
       width=25, height=8, unit="cm", dpi=300)
p.5 = makePlotEffectsErrorBars("10000", 5,
                               TITLE="")
saveFigure(p.5, "FigureS11_PTU-modelling-size-effect-10k-k5", 
       width=10, height=4)
ggsave(plot=p.5, file=paste0(figureDir, Sys.Date(), "FigureS10_PTU-modelling-size-effect-10k-k5.png"), 
       width=25, height=8, unit="cm", dpi=300)




# INVESTIGATE MTASES
MTase_hits = read.csv(paste0(dataDir, 'rmsFinder_plasmids_MT_results.csv'),
                      header=T, stringsAsFactors = F)
MTase_hits$plasmid = gsub(".[0-9]_.*", "\\1\\2", MTase_hits$qseqid)
MTase_hits$length.target = nchar(MTase_hits$target)
# Do some filtering out
# MTase: 55% (similarity)
# These are an e-value of 0.001 and a lateral coverage of >50%.
MTase_hits_filtered = MTase_hits %>% filter(coverage_threshold_met=="True", 
                                            pident>55)
# Get list of targets
summary.targets = MTase_hits_filtered %>% group_by(plasmid) %>%
  summarise(targets=paste(unique(target)))

# Also check if they have an RM system that recognises the same target
RM_systems = read.csv(paste0(dataDir, 'rmsFinder_plasmids_RMS.csv'),
                      header=T, stringsAsFactors = F)
  RM_systems$plasmid = gsub(".[0-9]_.*", "\\1\\2", RM_systems$prot_MT)
RM_systems$length.target = nchar(RM_systems$sequence)
# How many carry an R-M system?
print("Plasmids carrying a putative R-M system:")
length(table(RM_systems$plasmid))
# And how big are these plasmids?
median(plasmid_db[unique(RM_systems$plasmid),"Size"])


RM_systems$plasmid_and_target = paste(RM_systems$plasmid, RM_systems$sequence)
MTase_hits_filtered$plasmid_and_target = paste(MTase_hits_filtered$plasmid, MTase_hits_filtered$target)


print("How many plasmids carrying MTase also have a putative R-M system with the same predicted target?")
table(MTase_hits_filtered$plasmid_and_target %in% RM_systems$plasmid_and_target)

# We can exclude these further:
MTase_hits_filtered = MTase_hits_filtered[which(! MTase_hits_filtered$plasmid %in% RM_systems$plasmid),]

# Should we check that they recognise the same target?


# Just 4-6bp targets? not for now
#MTase_hits_filtered = MTase_hits_filtered[which(MTase_hits_filtered$length.target %in% c(4, 5, 6)),]
# Plasmids and number of hits
plasmids.with.hits = MTase_hits_filtered %>% group_by(plasmid) %>%
  summarise(N.mt=length(target),
            targets=paste(unique(target), collapse='-'))
# Normalise by length
plasmids.with.hits$size = plasmid_db[plasmids.with.hits$plasmid,"Size"]
plasmids.with.hits$N.mt.norm = plasmids.with.hits$N.mt/plasmids.with.hits$size
plasmids.with.hits = data.frame(plasmids.with.hits)
rownames(plasmids.with.hits) = plasmids.with.hits$plasmid


# Remove trailing numbers on rownames
rownames(plasmid_db) = gsub("\\..*", "", rownames(plasmid_db))

# Check overlap
table(MTase_hits$plasmid  %in% rownames(plasmid_db))

plasmid_db$MTase = sapply(rownames(plasmid_db), 
                          function(x) ifelse(x %in% MTase_hits_filtered$plasmid, 1, 0))
plasmid_db$N.mt.norm = sapply(rownames(plasmid_db), 
                          function(x) plasmids.with.hits[x, "N.mt.norm"])
plasmid_db$N.mt.norm[is.na(plasmid_db$N.mt.norm)] = 0
plasmid_db$targets = sapply(rownames(plasmid_db), 
                              function(x) plasmids.with.hits[x, "targets"])
#plasmid_db$MTase_target = sapply(rownames(plasmid_db), 
#                                 function(x) MTase_hits[x])


# Average by PTU
PTU_summary = plasmid_db %>% group_by(PTU, PTU.hostrange) %>%
  summarise(MT=sum(MTase)/length(MTase),
            MT.norm=mean(N.mt.norm),
            size=median(Size),
            n=length(MTase),
            targets=length(unique(na.omit(targets))),
            n.MTase=sum(MTase)) %>%
  mutate(targets.norm=targets/n.MTase)
PTU_summary$hostrange = ordered(PTU_summary$PTU.hostrange,
                                levels=c("I", "II", "III", "IV", "V", "VI"))
PTU_summary = PTU_summary[!is.na(PTU_summary$hostrange),]
PTU_summary = PTU_summary[which(PTU_summary$PTU!="-"),]
PTU_summary$hostrange.numeric = as.numeric(PTU_summary$hostrange)

plasmid_db_modelling = plasmid_db
plasmid_db_modelling$hostrange = ordered(plasmid_db_modelling$PTU.hostrange,
                               levels=c("I", "II", "III", "IV", "V", "VI"))

plasmid_db_modelling$hostrange.numeric = as.numeric(plasmid_db_modelling$hostrange)
plasmid_db_modelling = plasmid_db_modelling[!is.na(plasmid_db_modelling$hostrange.numeric),]

glm_model_all = glm(MTase ~ log10(Size)+hostrange.numeric, data=plasmid_db_modelling, family = "binomial")
glm_model_PTU = glm(MT ~ log10(size)+hostrange.numeric, data=PTU_summary)
glm_model_PTU_diversity = glm(targets ~ log10(size)+hostrange.numeric, data=PTU_summary)
af <- summary(aov(glm_model_PTU))[[1]]
afss = af$"Sum Sq"
print("Summary of aov of glm_model_PTU")
modelling.results.PTU = cbind(af,PctExp=afss/sum(afss)*100)
print(modelling.results.PTU)
write.csv(file=paste0(outputDir, "TableSX_results_for_GLM_PTU_all.csv"), modelling.results.PTU)

# Larger than 100kb
glm_model_all_100kb = glm(MTase ~ log10(Size)+hostrange.numeric, data=plasmid_db_modelling[which(plasmid_db_modelling$Size>100000),], family = "binomial")
glm_model_PTU_100kb = glm(MT ~ log10(size)+hostrange.numeric, data=PTU_summary[which(PTU_summary$size>100000),], family = "gaussian")
af <- summary(aov(glm_model_PTU_100kb))[[1]]
afss = af$"Sum Sq"
print("Summary of aov of glm_model_PTU_100kb")
modelling.results.PTU.100kb = cbind(af,PctExp=afss/sum(afss)*100)
print(modelling.results.PTU.100kb)
write.csv(file=paste0(outputDir, "TableSX_results_for_GLM_PTU_100kb.csv"), modelling.results.PTU.100kb)
#lm_model_MT_norm_all = lm(MT.norm ~ log10(size)+hostrange.numeric, data=PTU_summary)

# exponential to get the predictions 0-1
predictions = exp(predict(glm_model_all, data.frame(plasmid_db_modelling)))


size.values = seq(1e3, 2e6, 1e3)
#predictions.df = data.frame(size=size.values,
#                           MT=predictions)

#cor.test(log10(PTU_summary$size), PTU_summary$hostrange.numeric, method="spearman")

# #summary(glm(MT ~ log10(size)+hostrange.numeric, data=PTU_summary[which(PTU_summary$size>50000),]))
# summary(glm(MT ~ log10(size)+hostrange.numeric, data=PTU_summary[which(PTU_summary$size>100000),]))
# 
# 
# ggplot(PTU_summary, aes(x=log10(size), MT.norm))+
#   geom_point()
# 
# 
# 
# 
# summary(lm(MT.norm ~ log10(size)+hostrange.numeric, data=PTU_summary[which(PTU_summary$size>40000),]))
# summary(lm(MT.norm ~ log10(size)+hostrange.numeric, data=PTU_summary[which(PTU_summary$size>50000),]))
# summary(lm(MT.norm ~ log10(size)+hostrange.numeric, data=PTU_summary[which(PTU_summary$size>100000),]))



# ggplot(PTU_summary, aes(group=targets, y=targets, x=log10(size)))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_jitter(height=0.1, alpha=0.7)+
#   theme_bw()+
#   theme(panel.grid = element_blank())+
#   theme(axis.text=element_text(colour="black"))+
#   xlab("PTU size (log10)")+
#   ylab("Unique targets of MTases")+
#   theme(panel.border = element_blank(),
#         axis.line = element_line(colour = "black"))+
#   scale_y_continuous(breaks=seq(0,5))+
#   facet_wrap(~hostrange)

PTU_summary = PTU_summary %>%
  mutate(new_bin = cut(log10(size), 
                       breaks=c(2, 3, 4, 4.69897, 5, 10),
                       labels=c("<1kb", "1-10kb", "10-50kb", "50-100kb", ">100kb")))
PTU_summary$targets.norm 
p.unique.targets = ggplot(PTU_summary, aes(group=hostrange, y=targets, x=hostrange))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.1, height=0,alpha=0.7)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("Host range")+
  ylab("Unique MTase targets")+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_wrap(~new_bin)+
  stat_smooth(method="lm", aes(x=hostrange.numeric, group=1), 
              se=FALSE)
p.total.targets = ggplot(PTU_summary, aes(group=hostrange, y=, x=hostrange))+
  geom_boxplot(outlier.shape = NA, width=0.5)+
  geom_jitter(width=0.1, height=0,alpha=0.7)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("Host range")+
  ylab("Unique MTase targets")+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_wrap(~new_bin)+
  stat_smooth(method="lm", aes(x=hostrange.numeric, group=1), 
              se=FALSE)


p.norm = ggplot(PTU_summary[which(PTU_summary$size>100000),], aes(group=hostrange, x=hostrange, y=MT.norm))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("PTU host range")+
  ylab("Type II MTase density (mean per base)")+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"))

p.prop = ggplot(PTU_summary[which(PTU_summary$size>100000),], aes(group=hostrange, x=hostrange, y=MT))+
  geom_hline(yintercept = 1, linetype='dashed')+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.1, alpha=0.7)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("PTU host range")+
  ylab("Proportion with Type II MTase(s)")+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
  
p.combined = cowplot::plot_grid(p.prop+ggtitle("(b) PTUs >100kbp"), 
                                p.norm, 
                                align='h', axis='t')
cowplot::plot_grid(p.unique.targets, 
                   p.combined, nrow=2)


p.prop = ggplot(PTU_summary, aes(group=hostrange, x=hostrange, y=MT, colour=new_bin))+
  geom_hline(yintercept = 1, linetype='dashed')+
  geom_quasirandom(width=0.1)+
  stat_summary(fun=median, shape="-", colour="black", size=2)+
  scale_color_manual(values=size.colours)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("PTU host range")+
  ylab("Proportion carrying\nMTase(s)")+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_wrap(~new_bin, nrow=1)+
  theme(axis.title.y=element_text(size=12))+
  theme(legend.position = "none")
p.norm = ggplot(PTU_summary, aes(group=hostrange, x=hostrange, y=MT.norm+0.5e-6, colour=new_bin))+
  geom_quasirandom(width=0.25)+
  stat_summary(fun=median, shape="-", colour="black", size=3)+
  scale_color_manual(values=size.colours)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text=element_text(colour="black"))+
  xlab("PTU host range")+
  ylab("MTase density\n(per base)")+
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black"))+
  facet_wrap(~new_bin, nrow=1)+
  scale_y_log10()+
  theme(axis.title.y=element_text(size=12))+
    theme(legend.position = "none")
p.combined = cowplot::plot_grid(p.prop+ggtitle("(a)")+theme(legend.position = "none"), 
                                p.norm+ggtitle("(b)")+theme(legend.position = "none"), 
                                align='h', axis='t', nrow=2)
saveFigure(p.combined, "FigureS12_PTUs-MTase-combined",
       width=8, height=6)

#ggsave(plot=p.prop, filename=paste0(figureDir, Sys.Date(), "-fig-PTUs-host-range-density.pdf"),
#       width=6, height=4)

#ggsave(plot=p.norm, filename=paste0(figureDir, Sys.Date(), "-fig-PTUs-host-range-proportion.pdf"),
#       width=6, height=4)

summary.df = PTU_summary %>% group_by(hostrange, new_bin) %>%
  summarise(mean=mean(MT),
            n=length(MT.norm))

plasmid_db = plasmid_db %>% mutate(new_bin = cut(log10(Size), 
                       breaks=c(2, 3, 4, 4.69897, 5, 10),
                       labels=c("<1kb", "1-10kb", "10-50kb", "50-100kb", ">100kb")))

summary.df = plasmid_db %>% group_by(PTU.hostrange, new_bin) %>%
  summarise(mean=mean(MTase),
            n=length(MTase),
            nMT=sum(MTase))
summary.df$PTU.hostrange = ordered(summary.df$PTU.hostrange,
                                   levels=c("-", "I", "II", "III", "IV", "V", "VI"),
                                   labels=c("Unassigned\nplasmids", "I", "II", "III", "IV", "V", "VI"))
p.heatmap = ggplot(summary.df, aes(PTU.hostrange, new_bin, fill=mean))+
  geom_tile(colour="black", size=0.5)+
  scale_fill_continuous(low="white", high="red")+
  geom_text(aes(label=paste0(nMT, "/", n )), size=3.5)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(fill="Proportion with MTase")+
  theme(axis.text=element_text(colour="black"))+
  xlab("Host range \n(of parent PTU)")+
  ylab("Plasmid length")+
  theme(panel.border  = element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12))
 saveFigure(p.heatmap, "Figure6_MTase-heatmap.pdf",
       width=7.5, height=4)

