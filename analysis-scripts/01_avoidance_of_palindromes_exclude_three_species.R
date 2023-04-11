# Exclude species
source('setup.R')
source('01_modelling_functions.R')


EXCLUDE = c("Vibrio_parahaemolyticus", 
            "Vibrio_cholerae", 
            "Burkholderia_pseudomallei")
main.df.4 = main.df.4[which(! main.df.4$species %in% EXCLUDE),]
main.df.6 = main.df.6[which(! main.df.6$species %in% EXCLUDE),]

mean.50000.k6 = runModelMakePlot(subsampling_value = "50000", K=6, file_prefix="FigureS12_exclude-reanalysis")