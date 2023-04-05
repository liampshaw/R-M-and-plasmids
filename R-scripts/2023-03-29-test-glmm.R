# AFter running some of 01_avoidance_of_palindromes.Rmd
species.tree.subset = keep.tip(species.tree, tip=unique(summary.results.per.species$species))
tree.d = cophenetic(species.tree.subset)
chrom.scores = summary.results.per.species$value[which(summary.results.per.species$name=="chrom_score")]
names(chrom.scores) =  summary.results.per.species$species[which(summary.results.per.species$name=="chrom_score")]

# look for correlation
score.d = matrix(nrow=length(chrom.scores), ncol=length(chrom.scores))
for (i in 1:length(chrom.scores)){
  print(i)
  for (j in 1:length(chrom.scores)){
    score.d[i,j] = chrom.scores[i]-chrom.scores[j]
  }
}

rownames(score.d) = names(chrom.scores)
colnames(score.d) = names(chrom.scores)

tree.d = tree.d[rownames(score.d), colnames(score.d)]

plot(tree.d, score.d)


phylo.plot = midpoint(species.tree.subset)

p.species = ggtree(phylo.plot)+theme_tree2()
p.species +
  geom_fruit(geom=geom_tile,
             data=summary.results.per.species[which(summary.results.per.species$name=="chrom_score"),],
             mapping=aes(y=species,fill=value),
             offset=0.2, width=0.1)


# Fake phylogeny etc
# Packages
library(tidyverse)
library(MCMCglmm)
library(ape)
library(phylotools)
library(phytools)
library(phangorn)
library(ggtree)
library(ggtreeExtra)

# Really stupid data
set.seed(1)
N_tree = 30
true_tree_1 = midpoint(rtree(N_tree, tip.label = 1:N_tree))
true_tree_2 = midpoint(rtree(N_tree, tip.label = (N_tree+1):(N_tree*2)))
# Bind these two trees together
true_tree = bind.tree(true_tree_1, true_tree_2)
#true_tree = bind.tree(true_tree, rtree(1))
true_tree = midpoint(true_tree)
plot(true_tree)
# make true data that has deliberate phylogenetic structure from the two trees
true_means = c(rnorm(N_tree, mean=0, sd=0.1), rnorm(N_tree, mean=2, sd=0.1))
plot(true_means)
true_data = data.frame(species=1:(2*N_tree),
                       phenotype=true_means,
                       animal=true_tree$tip.label)
true_tree = force.ultrametric(true_tree)
true_tree = di2multi(true_tree)
true_tree.Node <- makeNodeLabel(true_tree, method = "number")
INtree.true <- inverseA(true_tree.Node, nodes="TIPS") ## Makes matrix of phylogeny
Ainv.true <- INtree.true$Ainv # true inverse to use in MCMC modelling

p.true = ggtree(true_tree)
p.true = p.true +
  geom_fruit(geom=geom_tile,
             data=true_data,
             mapping=aes(y=species,fill=phenotype),
             offset=0.2)


# Random tree - no link to the phenotype
random_tree = midpoint(rtree(2*N_tree, tip.label = 1:(2*N_tree)))
random_tree = force.ultrametric(random_tree)
random_tree = di2multi(random_tree)
random_tree.Node <- makeNodeLabel(random_tree, method = "number")
INtree.random<- inverseA(random_tree.Node, nodes="TIPS") ## Makes matrix of phylogeny
Ainv.random <- INtree.random$Ainv # true inverse to use in MCMC modelling

p.random = ggtree(random_tree) 
p.random = p.random +
  geom_fruit(geom=geom_tile,
             data=true_data,
             mapping=aes(y=species,fill=phenotype),
             offset=0.2)

# Look at correlation of cophenetic distance and score distance
tree.d = cophenetic(true_tree)
chrom.scores = true_data$phenotype
names(chrom.scores) =  true_data$species

# look for correlation
score.d = matrix(nrow=length(chrom.scores), ncol=length(chrom.scores))
for (i in 1:length(chrom.scores)){
  print(i)
  for (j in 1:length(chrom.scores)){
    score.d[i,j] = chrom.scores[i]-chrom.scores[j]
  }
}

rownames(score.d) = names(chrom.scores)
colnames(score.d) = names(chrom.scores)

tree.d = tree.d[rownames(score.d), colnames(score.d)]
plot(tree.d, score.d)
