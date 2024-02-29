#Fig5 Phylogenetic and developmental analysis of PNS microglia across the vertebrate tree of life.
library(parallel)
library(tidyverse)
library(ape)
library(geiger)
library(MCMCglmm)
library(coda)
library(ggtrendline)
library(treeplyr)
library(caper)
library(nlme)
library(picante)
library(GGally)
library(car)
setwd('E:/R/DRG/Fig5')


## Fig5h -----------------------------------------PIC-----------------------------------------------------------------------
traits <- read.csv('species_traits.csv', row.names = 1)
rownames(traits) <- traits$Species
glimpse(traits)



# Phylogenetic tree
oTree <- read.tree("vertebrate_tree_rooted.timed.tre")
str(oTree)
plot(oTree)
is.binary(oTree)
is.rooted(oTree)
is.ultrametric(oTree)

bCheck <- name.check(phy = oTree, data = traits, 
                     data.names = traits$Species)
bCheck

oTree <- drop.tip(oTree, bCheck$tree_not_data)
name.check(traits, oTree)

matches <- match(traits$Species, bCheck$data_not_tree, nomatch = 0)
traits <- subset(traits, matches == 0)
traits <- as.data.frame(traits)
# Remove zero length branches and replace with polytomies
oTree <- di2multi(oTree)
# Remove node labels 
oTree$node.label <- NULL




# calculate the PICs
pic_neu <- pic(log(traits$Neuron_size), oTree)
pic_p2ry12 <- pic(traits$P2RY12_counts, oTree)
pic_pro <- pic(traits$Proportion, oTree)



# pic_module_proportion
cor_pro <- cor.test(pic_neu, pic_pro, 
                    meth='spearman') 

# pic_module_counts
cor_counts <- cor.test(pic_neu, pic_p2ry12, 
                       meth='spearman')


## FigS5c -----------------------------------MCMCglmm, PGLS------------------------------------------------------------

# MCMCglmm
  datTraits <- read.csv('species_single_neuron_traits.csv')
  oTree <- read.tree("vertebrate_tree_rooted.timed.tre")
  str(oTree)
  is.binary(oTree)
  is.rooted(oTree)
  is.ultrametric(oTree)
  bCheck <- name.check(phy = oTree, data = datTraits, 
                       data.names = datTraits$Species)
  # Look at check
  bCheck
  oTree <- drop.tip(oTree, bCheck$tree_not_data)
  
  matches <- match(datTraits$Species, bCheck$data_not_tree, nomatch = 0)
  datTraits <- subset(datTraits, matches == 0)
  
  datTraits <- as.data.frame(datTraits)
  
  # Remove zero length branches and replace with polytomies
  oTree <- di2multi(oTree)
  # Remove node labels 
  oTree$node.label <- NULL
  
  # Get the inverse vcv matrix for the phylogeny
  #oTreeUltrametric <-  chronos(oTree, lambda=0)  
  #plot(oTreeUltrametric)
  inv.phylo <- inverseA(oTree, nodes = "TIPS", scale = T)$Ainv
  
  # Set up priors for MCMCglmm
  # Inverse Wishart with V = 1 and nu = 0.02
  # i.e. fairly uninformative priors
  prior <- list(G = list(G1 = list(V = 1, nu = 0.02)),
                R = list(V = 1, nu = 0.02))
  
  nitt <- 10000000
  thin <- 1000
  burnin <- 10000



for(nRep in 1:4){
model_mcmcglmm <- MCMCglmm(P2RY12_counts ~ log(Neuron_size), 
                             family = "poisson",
                             data = datTraits, 
                             random = ~ Species,
                             ginverse = list(Species = inv.phylo), 
                             prior = prior,
                             nitt = nitt, thin = thin, burnin = burnin,
                             verbose = TRUE)
    saveRDS(model_mcmcglmm, file = paste0("counts", nRep, ".rda"))
    }
summary(model_mcmcglmm)$solutions





# PGLS
np <- comparative.data(phy = oTree, data = as.data.frame(traits), 
                       names.col = Species, vcv = TRUE, 
                       na.omit = FALSE, warn.dropped = TRUE)
np$dropped$tips
np$dropped$unmatched.rows
# Fit a PGLS model
model.pgls <- pgls(Proportion ~ log(Neuron_size), 
                   data = np) 
anova(model.pgls)
summary(model.pgls)
coef(model.pgls)
lambda.profile <- pgls.profile(model.pgls, "lambda")
# Plot the likelihood profile
plot(lambda.profile)
