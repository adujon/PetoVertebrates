rm(list=ls()) #remove all objects into memmory
library(readxl)
library(ape)
library(nlme)
library(geiger)
library(ggeffects)
library(caper)
library(phylolm)
library(phytools)
library(MuMIn)
library(ggplot2)
library(rr2)
library(DHARMa)
library(dplyr)
library(tidyr)
library(ggpubr)
library(grid)
library(brms)
library(parallel)
library(tidybayes)
library(bayesplot)

################################################################################

#setup import and export folders
export_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/models"
import_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data"

#############################
### import phylogeny data ###
#############################

#import the phylogenetic tree for mammals
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data")
phylo_all <- read.nexus("MamPhy_fullPosterior_BDvr_DNAonly_4098sp_topoFree_NDexp_MCC_v2_target.tre") 

#rename the phylogenetic tree tips to only keep the species name
current_tip_labels <- phylo_all$tip.label
new_tip_labels <- c()

for(i in 1:length(current_tip_labels)){
  
  vec <- strsplit(current_tip_labels[i], "[_]")[[1]][1:3]
  
  if(vec[1] == ""){
    new_tip_labels[i] <- paste(vec[2], vec[3], sep = "_")
  } else {
    new_tip_labels[i] <- paste(vec[1], vec[2], sep = "_")
  }
  
  
}

df_label <- data.frame(current_tip_labels, new_tip_labels)
phylo_all$tip.label<-df_label[[2]][match(phylo_all$tip.label, df_label[[1]])]

##########################
### import cancer data ###
##########################

###############
#Compton's data
data_compton <-  read.csv("Compton et al Vertebrates.csv")
data_compton <- data.frame(data_compton) #convert to dataframe
colnames(data_compton)[11] <- "NNecropsies" 
colnames(data_compton)[12] <- "NNeoplasia" 
colnames(data_compton)[16] <- "NMalignant" 

#subset the datset to only keep mammals
data_compton <- subset(data_compton, data_compton$Class == "Mammalia")

#replace the -1 in the mass and maxl ongevity columns by NA
data_compton$adult_weight.g[data_compton$adult_weight.g == -1] <- NA
data_compton$max_longevity.months.[data_compton$max_longevity.months. == -1] <- NA

#convert the body mass to kg
data_compton$adult_mass_kg <- data_compton$adult_weight.g / 1000
data_compton$logMass <- log(data_compton$adult_mass_kg) #log transform the body mass
data_compton$logMaxLongevity <- log(data_compton$max_longevity.months.) #log transform the  Max Longevity in months

data_compton <- data_compton[, c("Species", "NNecropsies", "logMass", "logMaxLongevity", "NMalignant", "NNeoplasia")]
data_compton$study <- "Compton"

################
## Boddy's data
data_boddy <-  read.csv("Boddy cancer mammals.csv")
data_boddy <- data.frame(data_boddy) 
colnames(data_boddy)[2] <- "Species" 
colnames(data_boddy)[6] <- "NNecropsies" 
colnames(data_boddy)[7] <- "NNeoplasia" 
colnames(data_boddy)[8] <- "NMalignant" 

row.names(data_boddy) <- data_boddy$Species #change row names to species names
data_boddy$logMass <- log(data_boddy$adult_mass_kg) #log transform the body mass
data_boddy$logMaxLongevity <- log(data_boddy$max_lifespan_yr * 12) #log transform the life expectancy

#subset the data to only keep the columns we need
data_boddy <- data_boddy[, c("Species","NNecropsies", "logMass", "logMaxLongevity", "NMalignant", "NNeoplasia")]
data_boddy$study <- "Boddy"

################
#Vincze's data
data_vincze <-  read.csv("Vincze cancer mammals.csv")
data_vincze <- data.frame(data_vincze) 
data_vincze$logMass <- log(data_vincze$Mass) 
data_vincze$NMalignant <- data_vincze$NMalignant_CMR 
data_vincze$logMaxLongevity <- log(data_vincze$MaxLongevity) 
data_vincze <- data_vincze[, c("Species","NNecropsies", "logMass", "logMaxLongevity", "NMalignant")]
data_vincze$NNeoplasia <- NA
data_vincze$study <- "Vincze"

################
#combine the three datasets together
data <- rbind(data_compton, data_boddy, data_vincze)
data$SpeciesRepetitions  <- data$Species

#fix species names synonyms
sel_cal <- which(data$Species == "Cebuella_pygmaea") 
data$Species[sel_cal] <- "Callithrix_pygmaea"

sel_mac <- which(data$Species == "Macropus_rufus") 
data$Species[sel_mac ] <- "Callithrix_pygmaea"

## fix the name for the parma wallaby 
sel_wall <- which(data$Species == "Macropus_parma")
data$Species[sel_wall] <- "Notamacropus_parma"

## fix the name for the Common wallaroo 
sel_wall <- which(data$Species == "Macropus_robustus")
data$Species[sel_wall] <- "Osphranter_robustus"

## fix the name for the red-necked wallaby t
sel_wall <- which(data$Species == "Macropus_rufogriseus")
data$Species[sel_wall] <- "Notamacropus_rufogriseus"

##############################################################################################
############################# visualize the data ########################################
##############################################################################################

#histogram of the body mass
hist(data$logMass) #better when log transformed

#histogram of the number of necropsies
hist(data$NNecropsies)

#histovrame of log max longevity
hist(data$logMaxLongevity)

#compute the number of species with no individuals that died from malignant cancer
length(which(data$NMalignant == 0)) / nrow(data) 


##############################################################################################
#################### prepare the data for brms models ########################################
##############################################################################################

#create a dataframe that can be used to prune the tree of all species that we don't need
df_prune <- data.frame(Species = unique(data$Species))
rownames(df_prune) <- df_prune$Species

#fix the tree for two species 
phylo_all <- bind.tip(phylo_all, "Cervus_canadensis", where = which(phylo_all$tip.label=="Cervus_elaphus"),
                      edge.length=0.5, position = 0.5)
phylo_all <- bind.tip(phylo_all, "Gazella_marica", where = which(phylo_all$tip.label=="Gazella_subgutturosa"),
                      edge.length=0.5, position = 0.5)

## fix the name for the red-necked wallaby to match names in the csv
sel_wall <- which(phylo_all$tip.label == "Macropus_rufogriseus")
phylo_all$tip.label[sel_wall] <- "Notamacropus_rufogriseus"

## fix the name for the parma wallaby to match names in the csv
sel_wall <- which(phylo_all$tip.label == "Macropus_parma")
phylo_all$tip.label[sel_wall] <- "Notamacropus_parma"

## fix the name for the Common wallaroo to match names in the csv
sel_wall <- which(phylo_all$tip.label == "Macropus_robustus")
phylo_all$tip.label[sel_wall] <- "Osphranter_robustus"

## fix the name for the Red kangaroo to match names in the csv
sel_kan <- which(phylo_all$tip.label == "Macropus_rufus")
phylo_all$tip.label[sel_kan] <- "Osphranter_rufus"

## fix the name for the Giant eland to match names in the csv
sel_gia <- which(phylo_all$tip.label == "Tragelaphus_derbianus")
phylo_all$tip.label[sel_gia] <- "Taurotragus_derbianus"

## fix the name for the donkeys Equus africanus asinus to match names in the csv
sel_donkey <- which(phylo_all$tip.label == "Equus_africanus")
phylo_all$tip.label[sel_donkey] <- "Equus_asinus"

# fix the name for the  Pygmy marmoset to match names in the csv
sel_donkey <- which(phylo_all$tip.label == "Cebuella_pygmaea")
phylo_all$tip.label[sel_donkey] <- "Callithrix_pygmaea"

#check which species are both in the main dataframe and the phylogeny tree
obj_brms<- name.check(phylo_all,df_prune)
obj_brms

#drop the species from the tree for which we have no cancer data
phylo_brms <- drop.tip(phylo_all, obj_brms$tree_not_data)
phylo_brms

#compute the var-covar matrix from the tree to use in brms
vcv_brms <- vcv(phylo_brms, model = "Brownian")

##############################################################################################
#######################  Modelling the effect of logMass on neoplasia #######################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logmass_neoplasm <- brm(NNeoplasia | trials(NNecropsies) ~ logMass + study +  (1|gr(Species, cov = vcv_brms)) + (1|SpeciesRepetitions),
                             family =binomial(),
                             data = data,
                             data2 = list(vcv_brms = vcv_brms), 
                             chains = 4,
                             cores = 32,
                             iter = 10000, warmup = 5000,
                             thin = 10,
                             seed = 1,
                             control=list(max_treedepth = 12, adapt_delta = 0.99))

#inspect the chains for improper mixing
plot(fit_logmass_neoplasm)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_neoplasm, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_neoplasm, pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_neoplasm, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_neoplasm, ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_neoplasm, "logMass")
plot(cond_effects)

#export the model
setwd(export_dir)
saveRDS(fit_logmass_neoplasm, "Combined_data_mammals_neoplasm_logMass.rds")

##############################################################################################
#######################  Modelling the effect of logMass on malignancy #######################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logmass_malignant <- brm(NMalignant | trials(NNecropsies) ~ logMass + study +  (1|gr(Species, cov = vcv_brms)) + (1|SpeciesRepetitions),
                            family =binomial(),
                            data = data,
                            data2 = list(vcv_brms = vcv_brms), 
                            chains = 4,
                            cores = 32,
                            iter = 10000, warmup = 5000,
                            thin = 10,
                            seed = 1,
                            control=list(max_treedepth = 12, adapt_delta = 0.99))

#inspect the chains for improper mixing
plot(fit_logmass_malignant)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_malignant, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_malignant, pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_malignant, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_malignant, ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_malignant, "logMass")
plot(cond_effects)

#export the model
setwd(export_dir)
saveRDS(fit_logmass_malignant, "Combined_data_mammals_malignant_logMass.rds")


##############################################################################################
#######################  Modelling the effect of logMaxLongevity on neoplasia ################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_neoplasm <- brm(NNeoplasia | trials(NNecropsies) ~ logMaxLongevity + study +  (1|gr(Species, cov = vcv_brms)) + (1|SpeciesRepetitions),
                                     family =binomial(),
                                     data = data,
                                     data2 = list(vcv_brms = vcv_brms), 
                                     chains = 4,
                                     cores = 32,
                                     iter = 10000, warmup = 5000,
                                     thin = 10,
                                     seed = 1,
                                     control=list(max_treedepth = 12, adapt_delta = 0.99))

#inspect the chains for improper mixing
plot(fit_logMaxLongevity_neoplasm)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logMaxLongevity_neoplasm, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logMaxLongevity_neoplasm, pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logMaxLongevity_neoplasm, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logMaxLongevity_neoplasm, ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logMaxLongevity_neoplasm, "logMaxLongevity")
plot(cond_effects)

#export the model
setwd(export_dir)
saveRDS(fit_logMaxLongevity_neoplasm, "Combined_data_mammals_neoplasm_logMaxLongevity.rds")

##############################################################################################
#######################  Modelling the effect of logMaxLongevity on malignancy #######################
##############################################################################################

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_malignant <- brm(NMalignant | trials(NNecropsies) ~ logMaxLongevity + study +  (1|gr(Species, cov = vcv_brms)) + (1|SpeciesRepetitions),
                             family =binomial(),
                             data = data,
                             data2 = list(vcv_brms = vcv_brms), 
                             chains = 4,
                             cores = 32,
                             iter = 10000, warmup = 5000,
                             thin = 10,
                             seed = 1,
                             control=list(max_treedepth = 12, adapt_delta = 0.99))

#inspect the chains for improper mixing
plot(fit_logMaxLongevity_malignant)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logMaxLongevity_malignant, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logMaxLongevity_malignant, pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logMaxLongevity_malignant, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logMaxLongevity_malignant, ndraws = 100)


#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logMaxLongevity_malignant, "logMaxLongevity")
plot(cond_effects)

#export the model
setwd(export_dir)
saveRDS(fit_logMaxLongevity_malignant, "Combined_data_mammals_malignant_logMaxLongevity.rds")

