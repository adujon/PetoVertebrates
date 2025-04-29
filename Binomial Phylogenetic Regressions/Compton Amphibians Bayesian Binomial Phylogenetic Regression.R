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


#setup import and export folders
export_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/models"
import_dir <- "C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import phylogeny data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#import the phylogenetic tree for squamates
setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/data")
phylo_all <- read.tree("amph_shl_new_Consensus_7238.tre") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import cancer data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data <-  read.csv("Compton et al amphibians.csv")
data <- data.frame(data) #convert to dataframe
colnames(data)[11] <- "NNecropsies"  
colnames(data)[12] <- "NNeoplasia"  
colnames(data)[16] <- "NMalignant"  

#change the dataframe rown names to species names
row.names(data) <- data$Species

#convert the body mass to kg and maximum longevity to months
data$adult_mass_kg <- data$Body_Mass_Grams / 1000
data$max_longevity_months <- data$Maximum_Longevity_Years * 12

data$logMass <- log10(data$adult_mass_kg) #log transform the body mass
data$logMaxLongevity <- log10(data$max_longevity_months) #log transform the  Max Longevity in months

#remove species with no mass or longevity data
sel_rm <- c()
for(i in 1:nrow(data)){
  
  if(is.na(data$adult_mass_kg[i])&is.na(data$logMaxLongevity[i])){
    
    sel_rm <- c( sel_rm, i)
  }
  
}

data <- data[-sel_rm,]

#remove species not in tree
data <- data[-which(data$Species == "Lithobates_pipiens"),]
data <- data[-which(data$Species == "Lithobates_catesbeianus"),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Visualize the data ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#histogram of the body mass
hist(data$adult_mass_kg)
hist(data$logMass) #better when log transformed

#histogram of the number of necropsies
hist(data$NNecropsies)

#histovrame log longevity
hist(data$logMaxLongevity)

#compute the number of species with no individuals that died from malignant cancer
length(which(data$NNeoplasia == 0)) / nrow(data) 

#histogram of the number of individuals that died from malignant cancer
hist(data$NNeoplasia)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Prepare the data for brms models ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#check which species are both in the main dataframe and the phylogeny tree
obj_brms<- name.check(phylo_all,data)
obj_brms


## fix the name for the Mississippi gopher frog
sel_ran <- which(phylo_all$tip.label == "Rana_sevosa")
phylo_all$tip.label[sel_ran] <- "Lithobates_sevosus"


#recheck the conflicts have been revolved
obj_brms<- name.check(phylo_all,data)
obj_brms

#drop the species from the tree for which we have no cancer data
phylo_brms <- drop.tip(phylo_all, obj_brms$tree_not_data)


#compute the var-covar matrix from the tree to use in brms
vcv_brms_OU <- vcv(corMartins(1, phylo_brms, fixed=FALSE,form=~phylo_brms$tip.label)) #OU model
vcv_brms_Pagel <- vcv(corPagel(value = 1, phy = phylo_brms, fixed = FALSE, form=~phylo_brms$tip.label)) #Pagel model

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#######################  Create a table to store the effect sizes ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

effect_sizes_summary <- data.frame(Dataset = NULL,
                                   Tumor_Type = NULL,
                                   Variable = NULL,
                                   Effect_Size_OU = NULL,
                                   Effect_Size_Pagel = NULL,
                                   LOOCV_OU = NULL,
                                   LOOCV_Pagel = NULL,
                                   Diff_LOOCV_OU_Minus_Pagel = NULL,
                                   WAIC_OU = NULL,
                                   WAIC_Pagel = NULL,
                                   Diff_WAIC_OU_Minus_Pagel= NULL)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#######################  Modelling the effect of logMass on neoplasia ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#~~~~~~~~~~
## OU ----
#~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_neoplasia_OU <- brm(NNeoplasia | trials(NNecropsies) ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
                                family =binomial(),
                                data = data,
                                data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                chains = 4,
                                cores = 24,
                                iter = 10000, warmup = 5000,
                                thin = 10,
                                seed = 1,
                                control=list(max_treedepth = 12, adapt_delta = 0.99),
                                save_pars = save_pars(all = TRUE))

fit_logmass_neoplasia_OU

#inspect the chains for improper mixing 
plot(fit_logmass_neoplasia_OU)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_neoplasia_OU , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_neoplasia_OU , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_neoplasia_OU , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_neoplasia_OU , ndraws = 100)


#~~~~~~~~~~~
## Pagel ----
#~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_neoplasia_Pagel <- brm(NNeoplasia | trials(NNecropsies) ~ logMass + (1|gr(Species, cov = vcv_brms_Pagel)),
                                   family =binomial(),
                                   data = data,
                                   data2 = list(vcv_brms_Pagel = vcv_brms_Pagel), 
                                   chains = 4,
                                   cores = 24,
                                   iter = 10000, warmup = 5000,
                                   thin = 10,
                                   seed = 1,
                                   control=list(max_treedepth = 12, adapt_delta = 0.99),
                                   save_pars = save_pars(all = TRUE))

fit_logmass_neoplasia_Pagel

#inspect the chains for improper mixing 
plot(fit_logmass_neoplasia_Pagel)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_neoplasia_Pagel , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_neoplasia_Pagel , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_neoplasia_Pagel , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_neoplasia_Pagel , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_neoplasia_Pagel , "logMass")
plot(cond_effects)

#~~~~~~~~~~~~~~~~~
## Find the best model ----
#~~~~~~~~~~~~~~~~~

#compared models
model_OU <- fit_logmass_neoplasia_OU
model_Pagel <- fit_logmass_neoplasia_Pagel

#calculate the LOOCV values
loo_OU <- round(loo(model_OU)$elpd_loo, digits = 1) 
loo_Pagel <- round(loo(model_Pagel)$elpd_loo , digits = 1) 

#calculate the WAIC
WAIC_OU <- round(waic(model_OU)$waic, digits =  1)
WAIC_Pagel <- round(waic(model_Pagel)$waic, digits =  1)

#get the effect sizes for OU model
es_OU <- paste(round(fixef(model_OU)[2,1], digits = 2), " (95%CI:", 
               round(fixef(model_OU)[2,3], digits = 2), ",", 
               round(fixef(model_OU)[2,4], digits = 2), ")", sep ="")

#get the effect sizes Pagel model
es_Pagel <- paste(round(fixef(model_Pagel)[2,1], digits = 2), " (95%CI:", 
                  round(fixef(model_Pagel)[2,3], digits = 2), ",", 
                  round(fixef(model_Pagel)[2,4], digits = 2), ")", sep ="")

#combine the results into a vector
effect_sizes_tmp <- data.frame(Dataset = "Compton et al.",
                               Tumor_Type = "Neoplasia",
                               Variable = "logMass",
                               Effect_Size_OU = es_OU,
                               Effect_Size_Pagel = es_Pagel,
                               LOOCV_OU = loo_OU,
                               LOOCV_Pagel = loo_Pagel,
                               Diff_LOOCV_OU_Minus_Pagel = round((loo_OU - loo_Pagel), digits = 2),
                               WAIC_OU = WAIC_OU,
                               WAIC_Pagel = WAIC_Pagel,
                               Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU - WAIC_Pagel), digits = 2))

#add the effect size to main table
effect_sizes_summary <- rbind(effect_sizes_summary, effect_sizes_tmp)

#check the residuals of the best model
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_logmass_neoplasia_OU)),
  observedResponse = data$NNeoplasia,
  fittedPredictedResponse = apply(t(posterior_epred(fit_logmass_neoplasia_OU)), 1, mean),
  integerResponse = TRUE)

#plot the residuals
plot(model.check)

#export the best model
setwd(export_dir)
saveRDS(fit_logmass_neoplasia_OU , "Compton_amphibians_neoplasm_logMass.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Modelling the effect of logMass on malignancy ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#~~~~~~~~~~
## OU ----
#~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_malignant_OU <- brm(NMalignant | trials(NNecropsies) ~ logMass + (1|gr(Species, cov = vcv_brms_OU)),
                                family =binomial(),
                                data = data,
                                data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                chains = 4,
                                cores = 32,
                                iter = 10000, warmup = 5000,
                                thin = 10,
                                seed = 1,
                                control=list(max_treedepth = 12, adapt_delta = 0.99),
                                save_pars = save_pars(all = TRUE))

#inspect the chains for improper mixing
plot(fit_logmass_malignant_OU)

#check if the iteration in the chains are not auto correlated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_malignant_OU, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_malignant_OU, pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_malignant_OU, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_malignant_OU, ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_malignant_OU, "logMass")
plot(cond_effects)

#~~~~~~~~~~~
## Pagel ----
#~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_malignant_Pagel <- brm(NMalignant | trials(NNecropsies) ~ logMass + (1|gr(Species, cov = vcv_brms_Pagel)),
                                   family =binomial(),
                                   data = data,
                                   data2 = list(vcv_brms_Pagel = vcv_brms_Pagel), 
                                   chains = 4,
                                   cores = 24,
                                   iter = 10000, warmup = 5000,
                                   thin = 10,
                                   seed = 1,
                                   control=list(max_treedepth = 12, adapt_delta = 0.99),
                                   save_pars = save_pars(all = TRUE))

fit_logmass_malignant_Pagel

#inspect the chains for improper mixing 
plot(fit_logmass_malignant_Pagel)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_malignant_Pagel , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_malignant_Pagel , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_malignant_Pagel , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_malignant_Pagel , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_malignant_Pagel , "logMass")
plot(cond_effects)

#~~~~~~~~~~~~~~~~~
## Find the best model ----
#~~~~~~~~~~~~~~~~~

#compared models
model_OU <- fit_logmass_malignant_OU
model_Pagel <- fit_logmass_malignant_Pagel

#calculate the LOOCV values
loo_OU <- round(loo(model_OU)$elpd_loo, digits = 1) 
loo_Pagel <- round(loo(model_Pagel)$elpd_loo , digits = 1) 

#calculate the WAIC
WAIC_OU <- round(waic(model_OU)$waic, digits =  1)
WAIC_Pagel <- round(waic(model_Pagel)$waic, digits =  1)

#get the effect sizes for OU model
es_OU <- paste(round(fixef(model_OU)[2,1], digits = 2), " (95%CI:", 
               round(fixef(model_OU)[2,3], digits = 2), ",", 
               round(fixef(model_OU)[2,4], digits = 2), ")", sep ="")

#get the effect sizes Pagel model
es_Pagel <- paste(round(fixef(model_Pagel)[2,1], digits = 2), " (95%CI:", 
                  round(fixef(model_Pagel)[2,3], digits = 2), ",", 
                  round(fixef(model_Pagel)[2,4], digits = 2), ")", sep ="")

#combine the results into a vector
effect_sizes_tmp <- data.frame(Dataset = "Compton et al.",
                               Tumor_Type = "Malignant",
                               Variable = "logMass",
                               Effect_Size_OU = es_OU,
                               Effect_Size_Pagel = es_Pagel,
                               LOOCV_OU = loo_OU,
                               LOOCV_Pagel = loo_Pagel,
                               Diff_LOOCV_OU_Minus_Pagel = round((loo_OU - loo_Pagel), digits = 2),
                               WAIC_OU = WAIC_OU,
                               WAIC_Pagel = WAIC_Pagel,
                               Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU - WAIC_Pagel), digits = 2))

#add the effect size to main table
effect_sizes_summary <- rbind(effect_sizes_summary, effect_sizes_tmp)

#check the residuals of the best model
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_logmass_malignant_OU)),
  observedResponse = data$NMalignant,
  fittedPredictedResponse = apply(t(posterior_epred(fit_logmass_malignant_OU)), 1, mean),
  integerResponse = TRUE)

#plot the residuals
plot(model.check)

#export the model
setwd(export_dir)
saveRDS(fit_logmass_malignant_OU, "Compton_amphibians_malignant_logMass.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#  Modelling the effect of logMaxLongevity on neoplasia ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#~~~~~~~~~~~
## OU ----
#~~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_neoplasia_OU <- brm(NNeoplasia | trials(NNecropsies) ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                        chains = 4,
                                        cores = 32,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

#inspect the chains for improper mixing
plot(fit_logMaxLongevity_neoplasia_OU)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logMaxLongevity_neoplasia_OU , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logMaxLongevity_neoplasia_OU , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logMaxLongevity_neoplasia_OU , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logMaxLongevity_neoplasia_OU , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logMaxLongevity_neoplasia_OU , "logMaxLongevity")
plot(cond_effects)

#~~~~~~~~~~~
## Pagel----
#~~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_neoplasia_Pagel <- brm(NNeoplasia | trials(NNecropsies) ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_Pagel)),
                                           family =binomial(),
                                           data = data,
                                           data2 = list(vcv_brms_Pagel = vcv_brms_Pagel), 
                                           chains = 4,
                                           cores = 32,
                                           iter = 10000, warmup = 5000,
                                           thin = 10,
                                           seed = 1,
                                           control=list(max_treedepth = 12, adapt_delta = 0.99),
                                           save_pars = save_pars(all = TRUE))

#inspect the chains for improper mixing
plot(fit_logMaxLongevity_neoplasia_Pagel)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logMaxLongevity_neoplasia_Pagel , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logMaxLongevity_neoplasia_Pagel , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logMaxLongevity_neoplasia_Pagel , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logMaxLongevity_neoplasia_Pagel , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logMaxLongevity_neoplasia_Pagel , "logMaxLongevity")
plot(cond_effects)

#~~~~~~~~~~~~~~~~~
## Find the best model ----
#~~~~~~~~~~~~~~~~~

#compared models
model_OU <- fit_logMaxLongevity_neoplasia_OU
model_Pagel <- fit_logMaxLongevity_neoplasia_Pagel

#calculate the LOOCV values
loo_OU <- round(loo(model_OU)$elpd_loo, digits = 1) 
loo_Pagel <- round(loo(model_Pagel)$elpd_loo , digits = 1) 

#calculate the WAIC
WAIC_OU <- round(waic(model_OU)$waic, digits =  1)
WAIC_Pagel <- round(waic(model_Pagel)$waic, digits =  1)

#get the effect sizes for OU model
es_OU <- paste(round(fixef(model_OU)[2,1], digits = 2), " (95%CI:", 
               round(fixef(model_OU)[2,3], digits = 2), ",", 
               round(fixef(model_OU)[2,4], digits = 2), ")", sep ="")

#get the effect sizes Pagel model
es_Pagel <- paste(round(fixef(model_Pagel)[2,1], digits = 2), " (95%CI:", 
                  round(fixef(model_Pagel)[2,3], digits = 2), ",", 
                  round(fixef(model_Pagel)[2,4], digits = 2), ")", sep ="")

#combine the results into a vector
effect_sizes_tmp <- data.frame(Dataset = "Compton et al.",
                               Tumor_Type = "Neoplasia",
                               Variable = "logMaxLongevity",
                               Effect_Size_OU = es_OU,
                               Effect_Size_Pagel = es_Pagel,
                               LOOCV_OU = loo_OU,
                               LOOCV_Pagel = loo_Pagel,
                               Diff_LOOCV_OU_Minus_Pagel = round((loo_OU - loo_Pagel), digits = 2),
                               WAIC_OU = WAIC_OU,
                               WAIC_Pagel = WAIC_Pagel,
                               Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU - WAIC_Pagel), digits = 2))

#add the effect size to main table
effect_sizes_summary <- rbind(effect_sizes_summary, effect_sizes_tmp)

#check the residuals of the best model
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_logMaxLongevity_neoplasia_OU)),
  observedResponse = data$NNeoplasia,
  fittedPredictedResponse = apply(t(posterior_epred(fit_logMaxLongevity_neoplasia_OU)), 1, mean),
  integerResponse = TRUE)

#plot the residuals
plot(model.check)

#export the model
setwd(export_dir)
saveRDS(fit_logMaxLongevity_neoplasia_OU , "Compton_amphibians_neoplasm_logMaxLongevity.rds")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#  Modelling the effect of logMaxLongevity on malignancy ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#~~~~~~~~~~~
## OU ----
#~~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_malignant_OU <- brm(NMalignant | trials(NNecropsies) ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                        family =binomial(),
                                        data = data,
                                        data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                        chains = 4,
                                        cores = 32,
                                        iter = 10000, warmup = 5000,
                                        thin = 10,
                                        seed = 1,
                                        control=list(max_treedepth = 12, adapt_delta = 0.99),
                                        save_pars = save_pars(all = TRUE))

#inspect the chains for improper mixing
plot(fit_logMaxLongevity_malignant_OU)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logMaxLongevity_malignant_OU, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logMaxLongevity_malignant_OU, pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logMaxLongevity_malignant_OU, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logMaxLongevity_malignant_OU, ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logMaxLongevity_malignant_OU, "logMaxLongevity")
plot(cond_effects)

#~~~~~~~~~~~
## Pagel ----
#~~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logMaxLongevity_malignant_Pagel <- brm(NMalignant | trials(NNecropsies) ~ logMaxLongevity + (1|gr(Species, cov = vcv_brms_Pagel)),
                                           family =binomial(),
                                           data = data,
                                           data2 = list(vcv_brms_Pagel = vcv_brms_Pagel), 
                                           chains = 4,
                                           cores = 32,
                                           iter = 10000, warmup = 5000,
                                           thin = 10,
                                           seed = 1,
                                           control=list(max_treedepth = 12, adapt_delta = 0.99),
                                           save_pars = save_pars(all = TRUE))

#inspect the chains for improper mixing
plot(fit_logMaxLongevity_malignant_Pagel)

#check if the iteration in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logMaxLongevity_malignant_Pagel, pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logMaxLongevity_malignant_Pagel, pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logMaxLongevity_malignant_Pagel, pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logMaxLongevity_malignant_Pagel, ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logMaxLongevity_malignant_Pagel, "logMaxLongevity")
plot(cond_effects)

#~~~~~~~~~~~~~~~~~
## Find the best model ----
#~~~~~~~~~~~~~~~~~

#compared models
model_OU <- fit_logMaxLongevity_malignant_OU
model_Pagel <- fit_logMaxLongevity_malignant_Pagel

#calculate the LOOCV values
loo_OU <- round(loo(model_OU)$elpd_loo, digits = 1) 
loo_Pagel <- round(loo(model_Pagel)$elpd_loo , digits = 1) 

#calculate the WAIC
WAIC_OU <- round(waic(model_OU)$waic, digits =  1)
WAIC_Pagel <- round(waic(model_Pagel)$waic, digits =  1)

#get the effect sizes for OU model
es_OU <- paste(round(fixef(model_OU)[2,1], digits = 2), " (95%CI:", 
               round(fixef(model_OU)[2,3], digits = 2), ",", 
               round(fixef(model_OU)[2,4], digits = 2), ")", sep ="")

#get the effect sizes Pagel model
es_Pagel <- paste(round(fixef(model_Pagel)[2,1], digits = 2), " (95%CI:", 
                  round(fixef(model_Pagel)[2,3], digits = 2), ",", 
                  round(fixef(model_Pagel)[2,4], digits = 2), ")", sep ="")

#combine the results into a vector
effect_sizes_tmp <- data.frame(Dataset = "Compton et al.",
                               Tumor_Type = "Malignant",
                               Variable = "logMaxLongevity",
                               Effect_Size_OU = es_OU,
                               Effect_Size_Pagel = es_Pagel,
                               LOOCV_OU = loo_OU,
                               LOOCV_Pagel = loo_Pagel,
                               Diff_LOOCV_OU_Minus_Pagel = round((loo_OU - loo_Pagel), digits = 2),
                               WAIC_OU = WAIC_OU,
                               WAIC_Pagel = WAIC_Pagel,
                               Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU - WAIC_Pagel), digits = 2))

#add the effect size to main table
effect_sizes_summary <- rbind(effect_sizes_summary, effect_sizes_tmp)

#check the residuals of the best model
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_logMaxLongevity_malignant_OU)),
  observedResponse = data$NMalignant,
  fittedPredictedResponse = apply(t(posterior_epred(fit_logMaxLongevity_malignant_OU)), 1, mean),
  integerResponse = TRUE)

#plot the residuals
plot(model.check)

#export the model
setwd(export_dir)
saveRDS(fit_logMaxLongevity_malignant_OU, "Compton_amphibians_malignant_logMaxLongevity.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#  Modelling the effect of logMass and logMaxLongevity on neoplasia ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#~~~~~~~~~
# No interaction ----
#~~~~~~~~~

#~~~~~~~~~~
## OU ----
#~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_logmaxlongevity_neoplasia_OU <- brm(NNeoplasia | trials(NNecropsies) ~ logMass + logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                                family =binomial(),
                                                data = data,
                                                data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                                chains = 4,
                                                cores = 24,
                                                iter = 10000, warmup = 5000,
                                                thin = 10,
                                                seed = 1,
                                                control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                save_pars = save_pars(all = TRUE))

fit_logmass_logmaxlongevity_neoplasia_OU

#inspect the chains for improper mixing 
plot(fit_logmass_logmaxlongevity_neoplasia_OU)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_OU , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_OU , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_OU , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_OU , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_logmaxlongevity_neoplasia_OU , ndraws = 100)


#~~~~~~~~~~~
## Pagel ----
#~~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_logmaxlongevity_neoplasia_Pagel <- brm(NNeoplasia | trials(NNecropsies) ~ logMass + logMaxLongevity + (1|gr(Species, cov = vcv_brms_Pagel)),
                                                   family =binomial(),
                                                   data = data,
                                                   data2 = list(vcv_brms_Pagel = vcv_brms_Pagel), 
                                                   chains = 4,
                                                   cores = 24,
                                                   iter = 10000, warmup = 5000,
                                                   thin = 10,
                                                   seed = 1,
                                                   control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                   save_pars = save_pars(all = TRUE))

fit_logmass_logmaxlongevity_neoplasia_Pagel

#inspect the chains for improper mixing 
plot(fit_logmass_logmaxlongevity_neoplasia_Pagel)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_Pagel , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_Pagel , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_Pagel , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_Pagel , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_logmaxlongevity_neoplasia_Pagel , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_logmaxlongevity_neoplasia_Pagel , "logMass")
plot(cond_effects)

#~~~~~~~~~
# With interaction ----
#~~~~~~~~~

#~~~~~~~~~~
## OU ----
#~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_logmaxlongevity_neoplasia_OU_int <- brm(NNeoplasia | trials(NNecropsies) ~ logMass * logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                                    family =binomial(),
                                                    data = data,
                                                    data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                                    chains = 4,
                                                    cores = 24,
                                                    iter = 10000, warmup = 5000,
                                                    thin = 10,
                                                    seed = 1,
                                                    control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                    save_pars = save_pars(all = TRUE))

fit_logmass_logmaxlongevity_neoplasia_OU_int

#inspect the chains for improper mixing 
plot(fit_logmass_logmaxlongevity_neoplasia_OU_int)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_OU_int , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_OU_int , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_OU_int , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_OU_int , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_logmaxlongevity_neoplasia_OU_int , ndraws = 100)


#~~~~~~~~~~~
## Pagel ----
#~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_logmaxlongevity_neoplasia_Pagel_int <- brm(NNeoplasia | trials(NNecropsies) ~ logMass * logMaxLongevity + (1|gr(Species, cov = vcv_brms_Pagel)),
                                                       family =binomial(),
                                                       data = data,
                                                       data2 = list(vcv_brms_Pagel = vcv_brms_Pagel), 
                                                       chains = 4,
                                                       cores = 24,
                                                       iter = 10000, warmup = 5000,
                                                       thin = 10,
                                                       seed = 1,
                                                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                       save_pars = save_pars(all = TRUE))

fit_logmass_logmaxlongevity_neoplasia_Pagel_int

#inspect the chains for improper mixing 
plot(fit_logmass_logmaxlongevity_neoplasia_Pagel_int)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_Pagel_int , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_Pagel_int , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_Pagel_int , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_neoplasia_Pagel_int , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_logmaxlongevity_neoplasia_Pagel_int , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_logmaxlongevity_neoplasia_Pagel_int , "logMass")
plot(cond_effects)

#~~~~~~~~~~~~~~~~~
## Find the best model ----
#~~~~~~~~~~~~~~~~~
#compared models
model_OU <- fit_logmass_logmaxlongevity_neoplasia_OU
model_Pagel <- fit_logmass_logmaxlongevity_neoplasia_Pagel
model_OU_int <- fit_logmass_logmaxlongevity_neoplasia_OU_int
model_Pagel_int <- fit_logmass_logmaxlongevity_neoplasia_Pagel_int

#calculate the LOOCV values
loo_OU <- round(loo(model_OU)$elpd_loo, digits = 1) 
loo_Pagel <- round(loo(model_Pagel)$elpd_loo , digits = 1) 
loo_OU_int <- round(loo(model_OU_int)$elpd_loo, digits = 1) 
loo_Pagel_int <- round(loo(model_Pagel_int)$elpd_loo , digits = 1)

#calculate the WAIC
WAIC_OU <- round(waic(model_OU)$waic, digits =  1)
WAIC_Pagel <- round(waic(model_Pagel)$waic, digits =  1)
WAIC_OU_int <- round(waic(model_OU_int)$waic, digits =  1)
WAIC_Pagel_int <- round(waic(model_Pagel_int)$waic, digits =  1)

#get the logMass effect sizes for OU model (no interaction)
es_OU_logMass <- paste(round(fixef(model_OU)[2,1], digits = 2), " (95%CI:", 
                       round(fixef(model_OU)[2,3], digits = 2), ",", 
                       round(fixef(model_OU)[2,4], digits = 2), ")", sep ="")

#get the logMass effect sizes Pagel model (no interaction)
es_Pagel_logMass <- paste(round(fixef(model_Pagel)[2,1], digits = 2), " (95%CI:", 
                          round(fixef(model_Pagel)[2,3], digits = 2), ",", 
                          round(fixef(model_Pagel)[2,4], digits = 2), ")", sep ="")


#get the logMass effect sizes for OU model (with interaction)
es_OU_logMass_int <- paste(round(fixef(model_OU_int)[2,1], digits = 2), " (95%CI:", 
                           round(fixef(model_OU_int)[2,3], digits = 2), ",", 
                           round(fixef(model_OU_int)[2,4], digits = 2), ")", sep ="")

#get the logMass effect sizes Pagel model (with interaction)
es_Pagel_logMass_int <- paste(round(fixef(model_Pagel_int)[2,1], digits = 2), " (95%CI:", 
                              round(fixef(model_Pagel_int)[2,3], digits = 2), ",", 
                              round(fixef(model_Pagel_int)[2,4], digits = 2), ")", sep ="")

#get the logMaxLongevity effect sizes for OU model (no interaction)
es_OU_logMaxLongevity <- paste(round(fixef(model_OU)[3,1], digits = 2), " (95%CI:", 
                               round(fixef(model_OU)[3,3], digits = 2), ",", 
                               round(fixef(model_OU)[3,4], digits = 2), ")", sep ="")

#get the logMaxLongevity effect sizes Pagel model (no interaction)
es_Pagel_logMaxLongevity <- paste(round(fixef(model_Pagel)[3,1], digits = 2), " (95%CI:", 
                                  round(fixef(model_Pagel)[3,3], digits = 2), ",", 
                                  round(fixef(model_Pagel)[3,4], digits = 2), ")", sep ="")


#get the logMaxLongevity effect sizes for OU model (with interaction)
es_OU_logMaxLongevity_int <- paste(round(fixef(model_OU_int)[3,1], digits = 2), " (95%CI:", 
                                   round(fixef(model_OU_int)[3,3], digits = 2), ",", 
                                   round(fixef(model_OU_int)[3,4], digits = 2), ")", sep ="")

#get the logMaxLongevity effect sizes Pagel model (with interaction)
es_Pagel_logMaxLongevity_int <- paste(round(fixef(model_Pagel_int)[3,1], digits = 2), " (95%CI:", 
                                      round(fixef(model_Pagel_int)[3,3], digits = 2), ",", 
                                      round(fixef(model_Pagel_int)[3,4], digits = 2), ")", sep ="")


#combine the results for logMass into a vector (no interaction)
effect_sizes_tmp_logMass <- data.frame(Dataset = "Compton et al.",
                                       Tumor_Type = "Neoplasia",
                                       Variable = "logMass + logMaxLongevity (effect size for logMass)",
                                       Effect_Size_OU = es_OU_logMass,
                                       Effect_Size_Pagel = es_Pagel_logMass,
                                       LOOCV_OU = loo_OU,
                                       LOOCV_Pagel = loo_Pagel,
                                       Diff_LOOCV_OU_Minus_Pagel = round((loo_OU - loo_Pagel), digits = 2),
                                       WAIC_OU = WAIC_OU,
                                       WAIC_Pagel = WAIC_Pagel,
                                       Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU - WAIC_Pagel), digits = 2))

#combine the results for logMaxLongevity into a vector (no interaction)
effect_sizes_tmp_logMaxLongevity <- data.frame(Dataset = "Compton et al.",
                                               Tumor_Type = "Neoplasia",
                                               Variable = "logMass + logMaxLongevity (effect size for logMaxLongevity)",
                                               Effect_Size_OU = es_OU_logMaxLongevity,
                                               Effect_Size_Pagel = es_Pagel_logMaxLongevity,
                                               LOOCV_OU = loo_OU,
                                               LOOCV_Pagel = loo_Pagel,
                                               Diff_LOOCV_OU_Minus_Pagel = round((loo_OU - loo_Pagel), digits = 2),
                                               WAIC_OU = WAIC_OU,
                                               WAIC_Pagel = WAIC_Pagel,
                                               Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU - WAIC_Pagel), digits = 2))

#combine the results for logMass into a vector (with interaction)
effect_sizes_tmp_logMass_int <- data.frame(Dataset = "Compton et al.",
                                           Tumor_Type = "Neoplasia",
                                           Variable = "logMass * logMaxLongevity (effect size for logMass)",
                                           Effect_Size_OU = es_OU_logMass_int,
                                           Effect_Size_Pagel = es_Pagel_logMass_int,
                                           LOOCV_OU = loo_OU_int,
                                           LOOCV_Pagel = loo_Pagel_int,
                                           Diff_LOOCV_OU_Minus_Pagel = round((loo_OU_int - loo_Pagel_int), digits = 2),
                                           WAIC_OU = WAIC_OU_int,
                                           WAIC_Pagel = WAIC_Pagel_int,
                                           Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU_int - WAIC_Pagel_int), digits = 2))

#combine the results  for logMaxLongevity into a vector (with  interaction)
effect_sizes_tmp_logMaxLongevity_int <- data.frame(Dataset = "Compton et al.",
                                                   Tumor_Type = "Neoplasia",
                                                   Variable = "logMass * logMaxLongevity (effect sizes for logMaxLongevity)",
                                                   Effect_Size_OU = es_OU_logMaxLongevity_int,
                                                   Effect_Size_Pagel = es_Pagel_logMaxLongevity_int,
                                                   LOOCV_OU = loo_OU_int,
                                                   LOOCV_Pagel = loo_Pagel_int,
                                                   Diff_LOOCV_OU_Minus_Pagel = round((loo_OU_int - loo_Pagel_int), digits = 2),
                                                   WAIC_OU = WAIC_OU_int,
                                                   WAIC_Pagel = WAIC_Pagel_int,
                                                   Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU_int - WAIC_Pagel_int), digits = 2))

#add the effect size to main table
effect_sizes_summary <- rbind(effect_sizes_summary, 
                              effect_sizes_tmp_logMass,
                              effect_sizes_tmp_logMaxLongevity,
                              effect_sizes_tmp_logMass_int,
                              effect_sizes_tmp_logMaxLongevity_int)

#check the residuals of the best model
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_logmass_logmaxlongevity_neoplasia_OU_int)),
  observedResponse = data$NNeoplasia,
  fittedPredictedResponse = apply(t(posterior_epred(fit_logmass_logmaxlongevity_neoplasia_OU_int)), 1, mean),
  integerResponse = TRUE)

#plot the residuals
plot(model.check)

#export the best model
setwd(export_dir)
saveRDS(fit_logmass_logmaxlongevity_neoplasia_OU , "Compton_amphibians_neoplasm_logMass_logMaxLongevity.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
#  Modelling the effect of logMass and logMaxLongevity on malignancy ----
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

#~~~~~~~~~
# No interaction ----
#~~~~~~~~~

#~~~~~~~~~~
## OU ----
#~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_logmaxlongevity_malignant_OU <- brm(NMalignant | trials(NNecropsies) ~ logMass + logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                                family =binomial(),
                                                data = data,
                                                data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                                chains = 4,
                                                cores = 24,
                                                iter = 10000, warmup = 5000,
                                                thin = 10,
                                                seed = 1,
                                                control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                save_pars = save_pars(all = TRUE))

fit_logmass_logmaxlongevity_malignant_OU

#inspect the chains for improper mixing 
plot(fit_logmass_logmaxlongevity_malignant_OU)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_logmaxlongevity_malignant_OU , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_OU , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_OU , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_OU , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_logmaxlongevity_malignant_OU , ndraws = 100)


#~~~~~~~~~~~
## Pagel ----
#~~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_logmaxlongevity_malignant_Pagel <- brm(NMalignant | trials(NNecropsies) ~ logMass + logMaxLongevity + (1|gr(Species, cov = vcv_brms_Pagel)),
                                                   family =binomial(),
                                                   data = data,
                                                   data2 = list(vcv_brms_Pagel = vcv_brms_Pagel), 
                                                   chains = 4,
                                                   cores = 24,
                                                   iter = 10000, warmup = 5000,
                                                   thin = 10,
                                                   seed = 1,
                                                   control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                   save_pars = save_pars(all = TRUE))

fit_logmass_logmaxlongevity_malignant_Pagel

#inspect the chains for improper mixing 
plot(fit_logmass_logmaxlongevity_malignant_Pagel)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_logmaxlongevity_malignant_Pagel , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_Pagel , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_Pagel , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_Pagel , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_logmaxlongevity_malignant_Pagel , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_logmaxlongevity_malignant_Pagel , "logMass")
plot(cond_effects)

#~~~~~~~~~
# With interaction ----
#~~~~~~~~~

#~~~~~~~~~~
## OU ----
#~~~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_logmaxlongevity_malignant_OU_int <- brm(NMalignant | trials(NNecropsies) ~ logMass * logMaxLongevity + (1|gr(Species, cov = vcv_brms_OU)),
                                                    family =binomial(),
                                                    data = data,
                                                    data2 = list(vcv_brms_OU = vcv_brms_OU), 
                                                    chains = 4,
                                                    cores = 24,
                                                    iter = 10000, warmup = 5000,
                                                    thin = 10,
                                                    seed = 1,
                                                    control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                    save_pars = save_pars(all = TRUE))

fit_logmass_logmaxlongevity_malignant_OU_int

#inspect the chains for improper mixing 
plot(fit_logmass_logmaxlongevity_malignant_OU_int)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_logmaxlongevity_malignant_OU_int , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_OU_int , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_OU_int , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_OU_int , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_logmaxlongevity_malignant_OU_int , ndraws = 100)


#~~~~~~~~~~~
## Pagel ----
#~~~~~~~~

#fit the phylogenetic regression model using brms
fit_logmass_logmaxlongevity_malignant_Pagel_int <- brm(NMalignant | trials(NNecropsies) ~ logMass * logMaxLongevity + (1|gr(Species, cov = vcv_brms_Pagel)),
                                                       family =binomial(),
                                                       data = data,
                                                       data2 = list(vcv_brms_Pagel = vcv_brms_Pagel), 
                                                       chains = 4,
                                                       cores = 24,
                                                       iter = 10000, warmup = 5000,
                                                       thin = 10,
                                                       seed = 1,
                                                       control=list(max_treedepth = 12, adapt_delta = 0.99),
                                                       save_pars = save_pars(all = TRUE))

fit_logmass_logmaxlongevity_malignant_Pagel_int

#inspect the chains for improper mixing 
plot(fit_logmass_logmaxlongevity_malignant_Pagel_int)

#check if the iterations in the chains are not autocorrelated
#this will plot the autocorrelation function for each of the four chains
mcmc_acf(fit_logmass_logmaxlongevity_malignant_Pagel_int , pars = "b_Intercept", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_Pagel_int , pars = "b_logMass", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_Pagel_int , pars = "b_logMaxLongevity", lags = 10)
mcmc_acf(fit_logmass_logmaxlongevity_malignant_Pagel_int , pars = "sd_Species__Intercept", lags = 10)

#check if the posterior distribution matches the observed distribution
pp_check(fit_logmass_logmaxlongevity_malignant_Pagel_int , ndraws = 100)

#plot the conditional effect of mass for this model
cond_effects <- conditional_effects(fit_logmass_logmaxlongevity_malignant_Pagel_int , "logMass")
plot(cond_effects)

#~~~~~~~~~~~~~~~~~
## Find the best model ----
#~~~~~~~~~~~~~~~~~

#compared models
model_OU <- fit_logmass_logmaxlongevity_malignant_OU
model_Pagel <- fit_logmass_logmaxlongevity_malignant_Pagel
model_OU_int <- fit_logmass_logmaxlongevity_malignant_OU_int
model_Pagel_int <- fit_logmass_logmaxlongevity_malignant_Pagel_int

#calculate the LOOCV values
loo_OU <- round(loo(model_OU)$elpd_loo, digits = 1) 
loo_Pagel <- round(loo(model_Pagel)$elpd_loo , digits = 1) 
loo_OU_int <- round(loo(model_OU_int)$elpd_loo, digits = 1) 
loo_Pagel_int <- round(loo(model_Pagel_int)$elpd_loo , digits = 1)

#calculate the WAIC
WAIC_OU <- round(waic(model_OU)$waic, digits =  1)
WAIC_Pagel <- round(waic(model_Pagel)$waic, digits =  1)
WAIC_OU_int <- round(waic(model_OU_int)$waic, digits =  1)
WAIC_Pagel_int <- round(waic(model_Pagel_int)$waic, digits =  1)

#get the logMass effect sizes for OU model (no interaction)
es_OU_logMass <- paste(round(fixef(model_OU)[2,1], digits = 2), " (95%CI:", 
                       round(fixef(model_OU)[2,3], digits = 2), ",", 
                       round(fixef(model_OU)[2,4], digits = 2), ")", sep ="")

#get the logMass effect sizes Pagel model (no interaction)
es_Pagel_logMass <- paste(round(fixef(model_Pagel)[2,1], digits = 2), " (95%CI:", 
                          round(fixef(model_Pagel)[2,3], digits = 2), ",", 
                          round(fixef(model_Pagel)[2,4], digits = 2), ")", sep ="")


#get the logMass effect sizes for OU model (with interaction)
es_OU_logMass_int <- paste(round(fixef(model_OU_int)[2,1], digits = 2), " (95%CI:", 
                           round(fixef(model_OU_int)[2,3], digits = 2), ",", 
                           round(fixef(model_OU_int)[2,4], digits = 2), ")", sep ="")

#get the logMass effect sizes Pagel model (with interaction)
es_Pagel_logMass_int <- paste(round(fixef(model_Pagel_int)[2,1], digits = 2), " (95%CI:", 
                              round(fixef(model_Pagel_int)[2,3], digits = 2), ",", 
                              round(fixef(model_Pagel_int)[2,4], digits = 2), ")", sep ="")

#get the logMaxLongevity effect sizes for OU model (no interaction)
es_OU_logMaxLongevity <- paste(round(fixef(model_OU)[3,1], digits = 2), " (95%CI:", 
                               round(fixef(model_OU)[3,3], digits = 2), ",", 
                               round(fixef(model_OU)[3,4], digits = 2), ")", sep ="")

#get the logMaxLongevity effect sizes Pagel model (no interaction)
es_Pagel_logMaxLongevity <- paste(round(fixef(model_Pagel)[3,1], digits = 2), " (95%CI:", 
                                  round(fixef(model_Pagel)[3,3], digits = 2), ",", 
                                  round(fixef(model_Pagel)[3,4], digits = 2), ")", sep ="")


#get the logMaxLongevity effect sizes for OU model (with interaction)
es_OU_logMaxLongevity_int <- paste(round(fixef(model_OU_int)[3,1], digits = 2), " (95%CI:", 
                                   round(fixef(model_OU_int)[3,3], digits = 2), ",", 
                                   round(fixef(model_OU_int)[3,4], digits = 2), ")", sep ="")

#get the logMaxLongevity effect sizes Pagel model (with interaction)
es_Pagel_logMaxLongevity_int <- paste(round(fixef(model_Pagel_int)[3,1], digits = 2), " (95%CI:", 
                                      round(fixef(model_Pagel_int)[3,3], digits = 2), ",", 
                                      round(fixef(model_Pagel_int)[3,4], digits = 2), ")", sep ="")


#combine the results for logMass into a vector (no interaction)
effect_sizes_tmp_logMass <- data.frame(Dataset = "Compton et al.",
                                       Tumor_Type = "Malignancy",
                                       Variable = "logMass + logMaxLongevity (effect size for logMass)",
                                       Effect_Size_OU = es_OU_logMass,
                                       Effect_Size_Pagel = es_Pagel_logMass,
                                       LOOCV_OU = loo_OU,
                                       LOOCV_Pagel = loo_Pagel,
                                       Diff_LOOCV_OU_Minus_Pagel = round((loo_OU - loo_Pagel), digits = 2),
                                       WAIC_OU = WAIC_OU,
                                       WAIC_Pagel = WAIC_Pagel,
                                       Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU - WAIC_Pagel), digits = 2))

#combine the results for logMaxLongevity into a vector (no interaction)
effect_sizes_tmp_logMaxLongevity <- data.frame(Dataset = "Compton et al.",
                                               Tumor_Type = "Malignancy",
                                               Variable = "logMass + logMaxLongevity (effect size for logMaxLongevity)",
                                               Effect_Size_OU = es_OU_logMaxLongevity,
                                               Effect_Size_Pagel = es_Pagel_logMaxLongevity,
                                               LOOCV_OU = loo_OU,
                                               LOOCV_Pagel = loo_Pagel,
                                               Diff_LOOCV_OU_Minus_Pagel = round((loo_OU - loo_Pagel), digits = 2),
                                               WAIC_OU = WAIC_OU,
                                               WAIC_Pagel = WAIC_Pagel,
                                               Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU - WAIC_Pagel), digits = 2))

#combine the results for logMass into a vector (with interaction)
effect_sizes_tmp_logMass_int <- data.frame(Dataset = "Compton et al.",
                                           Tumor_Type = "Malignancy",
                                           Variable = "logMass * logMaxLongevity (effect size for logMass)",
                                           Effect_Size_OU = es_OU_logMass_int,
                                           Effect_Size_Pagel = es_Pagel_logMass_int,
                                           LOOCV_OU = loo_OU_int,
                                           LOOCV_Pagel = loo_Pagel_int,
                                           Diff_LOOCV_OU_Minus_Pagel = round((loo_OU_int - loo_Pagel_int), digits = 2),
                                           WAIC_OU = WAIC_OU_int,
                                           WAIC_Pagel = WAIC_Pagel_int,
                                           Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU_int - WAIC_Pagel_int), digits = 2))

#combine the results  for logMaxLongevity into a vector (with  interaction)
effect_sizes_tmp_logMaxLongevity_int <- data.frame(Dataset = "Compton et al.",
                                                   Tumor_Type = "Malignancy",
                                                   Variable = "logMass * logMaxLongevity (effect sizes for logMaxLongevity)",
                                                   Effect_Size_OU = es_OU_logMaxLongevity_int,
                                                   Effect_Size_Pagel = es_Pagel_logMaxLongevity_int,
                                                   LOOCV_OU = loo_OU_int,
                                                   LOOCV_Pagel = loo_Pagel_int,
                                                   Diff_LOOCV_OU_Minus_Pagel = round((loo_OU_int - loo_Pagel_int), digits = 2),
                                                   WAIC_OU = WAIC_OU_int,
                                                   WAIC_Pagel = WAIC_Pagel_int,
                                                   Diff_WAIC_OU_Minus_Pagel = round((WAIC_OU_int - WAIC_Pagel_int), digits = 2))

#add the effect size to main table
effect_sizes_summary <- rbind(effect_sizes_summary, 
                              effect_sizes_tmp_logMass,
                              effect_sizes_tmp_logMaxLongevity,
                              effect_sizes_tmp_logMass_int,
                              effect_sizes_tmp_logMaxLongevity_int)

#check the residuals of the best model
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(fit_logmass_logmaxlongevity_malignant_OU_int)),
  observedResponse = data$NMalignant,
  fittedPredictedResponse = apply(t(posterior_epred(fit_logmass_logmaxlongevity_malignant_OU_int)), 1, mean),
  integerResponse = TRUE)

#plot the residuals
plot(model.check)

#export the best model
setwd(export_dir)
saveRDS(fit_logmass_logmaxlongevity_malignant_OU  , "Compton_amphibians_malignant_logMass_logMaxLongevity.rds")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Export the effect sizes in a table
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

setwd("C:/Users/antoi/My Drive/Postdoc/Zoo data/Peto/manuscript/Royal Open/Model results tables")
write.csv(effect_sizes_summary, "Compton et al amphibians.csv",  row.names = F)
