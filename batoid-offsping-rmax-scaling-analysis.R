# Load required packages
library(ape) # to read phylogenetic tree
library(caper) # to plot pgls models
library(car) # for vif function
library(dplyr)

# load phylogenetic and model data------------------------------------------------------
phy <- read.tree("data/stein-et-al-single.tree") # 85 species phylogeny
dat <- read.csv("data/batoid_model_offspring.csv", stringsAsFactors = FALSE,header = TRUE)

# fit final 11 pgls models---------------
# outputs reported in main text

# combine phylogeny with dataset for use in pgls function
cd <- comparative.data(phy, dat, names.col = "BinomName") 

# fit pgls models with temp, depth, offspring size

# intercept only
zslope_0 <- pgls(log(rmax) ~ 1, data = cd, lambda = "ML") 

# rmax varies with adult mass only 
zslope_wt <- pgls(log(rmax) ~ log_wt, data = cd, lambda = "ML") 

# depth only
zslope_depth <- pgls(log(rmax) ~ depth_scaled, data = cd, lambda = "ML")

# temp only
zslope_invtemp <- pgls(log(rmax) ~ invtemp_scaled, data = cd, lambda = "ML") 

# adult mass + depth
zslope_wt_depth <- pgls(log(rmax) ~ log_wt + depth_scaled,
                        data = cd, lambda = "ML")

# adult mass + temp
zslope_wt_invtemp <- pgls(log(rmax) ~ log_wt + invtemp_scaled,
                          data = cd, lambda = "ML")

# adult mass * depth 
zslope_wt_x_depth <- pgls(log(rmax) ~ log_wt * depth_scaled, 
                          data = cd, lambda = "ML")

# adult mass * temp 
zslope_wt_x_invtemp <- pgls(log(rmax) ~ log_wt * invtemp_scaled, 
                            data = cd, lambda = "ML")

# add in 3 models with offspring weight (& temp / depth)----------------
# excluding models with adult and offspring size in the same model
# rmax varies with offspring size
zslope_off <- pgls(log(rmax) ~ log_offspring,
                   data = cd, lambda = "ML")

# rmax varies with offspring size and temperature
zslope_invtemp_off <- pgls(log(rmax) ~ log_offspring + invtemp_scaled,
                           data = cd, lambda = "ML")

# rmax varies with offspring size and depth
zslope_depth_off <- pgls(log(rmax) ~ log_offspring + depth_scaled,
                         data = cd, lambda = "ML")

#  saving all 11 models in list
capermodels <- list(zslope_0,
                    zslope_wt, 
                    zslope_depth, 
                    zslope_invtemp,
                    zslope_wt_depth,
                    zslope_wt_invtemp, 
                    zslope_wt_x_depth,
                    zslope_wt_x_invtemp, 
                    zslope_off,
                    zslope_invtemp_off,
                    zslope_depth_off
) 

# manually extracting information from each model object
aics <- sapply(capermodels, function (x) bbmle::AIC(x)) 
aiccs <- sapply(capermodels, function (x) x$aicc) 
formulas <- sapply(capermodels,   # tidying formulas for easier reading
                   function (x) deparse(formula(x), width.cutoff = 90L)) %>%
  sub("^log\\(\\w+\\)\\s\\~\\s", "", .) %>% # update to log10 if want it to just have log in the formula
  gsub("log", "ln", .)  %>% # to have ln instead of log in formula
  gsub("_scaled", "", .) %>%
  gsub("\\_offspring", "(Moffspring)", .)  %>%
  gsub("\\_wt", "(M)", .)  %>%
  gsub("sizeratio", "M/Moffspring", .)  %>%
  gsub("O", "o", .) # so order is a lower case
r2s <- sapply(capermodels, function (x) summary(x)$r.squared) 
ar2s <- sapply(capermodels, function (x) summary(x)$adj.r.squared)
LL <- sapply(capermodels, function (x) x$model$log.lik) 
ks <- sapply(capermodels, function (x) x$k) 

models_table <- data.frame(formulas, ks, LL, aics, aiccs, r2s, ar2s) %>%
  rename(Model = formulas, n = ks, AIC = aics,
         AICc = aiccs, R_sq = r2s, adj_R_sq = ar2s) %>%
  mutate(Model = as.character(formulas), 
         LL = round(LL, 1), 
         AIC = round(AIC, 1), AICc = round(AICc, 1), 
         dAIC = round(AIC - min(AIC), 2),
         dAICc = round(AICc - min(AICc), 2),
         R_sq = round(R_sq, 2), adj_R_sq = round(adj_R_sq, 2), 
         Weights = round(exp(-dAICc/2)/sum(exp(-dAICc/2)), 3)) 

models_table

### Top-ranked model diagnostic plots and analyses---------------------------------
topmod <- which(models_table$dAICc == 0)
bestmodel <- capermodels[[topmod]]

# Re-running the model using nlme:gls for diagnostics
# parameter estimates are obtained
# corPagel() is used as it's equivalent to `lambda = "ML"` in pgls
bestmodel_gls <- nlme::gls(formula(bestmodel), data = cd$data, 
                           correlation = corPagel(1, cd$phy)) 

# almost identical coefficent estimates
rbind(coef(bestmodel), coef(bestmodel_gls))

# residual plots almost identical as well, homoscedastic
plot(resid(bestmodel) ~ fitted(bestmodel)) # pgls
plot(bestmodel_gls) # gls

# estimate Variance inflation factors (on gls model) to test for collinearity
vif(bestmodel_gls)

# qqplot (on gls model) looks normal
qqnorm(bestmodel_gls)

# profile plot of lambda
plot(pgls.profile(bestmodel)) # lambda closer to 1 indicates strong phylogenetic signal

# fit full 32 candidate pgls models---------------
# outputs reported in supporting information

# combine phylogeny with dataset for use in pgls function
cd <- comparative.data(phy, dat, names.col = "BinomName") 

# 24 pgls models from Barrowclift et al., 2023 + 4 offspring size models

# intercept only
zslope_0 <- pgls(log(rmax) ~ 1, data = cd, lambda = "ML") 

# rmax varies with mass only 
zslope_wt <- pgls(log(rmax) ~ log_wt, data = cd, lambda = "ML") 

# depth only
zslope_depth <- pgls(log(rmax) ~ depth_scaled, data = cd, lambda = "ML")

# temp only
zslope_invtemp <- pgls(log(rmax) ~ invtemp_scaled, data = cd, lambda = "ML") 

# mass + depth
zslope_wt_depth <- pgls(log(rmax) ~ log_wt + depth_scaled,
                        data = cd, lambda = "ML")

# mass + temp
zslope_wt_invtemp <- pgls(log(rmax) ~ log_wt + invtemp_scaled,
                          data = cd, lambda = "ML")

# mass * depth 
zslope_wt_x_depth <- pgls(log(rmax) ~ log_wt * depth_scaled, 
                          data = cd, lambda = "ML")

# mass * temp 
zslope_wt_x_invtemp <- pgls(log(rmax) ~ log_wt * invtemp_scaled, 
                            data = cd, lambda = "ML")

# varying intercepts - rays (coded 0) and skates (coded 1)
zslope_2 <- pgls(log(rmax) ~ 1 + Order, data = cd, lambda = "ML") 

# mass + order
zslope_wt_o <- pgls(log(rmax) ~ log_wt + Order, data = cd, lambda = "ML") 

# depth + order
zslope_depth_o <- pgls(log(rmax) ~ depth_scaled + Order, data = cd, lambda = "ML")

# temp + order
zslope_temp_o <- pgls(log(rmax) ~ invtemp_scaled + Order, data = cd, lambda = "ML")

# mass + depth + order
zslope_wt_depth_o <- pgls(log(rmax) ~ log_wt + depth_scaled + Order,
                          data = cd, lambda = "ML")

# mass + temp + order
zslope_wt_invtemp_o <- pgls(log(rmax) ~ log_wt + invtemp_scaled + Order,
                            data = cd, lambda = "ML")

# mass * depth + order
zslope_wt_x_depth_o <- pgls(log(rmax) ~ log_wt * depth_scaled + Order, 
                            data = cd, lambda = "ML")

# mass * temp + order
zslope_wt_x_invtemp_o <- pgls(log(rmax) ~ log_wt  * invtemp_scaled + Order,
                              data = cd, lambda = "ML")

# mass + temp + depth
zslope_wt_invtemp_depth <- pgls(log(rmax) ~ log_wt + invtemp_scaled + depth_scaled,
                                data = cd, lambda = "ML")

# mass + temp * depth
zslope_wt_invtemp_x_depth <- pgls(log(rmax) ~ log_wt  + invtemp_scaled * depth_scaled,
                                  data = cd, lambda = "ML")

# rmax varies with temp-depth index only
zslope_PC1 <- pgls(log(rmax) ~ PC1, data = cd, lambda = "ML")

# rmax varies with mass and temp-depth index
zslope_wt_PC1 <- pgls(log(rmax) ~ log_wt + PC1,
                      data = cd, lambda = "ML")

# rmax varies with mass and temp-depth index, and the effect of mass scaling coefficient varies with temp-depth index
zslope_wt_x_PC1 <- pgls(log(rmax) ~ log_wt * PC1,
                        data = cd, lambda = "ML")

# rmax varies with temp-depth index, and order
zslope_PC1_o <- pgls(log(rmax) ~ PC1 + Order,
                     data = cd, lambda = "ML")

# rmax varies with mass, temp-depth index, and order
zslope_wt_PC1_o <- pgls(log(rmax) ~ log_wt + PC1 + Order,
                        data = cd, lambda = "ML")

# rmax varies with mass, temp-depth index, and order, and the effect of mass scaling coefficient varies with the temperature-depth index
zslope_wt_x_PC1_o <- pgls(log(rmax) ~ log_wt * PC1 + Order,
                          data = cd, lambda = "ML")

# add in 6 models with offspring mass----------------
# rmax varies with offspring size
zslope_off <- pgls(log(rmax) ~ log_offspring,
                   data = cd, lambda = "ML")

# rmax varies with adult mass and offspring mass
zslope_wt_off <- pgls(log(rmax) ~ log_wt  + log_offspring,
                      data = cd, lambda = "ML")

# would expect interaction where larger mother larger offspring size
# rmax varies with adult mass and offspring size, and the effect of mass scaling coefficient varies with offspring size
zslope_wt_x_off <- pgls(log(rmax) ~ log_wt  * log_offspring,
                        data = cd, lambda = "ML")

# rmax varies with adult mass, temperature, and offspring size
zslope_wt_invtemp_off <- pgls(log(rmax) ~ log_wt + log_offspring + invtemp_scaled,
                              data = cd, lambda = "ML")

# rmax varies with offspring size and temperature
zslope_invtemp_off <- pgls(log(rmax) ~ log_offspring + invtemp_scaled,
                           data = cd, lambda = "ML")

# rmax varies with offspring size and depth
zslope_depth_off <- pgls(log(rmax) ~ log_offspring + depth_scaled,
                         data = cd, lambda = "ML")

# 2 additional models using adult:offspring size ratio-----------
# rmax varies with size ratio only 
zslope_sr <- pgls(log(rmax) ~ log(sizeratio), data = cd, lambda = "ML") 

# rmax varies with size ratio and temperature
zslope_sr_invtemp <- pgls(log(rmax) ~ log(sizeratio) + invtemp_scaled,
                          data = cd, lambda = "ML")

#  saving all 32 models in list
capermodels <- list(zslope_0, # 24 original models
                    zslope_wt, 
                    zslope_depth, 
                    zslope_invtemp,
                    zslope_wt_depth,
                    zslope_wt_invtemp, 
                    zslope_wt_x_depth,
                    zslope_wt_x_invtemp, 
                    zslope_2,
                    zslope_wt_o, 
                    zslope_depth_o,
                    zslope_temp_o,
                    zslope_wt_depth_o,
                    zslope_wt_invtemp_o,
                    zslope_wt_x_depth_o,
                    zslope_wt_x_invtemp_o,
                    zslope_wt_invtemp_depth, 
                    zslope_wt_invtemp_x_depth,
                    zslope_PC1, 
                    zslope_wt_PC1,
                    zslope_wt_x_PC1,
                    zslope_PC1_o,
                    zslope_wt_PC1_o,
                    zslope_wt_x_PC1_o,
                    zslope_off, # 6 offspring size models
                    zslope_wt_off,
                    zslope_wt_x_off,
                    zslope_wt_invtemp_off,
                    zslope_invtemp_off,
                    zslope_depth_off,
                    zslope_sr, # 2 size ratio models
                    zslope_sr_invtemp
) 

# manually extracting information from each model object
aics <- sapply(capermodels, function (x) bbmle::AIC(x)) 
aiccs <- sapply(capermodels, function (x) x$aicc) 
formulas <- sapply(capermodels,   # tidying formulas for easier reading
                   function (x) deparse(formula(x), width.cutoff = 90L)) %>%
  sub("^log\\(\\w+\\)\\s\\~\\s", "", .) %>% 
  gsub("log", "ln", .)  %>% # to have ln instead of log in formula
  gsub("_scaled", "", .) %>%
  gsub("\\_offspring", "(Moffspring)", .)  %>%
  gsub("\\_wt", "(M)", .)  %>%
  gsub("sizeratio", "M/Moffspring", .)  %>%
  gsub("O", "o", .) # so order is a lower case
r2s <- sapply(capermodels, function (x) summary(x)$r.squared) 
ar2s <- sapply(capermodels, function (x) summary(x)$adj.r.squared)
LL <- sapply(capermodels, function (x) x$model$log.lik) 
ks <- sapply(capermodels, function (x) x$k) 

models_table <- data.frame(formulas, ks, LL, aics, aiccs, r2s, ar2s) %>%
  rename(Model = formulas, n = ks, AIC = aics,
         AICc = aiccs, R_sq = r2s, adj_R_sq = ar2s) %>%
  mutate(Model = as.character(formulas), 
         LL = round(LL, 1), 
         AIC = round(AIC, 1), AICc = round(AICc, 1), 
         dAIC = round(AIC - min(AIC), 2),
         dAICc = round(AICc - min(AICc), 2),
         R_sq = round(R_sq, 2), adj_R_sq = round(adj_R_sq, 2), 
         Weights = round(exp(-dAICc/2)/sum(exp(-dAICc/2)), 3)) 

models_table