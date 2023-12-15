# Mesocyclops phenology analysis
# 24 Sep 21 JLB
# updated JMH Aug 23

######################
# LIBRARIES
######################
library(tidyverse)
library(dplyr)
library(mgcv)
library(MuMIn)
# you also need brms, but we'll load that later because of its overlap with mgcv

######################
# LOAD DATA
######################
#
# Import data

ZP <-  read.csv(file.path(here::here("03_GeneratedData/"),"01_LEPASzpData.csv" ), row.names = 1)

# Tweaking logged retrocurva values so that zeros = NA's
# Calculating binary presence/absence of species
zpJLB3 <- ZP %>% 
  mutate(log_Mesocyclops = ifelse(Mesocyclops > 0, log(Mesocyclops), as.numeric(NA)),
         # presence = 1; absense = 0; NAs should be NAs
         pa_Mesocyclops= ifelse(Mesocyclops > 0, 1, 
                                ifelse(Mesocyclops == 0, 0, as.numeric(NA))),
         YearF = as.factor(Y),
         Sample_site = as.factor(Sample_site)) %>% 
  rename(Year = Y)




######################
# PLOT MESOCYCLOPS DATA
######################

# PRESENCE
ggplot(zpJLB3, aes(y = log_Mesocyclops, x= DOY, color = YearF)) +
  geom_line()+
  facet_wrap(vars(Sample_site))

ggplot(zpJLB3, aes(y = log_Mesocyclops, x= DOY, color = Sample_site)) +
  geom_line()+
  facet_wrap("YearF")

# PRESENCE/ABSENCE
ggplot(zpJLB3, aes(y = pa_Mesocyclops, x= DOY)) +
  geom_point()+
  facet_wrap(vars(Sample_site)) +
  stat_smooth(se = FALSE)


######################
# MODEL SELECTION
######################
# We'll do this in 2 steps. 1. do model selection using gamm for the presence and presence/absence models
# Then combine those models into a single model in brm

# PRESENCE MODEL SELECTION

# Build models

jm1 <- uGamm(log_Mesocyclops ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"), 
             REML = FALSE,
             data = zpJLB3,
             class = "gamm")


jm2 <- uGamm(log_Mesocyclops ~ s(DOY, by = Sample_site) + s(Year) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB3,
             class = "gamm")

jm3 <- uGamm(log_Mesocyclops ~ s(DOY) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB3,
             class = "gamm")

jm4 <- uGamm(log_Mesocyclops ~ s(DOY, by = Sample_site) + YearF + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB3,
             class = "gamm")

jm5 <- uGamm(log_Mesocyclops ~ s(DOY) + YearF + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB3,
             class = "gamm")

jm6 <- uGamm(log_Mesocyclops ~ te(DOY, Year) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB3,
             class = "gamm")

jm7 <- uGamm(log_Mesocyclops ~ te(DOY, Year, by = Sample_site) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB3,
             class = "gamm")

# dredge all the global models
dredge1m <- dredge(global.model = jm1, evaluate = TRUE, rank = AICc)
dredge2m <- dredge(global.model = jm2, evaluate = TRUE, rank = AICc)
dredge3m <- dredge(global.model = jm3, evaluate = TRUE, rank = AICc)
dredge4m <- dredge(global.model = jm4, evaluate = TRUE, rank = AICc)
dredge5m <- dredge(global.model = jm5, evaluate = TRUE, rank = AICc)
dredge6m <- dredge(global.model = jm6, evaluate = TRUE, rank = AICc)
dredge7m <- dredge(global.model = jm7, evaluate = TRUE, rank = AICc)

# combine all your dredging.
merge12m <- merge(dredge1m, dredge2m) # puts our model selection results into one table
merge34m <- merge(dredge3m, dredge4m)
merge56m <- merge(dredge5m, dredge6m)
merge1234m <- merge(merge12m,merge34m)
merge1t6m <- merge(merge1234m, merge56m)
mergeFm <- merge(merge1t6m, dredge7m)

#Best model for presence: ~s(DOY) + Sample_site + YearF
subset(mergeFm, delta <=2) 

# export for model sel table
write.csv(as.data.frame(mergeFm), "07_Tables/02c_MesocyclopsPresenseModSelTable.csv", row.names = FALSE)

# Best presence model: ~s(DOY) + Sample_site + YearF
mesoPres <- gamm(log_Mesocyclops ~ s(DOY) + YearF + s(Sample_site, bs = "re"),
                 REML = TRUE,
                 data = zpJLB3)


# Let's take a look at it
summary(mesoPres$gam)
plot(mesoPres$gam)

# PRESENCE/ABSENCE MODEL SELECTION
jm1pa <- uGamm(pa_Mesocyclops ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB3)
 
# won't converge 
#jm2pa <- uGamm(pa_Mesocyclops ~ s(DOY, by = Sample_site) + s(Year) + s(Sample_site, bs = "re"),
#                family = binomial("logit"),
#                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
#                data = zpJLB3)

jm3pa <- uGamm(pa_Mesocyclops ~ s(DOY) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB3)

jm4pa <- uGamm(pa_Mesocyclops ~ s(DOY, by = Sample_site) + YearF + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB3)

jm5pa <- uGamm(pa_Mesocyclops ~ s(DOY) + YearF + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB3)

jm6pa <- uGamm(pa_Mesocyclops ~ te(DOY, Year) + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB3)

#didn't coverge - Just won't do it.
# jm7pa <- uGamm(pa_Mesocyclops ~ te(DOY, Year, by = Sample_site),
#                family = binomial("logit"),
#                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
#                data = zpJLB3)


# Dredge models
dredge1pam <- dredge(global.model = jm1pa, evaluate = TRUE, rank = AICc)
dredge3pam <- dredge(global.model = jm3pa, evaluate = TRUE, rank = AICc)
dredge4pam <- dredge(global.model = jm4pa, evaluate = TRUE, rank = AICc)
dredge5pam <- dredge(global.model = jm5pa, evaluate = TRUE, rank = AICc)
dredge6pam <- dredge(global.model = jm6pa, evaluate = TRUE, rank = AICc)
 
# puts our model selection results into one table
merge13pam <- merge(dredge1pam, dredge3pam) 
merge134pam <- merge(merge13pam, dredge4pam) 
merge1345pam <- merge(merge134pam, dredge5pam) 
merge13456pam <- merge(merge1345pam, dredge6pam) 
merge13456pam

#Best model for presence: ~s(DOY) + Sample_site + YearF
subset(merge13456pam, delta <=10) 

# export for model sel table
write.csv(as.data.frame(mergeFinpam), "07_Tables/02c_MesocyclopsPresAbsModSelTable.csv", row.names = FALSE)

######################
# BUILD HURDLE MODEL
######################
# This requires the brms package
library(brms)

# HURDLE MODEL
mesoMod2 <- brm(
  bf(Mesocyclops ~  s(DOY) + YearF + s(Sample_site, bs = "re"), 
     hu ~ s(DOY)),
  data = zpJLB3, 
  family = hurdle_lognormal(),
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 30),
  iter = 3000)

# Let's look at the model
# All the Rhats are 1, indicating that the model successfully "converged"
# we'll talk about the coefs below
summary(mesoMod2)


# baysian version of R2 - 0.54
bayes_R2(mesoMod2)


# look at predicted v. observed
# this gives the estimates from the last 1000 draws from 4 chains - BIG DF
presMeso <- zpJLB3 %>%
  droplevels() %>% 
  tidybayes::add_predicted_draws(mesoMod2)

# summarise to look at medians
presM <- presMeso %>%
  group_by(.row, Sample_site, Year, DOY) %>%
  summarise(across(c(Mesocyclops, .prediction), median))

#plot predicted v. observed
# not a lot of bias issues = seems like it is underpredicting high values in 27 and 29
ggplot(presM, aes(y = log(.prediction + 0.0001), x = log(Mesocyclops + 0.0001), color = Sample_site)) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(vars(Sample_site))



# Traces
plot(mesoMod2)


# Posterior predicted checks
pp_check(mesoMod2, type = "dens_overlay", nsamples = 100)
pp_check(mesoMod2, type = "stat", stat = "median")
pp_check(mesoMod2, type = "error_scatter_avg")
# 
# # Leave one out cross validation
# # https://www.monicaalexander.com/posts/2020-28-02-bayes_viz/
library(loo)
loo1 <- loo::loo(mesoMod2, save_psis = TRUE)
# values over 0.7 are not good
plot(loo1)


# overall, I'm feeling pretty good about this model so let's look at the resutls

library(tidybayes)
library(bayesplot)


# PLOT THE SMOOTHS: How does retrocurva change with DOY and Y?
# Model tells us there is a DOY*Y interaction
# Get conditional smooths from model
rmCSM <- brms::conditional_smooths(mesoMod2)

# PRESENCE
plot(rmCSM)[1]


# let's just check that this looks like our presence gam model
plot(mesoPres$gam)

# PRESENCE/ABSENCE
# you have to page ahead one to see this one (label on legend should start with hu_")
# HIGHER values indicate a greater likelihood that meso is ABSENT
# This is telling us that meso is most likely absent before DOY 100
# There is a period b/w DOY 150 and 175 where it is most likely to be present.
# That period seems to be shifting a bit earlier, but not much
# there was a period between ~2003 and 2012 when meso was like present from DOY105 on
# before and after that it was slightly more likely absent
plot(rmCSM)[2]

#######################


# DIFFERENCES AMONG SITES: How does retrocurva change among sites
# Model tells us there is no site * Y interaction
# Get posteriors
rmMesoPos <- as.array(mesoMod2)
# what's in this
dimnames(rmMesoPos)

# PRESENCE: how much do the sites differ?
# circle is the median, wide bar is 50% CI and narrow bar is 90% CI
# read this as a site-specific intercept
# less retrocurva at the off shore sites 27 and 29 (and 8 is that more like the central basin or Sandusky?),
# slightly more at 36 and 3
# retro likes it productive!
# note these values are on a log scale
MesoPresentPlotSites <- bayesplot::mcmc_intervals(rmMesoPos, pars = c(
  "s_sSample_site_1[1]",
  "s_sSample_site_1[2]",
  "s_sSample_site_1[3]",
  "s_sSample_site_1[4]",
  "s_sSample_site_1[5]",
  "s_sSample_site_1[6]",
  "s_sSample_site_1[7]",
  "s_sSample_site_1[8]"))


save.image(file.path(here::here("04_SavedRImages"),"02c_Models_Mesocyclops_Rdat"))

saveRDS(mesoMod2, file.path(here::here("06_Rdat"),"02c_BRMSmodel_Mesocyclops.RDS"))

