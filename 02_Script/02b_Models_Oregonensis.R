# Oregonensis phenology analysis ------
# 24 Sep 21 JLB
# updated: Aug 23 - JMH

# libraries ----

library(tidyverse)
library(dplyr)
library(mgcv)
library(MuMIn)

# you also need brms, but we'll load that later because of its overlap with mgcv

# LOAD DATA ----
# Import data

ZP <-  read.csv(file.path(here::here("03_GeneratedData/"),"01_LEPASzpData.csv" ), row.names = 1)

# Tweaking logged retrocurva values so that zeros = NA's
# Calculating binary presence/absence of species
zpJLB2 <- ZP %>% 
  mutate(log_Soregonensis = ifelse(Soregonensis > 0, log(Soregonensis), as.numeric(NA)),
         # presence = 1; absense = 0; NAs should be NAs
         pa_Soregonensis = ifelse(Soregonensis > 0, 1, 
                                  ifelse(Soregonensis == 0, 0, as.numeric(NA))),
         YearF = as.factor(Y),
         Sample_site = as.factor(Sample_site)) %>% 
  rename(Year = Y) %>% 
  filter(Sample_site != "27-918", Sample_site != "29-905", Sample_site !=  "36-873",Sample_site !=  "37-890")



# Plot oregonensis data -------

# presence
ggplot(zpJLB2,
       aes(y = log_Soregonensis, x= DOY, color = YearF)) +
  geom_line()+
  facet_wrap(vars(Sample_site))

ggplot(zpJLB2,
      aes(y = log_Soregonensis, x= DOY, color = Sample_site)) +
  geom_point()+
  facet_wrap("YearF") 

# presence/absence
ggplot(zpJLB2, aes(y = pa_Soregonensis, x= DOY, color = Sample_site)) +
  geom_point()+
  #facet_wrap(vars(Sample_site)) +
  stat_smooth(se = FALSE)

# Model selection -----

# Presence model selection vv -------

jno1 <- uGamm(log_Soregonensis ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"), 
             REML = FALSE,
             data = zpJLB2,
             class = "gamm")

jno2 <- uGamm(log_Soregonensis ~ s(DOY, by = Sample_site) + s(Year) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB2,
             class = "gamm")

jno3 <- uGamm(log_Soregonensis ~ s(DOY) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB2,
             class = "gamm")

jno4 <- uGamm(log_Soregonensis ~ s(DOY, by = Sample_site) + YearF + s(Sample_site, bs = "re"),
             REML = log_Soregonensis,
             data = zpJLB2,
             class = "gamm")

jno5 <- uGamm(log_Soregonensis ~ s(DOY) + YearF +  + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB2,
             class = "gamm")

jno6 <- uGamm(log_Soregonensis ~ te(DOY, Year) +  s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB2,
             class = "gamm")

jno7 <- uGamm(log_Soregonensis ~ te(DOY, Year, by = Sample_site) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB2,
             class = "gamm")

# dredge all the global models
dredge1no <- dredge(global.model = jno1, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredge2no <- dredge(global.model = jno2, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredge3no <- dredge(global.model = jno3, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredge4no <- dredge(global.model = jno4, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredge5no <- dredge(global.model = jno5, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredge6no <- dredge(global.model = jno6, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredge7no <- dredge(global.model = jno7, evaluate = TRUE, rank = AICc,extra = list("R^2"))

# combine all your dredging.
merge12no <- merge(dredge1no, dredge2no) 
merge34no <- merge(dredge3no, dredge4no)
merge56no <- merge(dredge5no, dredge6no)
merge1234no <- merge(merge12no,merge34no)
merge1t6no <- merge(merge1234no, merge56no)
mergeFno <- merge(merge1t6no, dredge7no)


#Determine the best presence model based on AICc
subset(mergeFno, delta <=5) 

# export for model sel table
write.csv(as.data.frame(mergeFno), "07_Tables/02b_OregPresenseModSelTable.csv", row.names = FALSE)

#Best model
oregPresA <- gamm(log_Soregonensis ~ te(DOY, Year)  + s(Sample_site, bs = "re"),
                 REML = TRUE,
                 data = zpJLB2)

library(mgcViz)
summary(oregPresA$gam)
check(getViz(oregPresA$gam))
plot(getViz(oregPresA$gam))

# Presence/absence model selection vv -------

jo1pan <- uGamm(pa_Soregonensis ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB2)

jo2pan <- uGamm(pa_Soregonensis ~ s(DOY, by = Sample_site) + s(Year) +  s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB2)


# jo3pan <- uGamm(pa_Soregonensis ~ s(DOY) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
#                family = binomial("logit"),
#                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
#                data = zpJLB2)


jo4pan <- uGamm(pa_Soregonensis ~ s(DOY, by = Sample_site) + YearF + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB2)


jo5pan <- uGamm(pa_Soregonensis ~ s(DOY) + YearF + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB2)

jo6pan <- uGamm(pa_Soregonensis ~ te(DOY, Year) + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB2)

# jo7pan <- uGamm(pa_Soregonensis ~ te(DOY, Year, by = Sample_site) + s(Sample_site, bs = "re"),
#                family = binomial("logit"),
#                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
#                data = zpJLB2)

dredge1pano <- dredge(global.model = jo1pan, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))
dredge2pano <- dredge(global.model = jo2pan, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))
dredge4pano <- dredge(global.model = jo4pan, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))
dredge5pano <- dredge(global.model = jo5pan, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))
dredge6pano <- dredge(global.model = jo6pan, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))

merge12pano <- merge(dredge1pano, dredge2pano)
merge4_12pano <- merge(merge12pano, dredge4pano)
merge5_124pano <- merge(merge4_12pano, dredge5pano)
mergeFpano <- merge(merge5_124pano, dredge6pano)


# save csv comparing all models
write.csv(as.data.frame(mergeFpano), "07_Tables/02b_OregPresAbsModSelTable.csv", row.names = FALSE)

#Determnie best model based on AICc
subset(mergeFpano, delta <= 10) 

# Build hurdle model -----

#load required package
library(brms)

# HURDLE MODEL
# This runs a hurdle model which means that it contains a "regular" gam and a binomial gam
# this model uses a bayesian framework, much of the details of that are hidden

oregMod3 <- brm(
  bf(Soregonensis ~  t2(DOY, Year) + s(Sample_site, bs = "re"), 
     hu ~ s(Sample_site, bs = "re")),
  data = zpJLB2, 
  family = hurdle_lognormal(),
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.99),
  iter = 2000)


# Assess hurdle model -----
#look at the model
summary(oregMod3)
# All the Rhats are 1, indicating that the model successfully "converged"

# baysian version of R2
bayes_R2(oregMod3)

# look at predicted v. observed

presOreg3 <- zpJLB2 %>% 
  droplevels() %>% 
  tidybayes::add_predicted_draws(oregMod3) 

# summarise to look at medians
presSO3 <- presOreg3 %>% 
  group_by(.row, Sample_site, Year, DOY) %>% 
  summarise(across(c(Soregonensis, .prediction), median))

#plot predicted v. observed 

ggplot(presSO3, aes(y = log(.prediction + 0.0001), x = log(Soregonensis + 0.0001), color = Sample_site)) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(vars(Sample_site))


# look at the MCMC chains. They are only showing us the last 1000 runs, after a 1000 run "burnin" period
plot(oregMod3)

# load required libraries
library(tidybayes)
library(bayesplot)

# Get conditional smooths from model
rmCSO3 <- brms::conditional_smooths(oregMod3)

plot(rmCSO3)[1]
plot(oregPres3$gam)



# PRESENCE/ABSENCE
# you have to page ahead one to see this one (label on legend should start with hu_")
# HIGHER values indicate a greater likelihood that _____ is ABSENT
plot(rmCSO3)[2]

# Differences among sites-------
# Get posteriors
rmPosOreg3 <- as.array(oregMod3)





#Save your data -------
save.image(file.path(here::here("04_SavedRImages"),"02b_Models_Oregonensis_Rdat"))

saveRDS(oregMod3, file.path(here::here("06_Rdat"),"02b_BRMSmodel_Oregonensis.RDS"))

