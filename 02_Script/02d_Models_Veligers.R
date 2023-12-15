# Veligers phenology analysis --------
# 5 Oct 21 JLb
# updated JMH Sep 1
# updated JMH Oct 23

#Load required libraries
library(tidyverse)
library(mgcv)
library(MuMIn)
library(gamm4)
library(dplyr)


# Import data ----
ZP <-  read.csv(file.path(here::here("03_GeneratedData/"),"01_LEPASzpData.csv" ), row.names = 1)

# Tweaking logged  values so that zeros = NA's
# Calculating binary presence/absence of species
zpJLB5 <- ZP %>% 
  mutate(log_Veliger = ifelse(Veliger > 0, log(Veliger), as.numeric(NA)),
         # presence = 1; absence = 0; NAs should be NAs
         pa_Veliger = ifelse(Veliger > 0, 1, 
                             ifelse(Veliger == 0, 0, as.numeric(NA))),
         YearF = as.factor(Y),
         Sample_site = as.factor(Sample_site),
         Region = ifelse(Sample_site %in% c("36-873","37-890","29-905","27-918", "14-972"), "Unimodal","Bimodal"),
         Region = as.factor(Region)) %>% 
  rename(Year = Y)


# Plot veliger data ---------

# PRESENCE
ggplot(zpJLB5, aes(y = log_Veliger, x= DOY, color = YearF)) +
  geom_line()+
  facet_wrap(vars(Sample_site))

ggplot(zpJLB5, aes(y = log_Veliger, x= DOY, color = Sample_site)) +
  geom_line()+
  facet_wrap(vars(YearF))

# PRESENCE/ABSENCE
ggplot(zpJLB5, aes(y = pa_Veliger, x= DOY, color = YearF)) +
  geom_point()+
  facet_wrap(vars(YearF)) +
  stat_smooth(se = FALSE)

# Model Selection --------
# Presence Model Selection ------
jv1 <- uGamm(log_Veliger ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"), 
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv1_R <- uGamm(log_Veliger ~ s(DOY, by = Region) + s(Year, by = Region) + s(Sample_site, bs = "re"), 
               REML = FALSE,
               data = zpJLB5,
               class = "gamm")

jv1 <- uGamm(log_Veliger ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"), 
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv1_R <- uGamm(log_Veliger ~ s(DOY, by = Region) + s(Year, by = Region) + s(Sample_site, bs = "re"), 
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")


jv2 <- uGamm(log_Veliger ~ s(DOY, by = Sample_site) + s(Year) +  + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv2_R <- uGamm(log_Veliger ~ s(DOY, by = Region) + s(Year) +  + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv3 <- uGamm(log_Veliger~ s(DOY) + s(Year, by = Sample_site) +  + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv3_R <- uGamm(log_Veliger~ s(DOY) + s(Year, by = Region) +  + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

 
jv4 <- uGamm(log_Veliger ~ s(DOY, by = Sample_site) + YearF + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv4_R <- uGamm(log_Veliger ~ s(DOY, by = Region) + YearF  + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv5 <- uGamm(log_Veliger ~ s(DOY) + YearF +  + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv6 <- uGamm(log_Veliger ~ te(DOY, Year) +  s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv7 <- uGamm(log_Veliger ~ te(DOY, Year, by = Sample_site) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

jv7_R <- uGamm(log_Veliger ~ te(DOY, Year, by = Region) + s(Sample_site, bs = "re"),
             REML = FALSE,
             data = zpJLB5,
             class = "gamm")

# dredge all the global models
dredgev1 <- dredge(global.model = jv1, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev2 <- dredge(global.model = jv2, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev3 <- dredge(global.model = jv3, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev4 <- dredge(global.model = jv4, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev5 <- dredge(global.model = jv5, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev6 <- dredge(global.model = jv6, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev7 <- dredge(global.model = jv7, evaluate = TRUE, rank = AICc,extra = list("R^2"))

dredgev1_R <- dredge(global.model = jv1_R, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev2_R <- dredge(global.model = jv2_R, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev4_R <- dredge(global.model = jv3_R, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev5_R <- dredge(global.model = jv4_R, evaluate = TRUE, rank = AICc,extra = list("R^2"))
dredgev7_R <- dredge(global.model = jv7_R, evaluate = TRUE, rank = AICc,extra = list("R^2"))

# combine all your dredging.
merge12v <- merge(dredgev1, dredgev2) 
merge34v <- merge(dredgev3, dredgev4)
merge56v <- merge(dredgev5, dredgev6)
merge1234v <- merge(merge12v,merge34v)
merge1t6v <- merge(merge1234v, merge56v)
merge1t7v <- merge(merge1t6v, dredgev7)

merge1t7_1Rv <- merge(merge1t7v, dredgev1_R)
merge1t7_12Rv <- merge(merge1t7_1Rv, dredgev2_R)
merge1t7_124Rv <- merge(merge1t7_12Rv, dredgev4_R)
merge1t7_1245Rv <- merge(merge1t7_124Rv, dredgev5_R)
mergeFv <- merge(merge1t7_1245Rv, dredgev7_R)

# export for model sel table
write.csv(as.data.frame(mergeFv), "07_Tables/02c_VeligerPresenseModSelTable.csv", row.names = FALSE)


# determine best presence model based on AICc
subset(mergeFv, delta <=5) 

# Best presence model
velPres <- gamm(log_Veliger ~ s(DOY, by = Region) + YearF + s(Sample_site, bs = "re"),
                REML = TRUE,
                data = zpJLB5)

# Looking at the presence model
summary(velPres$gam)
plot(velPres$gam)


# Presence/absence model selection -------
# singular convergence
 # jv1pa <- uGamm(pa_Veliger ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
 #                family = binomial("logit"),
 #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
 #                data = zpJLB5)

jv1pa_R <- uGamm(pa_Veliger ~ s(DOY, by = Region) + s(Year, by = Region) + s(Sample_site, bs = "re"),
               family = binomial("logit"),
               control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
               data = zpJLB5)

#didn't work
 # jv2pa <- uGamm(pa_Veliger ~ s(DOY, by = Sample_site) + s(Year) + s(Sample_site, bs = "re"),
 #                family = binomial("logit"),
 #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
 #                data = zpJLB5)

#didn't work
# jv2pa_R <- uGamm(pa_Veliger ~ s(DOY, by = Region) + s(Year) + s(Sample_site, bs = "re"),
#                family = binomial("logit"),
#                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
#                data = zpJLB5)

#did work, doesn't now
 # jv3pa <- uGamm(pa_Veliger ~ s(DOY) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
 #                family = binomial("logit"),
 #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
 #                data = zpJLB5)

# didn't converge 
 # jv3pa_R <- uGamm(pa_Veliger ~ s(DOY) + s(Year, by = Region) + s(Sample_site, bs = "re"),
 #                family = binomial("logit"),
 #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 5000),
 #                data = zpJLB5)

# didn't converge 
 # jv4pa <- uGamm(pa_Veliger~ s(DOY, by = Sample_site) + YearF + s(Sample_site, bs = "re"),
 #                family = binomial("logit"),
 #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
 #                data = zpJLB5)
 
# didn't converge 
 # jv4pa_R <- uGamm(pa_Veliger~ s(DOY, by = Region) + YearF + s(Sample_site, bs = "re"),
 #                family = binomial("logit"),
 #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
 #                data = zpJLB5)

# didn't converge 
 # jv5pa <- uGamm(pa_Veliger ~ s(DOY) + YearF + s(Sample_site, bs = "re"),
 #                family = binomial("logit"),
 #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
 #                data = zpJLB5)

# jv6pa <- uGamm(pa_Veliger ~ te(DOY, Year) + s(Sample_site, bs = "re"),
#                family = binomial("logit"),
#                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
#                data = zpJLB5)

# didn't converge 
 # jv7pa <- uGamm(pa_Veliger ~ te(DOY, Year, by = Sample_site) + s(Sample_site, bs = "re"),
 #                family = binomial("logit"),
 #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
 #                data = zpJLB5)

# jv7pa_R <- uGamm(pa_Veliger ~ te(DOY, Year, by = Region) + s(Sample_site, bs = "re"),
#                family = binomial("logit"),
#                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
#                data = zpJLB5)


# Dredge models
dredge1vpa_R <- dredge(global.model = jv1pa_R, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))
# dredge3vpa_R <- dredge(global.model = jv2pa_R, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))
# dredge6vpa <- dredge(global.model = jv6pa, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))
# dredge7vpa_R <- dredge(global.model = jv7pa_R, evaluate = TRUE, rank = AICc,extra = c(R2 = function(x) summary(x$gam)$r.sq))


# STOPPED HERE, LOAD


# puts our model selection results into one table
# merge1R3R <- merge(dredge1vpa_R, dredge3vpa_R)
# merge1R3R6 <- merge(merge1R3R, dredge6vpa)
# merge1R3R67R <- merge(merge1R3R6, dredge7vpa_R)

# determine best model based on AICc
subset(dredge1vpa_R, delta <=2) 

# export for model sel table
write.csv(as.data.frame(merge1R3R67R), "07_Tables/02c_VeligerPresAbsModSelTable.csv", row.names = FALSE)




# Build hurdle model -------

# Load required library
library(brms)

# HURDLE MODEL
# This runs a hurdle model which means that it contains a "regular" gam and a binomial gam
# this model uses a bayesian framework, much of the details of that are hidden
velMod2 <- brm(
   bf(Veliger ~s(DOY, by = Region) + YearF + s(Sample_site, bs = "re"), 
      # pre data changes this included s(Y)
     hu ~ s(Sample_site, bs = "re")),
  data = zpJLB5, 
  family = hurdle_lognormal(),
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 30),
  iter = 4000)


# Assess hurdle model -------

# looking at the model
summary(velMod2)
# All the Rhats are 1, indicating that the model successfully "converged"

# baysian version of R2 
bayes_R2(velMod2)

# look at predicted v. observed
presVel <- zpJLB5 %>% 
  tidybayes::add_predicted_draws(velMod2) 

# summarise to look at medians
presVel2 <- presVel %>% 
  group_by(.row, Sample_site, Year, DOY) %>% 
  summarise(across(c(Veliger, .prediction), median))

#plot predicted v. observed 
ggplot(presVel2, aes(y = log(.prediction + 0.0001), x = log(Veliger + 0.0001), color = Sample_site)) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1) +
  facet_wrap(vars(Sample_site))



# look at the MCMC chains after a 1000 run "burnin" period
plot(velMod2)

# Plot the smooths -------

library(tidybayes)
library(bayesplot)

# Get conditional smooths from model
rmCSV <- brms::conditional_smooths(velMod2)

plot(rmCSV)


# save brms model
saveRDS(velMod2, file.path(here::here("06_Rdat"),"02d_BRMSmodel_Veligers.RDS"))

# Save/load
save.image(file.path(here::here("04_SavedRImages"),"02d_Models_Veligers_Rdat"))
# load(file.path(here::here("04_SavedRImages"),"02d_Models_Veligers_Rdat"))


