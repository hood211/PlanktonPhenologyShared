# Retrocurva phenology analysis ----
# 11 Aug 21 JMH
# updated Aug 23 JMH


# libraries ----
library(tidyverse)
library(mgcv)
library(MuMIn)
# you also need brms and related packages, but we'll load that later because of its overlap with mgcv

# Import data ----
ZP <-  read.csv(file.path(here::here("03_GeneratedData/"),"01_LEPASzpData.csv" ), row.names = 1)

# Tweaking logged retrocurva values so that zeros = NA's
# Calculating binary presence/absence of species
zpJMH <- ZP %>% 
  mutate(log_Dretrocurva = ifelse(Dretrocurva > 0, log(Dretrocurva), as.numeric(NA)),
         # presence = 1; absence = 0; NAs should be NAs
         pa_Dretrocurva = ifelse(Dretrocurva > 0, 1, 
                                ifelse(Dretrocurva == 0, 0, as.numeric(NA))),
         YearF = as.factor(Y),
         Sample_site = as.factor(Sample_site)) %>% 
  rename(Year = Y)

# Plot retrocurva data ------

# PRESENCE DATA
ggplot(zpJMH, aes(y = log_Dretrocurva, x= DOY, color = YearF)) +
  geom_line()+
  facet_wrap(vars(Sample_site))

# PRESENCE/ABSENCE DATA
ggplot(zpJMH, aes(y = pa_Dretrocurva, x= DOY)) +
  geom_point()+
  facet_wrap(vars(Sample_site)) +
  stat_smooth(se = FALSE)


# Model selection process -----
## Presence model selection-----

  jr1 <- uGamm(log_Dretrocurva ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"), 
               REML = FALSE,
               data = zpJMH,
               class = "gamm")
  
  
  jr2 <- uGamm(log_Dretrocurva ~ s(DOY, by = Sample_site) + s(Year) + s(Sample_site, bs = "re"),
               REML = FALSE,
               data = zpJMH,
               class = "gamm")
  
  # jr3 <- uGamm(log_Dretrocurva ~ s(DOY) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
  #              REML = FALSE,
  #              data = zpJMH,
  #              class = "gamm")

  
  jr4 <- uGamm(log_Dretrocurva ~ s(DOY, by = Sample_site) + YearF + s(Sample_site, bs = "re"),
               REML = FALSE,
               data = zpJMH,
               class = "gamm")
  
  jr5 <- uGamm(log_Dretrocurva ~ s(DOY) + YearF + s(Sample_site, bs = "re"),
               REML = FALSE,
               data = zpJMH,
               class = "gamm")

  jr6 <- uGamm(log_Dretrocurva ~ te(DOY, Year) + s(Sample_site, bs = "re"),
               REML = FALSE,
               data = zpJMH,
               class = "gamm")
  
  jr7 <- uGamm(log_Dretrocurva ~ te(DOY, Year, by = Sample_site) + s(Sample_site, bs = "re"),
               REML = FALSE,
               data = zpJMH,
               class = "gamm")

  # dredge all the global models
  dredge1 <- dredge(global.model = jr1, evaluate = TRUE, rank = AICc, extra = list("R^2"))
  dredge2 <- dredge(global.model = jr2, evaluate = TRUE, rank = AICc, extra = list("R^2"))
  # dredge3 <- dredge(global.model = jr3, evaluate = TRUE, rank = AICc, extra = list("R^2"))
  dredge4 <- dredge(global.model = jr4, evaluate = TRUE, rank = AICc, extra = list("R^2"))
  dredge5 <- dredge(global.model = jr5, evaluate = TRUE, rank = AICc, extra = list("R^2"))
  dredge6 <- dredge(global.model = jr6, evaluate = TRUE, rank = AICc, extra = list("R^2"))
  dredge7 <- dredge(global.model = jr7, evaluate = TRUE, rank = AICc, extra = list("R^2"))
  
  # combine dredging.
  merge12 <- merge(dredge1, dredge2) 
  merge124 <- merge(merge12, dredge4)
  merge1245 <- merge(merge124, dredge5)
  merge12456 <- merge(merge1245, dredge6)
  merge124567 <- merge(merge12456, dredge7)

  
  # Determine the most favorable model based on AICc
  subset(merge124567, delta <=2) 
  
  # export for model sel table
  write.csv(as.data.frame(merge124567), "07_Tables/02c_RetrocurvaPresenseModSelTable.csv", row.names = FALSE)
  
  # Best presence model
  retPres <- gamm(log_Dretrocurva ~ te(DOY, Year) + s(Sample_site, bs = "re"),
                  REML = TRUE,
                  data = zpJMH)
  
  # Let's take a look at it
  summary(retPres$gam)
  plot(retPres$gam)

  
# Presence/absence model selection -------
  
  # Build models
  #didn't coverge 
  # jr1pa <- uGamm(pa_Dretrocurva ~ s(DOY, by = Sample_site) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
  #                family = binomial("logit"),
  #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
  #                data = zpJMH)
  
  #didn't coverge 
  # jr2pa <- uGamm(pa_Dretrocurva ~ s(DOY, by = Sample_site) + s(Year) + s(Sample_site, bs = "re"),
  #                family = binomial("logit"),
  #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
  #                data = zpJMH)
  
  #didn't coverge 
  # jr3pa <- uGamm(pa_Dretrocurva ~ s(DOY) + s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
  #                family = binomial("logit"),
  #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
  #                data = zpJMH)
  
  #didn't coverge 
  # jr4pa <- uGamm(pa_Dretrocurva ~ s(DOY, by = Sample_site) + YearF * s(Sample_site, bs = "re"),
  #                family = binomial("logit"),
  #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
  #                data = zpJMH)
  
  #didn't coverge 
  # jr5pa <- uGamm(pa_Dretrocurva ~ s(DOY) + YearF + s(Sample_site, bs = "re"),
  #                family = binomial("logit"),
  #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
  #                data = zpJMH)
  
  jr6pa <- uGamm(pa_Dretrocurva ~ te(DOY, Year) + s(Sample_site, bs = "re"),
                 family = binomial("logit"),
                 control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
                 data = zpJMH)
  
  #didn't coverge - Just won't do it.
  # jr7pa <- uGamm(pa_Dretrocurva ~ te(DOY, Year, by = Sample_site) + s(Sample_site, bs = "re"),
  #                family = binomial("logit"),
  #                control = lmeControl(maxIter = 1e8, msMaxIter = 1e8, niterEM = 1000, niterPQL = 1000),
  #                data = zpJMH)
  
  
# by default the DOY * Y response surface model with site as a fixed factor is the most favorable model
Dredge_jr6pa <- dredge(global.model = jr6pa, evaluate = TRUE, rank = AICc)

subset(Dredge_jr6pa, delta <= 2)

# export for model sel table
write.csv(as.data.frame(Dredge_jr6pa), "07_Tables/02c_RetrocurvaPresAbsModSelTable.csv", row.names = FALSE)

# Build hurdle model -------

  # This requires the brms package
  library(brms)
  
  # HURDLE MODEL
  # This runs a hurdle model which means that it contains a "regular" gam and a binomial gam
  # this model uses a bayesian framework, much of the details of that are hidden

  retMod2 <- brm(
    # this is the regular gam with a tensor product smooth of DOY and Year - just interacting smooths.
    # we're giving it the "raw" data but it is modeling log retrocurva (after removing the zeros)
    bf(Dretrocurva ~  t2(DOY, Year) + s(Sample_site, bs = "re"), 
       # this is the binomial gam for just presence/absence
       hu ~ s(Sample_site, bs = "re")),
    data = zpJMH, 
    # this is where we prescribe the model
    family = hurdle_lognormal(),
    # this gives us 4 MCMC chains, which we distribute over the 4 cores on this computer (so they run simultaneously)
    control = list(adapt_delta = 0.99, max_treedepth = 30),
    chains = 4,
    cores = 4,
    # each chain has 2000 steps.
    iter = 3000)
  

# Assess hurdle model ----  
  
  # look at the model
  summary(retMod2)
  # All the Rhats are 1, indicating that the model successfully "converged"
  
  # baysian version of R2 - 0.53
  bayes_R2(retMod2)

  
  # look at predicted v. observed
  # this gives the estimates from the last 1000 draws from 4 chains - BIG DF
  pres <- zpJMH %>% 
    tidybayes::add_predicted_draws(retMod2) 

  
    # summarise to look at medians
  presS <- pres %>% 
    group_by(.row, Sample_site, Year, DOY) %>% 
    summarise(across(c(Dretrocurva, .prediction), median))
  
  #plot predicted v. observed - looks pretty good.
  # not a lot of bias issues = seems like it is underpredicting high values in 27 and 29
  ggplot(presS, aes(y = log(.prediction + 0.0001), x = log(Dretrocurva + 0.0001), color = Sample_site)) +
    geom_point()+
    geom_abline(intercept = 0, slope = 1) +
    facet_wrap(vars(Sample_site))

  
  # look at the MCMC chains; only showing us the last 1000 runs, after a 1000 run "burning" period
  plot(retMod2)
  

# Leave one out cross validation ----
  # https://www.monicaalexander.com/posts/2020-28-02-bayes_viz/
  library(loo)
  loo1 <- loo::loo(retMod2, save_psis = TRUE)
  # values over 0.7 are not good
  plot(loo1)

# Plot the smooths -------
  # Asking how retrocurva changes with DOY and Y
  
  # load required libraries
  library(tidybayes)
  library(bayesplot)
  
  # Model tells us there is a DOY*Y interaction--get conditional smooths from model
    rmCS <- brms::conditional_smooths(retMod2)
    
    # PRESENCE
    # this is saying that the peak biomass is occuring a bit earlier now than in the past
    # see how those highest isoclines (the yellow ones) are moving leftward as you go from 1995 to 2020?
    # Furthermore the period of high retro biomass is getting longer
    # see how the green lines are about DOY ~160-210 in 95 and are ~150-230 in 2019
    plot(rmCS)[1]

      # Check that this looks like our presence gam model - yup!
        plot(retPres$gam)
    
    # PRESENCE/ABSENCE
    # you have to page ahead one to see this one (label on legend should start with hu_")
    # HIGHER values indicate a greater likelihood that retro is ABSENT
    # This is telling us that retro is most likely absent before DOY 100
    # There is a period b/w DOY 150 and 175 where it is most likely to be present. 
    # That period seems to be shifting a bit earlier, but not much
    # there was a period between ~2003 and 2012 when retro was like present from DOY105 on
        # before and after that it was slightly more likely absent
    plot(rmCS)[2]
  
  
  
    
  
# Save your data -----  
# save.image(file.path(here::here("04_SavedRImages"),"02a_Models_Retrocurva_Rdat"))
    # load(file.path(here::here("04_SavedRImages"),"02a_Models_Retrocurva_Rdat"))
saveRDS(retMod2, file.path(here::here("06_Rdat"),"02a_BRMSmodel_Retrocurva.RDS"))
    
    