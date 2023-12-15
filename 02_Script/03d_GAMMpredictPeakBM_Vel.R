# Calcultate doy of peak biomass
# updated JMH Sept 2023
# updated JMH Oct 2023

# Library ----
library(tidyverse)
library(dplyr)
library(tidybayes)
library(modelr)

# Brms models ----
BrmsVel <- readRDS(file.path(here::here("06_Rdat"),"02d_BRMSmodel_Veligers.RDS"))

# Raw data ----
ZP <-  read.csv(file.path(here::here("03_GeneratedData/"),"01_LEPASzpData.csv" ), row.names = 1)

ZP2 <- ZP %>% 
  mutate(log_Dretrocurva = ifelse(Dretrocurva > 0, log(Dretrocurva), as.numeric(NA)),
         # presence = 1; absence = 0; NAs should be NAs
         pa_Dretrocurva = ifelse(Dretrocurva > 0, 1, 
                                 ifelse(Dretrocurva == 0, 0, as.numeric(NA))),
         log_Soregonensis = ifelse(Soregonensis > 0, log(Soregonensis), as.numeric(NA)),
         # presence = 1; absense = 0; NAs should be NAs
         pa_Soregonensis = ifelse(Soregonensis > 0, 1, 
                                  ifelse(Soregonensis == 0, 0, as.numeric(NA))),
         log_Mesocyclops = ifelse(Mesocyclops > 0, log(Mesocyclops), as.numeric(NA)),
         # presence = 1; absense = 0; NAs should be NAs
         pa_Mesocyclops= ifelse(Mesocyclops > 0, 1, 
                                ifelse(Mesocyclops == 0, 0, as.numeric(NA))),
         log_Veliger = ifelse(Veliger > 0, log(Veliger), as.numeric(NA)),
         # presence = 1; absense = 0; NAs should be NAs
         pa_Veliger = ifelse(Veliger > 0, 1, 
                             ifelse(Veliger == 0, 0, as.numeric(NA))),
         YearF = as.factor(Y),
         Sample_site = as.factor(Sample_site),
         Region = ifelse(Sample_site %in% c("36-873","37-890","29-905","27-918", "14-972"), "Unimodal","Bimodal"),
         Region = as.factor(Region)) %>% 
  rename(Year = Y)

# Predict responses ----

Pred_Vel95 <- ZP2 %>% 
  data_grid(Year = 1995, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "1995") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel96 <- ZP2 %>% 
  data_grid(Year = 1996, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "1996") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel97 <- ZP2 %>% 
  data_grid(Year = 1997, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "1997") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel98 <- ZP2 %>% 
  data_grid(Year = 1998, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "1998") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel99 <- ZP2 %>% 
  data_grid(Year = 1999, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "1999") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel00 <- ZP2 %>% 
  data_grid(Year = 2000, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2000") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel01 <- ZP2 %>% 
  data_grid(Year = 2001, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2001") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel02 <- ZP2 %>% 
  data_grid(Year = 2002, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2002") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel03 <- ZP2 %>% 
  data_grid(Year = 2003, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2003") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel04 <- ZP2 %>% 
  data_grid(Year = 2004, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2004") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel05 <- ZP2 %>% 
  data_grid(Year = 2005, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2005") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel06 <- ZP2 %>% 
  data_grid(Year = 2006, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2006") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel07 <- ZP2 %>% 
  data_grid(Year = 2007, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2007") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel08 <- ZP2 %>% 
  data_grid(Year = 2008, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2008") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel09 <- ZP2 %>% 
  data_grid(Year = 2009, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2009") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel10 <- ZP2 %>% 
  data_grid(Year = 2010, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2010") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel11 <- ZP2 %>% 
  data_grid(Year = 2011, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2011") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel12 <- ZP2 %>% 
  data_grid(Year = 2012, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2012") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel13 <- ZP2 %>% 
  data_grid(Year = 2013, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2013") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel14 <- ZP2 %>% 
  data_grid(Year = 2014, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2014") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel15 <- ZP2 %>% 
  data_grid(Year = 2015, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2015") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel16 <- ZP2 %>% 
  data_grid(Year = 2016, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2016") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel17 <- ZP2 %>% 
  data_grid(Year = 2017, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2017") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel18 <- ZP2 %>% 
  data_grid(Year = 2018, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2018") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel19 <- ZP2 %>% 
  data_grid(Year = 2019, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2019") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)

Pred_Vel20 <- ZP2 %>% 
  data_grid(Year = 2020, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            Region = c("Unimodal", "Bimodal"),
            YearF = "2020") %>% 
  add_epred_draws(BrmsVel, value = "pVelBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site, Region) %>% 
  median_qi(pVelBM)


Pred_Vel <- rbind(Pred_Vel95,
                  Pred_Vel96,
                  Pred_Vel97,
                  Pred_Vel98,
                  Pred_Vel99,
                  Pred_Vel00,
                  Pred_Vel01,
                  Pred_Vel02,
                  Pred_Vel03,
                  Pred_Vel04,
                  Pred_Vel05,
                  Pred_Vel06,
                  Pred_Vel07,
                  Pred_Vel08,
                  Pred_Vel09,
                  Pred_Vel10,
                  Pred_Vel11,
                  Pred_Vel12,
                  Pred_Vel13,
                  Pred_Vel14,
                  Pred_Vel15,
                  Pred_Vel16,
                  Pred_Vel17,
                  Pred_Vel18,
                  Pred_Vel19,
                  Pred_Vel20) %>% 
  filter(DOY > 125 & DOY < 275) %>% 
  # remove predictions for site x region comparisons that don't exist
  filter(!(Region == "Bimodal" & Sample_site %in% c("36-873","37-890","29-905","27-918", "14-972"))) %>% 
  filter(!(Region == "Unimodal" & Sample_site %in% c("16-970", "3-996", "8-994"))) 


ggplot(Pred_Vel,  aes(y = log10(pVelBM), x = DOY, color = as.factor(Year))) +
  geom_line() +
  facet_grid(.~Sample_site) +
  geom_hline(yintercept = log10(10))

### Predict timing of peak biomass ----
Pred_VelMax_Uni <- Pred_Vel %>% 
  filter(DOY > 125 & DOY < 275) %>% 
  filter(Region == "Unimodal") %>% 
  group_by(Year, Sample_site) %>% 
  slice_max(pVelBM, n = 1)  %>% 
  mutate(Moment = "Pk1")

Pred_VelMax_Bi1 <- Pred_Vel %>% 
  filter(DOY > 125 & DOY < 275) %>% 
  filter(Region == "Bimodal" & DOY < 180) %>% 
  group_by(Year, Sample_site) %>% 
  slice_max(pVelBM, n = 1) %>% 
  mutate(Moment = "Pk1")

Pred_VelMax_Bi2 <- Pred_Vel %>% 
  filter(DOY > 125 & DOY < 275) %>% 
  filter(Region == "Bimodal" & DOY > 180) %>% 
  group_by(Year, Sample_site) %>% 
  slice_max(pVelBM, n = 1)  %>% 
  mutate(Moment = "Pk2")

Pred_VelMax <- rbind(Pred_VelMax_Uni, Pred_VelMax_Bi1, Pred_VelMax_Bi2)

ggplot(Pred_VelMax, aes(y = DOY, x= Year, color = Sample_site, shape = Moment)) +
  geom_point() +
  geom_line()

# export
# write.csv(Pred_Vel, file.path(here::here("03_GeneratedData/"),"03c_GAMMpredictions_Vel.csv" ))
# write.csv(Pred_VelMax, file.path(here::here("03_GeneratedData/"),"03b_GAMMpredictPeakBM_Vel.csv" ))

# save.image(file.path(here::here("04_SavedRImages"),"03d_GAMMpredictions_Vel.RData"))
# load(file.path(here::here("04_SavedRImages"),"03d_GAMMpredictions_Vel.RData"))
