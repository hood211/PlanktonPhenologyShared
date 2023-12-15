# Calcultate doy of peak biomass
# updated JMH Sept 2023

# Library ----
library(tidyverse)
library(dplyr)
library(tidybayes)
library(modelr)

# Brms models ----
BrmsMeso <- readRDS(file.path(here::here("06_Rdat"),"02c_BRMSmodel_Mesocyclops.RDS"))

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
         Sample_site = as.factor(Sample_site)) %>% 
  rename(Year = Y)

# Predict responses ----
### Peak model results ----
Pred_Meso95 <- ZP2 %>% 
  data_grid(Year = 1995, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "1995") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso96 <- ZP2 %>% 
  data_grid(Year = 1996, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "1996") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso97 <- ZP2 %>% 
  data_grid(Year = 1997, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "1997") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso98 <- ZP2 %>% 
  data_grid(Year = 1998, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "1998") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso99 <- ZP2 %>% 
  data_grid(Year = 1999, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "1999") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso00 <- ZP2 %>% 
  data_grid(Year = 2000, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2000") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso01 <- ZP2 %>% 
  data_grid(Year = 2001, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2001") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso02 <- ZP2 %>% 
  data_grid(Year = 2002, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2002") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso03 <- ZP2 %>% 
  data_grid(Year = 2003, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2003") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso04 <- ZP2 %>% 
  data_grid(Year = 2004, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2004") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso05 <- ZP2 %>% 
  data_grid(Year = 2005, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2005") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso06 <- ZP2 %>% 
  data_grid(Year = 2006, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2006") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso07 <- ZP2 %>% 
  data_grid(Year = 2007, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2007") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso08 <- ZP2 %>% 
  data_grid(Year = 2008, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2008") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso09 <- ZP2 %>% 
  data_grid(Year = 2009, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2009") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso10 <- ZP2 %>% 
  data_grid(Year = 2010, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2010") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso11 <- ZP2 %>% 
  data_grid(Year = 2011, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2011") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso12 <- ZP2 %>% 
  data_grid(Year = 2012, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2012") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso13 <- ZP2 %>% 
  data_grid(Year = 2013, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2013") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso14 <- ZP2 %>% 
  data_grid(Year = 2014, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2014") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso15 <- ZP2 %>% 
  data_grid(Year = 2015, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2015") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso16 <- ZP2 %>% 
  data_grid(Year = 2016, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2016") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso17 <- ZP2 %>% 
  data_grid(Year = 2017, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2017") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso18 <- ZP2 %>% 
  data_grid(Year = 2018, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2018") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso19 <- ZP2 %>% 
  data_grid(Year = 2019, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2019") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)

Pred_Meso20 <- ZP2 %>% 
  data_grid(Year = 2020, 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("36-873", "37-890","27-918", "29-905" ,"14-972", "16-970", "3-996", "8-994"),
            YearF = "2020") %>% 
  add_epred_draws(BrmsMeso, value = "pMesoBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pMesoBM)


Pred_Meso <- rbind(Pred_Meso95,
                   Pred_Meso96,
                   Pred_Meso97,
                   Pred_Meso98,
                   Pred_Meso99,
                   Pred_Meso00,
                   Pred_Meso01,
                   Pred_Meso02,
                   Pred_Meso03,
                   Pred_Meso04,
                   Pred_Meso05,
                   Pred_Meso06,
                   Pred_Meso07,
                   Pred_Meso08,
                   Pred_Meso09,
                   Pred_Meso10,
                   Pred_Meso11,
                   Pred_Meso12,
                   Pred_Meso13,
                   Pred_Meso14,
                   Pred_Meso15,
                   Pred_Meso16,
                   Pred_Meso17,
                   Pred_Meso18,
                   Pred_Meso19,
                   Pred_Meso20) %>% 
  filter(DOY > 125 & DOY < 275)


ggplot(Pred_Meso,  aes(y = log10(pMesoBM), x = DOY, color = as.factor(Year))) +
  geom_line() +
  facet_wrap(vars(Sample_site)) +
  geom_hline(yintercept = log10(2.5))

### Predict timing of peak biomass ----
Pred_MesoMax <- Pred_Meso %>% 
  filter(DOY > 125 & DOY < 275) %>% 
  group_by(Year, Sample_site) %>% 
  slice_max(pMesoBM, n = 1)

ggplot(Pred_MesoMax, aes(y = DOY, x= Year, color = Sample_site)) +
  geom_point()

# export
# write.csv(Pred_Meso, file.path(here::here("03_GeneratedData/"),"03c_GAMMpredictions_Meso.csv" ))
# write.csv(Pred_MesoMax, file.path(here::here("03_GeneratedData/"),"03b_GAMMpredictPeakBM_Meso.csv" ))


# save.image(file.path(here::here("04_SavedRImages"),"03c_GAMMpredictions_Meso_Rdat"))
# load(file.path(here::here("04_SavedRImages"),"03c_GAMMpredictions_Meso_Rdat"))
