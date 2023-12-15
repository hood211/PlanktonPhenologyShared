# Calcultate doy of peak biomass
# updated JMH Sept 2023

# Library ----
library(tidyverse)
library(dplyr)
library(tidybayes)
library(modelr)

# Brms models ----
BrmsOreg <- readRDS(file.path(here::here("06_Rdat"),"02b_BRMSmodel_Oregonensis.RDS"))

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
## Retro ----
Pred_Oreg <- ZP2 %>% 
  data_grid(Year = seq_range(Year, by = 1), 
            DOY = seq_range(DOY, by = 1), 
            Sample_site = c("14-972", "16-970", "3-996", "8-994")) %>% 
  add_epred_draws(BrmsOreg, value = "pOregBM", dpar = T) %>% 
  group_by(Year, DOY, Sample_site) %>% 
  median_qi(pOregBM)


### Predict timing of peak biomass ----
Pred_OregMax <- Pred_Oreg %>% 
  filter(DOY > 125 & DOY < 275) %>% 
  group_by(Year, Sample_site) %>% 
  slice_max(pOregBM, n = 1)

ggplot(Pred_OregMax, aes(y = DOY, x= Year, color = Sample_site)) +
  geom_point()

# export
# write.csv(Pred_OregMax, file.path(here::here("03_GeneratedData/"),"03b_GAMMpredictPeakBM_Oreg.csv" ))
# write.csv(Pred_Oreg, file.path(here::here("03_GeneratedData/"),"03b_GAMMpredictions_Oreg.csv" ))

# save.image(file.path(here::here("04_SavedRImages"),"03b_GAMMpredictPeakBM_Oreg_Rdat"))
# load(file.path(here::here("04_SavedRImages"),"03b_GAMMpredictPeakBM_Oreg_Rdat"))
