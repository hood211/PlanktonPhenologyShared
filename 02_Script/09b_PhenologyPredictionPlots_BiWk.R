# Predict phenol plots
# JMH Nov 2023

# libraries ----
library(tidyverse)
library(tidygam)
library(ggeffects)
library(viridis)
library(RColorBrewer)
library(ggrepel)

# GLMER RDS files ----

## Ret ----
RDS_Ret_15 <- readRDS(file.path(here::here("06_Rdat"),"08ab_GAMPhenRet_15per_BiWk_glmer.RDS"))
RDS_Ret_50 <- readRDS(file.path(here::here("06_Rdat"),"08ab_GAMPhenRet_50per_BiWk_glmer.RDS"))
RDS_Ret_85 <- readRDS(file.path(here::here("06_Rdat"),"08ab_GAMPhenRet_85per_BiWk_glmer.RDS"))
RDS_Ret_PB <- readRDS(file.path(here::here("06_Rdat"),"08ab_GAMPhenRet_PB_BiWk_glmer.RDS"))
 
## Oreg ----
RDS_Oreg_15 <- readRDS(file.path(here::here("06_Rdat"),"08bb_GAMPhenOreg_15per_BiWk_glmer.RDS"))
RDS_Oreg_50 <- readRDS(file.path(here::here("06_Rdat"),"08bb_GAMPhenOreg_50per_BiWk_glmer.RDS"))
RDS_Oreg_85 <- readRDS(file.path(here::here("06_Rdat"),"08bb_GAMPhenOreg_85per_BiWk_glmer.RDS"))
RDS_Oreg_PB <- readRDS(file.path(here::here("06_Rdat"),"08bb_GAMPhenOreg_PB_BiWk_glmer.RDS"))
     
## Meso ----
# no peak biomass
RDS_Meso_15 <- readRDS(file.path(here::here("06_Rdat"),"08cb_GAMPhenMeso_15per_BiWk_glmer.RDS"))
RDS_Meso_50 <- readRDS(file.path(here::here("06_Rdat"),"08cb_GAMPhenMeso_50per_BiWk_glmer.RDS"))
RDS_Meso_85 <- readRDS(file.path(here::here("06_Rdat"),"08cb_GAMPhenMeso_85per_BiWk_glmer.RDS"))
    
## Veligers ----
# No peak biomass
RDS_Vel_15 <- readRDS(file.path(here::here("06_Rdat"),"08db_GAMPhenVel_15per_BiWk_glmer.RDS"))
RDS_Vel_50 <- readRDS(file.path(here::here("06_Rdat"),"08db_GAMPhenVel_50per_BiWk_glmer.RDS"))
RDS_Vel_85 <- readRDS(file.path(here::here("06_Rdat"),"08db_GAMPhenVel_85per_BiWk_glmer.RDS"))

# Data files ----

## Ret ----
PhenPred_15per_Ret <- read.csv(here::here("03_GeneratedData", "08ab_GAMPhenRet_15per_BiWk_glmer.csv")) %>% 
                      mutate(plank_per_ha = as.numeric("NA"))
PhenPred_50per_Ret <- read.csv(here::here("03_GeneratedData", "08ab_GAMPhenRet_50per_BiWk_glmer.csv"))%>% 
  mutate(plank_per_ha = as.numeric("NA"))
PhenPred_85per_Ret <- read.csv(here::here("03_GeneratedData", "08ab_GAMPhenRet_85per_BiWk_glmer.csv"))%>% 
  mutate(plank_per_ha = as.numeric("NA"))
PhenPred_PB_Ret <- read.csv(here::here("03_GeneratedData", "08ab_GAMPhenRet_PB_BiWk_glmer.csv")) %>% 
  mutate(M = "Year",
    Sample_site = "WB") %>% 
  select("X","Taxa","Year","Sample_site","Moment","M","DOY","DOYwinStart", "AvgTemp","Byth","EdiblePP","InEdiblePP","Lepto","plank_per_ha")

PhenPred_ALL_Ret <- rbind(PhenPred_15per_Ret, PhenPred_50per_Ret, PhenPred_85per_Ret, PhenPred_PB_Ret)

## Oreg ----
PhenPred_15per_Oreg <- read.csv(here::here("03_GeneratedData", "08bb_GAMPhenOreg_15per_BiWk_glmer.csv")) %>% 
  mutate(plank_per_ha = as.numeric("NA"))
PhenPred_50per_Oreg <- read.csv(here::here("03_GeneratedData", "08bb_GAMPhenOreg_50per_BiWk_glmer.csv")) %>% 
  mutate(plank_per_ha = as.numeric("NA"))
PhenPred_85per_Oreg <- read.csv(here::here("03_GeneratedData", "08bb_GAMPhenOreg_85per_BiWk_glmer.csv")) %>% 
  mutate(plank_per_ha = as.numeric("NA"))
PhenPred_PB_Oreg <- read.csv(here::here("03_GeneratedData", "08bb_GAMPhenOreg_PB_BiWk_glmer.csv"))%>% 
  mutate(M = "Year",
         Sample_site = "WB") %>% 
  select("X","Taxa","Year","Sample_site","Moment","M","DOY","DOYwinStart", "AvgTemp","Byth","EdiblePP","InEdiblePP","Lepto","plank_per_ha")

PhenPred_ALL_Oreg <- rbind(PhenPred_15per_Oreg, PhenPred_50per_Oreg, PhenPred_85per_Oreg, PhenPred_PB_Oreg)

## Meso ----
# no peak biomass
PhenPred_15per_Meso <- read.csv(here::here("03_GeneratedData", "08cb_GAMPhenMeso_15per_BiWk_glmer.csv")) %>% 
  mutate(plank_per_ha = as.numeric("NA"))
PhenPred_50per_Meso <- read.csv(here::here("03_GeneratedData", "08cb_GAMPhenMeso_50per_BiWk_glmer.csv")) %>% 
  mutate(plank_per_ha = as.numeric("NA"))
PhenPred_85per_Meso <- read.csv(here::here("03_GeneratedData", "08cb_GAMPhenMeso_85per_BiWk_glmer.csv"))

PhenPred_ALL_Meso <- rbind(PhenPred_15per_Meso, PhenPred_50per_Meso, PhenPred_85per_Meso)

## Veligers ----
# No peak biomass
PhenPred_15per_Vel <- read.csv(here::here("03_GeneratedData", "08db_GAMPhenVel_15per_BiWk_glmer.csv"))
PhenPred_50per_Vel <- read.csv(here::here("03_GeneratedData", "08db_GAMPhenVel_50per_BiWk_glmer.csv")) %>% 
  mutate(plank_per_ha = as.numeric("NA"))
PhenPred_85per_Vel <- read.csv(here::here("03_GeneratedData", "08db_GAMPhenVel_85per_BiWk_glmer.csv")) %>% 
  mutate(plank_per_ha = as.numeric("NA"))

PhenPred_ALL_Vel <- rbind(PhenPred_15per_Vel, PhenPred_50per_Vel, PhenPred_85per_Vel)

## Combine all ----

PhenPred_ALL <-  rbind(PhenPred_ALL_Ret %>% 
                         mutate(MesoZPfood = as.numeric("NA"),
                                PPLthen10 = as.numeric("NA")) %>% 
                         select("X","Taxa","Year","Sample_site","Moment","M","DOY","DOYwinStart", 
                                "AvgTemp","Byth","EdiblePP","InEdiblePP","Lepto","plank_per_ha", "MesoZPfood", "PPLthen10"), 
                       PhenPred_ALL_Oreg %>% 
                         mutate(MesoZPfood = as.numeric("NA"),
                                PPLthen10 = as.numeric("NA")) %>% 
                         select("X","Taxa","Year","Sample_site","Moment","M","DOY","DOYwinStart", 
                                "AvgTemp","Byth","EdiblePP","InEdiblePP","Lepto","plank_per_ha", "MesoZPfood", "PPLthen10"),
                       PhenPred_ALL_Meso %>% 
                         mutate(PPLthen10 = as.numeric("NA")) %>% 
                         select("X","Taxa","Year","Sample_site","Moment","M","DOY","DOYwinStart", 
                                "AvgTemp","Byth","EdiblePP","InEdiblePP","Lepto","plank_per_ha", "MesoZPfood", "PPLthen10"),
                       PhenPred_ALL_Vel %>% 
                         mutate(MesoZPfood = as.numeric("NA")) %>% 
                         select("X","Taxa","Year","Sample_site","Moment","M","DOY","DOYwinStart", 
                                "AvgTemp","Byth","EdiblePP","InEdiblePP","Lepto","plank_per_ha", "MesoZPfood", "PPLthen10")) %>% 
                  select(-X) %>% 
                  mutate(Taxa = case_when(Taxa == "Dretrocurva" ~ "Ret",
                                          Taxa == "Soregonensis" ~ "Oreg",
                                          Taxa == "Mesocyclops" ~ "Meso",
                                          Taxa == "Veliger" ~ "Vel"),
                         Moment = ifelse(Moment == "PB1", "PB", Moment)) %>% 
                  rename(Season = Moment)

# Generate model pred ----
# ggpredict generates conditional predictions

## Ret ----
### 15% ----

gge_RDS_Ret_15_Byth <- ggpredict(RDS_Ret_15, 
                                 terms = "Byth",
                                 condition = c(EdiblePP = 0, Lepto = 0, InEdiblePP = 0, Sample_site = 0)) %>% 
                        as_tibble() %>% 
                        mutate(Predictor = "Byth",
                               group = "Null")

gge_RDS_Ret_15_InEdiblePP <- ggpredict(RDS_Ret_15, 
                                 terms = "InEdiblePP",
                                 condition = c(EdiblePP = 0, Lepto = 0, Byth = 0, Sample_site = 0))%>% 
                                as_tibble() %>% 
                                mutate(Predictor = "InEdiblePP",
                                       group = "Null")



gge_RDS_Ret_15_EdibleLepto <- ggpredict(RDS_Ret_15, 
                                       terms = c("EdiblePP", "Lepto [-0.68, 0.82]"),
                                       condition = c(InEdiblePP = 0, Byth = 0, Sample_site = 0)) %>% 
                                        as_tibble() %>% 
                                        mutate(Predictor = "EdiblePP_b_Lepto",
                                               group = ifelse(group == -0.68, "Low", "High"))


gge_RDS_Ret_15_df <- rbind(gge_RDS_Ret_15_Byth, gge_RDS_Ret_15_InEdiblePP, gge_RDS_Ret_15_EdibleLepto) %>% 
                        mutate(Taxa = "Dretrocurva",
                               Season = "per15",
                               Pval = case_when(Predictor == "Byth" ~ summary(RDS_Ret_15)$coefficients[2,4],
                                                Predictor == "InEdiblePP" ~ summary(RDS_Ret_15)$coefficients[5,4],
                                                Predictor == "EdiblePP_b_Lepto" ~ summary(RDS_Ret_15)$coefficients[6,4]),
                               PseudoR2 = cor(model.response(model.frame(RDS_Ret_15)),predict(RDS_Ret_15,type="response"))^2)


### 50% ----

gge_RDS_Ret_50_InEdiblePP <- ggpredict(RDS_Ret_50, 
                                 terms = "InEdiblePP",
                                 condition = c(Byth = 0, Sample_site = 0)) %>% 
                        as_tibble() %>% 
                        mutate(Predictor = "InEdiblePP",
                               group = "Null")

gge_RDS_Ret_50_Byth <- ggpredict(RDS_Ret_50, 
                                       terms = "Byth",
                                       condition = c(InEdiblePP = 0, Sample_site = 0)) %>% 
                        as_tibble() %>% 
                        mutate(Predictor = "Byth",
                               group = "Null")

gge_RDS_Ret_50_df <- rbind(gge_RDS_Ret_50_InEdiblePP, gge_RDS_Ret_50_Byth) %>% 
                      mutate(Taxa = "Dretrocurva",
                             Season = "per50",
                             Pval = case_when(Predictor == "Byth" ~ summary(RDS_Ret_50)$coefficients[2,4],
                                              Predictor == "InEdiblePP" ~ summary(RDS_Ret_50)$coefficients[3,4]),
                             PseudoR2 = cor(model.response(model.frame(RDS_Ret_50)),predict(RDS_Ret_50,type="response"))^2)

### 85% ----

gge_RDS_Ret_85_Byth <- ggpredict(RDS_Ret_85, 
                                       terms = "Byth",
                                       condition = c(InEdiblePP = 0, Lepto = 0, Sample_site = 0)) %>% 
                            as_tibble() %>% 
                            mutate(Predictor = "Byth",
                                   group = "Null")

gge_RDS_Ret_85_InEdiblePPLepto <- ggpredict(RDS_Ret_85, 
                                  terms = c("InEdiblePP", "Lepto [-0.71, 0.84]"),
                                 condition = c(Byth = 0, Sample_site = 0)) %>% 
                          as_tibble() %>% 
                          mutate(Predictor = "InEdiblePP_b_Lepto",
                                 group = ifelse(group == -0.71, "Low", "High"))

gge_RDS_Ret_85_df <- rbind(gge_RDS_Ret_85_Byth, gge_RDS_Ret_85_InEdiblePPLepto) %>% 
                      mutate(Taxa = "Dretrocurva",
                             Season = "per85",
                             Pval = case_when(Predictor == "Byth" ~ summary(RDS_Ret_85)$coefficients[2,4],
                             Predictor == "InEdiblePP_b_Lepto" ~ summary(RDS_Ret_85)$coefficients[5,4]),
                             PseudoR2 = cor(model.response(model.frame(RDS_Ret_85)),predict(RDS_Ret_85,type="response"))^2)

### PB -----

gge_RDS_Ret_PB_Byth <- ggpredict(RDS_Ret_PB, 
                                 terms = "Byth",
                                 condition = c(InEdiblePP = 0, plank_per_ha = 0, Sample_site = 0)) %>% 
                      as_tibble() %>% 
                      mutate(Predictor = "Byth",
                             group = "Null")

gge_RDS_Ret_PB_InEdPPPlank <- ggpredict(RDS_Ret_PB, 
                                 terms = c("InEdiblePP", "plank_per_ha [-0.74, 0.62]"),
                                 condition = c(Byth = 0, Sample_site = 0)) %>% 
                      as_tibble() %>% 
                      mutate(Predictor = "InEdiblePP_b_Plank",
                             group = ifelse(group == -0.74, "Low", "High"))

gge_RDS_Ret_PB_df <- rbind(gge_RDS_Ret_PB_Byth, gge_RDS_Ret_PB_InEdPPPlank) %>% 
                      mutate(Taxa = "Dretrocurva",
                             Season = "PB",
                             Pval = case_when(Predictor == "Byth" ~ summary(RDS_Ret_PB)$coefficients[2,4],
                                              Predictor == "InEdiblePP_b_Plank" ~ summary(RDS_Ret_PB)$coefficients[5,4]),
                             PseudoR2 = cor(model.response(model.frame(RDS_Ret_PB)),predict(RDS_Ret_PB,type="response"))^2)

### Combine ----
gge_RDS_Ret_df <- rbind(gge_RDS_Ret_15_df, gge_RDS_Ret_50_df, gge_RDS_Ret_85_df, gge_RDS_Ret_PB_df)

## Oreg ----

### 15% ----

gge_RDS_Oreg_15_AvgTemp <- ggpredict(RDS_Oreg_15, 
                                 terms = "AvgTemp",
                                 condition = c(Byth = 0, InEdiblePP = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "AvgTemp",
         group = "Null")

gge_RDS_Oreg_15_InEdiblePPByth <- ggpredict(RDS_Oreg_15, 
                                     terms = c("InEdiblePP", "Byth [-0.96, 0.79]"),
                                     condition = c(AvgTemp = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP_b_Byth",
         group = ifelse(group == -0.96, "Low", "High"))

gge_RDS_Oreg_15_df <- rbind(gge_RDS_Oreg_15_AvgTemp, gge_RDS_Oreg_15_InEdiblePPByth) %>% 
  mutate(Taxa = "Soregonensis",
         Season = "per15",
         Pval = case_when(Predictor == "AvgTemp" ~ summary(RDS_Oreg_15)$coefficients[2,4],
                          Predictor == "InEdiblePP_b_Byth" ~ summary(RDS_Oreg_15)$coefficients[5,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Oreg_15)),predict(RDS_Oreg_15,type="response"))^2) 

### 50% ----

gge_RDS_Oreg_50_Lepto <- ggpredict(RDS_Oreg_50, 
                                     terms = "Lepto",
                                     condition = c(Byth = 0, InEdiblePP = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "Lepto",
         group = "Null")

gge_RDS_Oreg_50_InEdiblePPByth <- ggpredict(RDS_Oreg_50, 
                                            terms = c("InEdiblePP", "Byth [-0.98, 0.84]"),
                                            condition = c(Lepto = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP_b_Byth",
         group = ifelse(group == -0.98, "Low", "High"))

gge_RDS_Oreg_50_df <- rbind(gge_RDS_Oreg_50_Lepto, gge_RDS_Oreg_50_InEdiblePPByth) %>% 
  mutate(Taxa = "Soregonensis",
         Season = "per50",
         Pval = case_when(Predictor == "Lepto" ~ summary(RDS_Oreg_50)$coefficients[4,4],
                          Predictor == "InEdiblePP_b_Byth" ~ summary(RDS_Oreg_50)$coefficients[5,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Oreg_50)),predict(RDS_Oreg_50,type="response"))^2) 


### 85% ----


gge_RDS_Oreg_85_TempByth <- ggpredict(RDS_Oreg_85, 
                                   terms = c("AvgTemp", "Byth [-0.96, 0.8]")) %>% 
  as_tibble() %>% 
  mutate(Predictor = "AvgTemp_b_Byth",
         group = ifelse(group == -0.96, "Low", "High"))

gge_RDS_Oreg_85_df <- rbind(gge_RDS_Oreg_85_TempByth) %>% 
  mutate(Taxa = "Soregonensis",
         Season = "per85",
         Pval = case_when(Predictor == "AvgTemp_b_Byth" ~ summary(RDS_Oreg_85)$coefficients[4,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Oreg_85)),predict(RDS_Oreg_85,type="response"))^2) 


### PB ----

gge_RDS_Oreg_PB_Byth <- ggpredict(RDS_Oreg_PB, 
                                  terms = "Byth",
                                  condition = c(InEdiblePP = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "Byth",
         group = "Null")

gge_RDS_Oreg_PB_InEdiblePP <- ggpredict(RDS_Oreg_PB, 
                                  terms = "InEdiblePP",
                                  condition = c(Byth = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP",
         group = "Null")

gge_RDS_Oreg_PB_df <- rbind(gge_RDS_Oreg_PB_Byth, gge_RDS_Oreg_PB_InEdiblePP) %>% 
  mutate(Taxa = "Soregonensis",
         Season = "PB",
         Pval = case_when(Predictor == "Byth" ~ summary(RDS_Oreg_PB)$coefficients[2,4],
                          Predictor == "InEdiblePP" ~ summary(RDS_Oreg_PB)$coefficients[3,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Oreg_PB)),predict(RDS_Oreg_PB,type="response"))^2) 

### Combine ----
gge_RDS_Oreg_df <- rbind(gge_RDS_Oreg_15_df, gge_RDS_Oreg_50_df, gge_RDS_Oreg_85_df, gge_RDS_Oreg_PB_df)

## Meso ----
# no peak biomass

### 15% ----

gge_RDS_Meso_15_AvgTemp <- ggpredict(RDS_Meso_15, 
                                         terms = "AvgTemp",
                                         condition = c(Byth = 0, InEdiblePP = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "AvgTemp",
         group = "Null")

gge_RDS_Meso_15_InEdiblePPByth <- ggpredict(RDS_Meso_15, 
                                     terms = c("InEdiblePP", "Byth [-0.8, 0.75]"),
                                     condition = c(AvgTemp = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP_b_Byth",
         group = ifelse(group == -0.8, "Low", "High"))

gge_RDS_Meso_15_df <- rbind(gge_RDS_Meso_15_AvgTemp, gge_RDS_Meso_15_InEdiblePPByth) %>% 
  mutate(Taxa = "Mesocyclops",
         Season = "per15",
         Pval = case_when(Predictor == "AvgTemp" ~ summary(RDS_Meso_15)$coefficients[2,4],
                          Predictor == "InEdiblePP_b_Byth" ~ summary(RDS_Meso_15)$coefficients[5,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Meso_15)),predict(RDS_Meso_15,type="response"))^2) 


### 50% ----

gge_RDS_Meso_50_InEdiblePPByth <- ggpredict(RDS_Meso_50, 
                                            terms = c("InEdiblePP", "Byth [-0.93, 0.81]"),
                                            condition = c(Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP_b_Byth",
         group = ifelse(group == -0.93, "Low", "High"))

gge_RDS_Meso_50_df <- rbind(gge_RDS_Meso_50_InEdiblePPByth) %>% 
  mutate(Taxa = "Mesocyclops",
         Season = "per50",
         Pval = case_when(Predictor == "InEdiblePP_b_Byth" ~ summary(RDS_Meso_50)$coefficients[4,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Meso_50)),predict(RDS_Meso_50,type="response"))^2) 

### 85% ----

gge_RDS_Meso_85_InEdiblePP <- ggpredict(RDS_Meso_85, 
                                            terms = c("InEdiblePP"),
                                            condition = c(plank_per_ha = 0, Byth = 0, MesoZPfood = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP",
         group = "Null")

gge_RDS_Meso_85_plank_per_ha <- ggpredict(RDS_Meso_85, 
                                        terms = c("plank_per_ha"),
                                        condition = c(InEdiblePP = 0, Byth = 0, MesoZPfood = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "plank_per_ha",
         group = "Null")

gge_RDS_Meso_85_MesoZPfoodByth <- ggpredict(RDS_Meso_85, 
                                          terms = c("MesoZPfood", "Byth [-0.8, 0.86]"),
                                          condition = c(InEdiblePP = 0, plank_per_ha = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "MesoZPfood_b_Byth",
         group = ifelse(group == -0.8, "Low", "High"))

gge_RDS_Meso_85_df <- rbind(gge_RDS_Meso_85_InEdiblePP, gge_RDS_Meso_85_plank_per_ha, gge_RDS_Meso_85_MesoZPfoodByth) %>% 
  mutate(Taxa = "Mesocyclops",
         Season = "per85",
         Pval = case_when(Predictor == "InEdiblePP" ~ summary(RDS_Meso_85)$coefficients[2,4],
                          Predictor == "plank_per_ha" ~ summary(RDS_Meso_85)$coefficients[3,4],
                          Predictor == "MesoZPfood_b_Byth" ~ summary(RDS_Meso_85)$coefficients[6,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Meso_85)),predict(RDS_Meso_85,type="response"))^2) 

### Combine ----
gge_RDS_Meso_df <- rbind(gge_RDS_Meso_15_df, gge_RDS_Meso_50_df, gge_RDS_Meso_85_df)

## Veligers ----
# No peak biomass
### 15% ----

gge_RDS_Vel_15_plank_per_ha <- ggpredict(RDS_Vel_15, 
                                         terms = "plank_per_ha",
                                         condition = c(InEdiblePP = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "plank_per_ha",
         group = "Null")

gge_RDS_Vel_15_InEdiblePP <- ggpredict(RDS_Vel_15, 
                                         terms = "InEdiblePP",
                                         condition = c(plank_per_ha = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP",
         group = "Null")

gge_RDS_Vel_15_df <- rbind(gge_RDS_Vel_15_plank_per_ha, gge_RDS_Vel_15_InEdiblePP) %>% 
  mutate(Taxa = "Veliger",
         Season = "per15",
         Pval = case_when(Predictor == "plank_per_ha" ~ summary(RDS_Vel_15)$coefficients[2,4],
                          Predictor == "InEdiblePP" ~ summary(RDS_Vel_15)$coefficients[3,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Vel_15)),predict(RDS_Vel_15,type="response"))^2) 


### 50% ----

gge_RDS_Vel_50_InEdiblePPByth <- ggpredict(RDS_Vel_50, 
                                            terms = c("InEdiblePP", "Byth [-0.88, 0.84]"),
                                            condition = c(Lepto = 0, PPLthen10 = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP_b_Byth",
         group = ifelse(group == -0.88, "Low", "High"))

gge_RDS_Vel_50_PPLthen10Lepto <- ggpredict(RDS_Vel_50, 
                                           terms = c("PPLthen10", "Lepto [-0.68, 0.83]"),
                                           condition = c(InEdiblePP = 0, Byth = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "PPLthen10_b_Lepto",
         group = ifelse(group == -0.68, "Low", "High"))

gge_RDS_Vel_50_df <- rbind(gge_RDS_Vel_50_InEdiblePPByth, gge_RDS_Vel_50_PPLthen10Lepto) %>% 
  mutate(Taxa = "Veliger",
         Season = "per50",
         Pval = case_when(Predictor == "InEdiblePP_b_Byth" ~ summary(RDS_Vel_50)$coefficients[6,4],
                          Predictor == "PPLthen10_b_Lepto" ~ summary(RDS_Vel_50)$coefficients[7,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Vel_50)),predict(RDS_Vel_50,type="response"))^2) 

### 85% ----

gge_RDS_Vel_85_InEdiblePP <- ggpredict(RDS_Vel_85, 
                                         terms = "InEdiblePP",
                                         condition = c(Byth = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "InEdiblePP",
         group = "Null")

gge_RDS_Vel_85_Byth <- ggpredict(RDS_Vel_85, 
                                       terms = "Byth",
                                       condition = c(InEdiblePP = 0, Sample_site = 0)) %>% 
  as_tibble() %>% 
  mutate(Predictor = "Byth",
         group = "Null")

gge_RDS_Vel_85_df <- rbind(gge_RDS_Vel_85_InEdiblePP, gge_RDS_Vel_85_Byth) %>% 
  mutate(Taxa = "Veliger",
         Season = "per85",
         Pval = case_when(Predictor == "InEdiblePP" ~ summary(RDS_Vel_85)$coefficients[2,4],
                          Predictor == "Byth" ~ summary(RDS_Vel_85)$coefficients[3,4]),
         PseudoR2 = cor(model.response(model.frame(RDS_Vel_85)),predict(RDS_Vel_85,type="response"))^2) 

### Combine ----
gge_RDS_Vel_df <- rbind(gge_RDS_Vel_15_df, gge_RDS_Vel_50_df, gge_RDS_Vel_85_df)

## Combine all ----
gge_RDS_ALL_df <- rbind(gge_RDS_Ret_df, gge_RDS_Oreg_df, gge_RDS_Meso_df, gge_RDS_Vel_df) %>% 
                  as_tibble() %>% 
                  mutate(across(group:Season, factor))  %>% 
                mutate(Taxa = case_when(Taxa == "Dretrocurva" ~ "Ret",
                                        Taxa == "Soregonensis" ~ "Oreg",
                                        Taxa == "Mesocyclops" ~ "Meso",
                                        Taxa == "Veliger" ~ "Vel"))
                  

# Median of Seasonal moments
PhenPred_ALL_medianSeas <- PhenPred_ALL %>% 
                            group_by(Taxa, Season) %>% 
                            summarise(medianDOY = median(DOY, na.rm = FALSE))

# Center on median of season
gge_RDS_ALL_df2 <- gge_RDS_ALL_df %>% 
                    left_join(PhenPred_ALL_medianSeas, by = c("Taxa", "Season")) %>% 
                    mutate(predictedC = predicted - medianDOY,
                           group = fct_relevel(as.factor(group), "Null", "Low", "High")) %>% 
                    mutate(Predictor = as.factor(Predictor),
                           Predictor = fct_recode(Predictor,
                                                  Temp = "AvgTemp",
                                                  "Temp x Byth" = "AvgTemp_b_Byth",
                                                  # PPngd = "EdiblePP",
                                                  "PPngd x Lepto" = "EdiblePP_b_Lepto",
                                                  PPgd = "InEdiblePP",
                                                  "PPgd x Lepto" = "InEdiblePP_b_Lepto",
                                                  "PPgd x Byth" = "InEdiblePP_b_Byth",
                                                  "PPgd x Plank" = "InEdiblePP_b_Plank",
                                                  # PPL10 = "PPLthen10",
                                                  "PPL10 x Lepto" = "PPLthen10_b_Lepto",
                                                  # ZPbms = "MesoZPfood", 
                                                  "ZPbms x Byth" = "MesoZPfood_b_Byth",
                                                  Lepto = "Lepto",
                                                  Byth = "Byth",
                                                  Plank = "plank_per_ha"),
                           Predictor = fct_relevel(Predictor, "Temp", "Temp x Byth",
                                                   "PPngd x Lepto",
                                                   "PPgd", "PPgd x Lepto", "PPgd x Byth", "PPgd x Plank",
                                                   "PPL10 x Lepto", "ZPbms x Byth",
                                                   "Lepto", "Byth", "Plank"),
                           Season = as.factor(Season),
                           Season = fct_relevel(Season, "per15", "per50", "PB", "per85"),
                           Taxa = as.factor(Taxa),
                           Taxa = fct_relevel(Taxa, "Ret", "Vel", "Oreg", "Meso")) 

# distinct values to include R2
gge_RDS_ALL_df2_distinctR2 <- gge_RDS_ALL_df2 %>% 
                              select(Taxa, Season, PseudoR2) %>% 
                              distinct()

gge_RDS_ALL_df2_distinctPval <- gge_RDS_ALL_df2 %>% 
  select(Taxa, Season, Predictor, Pval, PseudoR2) %>% 
  distinct() %>% 
  mutate(PvalStar = case_when(Pval <= 0.001 ~ "***",
                              Pval > 0.001 & Pval <= 0.01 ~ "**",
                              Pval > 0.01 & Pval <= 0.05 ~ "*",
                              Pval > 0.05 & Pval <= 0.1 ~ ".",
                              Pval > 0.1 ~ ""))

## prepare rug ----


# need to align names with gge_RDS_ALL_df2_distinctPval
PhenPred_ALL_rug <- PhenPred_ALL %>% 
                    pivot_longer(cols = AvgTemp:PPLthen10, names_to = "Predictor", values_to = "Values") %>% 
                    left_join(gge_RDS_ALL_df2_distinctPval, by = c("Taxa", "Season", "Predictor")) %>% 
                    filter(!is.na(Pval))


## Predictor names ----
gge_RDS_ALL_df2Names <- gge_RDS_ALL_df2 %>% 
                          select(Season, Taxa, Predictor) %>% 
                          distinct() %>% 
                                group_by(Taxa, Season) %>% 
                                # tweak taxa name placement
                                mutate(Predictor2 = as.numeric(as.factor(Predictor)),
                                       Predictor2r = rank(Predictor2),
                                       # numbers based on eyeballing plot
                                       DOYcMax2 = (-25 - Predictor2r*7)) %>% 
                                left_join(gge_RDS_ALL_df2_distinctPval, by = c("Taxa", "Season", "Predictor"))

# Colors ---- 
phencolors_Pred <- c("Temp", "Temp x Byth",
                     "PPngd x Lepto", 
                     "PPgd",  "PPgd x Lepto", "PPgd x Byth","PPL10 x Lepto", "PPgd x Plank",
                     "ZPbms x Byth",
                     "Lepto", "Byth", 
                     "Plank")
phencolors <- c("#DD3497","#DD3497",
                "#66FF00", 
                "#66A61E", "#66A61E", "#66A61E", "#66A61E", 
                "#005A32", 
                "#A6761D", 
                "#E6AB02", "#D95F02", 
                "#666666")


# Plot ----

# Get facet labels correct
Taxa.labels <- c("D. retrocurva", "Veliger", "S. oregonensis", "Mesocyclops")
names(Taxa.labels) <- c("Ret", "Vel", "Oreg", "Meso")


levels(gge_RDS_ALL_df2$Season) <- c("DOY[start]",
                                    "DOY[middle]",
                                    "DOY[PB]",
                                    "DOY[end]")
levels(gge_RDS_ALL_df2Names$Season) <- c("DOY[start]",
                                    "DOY[middle]",
                                    "DOY[PB]",
                                    "DOY[end]")
levels(gge_RDS_ALL_df2_distinctR2$Season) <- c("DOY[start]",
                                         "DOY[middle]",
                                         "DOY[PB]",
                                         "DOY[end]")

# Get line labels correct
PhenLegendLabels <- c("Temp", expression(paste("Temp x ", italic("B. longimanus"))),
                      expression(paste(PP[NGD], " x ",italic("L. kindtii"))), 
                      expression(PP[GD]), 
                      expression(paste(PP[GD], " x ",italic("L. kindtii"))), 
                      expression(paste(PP[GD], " x ",italic("B. longimanus"))), 
                      expression(paste(PP[GD], " x Planktivorous fish")),
                      expression(paste(PP[L10], " x ",italic("L. kindtii"))),
                      expression(paste(ZP[BM], " x ",italic("B. longimanus"))),
                      expression(italic("L. kindtii")), 
                      expression(italic("B. longimanus")) ,
                      "Planktivorous fish")
PhenLegendGuideParam <- c(label.hjust = 0)

png("05_Figures/09b_PhenologyPredictions3_BiWk20231202_glmer.png",
    units = "in", height = 9, width = 12, res = 300)
ggplot() +
  geom_line(data = gge_RDS_ALL_df2, 
            aes(y = predictedC, x = x, color = Predictor, linetype = group),
            size = 1.1) +
  # NEED TO REMOVE THINGS NOT IN MODEL FROM RUG
  # geom_rug(data = PhenPred_ALL_rug ,
  #          aes(x = Values, color = Predictor), sides = "b") +
  geom_text(data = gge_RDS_ALL_df2Names,
            aes(x = -7.4, y = DOYcMax2-10, label = Predictor, color = Predictor),
            hjust = 0, fontface = "bold") +
  geom_text(data = gge_RDS_ALL_df2Names,
            aes(x = -7.5, y = DOYcMax2-10, label = PvalStar, color = Predictor),
            hjust = 1, vjust = 0.5, fontface = "bold") +
  geom_text(data = gge_RDS_ALL_df2_distinctR2, aes(y = -40, x = 1, label = round(PseudoR2,2)),
            color = "black", hjust = 0, fontface = "bold") +
  scale_linetype_manual(values = c("solid", "dotted", "dashed"),
                        name = c("")) +
  facet_grid(Taxa ~ Season,
             labeller = labeller(Taxa = Taxa.labels,
                                 Season = Seas.labels)) +
ylab("Centered DOY") +
  xlab("Centered and standardized predictor values") +
  # ylim(-22,22) +
  scale_color_manual(values = phencolors, name = "",
                     labels = PhenLegendLabels,
                     guide = guide_legend(label.hjust = 0)) +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        # legend.position = "none",
        strip.text.y.right = element_text(face = "italic", size = 16),
        strip.text.x.top =  element_text(face = "bold", size = 16),
        strip.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "white", color = "black")) 
dev.off()


png("05_Figures/09b_PhenologyPredictions3_groupedbypredictor_BiWk20231112_glmer.png",
    units = "in", height = 9, width = 18, res = 300)
ggplot() +
  geom_line(data = gge_RDS_ALL_df2, 
            aes(y = predictedC, x = x, color = Taxa, linetype = group),
            size = 1.1) +
  geom_text(data = gge_RDS_ALL_df2_distinctR2, aes(y = -30, x = 1, label = round(PseudoR2,2)),
            color = "black", hjust = 0, fontface = "bold") +
  scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  facet_grid(Season ~ Predictor) +
  ylab("Centered DOY") +
  xlab("Centered and standardized predictor values") +
  # ylim(-22,22) +
  # scale_color_manual(values = phencolors, name = "",
  #                    # labels = PhenLegendLabels,
  #                    guide = guide_legend(label.hjust = 0)) +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        # legend.position = "none",
        strip.text = element_text(face = "bold.italic", size = 16),
        strip.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "white", color = "black")) 
dev.off()

# save/load ----
# save.image(file.path(here::here("04_SavedRImages"),
#                      "09b_PhenologyPreedictionPlots_BiWk_glmer.RData"))
# load(file.path(here::here("04_SavedRImages"),
#                      "09b_PhenologyPreedictionPlots_BiWk_glmer.RData"))



