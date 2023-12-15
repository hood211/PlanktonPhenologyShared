
# Predict Meso phenol
# JMH Nov 23

# libaries ----
library(tidyverse)
library(ggcorrplot)
library(vegan)
library(lme4)
library(MuMIn)
library(sjPlot)
library(ggeffects)

# data ----
## phenology ----
Phen <- read.csv(file.path(here::here("03_GeneratedData"),"05_ZPPhenologyMetrics.csv"),
                 row.names = 1) %>% 
  filter(Taxa == "Mesocyclops") %>% 
  mutate(DOY_moment = as.POSIXct(DOY_moment, format = "%Y-%m-%d"),
         M = as.character(strftime(DOY_moment, format = "%b")),
         DOY = as.numeric(strftime(DOY_moment, format = "%j")),
         DOYwinStart = DOY - 30) %>% 
  select(-DOY_moment)

## predictors ----
### basin wide avg ----
Pred_BW <- read.csv(file.path(here::here("03_GeneratedData"),"07b_Predictors_LogedAnoms_BasinAvg_BiWk.csv"),
                    row.names = 1) %>% 
  pivot_longer(cols = bw9:bw20, names_to = "BiWk", values_to = "BiWkAnom") %>% 
  pivot_wider(id_cols = c(Year, BiWk), names_from = Response, values_from = BiWkAnom) %>% 
  mutate(BiWk2 = as.numeric(str_sub(BiWk, 3, 4))*2,
         # this is the date of the biweek average
         BiWkDate = paste0(Year,"-", BiWk2, "-1  "),
         BiWkDate2 = as.POSIXct(BiWkDate, format = "%Y-%U-%u"),
         DOY = strftime(BiWkDate2, format = "%j"))  %>% 
  # drop predictors I don't want
  select(-c(MeanIce, Chla, PPLthen10))


### site group -----
Pred_SG <- read.csv(file.path(here::here("03_GeneratedData"),"07b_Predictors_LogedAnoms_SiteGroup_BiWk.csv"),
                    row.names = 1) %>% 
  pivot_longer(cols = bw9:bw20, names_to = "BiWk", values_to = "BiWkAnom") %>% 
  pivot_wider(id_cols = c(Sample_site, Year, BiWk), names_from = Response, values_from = BiWkAnom) %>% 
  mutate(BiWk2 = as.numeric(str_sub(BiWk, 3, 4))*2,
         # this is the date of the biweek average
         BiWkDate = paste0(Year,"-", BiWk2, "-1  "),
         BiWkDate2 = as.POSIXct(BiWkDate, format = "%Y-%U-%u"),
         DOY = strftime(BiWkDate2, format = "%j"))   %>% 
  # drop predictors I don't want
  select(-c(MeanIce, Chla, PPLthen10))


# Join Pred & Phen ----
## Basin wide ----
Phen_BW <- Phen %>% 
  filter(Moment == "PB1") %>% 
  group_by(Taxa, Moment, Year) %>%
  summarise(across(DOY:DOYwinStart, \(x) mean(x, na.rm = T)))

PhenPred_BW2 <- as.data.frame(matrix(nrow = 1, ncol = 12))
names(PhenPred_BW2) <- c("Taxa","Moment","Year","DOY","DOYwinStart","AvgTemp",
                         "Byth","EdiblePP","InEdiblePP","Lepto", "MesoZPfood","plank_per_ha")

# Get avg of predictors
for(i in 1:dim(Phen_BW)[1]){
  # i = 1
  Phen_BW_i <- Phen_BW[i,]
  
  Pred_BW_i <- Pred_BW %>% 
    filter(Year == Phen_BW_i$Year & 
             (DOY >= Phen_BW_i$DOYwinStart & DOY <= Phen_BW_i$DOY)) %>% 
    ungroup() %>% 
    summarise(across(AvgTemp:plank_per_ha, \(x) mean(x, na.rm = T)))
  
  PhenPred_BW2_i <- cbind(Phen_BW_i, Pred_BW_i)
  
  PhenPred_BW2[i,] <- PhenPred_BW2_i
}

## sample site ----
Phen_SG <- Phen %>% 
  filter(Moment != "PB1") 


PhenPred_SG2 <- as.data.frame(matrix(nrow = 1, ncol = 14))
names(PhenPred_SG2) <- c("Taxa", "Year", "Sample_site", "Moment","M", "DOY", "DOYwinStart", "AvgTemp", 
                         "Byth", "EdiblePP","InEdiblePP","Lepto", "MesoZPfood", "plank_per_ha")

# Get avg of predictors
for(i in 1:dim(Phen_SG)[1]){
  # i = 1
  Phen_SG_i <-Phen_SG[i,]
  
  Pred_SG_i <- Pred_SG %>% 
    filter(Year == Phen_SG_i$Year & 
             Sample_site == Phen_SG_i$Sample_site & 
             (DOY >= Phen_SG_i$DOYwinStart & DOY <= Phen_SG_i$DOY)) %>% 
    ungroup() %>% 
    summarise(across(AvgTemp:plank_per_ha, \(x) mean(x, na.rm = T)))
  
  PhenPred_SG_i <- cbind(Phen_SG_i, Pred_SG_i)
  
  PhenPred_SG2[i,] <- PhenPred_SG_i
}


# TS plot forPhen and metrics
ggplot(PhenPred_BW2 %>% 
         pivot_longer(cols = AvgTemp:plank_per_ha, names_to = "Pred",
                      values_to = "BiWkAnom"), 
       aes(y = BiWkAnom, x = Year)) +
  geom_point() +
  facet_grid(Pred ~ .)

ggplot(PhenPred_SG2 %>% 
         pivot_longer(cols = AvgTemp:plank_per_ha, names_to = "Pred",
                      values_to = "BiWkAnom"), 
       aes(y = BiWkAnom, x = Year, color = Sample_site)) +
  geom_point() +
  facet_grid(Pred ~ Moment)

# cor plot for each taxa x phen metric

pdf("05_Figures/08cb_PredictorsCorr4PhenMetrics_Meso_BW_BiWk.pdf")
for(i in 1:(length(unique(PhenPred_BW2$Moment)))){
  Response_i <- unique(PhenPred_BW2$Moment)[i]
  
  cor_i <- cor(PhenPred_BW2[PhenPred_BW2$Moment == Response_i,-c(1:3,5)], 
               use = "pairwise.complete.obs", method = "spearman")
  pval_i <- cor_pmat(PhenPred_BW2[PhenPred_BW2$Moment == Response_i,-c(1:3,5)], 
                     use = "pairwise.complete.obs")
  corrPlot_i <-  ggcorrplot(cor_i, p.mat = pval_i, 
                            insig = "blank", outline.col = "white", type = "lower",
                            lab = T, show.diag = T) + #
    ggtitle(Response_i)
  print(corrPlot_i)
}
dev.off()

pdf("05_Figures/08cb_PredictorsCorr4PhenMetrics_Meso_SiteGroup_BiWk.pdf")
for(i in 1:(length(unique(PhenPred_SG2$Moment)))){
  Response_i <- unique(PhenPred_SG2$Moment)[i]
  
  cor_i <- cor(PhenPred_SG2[PhenPred_SG2$Moment == Response_i,-c(1:5,6)], 
               use = "pairwise.complete.obs", method = "spearman")
  pval_i <- cor_pmat(PhenPred_SG2[PhenPred_SG2$Moment == Response_i,-c(1:5,6)], 
                     use = "pairwise.complete.obs")
  corrPlot_i <-  ggcorrplot(cor_i, p.mat = pval_i, 
                            insig = "blank", outline.col = "white", type = "lower",
                            lab = T, show.diag = T) + #
    ggtitle(Response_i)
  print(corrPlot_i)
}
dev.off()



# GAM of phen metrics ----
options(na.action = "na.fail")

## 15%  ----
PhenPred_15per <- PhenPred_SG2 %>% 
                  filter(Moment == "per15") %>% 
                  # loosing a lot of data bc fish
                  # fish  < 2, but NS
                  select(-plank_per_ha) %>%
                  filter(complete.cases(.)) %>% 
                  mutate(Sample_site = as.factor(Sample_site))


### model sel ----
lmm.15per.global <- glmer(DOY ~ AvgTemp + 
                            Byth * MesoZPfood + 
                            Byth * InEdiblePP + 
                            Lepto * MesoZPfood + 
                            Lepto * InEdiblePP +
                            # plank_per_ha * MesoZPfood +
                            # plank_per_ha * InEdiblePP +
                            (1|Sample_site),
                          family = Gamma(link = "identity"),
                          glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                          data = PhenPred_15per)

qqnorm(residuals(lmm.15per.global)); qqline(residuals(lmm.15per.global))
# lots of singular fits for the more complex models, but not an issues for top models
lmm.15per.global.Dredge <- dredge(lmm.15per.global, rank = "AICc")
subset(lmm.15per.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.15per.global.Dredge), 
          file.path(here::here("07_Tables"),"08cb_ModelSelTable_MesoPhen_15Per_BiWK_glmer.csv"))

### most likely ----
lmm.15per.MostLik <- glmer(DOY ~ AvgTemp + Byth*InEdiblePP +  (1|Sample_site),
                           family = Gamma(link = "identity"),
                           glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                           data = PhenPred_15per)

qqnorm(residuals(lmm.15per.MostLik)); qqline(residuals(lmm.15per.MostLik))
summary(lmm.15per.MostLik)
performance::r2(lmm.15per.MostLik)
cor(model.response(model.frame(lmm.15per.MostLik)),predict(lmm.15per.MostLik,type="response"))^2

plot(ggpredict(lmm.15per.MostLik, "AvgTemp"), residuals = TRUE)

plot(ggpredict(lmm.15per.MostLik, c("InEdiblePP", "Byth [-0.8,0.75]")), residuals = TRUE)


saveRDS(lmm.15per.MostLik, here::here("06_Rdat", "08cb_GAMPhenMeso_15per_BiWk_glmer.RDS"))
write.csv(PhenPred_15per, here::here("03_GeneratedData", "08cb_GAMPhenMeso_15per_BiWk_glmer.csv"))



## 50%  ----

PhenPred_50per <- PhenPred_SG2 %>% 
  filter(Moment == "per50") %>% 
  # loosing a lot of data bc fish
  # fish in <2, but NS
  select(-plank_per_ha) %>%
  filter(complete.cases(.)) %>% 
  mutate(Sample_site = as.factor(Sample_site))


### model sel ----
lmm.50per.global <- glmer(DOY ~ AvgTemp + 
                            Byth * MesoZPfood + 
                            Byth * InEdiblePP + 
                            Lepto * MesoZPfood + 
                            Lepto * InEdiblePP +
                            # plank_per_ha * MesoZPfood +
                            # plank_per_ha * InEdiblePP +
                            (1|Sample_site),
                          family = Gamma(link = "identity"),
                          glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                          data = PhenPred_50per)

qqnorm(residuals(lmm.50per.global)); qqline(residuals(lmm.50per.global))
# lots of singular fits for the more complex models, but not an issues for top models
lmm.50per.global.Dredge <- dredge(lmm.50per.global, rank = "AICc")
subset(lmm.50per.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.50per.global.Dredge), 
          file.path(here::here("07_Tables"),"08cb_ModelSelTable_MesoPhen_50per_BiWK_glmer.csv"))

### most likely ----
lmm.50per.MostLik <- glmer(DOY ~ Byth*InEdiblePP + (1|Sample_site),
                           family = Gamma(link = "identity"),
                           glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                           data = PhenPred_50per)

qqnorm(residuals(lmm.50per.MostLik)); qqline(residuals(lmm.50per.MostLik))
summary(lmm.50per.MostLik)
performance::r2(lmm.50per.MostLik)
cor(model.response(model.frame(lmm.50per.MostLik)),predict(lmm.50per.MostLik,type="response"))^2

plot(ggpredict(lmm.50per.MostLik, c("InEdiblePP", "Byth [-0.93,0.81]")), residuals = TRUE)


saveRDS(lmm.50per.MostLik, here::here("06_Rdat", "08cb_GAMPhenMeso_50per_BiWk_glmer.RDS"))
write.csv(PhenPred_50per, here::here("03_GeneratedData", "08cb_GAMPhenMeso_50per_BiWk_glmer.csv"))



## 85%  ----

PhenPred_85per <- PhenPred_SG2 %>% 
  filter(Moment == "per85") %>% 
  # loosing a lot of data bc fish
  # select(-plank_per_ha) %>%
  filter(complete.cases(.)) %>% 
  mutate(Sample_site = as.factor(Sample_site))


### model sel ----
lmm.85per.global <- glmer(DOY ~ AvgTemp + 
                            Byth * MesoZPfood + 
                            Byth * InEdiblePP + 
                            Lepto * MesoZPfood + 
                            Lepto * InEdiblePP +
                            plank_per_ha * MesoZPfood +
                            plank_per_ha * InEdiblePP +
                            (1|Sample_site),
                          family = Gamma(link = "identity"),
                          glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                          data = PhenPred_85per)

qqnorm(residuals(lmm.85per.global)); qqline(residuals(lmm.85per.global))
# lots of singular fits for the more complex models, but not an issues for top models
lmm.85per.global.Dredge <- dredge(lmm.85per.global, rank = "AICc")
subset(lmm.85per.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.85per.global.Dredge), 
          file.path(here::here("07_Tables"),"08cb_ModelSelTable_MesoPhen_85per_BiWK_glmer.csv"))

### most likely ----
lmm.85per.MostLik <- glmer(DOY ~ InEdiblePP + plank_per_ha + Byth*MesoZPfood + (1|Sample_site),
                           family = Gamma(link = "identity"),
                           glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                           data = PhenPred_85per)

qqnorm(residuals(lmm.85per.MostLik)); qqline(residuals(lmm.85per.MostLik))
# there are two very low MesoZP values. Removing doesn't influence Sig of Interaction term
summary(lmm.85per.MostLik)
performance::r2(lmm.85per.MostLik)
cor(model.response(model.frame(lmm.85per.MostLik)),predict(lmm.85per.MostLik,type="response"))^2

plot(ggpredict(lmm.85per.MostLik, "InEdiblePP"), residuals = TRUE)
plot(ggpredict(lmm.85per.MostLik, "plank_per_ha"), residuals = TRUE)
plot(ggpredict(lmm.85per.MostLik, c("MesoZPfood", "Byth [-0.8,0.86]")), residuals = TRUE)

saveRDS(lmm.85per.MostLik, here::here("06_Rdat", "08cb_GAMPhenMeso_85per_BiWk_glmer.RDS"))
write.csv(PhenPred_85per, here::here("03_GeneratedData", "08cb_GAMPhenMeso_85per_BiWk_glmer.csv"))


## Peak Biomass  ----

# No variation in PB


# save/load
# save.image(file.path(here::here("04_SavedRImages"),"08cb_PredictPhenol_Meso_BiWk_glmer.RData"))
# load(file.path(here::here("04_SavedRImages"),"08cb_PredictPhenol_Meso_BiWk_glmer.RData"))
