# Predict Oreg phenol
# JMH Nov 23, Dec 23

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
  filter(Taxa == "Soregonensis") %>% 
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
         BiWkDate = paste0(Year,"-", BiWk2, "-1"),
         BiWkDate2 = as.POSIXct(BiWkDate, format = "%Y-%U-%u"),
         DOY = strftime(BiWkDate2, format = "%j")) %>% 
  select(-c("Chla","MeanIce","MesoZPfood", "PPLthen10"))

### site group -----
# note: no longer site group, retaining object name for simplicity
Pred_SG <- read.csv(file.path(here::here("03_GeneratedData"),"07b_Predictors_LogedAnoms_Site_BiWk.csv"),
                    row.names = 1) %>% 
  pivot_longer(cols = bw9:bw20, names_to = "BiWk", values_to = "BiWkAnom") %>% 
  pivot_wider(id_cols = c(Sample_site, Year, BiWk), names_from = Response, values_from = BiWkAnom) %>% 
  mutate(BiWk2 = as.numeric(str_sub(BiWk, 3, 4))*2,
         # this is the date of the biweek average
         BiWkDate = paste0(Year,"-", BiWk2, "-1"),
         BiWkDate2 = as.POSIXct(BiWkDate, format = "%Y-%U-%u"),
         DOY = strftime(BiWkDate2, format = "%j"))  %>% 
  select(-c("Chla","MeanIce","MesoZPfood", "PPLthen10")) %>% 
  filter(Sample_site %in% c("14-972", "16-970", "3-996",  "8-994"))


# Join Pred & Phen ----
## Basin wide ----
Phen_BW <- Phen %>% 
  filter(Moment == "PB1") %>% 
  group_by(Taxa, Moment, Year) %>%
  summarise(across(DOY:DOYwinStart, \(x) mean(x, na.rm = T)))

PhenPred_BW2 <- as.data.frame(matrix(nrow = 1, ncol = 11))
names(PhenPred_BW2) <- c("Taxa","Moment","Year","DOY","DOYwinStart",
                         "AvgTemp","Byth","EdiblePP","InEdiblePP","Lepto","plank_per_ha")

# Get avg of predictors
for(i in 1:dim(Phen_BW)[1]){
  # i = 1
  Phen_BW_i <-Phen_BW[i,]
  
  Pred_BW_i <- Pred_BW %>% 
    filter(Year == Phen_BW_i$Year & 
             (DOY >= Phen_BW_i$DOYwinStart & DOY <= Phen_BW_i$DOY)) %>% 
    ungroup() %>% 
    summarise(across(AvgTemp:plank_per_ha, \(x) mean(x, na.rm = T)))
  
  PhenPred_BW2_i <- cbind(Phen_BW_i, Pred_BW_i)
  
  PhenPred_BW2[i,] <- PhenPred_BW2_i
}


# TS plot forPhen and metrics
ggplot(PhenPred_BW2 %>% 
         pivot_longer(cols = AvgTemp:plank_per_ha, names_to = "Pred",
                      values_to = "BiWkAnom"), 
       aes(y = BiWkAnom, x = Year)) +
  geom_point() +
  facet_grid(Pred ~ Moment)

# cor plot for each taxa x phen metric

pdf("05_Figures/08bb_PredictorsCorr4PhenMetrics_Oreg_BW_BiWk.pdf")
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


## Sample site ----
# using site, not site group, but retaining SG for simplicity
Phen_SG <- Phen %>% 
  filter(Moment != "PB1") 


PhenPred_SG2 <- as.data.frame(matrix(nrow = 1, ncol = 13))
names(PhenPred_SG2) <- c("Taxa", "Year", "Sample_site", "Moment","M", "DOY", "DOYwinStart", "AvgTemp", 
                         "Byth", "EdiblePP","InEdiblePP","Lepto", "plank_per_ha")


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
ggplot(PhenPred_SG2 %>% 
         pivot_longer(cols = AvgTemp:plank_per_ha, names_to = "Pred",
                      values_to = "BiWkAnom"), 
       aes(y = BiWkAnom, x = Year, color = Sample_site)) +
  geom_point() +
  facet_grid(Pred ~ Moment)

# cor plot for each taxa x phen metric

pdf("05_Figures/08bb_PredictorsCorr4PhenMetrics_Oreg_SiteGroup_BiWk.pdf")
for(i in 1:(length(unique(PhenPred_SG2$Moment)))){
  Response_i <- unique(PhenPred_SG2$Moment)[i]
  
  cor_i <- cor(PhenPred_SG2[PhenPred_SG2$Moment == Response_i,-c(1:5,7)], 
               use = "pairwise.complete.obs", method = "spearman")
  pval_i <- cor_pmat(PhenPred_SG2[PhenPred_SG2$Moment == Response_i,-c(1:5,7)], 
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
                  mutate(Sample_site = as.factor(Sample_site)) %>% 
                  # loosing a lot b/c of fish
                  # fish not in models <2
                  select(-plank_per_ha) %>%
                  filter(complete.cases(.)) 

### model sel ----
lmm.15per.global <- glmer(DOY ~ AvgTemp + 
                            Byth * EdiblePP + 
                            Byth * InEdiblePP + 
                            Lepto * EdiblePP + 
                            Lepto * InEdiblePP +
                            # plank_per_ha * EdiblePP +
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
          file.path(here::here("07_Tables"),"08bb_ModelSelTable_OregPhen_15Per_BiWK_glmer.csv"))

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
plot(ggpredict(lmm.15per.MostLik, c("InEdiblePP", "Byth [-0.96,0.79]")), residuals = TRUE)
plot_model(lmm.15per.MostLik, type = "re")
# https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

plot(ggpredict(lmm.15per.MostLik, terms = c("AvgTemp", "Sample_site")), add.data = TRUE, facets = TRUE)
plot(ggpredict(lmm.15per.MostLik, terms = c("InEdiblePP", "Byth [-0.96,0.79]", "Sample_site")), add.data = TRUE, facets = TRUE)

saveRDS(lmm.15per.MostLik, here::here("06_Rdat", "08bb_GAMPhenOreg_15per_BiWk_glmer.RDS"))
write.csv(PhenPred_15per, here::here("03_GeneratedData", "08bb_GAMPhenOreg_15per_BiWk_glmer.csv"))



## 50%  ----

PhenPred_50per <- PhenPred_SG2 %>% 
  filter(Moment == "per50")  %>%
  mutate(Sample_site = as.factor(Sample_site)) %>% 
  # fish not <2
  select(-plank_per_ha) %>%
  filter(complete.cases(.)) 


### model sel ----
# sample site isn't doing anything, removing RE
lmm.50per.global <- glm(DOY ~ AvgTemp + 
                            Byth * EdiblePP + 
                            Byth * InEdiblePP + 
                            Lepto * EdiblePP + 
                            Lepto * InEdiblePP, # +
                            # plank_per_ha * EdiblePP +
                            # plank_per_ha * InEdiblePP +
                            # (1|Sample_site),
                          family = Gamma(link = "identity"),
                          # glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                          data = PhenPred_50per)

qqnorm(residuals(lmm.50per.global)); qqline(residuals(lmm.50per.global))
lmm.50per.global.Dredge <- dredge(lmm.50per.global, rank = "AICc")
subset(lmm.50per.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.50per.global.Dredge), 
          file.path(here::here("07_Tables"),"08bb_ModelSelTable_OregPhen_50Per_BiWK_glmer.csv"))

### most likely ----
lmm.50per.MostLik <- glm(DOY ~ Byth * InEdiblePP + Lepto, # + (1|Sample_site),
                           family = Gamma(link = "identity"),
                           # glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                           data = PhenPred_50per)

qqnorm(residuals(lmm.50per.MostLik)); qqline(residuals(lmm.50per.MostLik))
summary(lmm.50per.MostLik)
performance::r2(lmm.50per.MostLik)
cor(model.response(model.frame(lmm.50per.MostLik)),predict(lmm.50per.MostLik,type="response"))^2

plot(ggpredict(lmm.50per.MostLik, "Lepto"), residuals = TRUE)
plot(ggpredict(lmm.50per.MostLik, c("InEdiblePP", "Byth [-0.98,0.84]")), residuals = TRUE)


saveRDS(lmm.50per.MostLik, here::here("06_Rdat", "08bb_GAMPhenOreg_50per_BiWk_glmer.RDS"))
write.csv(PhenPred_50per, here::here("03_GeneratedData", "08bb_GAMPhenOreg_50per_BiWk_glmer.csv"))



## 85%  ----

PhenPred_85per <- PhenPred_SG2 %>% 
  filter(Moment == "per85")  %>%
  mutate(Sample_site = as.factor(Sample_site)) %>% 
  # fish not <2
  select(-plank_per_ha) %>%
  # removes 2 outliers having big effect on model
  filter(DOY > 210) %>% 
  filter(complete.cases(.)) 

### model sel ----
# added temp X Byth after examining data
lmm.85per.global <- glm(DOY ~ AvgTemp * Byth + 
                            Byth * EdiblePP + 
                            Byth * InEdiblePP + 
                            Lepto * EdiblePP + 
                            Lepto * InEdiblePP,# +
                            # plank_per_ha * EdiblePP +
                            # plank_per_ha * InEdiblePP,# +
                            # (1|Sample_site),
                          family = Gamma(link = "identity"),
                          # glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                          data = PhenPred_85per)

qqnorm(residuals(lmm.85per.global)); qqline(residuals(lmm.85per.global))
lmm.85per.global.Dredge <- dredge(lmm.85per.global, rank = "AICc")
subset(lmm.85per.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.85per.global.Dredge), 
          file.path(here::here("07_Tables"),"08bb_ModelSelTable_OregPhen_85Per_BiWK_glmer.csv"))

### most likely ----
lmm.85per.MostLik <- glm(DOY ~ Byth * AvgTemp,# + (1|Sample_site),
                           family = Gamma(link = "identity"),
                           # glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                           data = PhenPred_85per)

qqnorm(residuals(lmm.85per.MostLik)); qqline(residuals(lmm.85per.MostLik))
summary(lmm.85per.MostLik)
performance::r2(lmm.85per.MostLik)
cor(model.response(model.frame(lmm.85per.MostLik)),predict(lmm.85per.MostLik,type="response"))^2

 plot(ggpredict(lmm.85per.MostLik, c("AvgTemp", "Byth [-0.96,0.8]")), residuals = TRUE)


saveRDS(lmm.85per.MostLik, here::here("06_Rdat", "08bb_GAMPhenOreg_85per_BiWk_glmer.RDS"))
write.csv(PhenPred_85per, here::here("03_GeneratedData", "08bb_GAMPhenOreg_85per_BiWk_glmer.csv"))




## Peak Biomass  ----

PhenPred_PB <- PhenPred_BW2 %>%
  filter(Moment == "PB1")%>%
  # loosing a lot b/c of fish
  # select(-plank_per_ha) %>% # not < 2
  filter(complete.cases(.)) %>% 
  filter(Year > 2002)
  # removing 1 big outlier, changes imp of edible PP
  # filter(EdiblePP >-2)

plot(PhenPred_BW2$DOY ~ PhenPred_BW2$Year)

### model sel ----
lmm.PB.global <- glm(DOY ~ AvgTemp * Byth + Lepto + EdiblePP + InEdiblePP, # + plank_per_ha,
                     family = Gamma(link = "identity"), # no family/link improves
                     data = PhenPred_PB)

qqnorm(residuals(lmm.PB.global)); qqline(residuals(lmm.PB.global))
lmm.PB.global.Dredge <- dredge(lmm.PB.global, rank = "AICc")
subset(lmm.PB.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.PB.global.Dredge), 
          file.path(here::here("07_Tables"),"08bb_ModelSelTable_OregPhen_PB_BiWK_glmer.csv"))

### most likely ----
lmm.PB.MostLik <- glm(DOY ~ Byth + InEdiblePP,
                      family = Gamma(link = "identity"),
                      data = PhenPred_PB)

qqnorm(residuals(lmm.PB.MostLik)); qqline(residuals(lmm.PB.MostLik))
summary(lmm.PB.MostLik)
performance::r2(lmm.PB.MostLik)
cor(model.response(model.frame(lmm.PB.MostLik)),predict(lmm.PB.MostLik,type="response"))^2

plot(ggpredict(lmm.PB.MostLik, "Byth"), residuals = TRUE)
plot(ggpredict(lmm.PB.MostLik, "InEdiblePP"), residuals = TRUE)


saveRDS(lmm.PB.MostLik, here::here("06_Rdat", "08bb_GAMPhenOreg_PB_BiWk_glmer.RDS"))
write.csv(PhenPred_PB, here::here("03_GeneratedData", "08bb_GAMPhenOreg_PB_BiWk_glmer.csv"))

# save/load ----
# save.image(file.path(here::here("04_SavedRImages"),"08bb_PredictPhenol_Oreg_BiWk_glmer.RData"))
# load(file.path(here::here("04_SavedRImages"),"08bb_PredictPhenol_Oreg_BiWk_glmer.RData"))
