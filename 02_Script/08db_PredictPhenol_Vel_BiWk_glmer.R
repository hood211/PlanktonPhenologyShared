# Predict Vel phenol
# JMH Nov 23, dec 23

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
  filter(Taxa == "Veliger") %>% 
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
  select(-c(MeanIce, Chla, MesoZPfood)) 

### site group -----
Pred_SG <- read.csv(file.path(here::here("03_GeneratedData"),"07b_Predictors_LogedAnoms_SiteGroup_BiWk.csv"),
                    row.names = 1) %>% 
  pivot_longer(cols = bw9:bw20, names_to = "BiWk", values_to = "BiWkAnom") %>% 
  pivot_wider(id_cols = c(Sample_site, Year, BiWk), names_from = Response, values_from = BiWkAnom) %>% 
  mutate(BiWk2 = as.numeric(str_sub(BiWk, 3, 4))*2,
         # this is the date of the biweek average
         BiWkDate = paste0(Year,"-", BiWk2, "-1  "),
         BiWkDate2 = as.POSIXct(BiWkDate, format = "%Y-%U-%u"),
         DOY = strftime(BiWkDate2, format = "%j"),
         Sample_site = as.factor(Sample_site))  %>% 
  # drop predictors I don't want
  select(-c(MeanIce, Chla, MesoZPfood)) 


# Join Pred & Phen ----

## Sample site ----
Phen_SG <- Phen %>% 
  filter(!(Moment %in% c("PB1", "PB2"))) 

PhenPred_SG2 <- as.data.frame(matrix(nrow = 1, ncol = 14))
names(PhenPred_SG2) <- c("Taxa", "Year", "Sample_site", "Moment","M", "DOY", "DOYwinStart", "AvgTemp", 
                         "Byth", "EdiblePP","InEdiblePP","Lepto", "PPLthen10", "plank_per_ha")

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

pdf("05_Figures/08db_PredictorsCorr4PhenMetrics_Vel_SiteGroup_BiWk.pdf")
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
                  # loosing a lot b/c of fish (167 v. 83)
                  # fish not in top models
                  # select(-plank_per_ha) %>%
                  filter(complete.cases(.)) 
                

### model sel ----
lmm.15per.global <- glmer(DOY ~ AvgTemp + 
                            Byth * PPLthen10 + 
                            Byth * InEdiblePP + 
                            Lepto * PPLthen10 + 
                            Lepto * InEdiblePP +
                            plank_per_ha * PPLthen10 +
                            plank_per_ha * InEdiblePP +
                            (1|Sample_site),
                          family = Gamma(link = "identity"),
                          glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                          data = PhenPred_15per)

qqnorm(residuals(lmm.15per.global)); qqline(residuals(lmm.15per.global))
# lots of singular fits for the more complex models, but not an issues for top models
lmm.15per.global.Dredge <- dredge(lmm.15per.global, rank = "AICc", m.max = 6)
subset(lmm.15per.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.15per.global.Dredge), 
          file.path(here::here("07_Tables"),"08db_ModelSelTable_VelPhen_15per_BiWK_glmer.csv"))

### most likely ----
lmm.15per.MostLik <- glmer(DOY ~ plank_per_ha + InEdiblePP + (1|Sample_site),
                           family = Gamma(link = "identity"),
                           # glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                           data = PhenPred_15per)


qqnorm(residuals(lmm.15per.MostLik)); qqline(residuals(lmm.15per.MostLik))
summary(lmm.15per.MostLik)
performance::r2(lmm.15per.MostLik)
cor(model.response(model.frame(lmm.15per.MostLik)),predict(lmm.15per.MostLik,type="response"))^2

plot_model(lmm.15per.MostLik, type = "re")
plot(ggpredict(lmm.15per.MostLik, "plank_per_ha"), residuals = TRUE)
plot(ggpredict(lmm.15per.MostLik, "InEdiblePP"), residuals = TRUE)
# https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html


saveRDS(lmm.15per.MostLik, here::here("06_Rdat", "08db_GAMPhenVel_15per_BiWk_glmer.RDS"))
write.csv(PhenPred_15per, here::here("03_GeneratedData", "08db_GAMPhenVel_15per_BiWk_glmer.csv"))


## 50%  ----

PhenPred_50per <- PhenPred_SG2 %>% 
  filter(Moment == "per50") %>%
  mutate(Sample_site = as.factor(Sample_site)) %>% 
  # loosing a lot b/c of fish (187 v 115)
  # fish <2, but NS
  select(-plank_per_ha) %>%
  filter(complete.cases(.)) 


### model sel ----
lmm.50per.global <- glmer(DOY ~ AvgTemp + 
                            Byth * PPLthen10 + 
                            Byth * InEdiblePP + 
                            Lepto * PPLthen10 + 
                            Lepto * InEdiblePP +
                            # plank_per_ha * PPLthen10 +
                            # plank_per_ha * InEdiblePP +
                            (1|Sample_site),
                          family = Gamma(link = "identity"),
                          glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                          data = PhenPred_50per)

qqnorm(residuals(lmm.50per.global)); qqline(residuals(lmm.50per.global))
# lots of singular fits for the more complex models, but not an issues for top models
lmm.50per.global.Dredge <- dredge(lmm.50per.global, rank = "AICc", m.max = 6)
subset(lmm.50per.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.50per.global.Dredge), 
          file.path(here::here("07_Tables"),"08db_ModelSelTable_VelPhen_50per_BiWK_glmer.csv"))

### most likely ----
lmm.50per.MostLik <- glmer(DOY ~ Byth * InEdiblePP + Lepto * PPLthen10  + (1|Sample_site),
                           family = Gamma(link = "identity"),
                           glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                           data = PhenPred_50per)


qqnorm(residuals(lmm.50per.MostLik)); qqline(residuals(lmm.50per.MostLik))
summary(lmm.50per.MostLik)
performance::r2(lmm.50per.MostLik)
cor(model.response(model.frame(lmm.50per.MostLik)),predict(lmm.50per.MostLik,type="response"))^2

plot_model(lmm.50per.MostLik, type = "re")
plot(ggpredict(lmm.50per.MostLik, c("InEdiblePP", "Byth [-0.88,0.84]")), residuals = TRUE)
plot(ggpredict(lmm.50per.MostLik, c("PPLthen10", "Lepto [-0.68,0.83]")), residuals = TRUE)

# https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html


saveRDS(lmm.50per.MostLik, here::here("06_Rdat", "08db_GAMPhenVel_50per_BiWk_glmer.RDS"))
write.csv(PhenPred_50per, here::here("03_GeneratedData", "08db_GAMPhenVel_50per_BiWk_glmer.csv"))

## 85%  ----

PhenPred_85per <- PhenPred_SG2 %>% 
  filter(Moment == "per85") %>%
  mutate(Sample_site = as.factor(Sample_site)) %>% 
  # loosing a lot b/c of fish (172 v 126)
  # fish <2, but NS
  select(-plank_per_ha) %>%
  filter(complete.cases(.)) 

### model sel ----
lmm.85per.global <- glmer(DOY ~ AvgTemp + 
                            Byth * PPLthen10 + 
                            Byth * InEdiblePP + 
                            Lepto * PPLthen10 + 
                            Lepto * InEdiblePP +
                            # plank_per_ha * PPLthen10 +
                            # plank_per_ha * InEdiblePP +
                            (1|Sample_site),
                          family = Gamma(link = "identity"),
                          glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                          data = PhenPred_85per)

qqnorm(residuals(lmm.85per.global)); qqline(residuals(lmm.85per.global))
# lots of singular fits for the more complex models, but not an issues for top models
lmm.85per.global.Dredge <- dredge(lmm.85per.global, rank = "AICc", m.max = 6)
subset(lmm.85per.global.Dredge, delta <2)


write.csv(as.data.frame(lmm.85per.global.Dredge), 
          file.path(here::here("07_Tables"),"08db_ModelSelTable_VelPhen_85per_BiWK_glmer.csv"))

### most likely ----
lmm.85per.MostLik <- glmer(DOY ~ InEdiblePP + Byth  + (1|Sample_site),
                           family = Gamma(link = "identity"),
                           glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),
                           data = PhenPred_85per)


qqnorm(residuals(lmm.85per.MostLik)); qqline(residuals(lmm.85per.MostLik))
summary(lmm.85per.MostLik)
performance::r2(lmm.85per.MostLik)
cor(model.response(model.frame(lmm.85per.MostLik)),predict(lmm.85per.MostLik,type="response"))^2

plot_model(lmm.85per.MostLik, type = "re")
plot(ggpredict(lmm.85per.MostLik, "Byth"), residuals = TRUE)
plot(ggpredict(lmm.85per.MostLik, "InEdiblePP"), residuals = TRUE)
# https://sscc.wisc.edu/sscc/pubs/MM/MM_DiagInfer.html
# https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html


saveRDS(lmm.85per.MostLik, here::here("06_Rdat", "08db_GAMPhenVel_85per_BiWk_glmer.RDS"))
write.csv(PhenPred_85per, here::here("03_GeneratedData", "08db_GAMPhenVel_85per_BiWk_glmer.csv"))


## Peak Biomass  ----

# No variation in PB


# save/load
# save.image(file.path(here::here("04_SavedRImages"),"08db_PredictPhenol_Vel_BiWk_glmer.RData"))
# load(file.path(here::here("04_SavedRImages"),"08db_PredictPhenol_Vel_BiWk_glmer.RData"))
