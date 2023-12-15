# Fits GAMS to start, middle, end of season
# JMH Dec 2023


library(tidyverse)
library(mgcv)
library(mgcViz)
library(MuMIn)
library(ggeffects)

# Data files ----
# phenology metrics
ZPphen <- read.csv("03_GeneratedData/05_ZPPhenologyMetrics.csv", row.names = 1) %>% 
              mutate(DOY = as.numeric(as.character(strftime(as.POSIXct(DOY_moment, format = "%Y-%m-%d"), format= "%j"))),
                     Sample_site = as.factor(Sample_site))


# how do things vary over time
## ret ----
### 15% ----
ZPphen_ret_15 <- ZPphen %>% 
                        filter(Taxa == "Dretrocurva" & Moment == "per15") %>% 
                        mutate(Sample_site2 = ifelse(Sample_site %in% c("14-972", "16-970"), "s1416",
                                                        ifelse(Sample_site %in% c("37-890", "8-994"), "s37_8", "s27_29_3_36")),
                               Sample_site2 = as.factor(Sample_site2))

gm_ret_15_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                      family = gaussian(link = "identity"),
                      method = "REML",
                 data = ZPphen_ret_15)

gratia::draw(gm_ret_15_bySS, residuals = TRUE)

gm_ret_15_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                      family = gaussian(link = "identity"),
                      method = "REML",
                      data = ZPphen_ret_15)

gratia::draw(gm_ret_15_bySS2, residuals = TRUE)

gm_ret_15_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                   family = gaussian(link = "identity"),
                   method = "REML",
                 data = ZPphen_ret_15)

gratia::draw(gm_ret_15_0, residuals = TRUE)

model.sel(gm_ret_15_bySS, gm_ret_15_0, gm_ret_15_bySS2)

gm_ret_15_ML <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                   family = gaussian(link = "identity"),
                   method = "ML",
                   data = ZPphen_ret_15)

summary(gm_ret_15_ML)
check(getViz(gm_ret_15_ML))
print(plot(getViz(gm_ret_15_ML), allTerms = T), pages = 1)
gratia::draw(gm_ret_15_ML, residuals = TRUE)

saveRDS(gm_ret_15_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_ret15.RDS"))

### 50% ----
ZPphen_ret_50 <- ZPphen %>% 
  filter(Taxa == "Dretrocurva" & Moment == "per50") %>% 
  mutate(Sample_site2 = case_when(Sample_site %in% c("3-996","8-994") ~ "s38",
                                  Sample_site %in% c("14-972","16-970") ~ "s1416",
                                  Sample_site %in% c("27-918","29-905") ~ "s2729",
                                  Sample_site %in% c("36-873","37-890") ~ "s3637"),
         Sample_site2 = as.factor(Sample_site2))

gm_ret_50_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                      family = Gamma(link = "identity"),
                      method = "REML",
                      data = ZPphen_ret_50)

gratia::draw(gm_ret_50_bySS, residuals = TRUE)

gm_ret_50_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                      family = Gamma(link = "identity"),
                      method = "REML",
                      data = ZPphen_ret_50)

gratia::draw(gm_ret_50_bySS2, residuals = TRUE)

gm_ret_50_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                   family = Gamma(link = "identity"),
                   method = "REML",
                   data = ZPphen_ret_50)

model.sel(gm_ret_50_bySS, gm_ret_50_0, gm_ret_50_bySS2)

gm_ret_50_ML <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                    family = Gamma(link = "identity"),
                    method = "ML",
                    data = ZPphen_ret_50)

summary(gm_ret_50_ML)
check(getViz(gm_ret_50_ML))
print(plot(getViz(gm_ret_50_ML), allTerms = T), pages = 1)
gratia::draw(gm_ret_50_ML, residuals = TRUE)
gm_ret_50_ML_df <- ggpredict(gm_ret_50_ML, c("Year", "Sample_site"))
plot(gm_ret_50_ML_df) 

saveRDS(gm_ret_50_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_ret50.RDS"))

### 85% ----
ZPphen_ret_85 <- ZPphen %>% 
  filter(Taxa == "Dretrocurva" & Moment == "per85") %>% 
  mutate(Sample_site2 = ifelse(Sample_site %in% c("3-996", "8-994", "16-970", "14-972"), "East", "West"),
         Sample_site2 = as.factor(Sample_site2))

gm_ret_85_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                      family = Gamma(link = "identity"),
                      method = "REML",
                      data = ZPphen_ret_85)

gm_ret_85_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                      family = Gamma(link = "identity"),
                      method = "REML",
                      data = ZPphen_ret_85)
gratia::draw(gm_ret_85_bySS2, residuals = TRUE)

gm_ret_85_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                   family = Gamma(link = "identity"),
                   method = "REML",
                   data = ZPphen_ret_85)

model.sel(gm_ret_85_bySS, gm_ret_85_0, gm_ret_85_bySS2)

gm_ret_85_ML <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                    family = Gamma(link = "identity"),
                    method = "ML",
                    data = ZPphen_ret_85)

gratia::draw(gm_ret_85_ML, residuals = TRUE)

summary(gm_ret_85_ML)
check(getViz(gm_ret_85_ML))
print(plot(getViz(gm_ret_85_ML), allTerms = T), pages = 1)

plot(ggpredict(gm_ret_85_ML, c("Year", "Sample_site2"))) 

saveRDS(gm_ret_85_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_ret85.RDS"))

## meso ----
### 15% ----
ZPphen_Meso_15 <- ZPphen %>% 
  filter(Taxa == "Mesocyclops" & Moment == "per15") %>% 
  mutate(Sample_site2 = case_when(Sample_site %in% c("3-996") ~ "s3",
                                  Sample_site %in% c("14-972","16-970", "27-918", "29-905") ~ "s14162729",
                                  # Sample_site %in% c("29-905") ~ "s29",
                                  Sample_site %in% c("36-873","37-890", "8-994") ~ "s36378"),
         Sample_site2 = as.factor(Sample_site2))

gm_Meso_15_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                      family = gaussian(link = "identity"),
                      method = "REML",
                      data = ZPphen_Meso_15)

gratia::draw(gm_Meso_15_bySS, residuals = TRUE)


gm_Meso_15_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                       family = gaussian(link = "identity"),
                       method = "REML",
                       data = ZPphen_Meso_15)

gratia::draw(gm_Meso_15_bySS2, residuals = TRUE)


gm_Meso_15_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                   family = gaussian(link = "identity"),
                   method = "REML",
                   data = ZPphen_Meso_15)

gratia::draw(gm_Meso_15_0, residuals = TRUE)

model.sel(gm_Meso_15_bySS, gm_Meso_15_0, gm_Meso_15_bySS2)

gm_Meso_15_ML <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                    family = gaussian(link = "identity"),
                    method = "ML",
                    data = ZPphen_Meso_15)

summary(gm_Meso_15_ML)
check(getViz(gm_Meso_15_ML))
print(plot(getViz(gm_Meso_15_ML), allTerms = T), pages = 1)

saveRDS(gm_Meso_15_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Meso15.RDS"))

### 50% ----
ZPphen_Meso_50 <- ZPphen %>% 
  filter(Taxa == "Mesocyclops" & Moment == "per50") %>% 
  mutate(Sample_site2 = case_when(#Sample_site %in% c("3-996") ~ "s3",
                                  Sample_site %in% c("14-972","16-970") ~ "s1416",
                                  Sample_site %in% c("27-918", "29-905", "3-996") ~ "s27293",
                                  Sample_site %in% c("36-873","37-890", "8-994") ~ "s36378"),
         Sample_site2 = as.factor(Sample_site2))

gm_Meso_50_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                       family = Gamma(link = "identity"),
                       method = "REML",
                       data = ZPphen_Meso_50)

gratia::draw(gm_Meso_50_bySS, residuals = TRUE)


gm_Meso_50_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                        family = Gamma(link = "identity"),
                        method = "REML",
                        data = ZPphen_Meso_50)

gratia::draw(gm_Meso_50_bySS2, residuals = TRUE)


gm_Meso_50_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                    family = Gamma(link = "identity"),
                    method = "REML",
                    data = ZPphen_Meso_50)

gratia::draw(gm_Meso_50_0, residuals = TRUE)

model.sel(gm_Meso_50_bySS, gm_Meso_50_0, gm_Meso_50_bySS2)

gm_Meso_50_ML <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                     family = Gamma(link = "identity"),
                     method = "ML",
                     data = ZPphen_Meso_50)

summary(gm_Meso_50_ML)
check(getViz(gm_Meso_50_ML))
print(plot(getViz(gm_Meso_50_ML), allTerms = T), pages = 1)
gratia::draw(gm_Meso_50_ML, residuals = TRUE)
gratia::appraise(gm_Meso_50_ML)

saveRDS(gm_Meso_50_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Meso50.RDS"))

### 85% ----
ZPphen_Meso_85 <- ZPphen %>% 
  filter(Taxa == "Mesocyclops" & Moment == "per85") %>% 
  mutate(Sample_site2 = case_when(#Sample_site %in% c("3-996") ~ "s3",
    Sample_site %in% c("14-972","16-970", "27-918", "29-905") ~ "s14162729",
    Sample_site %in% c("36-873","37-890", "8-994", "3-996") ~ "s363738"),
    Sample_site2 = as.factor(Sample_site2))

gm_Meso_85_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                       family = Gamma(link = "identity"),
                       method = "REML",
                       data = ZPphen_Meso_85)

gratia::draw(gm_Meso_85_bySS, residuals = TRUE)


gm_Meso_85_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                        family = Gamma(link = "identity"),
                        method = "REML",
                        data = ZPphen_Meso_85)

gratia::draw(gm_Meso_85_bySS2, residuals = TRUE)


gm_Meso_85_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                    family = Gamma(link = "identity"),
                    method = "REML",
                    data = ZPphen_Meso_85)

gratia::draw(gm_Meso_85_0, residuals = TRUE)

model.sel(gm_Meso_85_bySS, gm_Meso_85_0, gm_Meso_85_bySS2)

gm_Meso_85_ML <- gam(DOY ~ s(Year, by = Sample_site2, bs = "tp") + s(Sample_site, bs = "re"),
                     family = Gamma(link = "identity"),
                     method = "ML",
                     data = ZPphen_Meso_85)

summary(gm_Meso_85_ML)
check(getViz(gm_Meso_85_ML))
print(plot(getViz(gm_Meso_85_ML), allTerms = T), pages = 1)

saveRDS(gm_Meso_85_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Meso85.RDS"))


## Veligers ----
### 15% ----
ZPphen_Vel_15 <- ZPphen %>% 
  filter(Taxa == "Veliger" & Moment == "per15") %>% 
  mutate(Sample_site2 = case_when(#Sample_site %in% c("3-996") ~ "s3",
    Sample_site %in% c("14-972","29-905") ~ "s1429",
    Sample_site %in% c("27-918","8-994") ~ "s278",
    Sample_site %in% c("36-873","37-890", "3-996", "16-970") ~ "s3637163"),
    Sample_site2 = as.factor(Sample_site2),
    Sample_site1 = case_when(Sample_site %in% c("16-970", "3-996", "8-994") ~ "Unimodal",
                             Sample_site %in% c("36-873","37-890","29-905","27-918", "14-972") ~ "Bimodal"),
    Sample_site1 = as.factor(Sample_site1))

gm_Vel_15_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                       family = Gamma(link = "identity"),
                       method = "REML",
                       data = ZPphen_Vel_15)

summary(gm_Vel_15_bySS)
gratia::draw(gm_Vel_15_bySS, residuals = TRUE)
check(getViz(gm_Vel_15_bySS))

gm_Vel_15_bySS1 <- gam(DOY ~ s(Year, by = Sample_site1) + s(Sample_site, bs = "re"),
                       family = Gamma(link = "identity"),
                       method = "REML",
                       data = ZPphen_Vel_15)

summary(gm_Vel_15_bySS1)
gratia::draw(gm_Vel_15_bySS1, residuals = TRUE)
check(getViz(gm_Vel_15_bySS1))

gm_Vel_15_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                      family = Gamma(link = "identity"),
                      method = "REML",
                      data = ZPphen_Vel_15)

summary(gm_Vel_15_bySS2)
gratia::draw(gm_Vel_15_bySS2, residuals = TRUE)
check(getViz(gm_Vel_15_bySS2))

gm_Vel_15_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                    family = Gamma(link = "identity"),
                    method = "REML",
                    data = ZPphen_Vel_15)

gratia::draw(gm_Vel_15_0, residuals = TRUE)

model.sel(gm_Vel_15_bySS, gm_Vel_15_0, gm_Vel_15_bySS1, gm_Vel_15_bySS2)

gm_Vel_15_ML <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                     family = Gamma(link = "identity"),
                     method = "ML",
                     data = ZPphen_Vel_15)

summary(gm_Vel_15_ML)
check(getViz(gm_Vel_15_ML))
print(plot(getViz(gm_Vel_15_ML), allTerms = T), pages = 1)

saveRDS(gm_Vel_15_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Vel15.RDS"))

### 50% ----
ZPphen_Vel_50 <- ZPphen %>% 
  filter(Taxa == "Veliger" & Moment == "per50") %>% 
  mutate(Sample_site2 = case_when(#Sample_site %in% c("3-996") ~ "s3",
                        Sample_site %in% c("14-972","16-970") ~ "s1416",
                        Sample_site %in% c("27-918","29-905") ~ "s2729",
                        Sample_site %in% c("3-996","8-994") ~ "s38",
                        Sample_site %in% c("36-873","37-890") ~ "s3637"),
    Sample_site2 = as.factor(Sample_site2),
    Sample_site1 = case_when(Sample_site %in% c("16-970", "3-996", "8-994") ~ "Unimodal",
                             Sample_site %in% c("36-873","37-890","29-905","27-918", "14-972") ~ "Bimodal"),
    Sample_site1 = as.factor(Sample_site1))

gm_Vel_50_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                      family = gaussian(link = "identity"),
                      method = "REML",
                      data = ZPphen_Vel_50)

summary(gm_Vel_50_bySS)
gratia::draw(gm_Vel_50_bySS, residuals = TRUE)
check(getViz(gm_Vel_50_bySS))

gm_Vel_50_bySS1 <- gam(DOY ~ s(Year, by = Sample_site1) + s(Sample_site, bs = "re"),
                       family = gaussian(link = "identity"),
                       method = "REML",
                       data = ZPphen_Vel_50)

summary(gm_Vel_50_bySS1)
gratia::draw(gm_Vel_50_bySS1, residuals = TRUE)
check(getViz(gm_Vel_50_bySS1))

gm_Vel_50_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                       family = gaussian(link = "identity"),
                       method = "REML",
                       data = ZPphen_Vel_50)

summary(gm_Vel_50_bySS2)
gratia::draw(gm_Vel_50_bySS2, residuals = TRUE)
check(getViz(gm_Vel_50_bySS2))

gm_Vel_50_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                   family = gaussian(link = "identity"),
                   method = "REML",
                   data = ZPphen_Vel_50)

summary(gm_Vel_50_0)
gratia::draw(gm_Vel_50_0, residuals = TRUE)
check(getViz(gm_Vel_50_0))

model.sel(gm_Vel_50_bySS, gm_Vel_50_0, gm_Vel_50_bySS1, gm_Vel_50_bySS2)

gm_Vel_50_ML <- gam(DOY ~ s(Year, by = Sample_site1) + s(Sample_site, bs = "re"),
                    family = gaussian(link = "identity"),
                    method = "ML",
                    data = ZPphen_Vel_50)

summary(gm_Vel_50_ML)
check(getViz(gm_Vel_50_ML))
gratia::draw(gm_Vel_50_ML, residuals = TRUE)
print(plot(getViz(gm_Vel_50_ML), allTerms = T), pages = 1)

saveRDS(gm_Vel_50_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Vel50.RDS"))


### 85% ----
ZPphen_Vel_85 <- ZPphen %>% 
  filter(Taxa == "Veliger" & Moment == "per85") %>% 
  mutate(Sample_site2 = case_when(#Sample_site %in% c("3-996") ~ "s3",
    Sample_site %in% c("14-972","16-970") ~ "s1416",
    Sample_site %in% c("27-918","29-905") ~ "s2729",
    Sample_site %in% c("3-996","8-994") ~ "s38",
    Sample_site %in% c("36-873","37-890") ~ "s3637"),
    Sample_site2 = as.factor(Sample_site2),
    Sample_site1 = case_when(Sample_site %in% c("16-970", "3-996", "8-994") ~ "Unimodal",
                             Sample_site %in% c("36-873","37-890","29-905","27-918", "14-972") ~ "Bimodal"),
    Sample_site1 = as.factor(Sample_site1))

gm_Vel_85_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                      family = gaussian(link = "identity"),
                      method = "REML",
                      data = ZPphen_Vel_85)

summary(gm_Vel_85_bySS)
gratia::draw(gm_Vel_85_bySS, residuals = TRUE)
check(getViz(gm_Vel_85_bySS))

gm_Vel_85_bySS1 <- gam(DOY ~ s(Year, by = Sample_site1) + s(Sample_site, bs = "re"),
                       family = gaussian(link = "identity"),
                       method = "REML",
                       data = ZPphen_Vel_85)

summary(gm_Vel_85_bySS1)
gratia::draw(gm_Vel_85_bySS1, residuals = TRUE)
check(getViz(gm_Vel_85_bySS1))

gm_Vel_85_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                       family = gaussian(link = "identity"),
                       method = "REML",
                       data = ZPphen_Vel_85)

summary(gm_Vel_85_bySS2)
gratia::draw(gm_Vel_85_bySS2, residuals = TRUE)
check(getViz(gm_Vel_85_bySS2))

gm_Vel_85_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                   family = gaussian(link = "identity"),
                   method = "REML",
                   data = ZPphen_Vel_85)

summary(gm_Vel_85_0)
gratia::draw(gm_Vel_85_0, residuals = TRUE)
check(getViz(gm_Vel_85_0))

model.sel(gm_Vel_85_bySS, gm_Vel_85_0, gm_Vel_85_bySS1, gm_Vel_85_bySS2)

gm_Vel_85_ML <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                    family = gaussian(link = "identity"),
                    method = "ML",
                    data = ZPphen_Vel_85)

summary(gm_Vel_85_ML)
check(getViz(gm_Vel_85_ML))
gratia::draw(gm_Vel_85_ML, residuals = TRUE)
print(plot(getViz(gm_Vel_85_ML), allTerms = T), pages = 1)

saveRDS(gm_Vel_85_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Vel85.RDS"))

## Orego ----
### 15% ----
ZPphen_Orego_15 <- ZPphen %>% 
  filter(Taxa == "Soregonensis" & Moment == "per15") %>% 
  mutate(Sample_site2 = case_when(#Sample_site %in% c("3-996") ~ "s3",
    Sample_site %in% c("14-972","16-970") ~ "s1416",
    Sample_site %in% c("3-996","8-994") ~ "s38"),
    Sample_site2 = as.factor(Sample_site2))

gm_Orego_15_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                      family = gaussian(link = "identity"),
                      method = "REML",
                      data = ZPphen_Orego_15)

summary(gm_Orego_15_bySS)
gratia::draw(gm_Orego_15_bySS, residuals = TRUE)
check(getViz(gm_Orego_15_bySS))


gm_Orego_15_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                       family = gaussian(link = "identity"),
                       method = "REML",
                       data = ZPphen_Orego_15)

summary(gm_Orego_15_bySS2)
gratia::draw(gm_Orego_15_bySS2, residuals = TRUE)
check(getViz(gm_Orego_15_bySS2))

gm_Orego_15_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                   family = gaussian(link = "identity"),
                   method = "REML",
                   data = ZPphen_Orego_15)

summary(gm_Orego_15_0)
gratia::draw(gm_Orego_15_0, residuals = TRUE)
check(getViz(gm_Orego_15_0))

model.sel(gm_Orego_15_bySS, gm_Orego_15_0, gm_Orego_15_bySS2)

gm_Orego_15_ML <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                    family = gaussian(link = "identity"),
                    method = "ML",
                    data = ZPphen_Orego_15)

summary(gm_Orego_15_ML)
check(getViz(gm_Orego_15_ML))
gratia::draw(gm_Orego_15_ML, residuals = TRUE)
print(plot(getViz(gm_Orego_15_ML), allTerms = T), pages = 1)

saveRDS(gm_Orego_15_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Oreg15.RDS"))

### 50% ----
ZPphen_Orego_50 <- ZPphen %>% 
  filter(Taxa == "Soregonensis" & Moment == "per50") %>% 
  mutate(Sample_site2 = case_when(#Sample_site %in% c("3-996") ~ "s3",
    Sample_site %in% c("14-972","16-970") ~ "s1416",
    Sample_site %in% c("3-996","8-994") ~ "s38"),
    Sample_site2 = as.factor(Sample_site2))

gm_Orego_50_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                        family = gaussian(link = "identity"),
                        method = "REML",
                        data = ZPphen_Orego_50)

summary(gm_Orego_50_bySS)
gratia::draw(gm_Orego_50_bySS, residuals = TRUE)
check(getViz(gm_Orego_50_bySS))


gm_Orego_50_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                         family = gaussian(link = "identity"),
                         method = "REML",
                         data = ZPphen_Orego_50)

summary(gm_Orego_50_bySS2)
gratia::draw(gm_Orego_50_bySS2, residuals = TRUE)
check(getViz(gm_Orego_50_bySS2))

gm_Orego_50_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                     family = gaussian(link = "identity"),
                     method = "REML",
                     data = ZPphen_Orego_50)

summary(gm_Orego_50_0)
gratia::draw(gm_Orego_50_0, residuals = TRUE)
check(getViz(gm_Orego_50_0))

model.sel(gm_Orego_50_bySS, gm_Orego_50_0, gm_Orego_50_bySS2)

gm_Orego_50_ML <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                      family = gaussian(link = "identity"),
                      method = "ML",
                      data = ZPphen_Orego_50)

summary(gm_Orego_50_ML)
check(getViz(gm_Orego_50_ML))
gratia::draw(gm_Orego_50_ML, residuals = TRUE)
print(plot(getViz(gm_Orego_50_ML), allTerms = T), pages = 1)

saveRDS(gm_Orego_50_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Oreg50.RDS"))

### 85% ----
ZPphen_Orego_85 <- ZPphen %>% 
  filter(Taxa == "Soregonensis" & Moment == "per85") %>% 
  mutate(Sample_site2 = case_when(#Sample_site %in% c("3-996") ~ "s3",
    Sample_site %in% c("14-972","16-970") ~ "s1416",
    Sample_site %in% c("3-996","8-994") ~ "s38"),
    Sample_site2 = as.factor(Sample_site2))

gm_Orego_85_bySS <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                        family = gaussian(link = "identity"),
                        method = "REML",
                        data = ZPphen_Orego_85)

summary(gm_Orego_85_bySS)
gratia::draw(gm_Orego_85_bySS, residuals = TRUE)
check(getViz(gm_Orego_85_bySS))


gm_Orego_85_bySS2 <- gam(DOY ~ s(Year, by = Sample_site2) + s(Sample_site, bs = "re"),
                         family = gaussian(link = "identity"),
                         method = "REML",
                         data = ZPphen_Orego_85)

summary(gm_Orego_85_bySS2)
gratia::draw(gm_Orego_85_bySS2, residuals = TRUE)
check(getViz(gm_Orego_85_bySS2))

gm_Orego_85_0 <- gam(DOY ~ s(Year) + s(Sample_site, bs = "re"),
                     family = gaussian(link = "identity"),
                     method = "REML",
                     data = ZPphen_Orego_85)

summary(gm_Orego_85_0)
gratia::draw(gm_Orego_85_0, residuals = TRUE)
check(getViz(gm_Orego_85_0))

model.sel(gm_Orego_85_bySS, gm_Orego_85_0, gm_Orego_85_bySS2)

gm_Orego_85_ML <- gam(DOY ~ s(Year, by = Sample_site) + s(Sample_site, bs = "re"),
                      family = gaussian(link = "identity"),
                      method = "ML",
                      data = ZPphen_Orego_85)

summary(gm_Orego_85_ML)
check(getViz(gm_Orego_85_ML))
gratia::draw(gm_Orego_85_ML, residuals = TRUE)
print(plot(getViz(gm_Orego_85_ML), allTerms = T), pages = 1)

saveRDS(gm_Orego_85_ML, here::here("06_Rdat", "10_GAMS_seasMoments_ML_Oreg85.RDS"))


# save/load ----
# save.image(here::here("04_SavedRImages", "10_GAMS_seasonalMoments.Rdata"))
# load(here::here("04_SavedRImages", "10_GAMS_seasonalMoments.Rdata"))
