# - Cleans up phenology metrics data and prepares for modeling
# - Generates figures 2, 3, S5, S6, S7
# - Models relationship between Biomass and start, middle, and end of season
# JMH, Dec 23

# Libraries ----
library(tidyverse)
library(tidybayes)
library(cowplot)
library(grid)
library(gridExtra)
library(MuMIn)

# Data ----
## Percentile mass data ----

PerMass <- read.csv("03_GeneratedData/04_ZPCumSumPercentiles_LOESSEst.csv", row.names = 1) %>% 
        mutate(across(c(Taxa, Site, CumSumPer = Moment), as.factor)) %>% 
        mutate(DOY.cs = round(DOY.cs,0),
               DOY_CumSumPer_date = as.Date(paste("2010-",DOY.cs), format = "%Y-%j"),
               Taxa = fct_relevel(Taxa, "Veliger","Mesocyclops","Dretrocurva", "Soregonensis")) %>% 
        rename(Year = "Y", Sample_site = Site) %>% 
        filter(Year >= 1995 & Year <= 2022)


 

## Seasonality ----
# see 03's GAMM predict Peak BM
Seas_ret <- read.csv(file.path(here::here("03_GeneratedData/"),"03a_GAMMpredictions_Ret.csv" ), row.names = 1) %>% 
            mutate(Taxa = "Dretrocurva") %>% 
            select("Taxa","Sample_site", "Year", "DOY", "ppmedianBM" = "pRetBM",
                   ".lower", ".upper", ".width", ".point", ".interval") %>% 
  mutate(Region = "ALL")
Seas_oreg <- read.csv(file.path(here::here("03_GeneratedData/"),"03b_GAMMpredictions_Oreg.csv" ), row.names = 1)%>% 
            mutate(Taxa = "Soregonensis") %>% 
            select("Taxa","Sample_site", "Year", "DOY", "ppmedianBM" = "pOregBM",
                   ".lower", ".upper", ".width", ".point", ".interval") %>% 
  mutate(Region = "ALL")
Seas_meso <- read.csv(file.path(here::here("03_GeneratedData/"),"03c_GAMMpredictions_Meso.csv" ), row.names = 1)%>% 
          mutate(Taxa = "Mesocyclops") %>% 
          select("Taxa","Sample_site", "Year", "DOY", "ppmedianBM" = "pMesoBM",
                 ".lower", ".upper", ".width", ".point", ".interval") %>% 
          mutate(Region = "ALL")
Seas_vel <- read.csv(file.path(here::here("03_GeneratedData/"),"03c_GAMMpredictions_Vel.csv" ), row.names = 1)%>% 
          mutate(Taxa = "Veliger") %>% 
          select("Taxa","Sample_site", "Year", "DOY", "ppmedianBM" = "pVelBM",
                 ".lower", ".upper", ".width", ".point", ".interval",  "Region") 

Seas_pred <- rbind(Seas_ret, Seas_oreg, Seas_meso, Seas_vel) %>% 
        mutate(across(Taxa:Sample_site, as.factor)) %>% 
        filter(DOY >= 125 & DOY <= 275) %>% 
        mutate(Date = as.Date(paste("2010-",DOY), format = "%Y-%j"),
               Taxa = fct_relevel(Taxa, "Dretrocurva", "Veliger", "Soregonensis","Mesocyclops")) %>% 
        filter(Year >= 1995 & Year <= 2022) 


## Peak biomass ----
PB_ret <- read.csv(file.path(here::here("03_GeneratedData/"),"03a_GAMMpredictPeakBM_Ret.csv" ), row.names = 1) %>% 
        mutate(Taxa = "Dretrocurva") %>% 
        select("Taxa","Sample_site", "Year", "DOYpeakBM" = "DOY", "ppmedianBM" = "pRetBM",
               ".lower", ".upper", ".width", ".point", ".interval") %>% 
  mutate(Moment = "Pk1")

PB_oreg <- read.csv(file.path(here::here("03_GeneratedData/"),"03b_GAMMpredictPeakBM_Oreg.csv" ), row.names = 1) %>% 
        mutate(Taxa = "Soregonensis") %>% 
        select("Taxa","Sample_site", "Year", "DOYpeakBM" = "DOY", "ppmedianBM" = "pOregBM",
               ".lower", ".upper", ".width", ".point", ".interval") %>% 
  mutate(Moment = "Pk1")

PB_meso <- read.csv(file.path(here::here("03_GeneratedData/"),"03b_GAMMpredictPeakBM_Meso.csv" ), row.names = 1) %>% 
        mutate(Taxa = "Mesocyclops") %>% 
        select("Taxa","Sample_site", "Year", "DOYpeakBM" = "DOY", "ppmedianBM" = "pMesoBM",
               ".lower", ".upper", ".width", ".point", ".interval") %>% 
        mutate(Moment = "Pk1")

PB_vel <- read.csv(file.path(here::here("03_GeneratedData/"),"03b_GAMMpredictPeakBM_Vel.csv" ), row.names = 1) %>% 
        mutate(Taxa = "Veliger") %>% 
        select("Taxa","Sample_site", "Year", "DOYpeakBM" = "DOY", "ppmedianBM" = "pVelBM",
               ".lower", ".upper", ".width", ".point", ".interval", "Moment") 


PB_pred <- rbind(PB_ret, PB_oreg, PB_meso, PB_vel) %>% 
        mutate(across(Taxa:Sample_site, as.factor)) %>% 
        mutate(DOYpeakBM_date = as.Date(paste("2010-",DOYpeakBM), format = "%Y-%j"),
               Taxa = fct_relevel(Taxa, "Dretrocurva", "Veliger","Soregonensis","Mesocyclops")) %>% 
        filter(Year >= 1995 & Year <= 2022)

## combine CumSum & PB $ plot ----

ZP_moments <- rbind(PerMass %>%
                      select(Taxa, Year, Sample_site, "Moment" = "CumSumPer", "DOY_moment" = "DOY_CumSumPer_date"),
                    PB_pred %>%
                      mutate(Moment = ifelse(Moment == "Pk1", "PB1", "PB2")) %>% 
                      select(Taxa, Year, Sample_site, Moment, "DOY_moment" = "DOYpeakBM_date")) %>% 
  mutate(Moment = fct_relevel(Moment, "PB1", "PB2", "per15", "per50", "per85"))


ZP_moments2 <- ZP_moments %>% 
                  mutate(Sample_site = as.factor(Sample_site),
                         Sample_site = fct_recode(Sample_site, 
                                                  "14" = "14-972","16" = "16-970",
                                                  "27" = "27-918", "29" = "29-905",
                                                  "3" = "3-996", "36" = "36-873",
                                                  "37" = "37-890", "8" = "8-994"),
                         Sample_site = fct_relevel(Sample_site, "36", "37", "27", "29", "14", "16", "3", "8"),
                         Taxa = as.factor(Taxa),
                         Taxa = fct_relevel(Taxa, "Dretrocurva", "Veliger","Soregonensis","Mesocyclops"),
                         Taxa = fct_recode(Taxa,"Retrocurva" = "Dretrocurva", "Oregonensis" = "Soregonensis"))


# Fig 2 - Ret, Oreg, Meso PB ----
TaxaLables <- c("D. retrocurva", "Veligers",  "S. oregonensis","Mesocyclops")
names(TaxaLables) <- c("Dretrocurva","Veliger", "Soregonensis","Mesocyclops")


png(file.path(here::here("05_Figures"), "05_Fig2_SeasonalityPhenology.png"), 
    units = "in", height = 10, width = 8, res = 300)
ggplot() +
    geom_raster(data = Seas_pred %>% 
                    filter(Sample_site == "14-972" & Taxa != "Veliger"), 
                aes(y = Year, x = Date, fill = ppmedianBM)) +
    geom_contour(data = Seas_pred %>% 
                     filter(Sample_site == "14-972" & Taxa != "Veliger"),  
                 aes(y = Year, x = Date, z = ppmedianBM),
                 color = "white", binwidth = 50) +
    geom_point(data = ZP_moments %>% 
                   filter(Sample_site == "14-972" & Moment == "PB1" & Taxa != "Veliger"), #
               aes(y = Year, x = DOY_moment, shape = Moment, color = Moment),
               size = 2) +
    xlab("Month") +
    scale_fill_distiller(palette = "Spectral", name = expression(paste("Biomass (µg ",L^-1,")"))) +
    scale_shape_manual(values = c(3, 15, 16, 17),
                       labels = c("Peak biomass", "15% CB", "50% CB", "85% CB"),
                       name = "Phenology") +
    scale_color_manual(values = c("black", "purple", "purple", "purple"),
                       labels = c("Peak biomass", "15% CB", "50% CB", "85% CB"),
                       name = "Phenology") +
    guides(colour = guide_legend(override.aes = list(size=4.5))) +
    facet_wrap(~Taxa, nrow = 3, ncol = 1,
               labeller = labeller(Taxa = TaxaLables)) +
    theme_bw() +
    theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.position = "right",
          strip.text = element_text(face = "bold.italic", size = 16),
          strip.background = element_rect(color = "transparent"),
          panel.background = element_rect(fill = "white", color = "black"))
dev.off()    

# Fig. S7 ----
png(file.path(here::here("05_Figures"), "05_FigS7_SeasonalityPhenology_AllSites.png"), 
    units = "in", height = 10, width = 16, res = 300)
ggplot() +
   geom_line(data = ZP_moments2, 
             aes(x = Year, y = DOY_moment, color = Moment), 
             size = 0.5, alpha = 0.35) +
  geom_point(data = ZP_moments2, 
             aes(x = Year, y = DOY_moment, shape = Moment, color = Moment),
             size = 2) +
  xlab("Year") +
  ylab("Month") +
  scale_shape_manual(values = c(3, 6, 15, 16, 17),
                     labels = c(expression(DOY[PB1]), expression(DOY[PB2]),
                                expression(DOY[start]), expression(DOY[middle]), expression(DOY[end])),
                     name = "Phenology") +
  scale_color_manual(values = c("#2b83ba", "#d7191c", "#fc8d59", "gold", "#99d594"),
                     labels = c(expression(DOY[PB1]), expression(DOY[PB2]),
                                expression(DOY[start]), expression(DOY[middle]), expression(DOY[end])),
                     name = "Phenology") +
  facet_grid(Taxa ~ Sample_site) +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 16),
        panel.grid.major = element_line(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "right",
        strip.text = element_text(face = "bold.italic", size = 16),
        strip.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "white", color = "black"))
dev.off()

# Writing about moments ----
## Summary ----
PerMass %>% 
    mutate(TaxaMoment = paste0(Taxa,"_", Moment)) %>% 
    split(.$TaxaMoment) %>% 
    map(summary)

PerMass2 <- PerMass %>% 
  select(Taxa, Year, Sample_site, Moment, DOY.cs) %>% 
  pivot_wider(id_cols = c(Taxa, Year, Sample_site), names_from = Moment, values_from = DOY.cs) %>% 
  mutate(SeasDur = per85 - per15) %>% 
  pivot_longer(cols = per15:SeasDur, names_to = "Metric", values_to = "Value")

## Duration ----
PerMassw <- PerMass %>% 
  select(Taxa, Year, Sample_site, Moment, DOY.cs) %>% 
  pivot_wider(id_cols = c(Taxa, Year, Sample_site), names_from = Moment, values_from = DOY.cs) %>% 
  mutate(SeasDur = per85 - per15)

ggplot(PerMassw, aes(y = SeasDur, x = Year, color = Sample_site)) +
  geom_point() +
  facet_grid(Taxa ~ Sample_site) +
  stat_smooth()

ggplot(PerMassw, aes(y = per85, x = per15, color = Sample_site)) +
  geom_point() +
  facet_grid(Taxa ~ Sample_site) +
  geom_abline(intercept = 0, slope = 1) +
  geom_abline(intercept = 50, slope = 1)

PerMassS <- PerMass2 %>% 
  group_by(Taxa, Metric) %>% 
  summarise(min = min(Value, na.rm = T),
            per10 = quantile(Value, 0.10, na.rm = T),
            median = median(Value, na.rm = T),
            mean = mean(Value, na.rm = T),
            per90 = quantile(Value, 0.90, na.rm = T),
            max = max(Value, na.rm = T)) %>% 
  arrange(Taxa, Metric) %>% 
  select(Taxa, Metric, mean, median, per10, per90, min, max)

# Phenology v BM reg ----
# CumSum = cummulative sum at moment
# Biomass = biomass at moment
# MaxYearCumSum = is maximum cum some for year
ggplot(PerMass, aes(y = Biomass.bm, x = CumSum.bm, color = Sample_site)) +
  geom_point() +
  facet_grid(Taxa ~ Moment) +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_log10() +
  scale_y_log10()


# biomass
ggplot(PerMass %>% 
         filter(Biomass.bm > 0), aes(y = DOY.cs, x = Biomass.bm, color = Sample_site)) + #
  geom_point() +
  facet_grid(Taxa ~ Moment, scales = "free_x") +
  scale_x_log10() +
  stat_smooth(se = FALSE, method = "lm")


# CumSum
ggplot(PerMass %>% 
         filter(CumSum.bm > 0), aes(y = DOY.cs, x = CumSum.bm)) + #, color = Sample_site
  geom_point() +
  facet_grid(Taxa ~ Moment, scales = "free_x") +
  scale_x_log10() +
  stat_smooth(se = FALSE, method = "lm")

# Max CumSum biomass
ggplot(PerMass %>% 
         filter(MaxYearCumSum.bm > 0), aes(y = DOY.cs, x = MaxYearCumSum.bm, color = Sample_site)) + #, color = Sample_site
  geom_point() +
  facet_grid(Taxa ~ Moment, scales = "free_x") +
  scale_x_log10() +
  stat_smooth(se = FALSE, method = "lm")

# checking to see that data is the same
ggplot(PerMass %>% 
         filter(MaxYearCumSum.bm > 0), aes(x = MaxYearCumSum.bm)) + #, color = Sample_site
  geom_histogram() +
  facet_grid(Taxa ~ Moment, scales = "free_x") +
  scale_x_log10() 



## scale for modeling ---- 

PerMass_scales <- PerMass %>% 
                    mutate(Taxa_moment = paste0(Taxa,"_", Moment)) %>% 
                    select(Taxa, Sample_site, Year, Moment, Taxa_moment, DOY.cs, Biomass.bm) %>% 
                    mutate(NonZeroMin = case_when(Taxa == "Veliger" ~ 0.0331762,
                                           Taxa == "Mesocyclops" ~ 0.01899621,
                                           Taxa == "Dretrocurva" ~ 0.02804686,
                                           Taxa == "Soregonensis" ~ 0.09223177)) %>% 
                    mutate(Biomass.bm = log10(Biomass.bm + NonZeroMin)) %>% 
                    group_by(Taxa, Sample_site, Moment, Taxa_moment) %>% 
                    mutate(across(DOY.cs:Biomass.bm, \(x) scale(x, center = T, scale = T))) %>% 
                    # mutate(Yearf = case_when(Year >= 1995 & Year < 2000 ~ "1995_1999",
                    #                          Year >= 2000 & Year < 2005 ~ "2000_2004",
                    #                          Year >= 2005 & Year < 2010 ~ "2005_2009",
                    #                          Year >= 2010 & Year < 2015 ~ "2010_2014",
                    #                          Year >= 2015  ~ "2015_2022")) %>% 
                    mutate(Yearf = case_when(Year >= 1995 & Year < 2003 ~ "1995_2003",
                                             Year >= 2003 & Year < 2015 ~ "2000_2014",
                                             Year >= 2015  ~ "2015_2022"),
                           Yearf = fct_relevel(as.factor(Yearf), "1995_2003", "2000_2014", "2015_2022"))
                    




PerMass_scales2 <- cbind(PerMass_scales[,c(1:5,9)], PerMass_scales$DOY.cs[,1], PerMass_scales$Biomass.bm[,1])
names(PerMass_scales2)[7] <- "DOY.csS"
names(PerMass_scales2)[8] <- "Biomass.bmS"

ggplot(PerMass_scales2, aes(y = Biomass.bmS, x = DOY.csS, color = Yearf)) + #, color = Sample_site
  geom_point() +
  facet_grid(Taxa ~ Moment, scales = "free_x") +
  # scale_x_log10() +
  stat_smooth(se = FALSE, method = "lm")


options(na.action = "na.fail")


Taxa_momentList_i <- unique(PerMass_scales2$Taxa_moment)

CumSumBioLMfun <- function(Taxa_momentList_i){
  # Taxa_momentList_i <- "Dretrocurva_per15"
  PerMass_scales2_i <- PerMass_scales2 %>% 
    filter(Taxa_moment == Taxa_momentList_i)
  
  lm_Taxa_momentList_iG <- lm(Biomass.bmS ~ DOY.csS*Sample_site + DOY.csS*Yearf, data = PerMass_scales2_i)
  lm_Taxa_momentList_D <- dredge(lm_Taxa_momentList_iG, rank = "AICc",
                                 extra = list("R2" = function(x) {s <-  summary(x)
                                 c(Rsq = s$r.squared, adjRsq = s$adj.r.squared)}))
  model.sel(lm_Taxa_momentList_D)
}


## biomass - DOY reg ----
# vel
  # 15 - biomass only, -
  Veliger_per15ms <- CumSumBioLMfun("Veliger_per15")
  model.sel(Veliger_per15ms)
  Veliger_per15ms_lm <- get.models(Veliger_per15ms, delta == 0)[[1]]
  summary(Veliger_per15ms_lm)
  
  # 50 - biomass only, -
  Veliger_per50ms <- CumSumBioLMfun("Veliger_per50")
  model.sel(Veliger_per50ms)
  Veliger_per50ms_lm <- get.models(Veliger_per50ms, delta == 0)[[1]]
  summary(Veliger_per50ms_lm)
  
  # 85 - biomass only, -
  Veliger_per85ms <- CumSumBioLMfun("Veliger_per85")
  model.sel(Veliger_per85ms)
  Veliger_per85ms_lm <- get.models(Veliger_per85ms, delta == 0)[[1]]
  summary(Veliger_per85ms_lm)

  # meso
  # 15 - biomass only, -
  Mesocyclops_per15ms <- CumSumBioLMfun("Mesocyclops_per15")
  model.sel(Mesocyclops_per15ms)
  Mesocyclops_per15ms_lm <- get.models(Mesocyclops_per15ms, delta == 0)[[1]]
  summary(Mesocyclops_per15ms_lm)
  
  # 50 - biomass only, -
  Mesocyclops_per50ms <- CumSumBioLMfun("Mesocyclops_per50")
  model.sel(Mesocyclops_per50ms)
  Mesocyclops_per50ms_lm <- get.models(Mesocyclops_per50ms, delta == 0)[[1]]
  summary(Mesocyclops_per50ms_lm)
  
  # 85 - biomass only, -
  Mesocyclops_per85ms <- CumSumBioLMfun("Mesocyclops_per85")
  model.sel(Mesocyclops_per85ms)
  Mesocyclops_per85ms_lm <- get.models(Mesocyclops_per85ms, delta == 0)[[1]]
  summary(Mesocyclops_per85ms_lm)
  
# ret
  # 15 - biomass only, +
  Dretrocurva_per15ms <- CumSumBioLMfun("Dretrocurva_per15")
  model.sel(Dretrocurva_per15ms)
  Dretrocurva_per15ms_lm <- get.models(Dretrocurva_per15ms, delta == 0)[[1]]
  summary(Dretrocurva_per15ms_lm)
  
  # 50 - biomass only, -
  Dretrocurva_per50ms <- CumSumBioLMfun("Dretrocurva_per50")
  model.sel(Dretrocurva_per50ms)
  Dretrocurva_per50ms_lm <- get.models(Dretrocurva_per50ms, delta == 0)[[1]]
  summary(Dretrocurva_per50ms_lm)
  
  # 85 - biomass only, 0
  Dretrocurva_per85ms <- CumSumBioLMfun("Dretrocurva_per85")
  model.sel(Dretrocurva_per85ms)
  Dretrocurva_per85ms_lm <- get.models(Dretrocurva_per85ms, delta == 0)[[1]]
  summary(Dretrocurva_per85ms_lm)
  
# oreg
  # 15 - biomass only, -
  Soregonensis_per15ms <- CumSumBioLMfun("Soregonensis_per15")
  model.sel(Soregonensis_per15ms)
  Soregonensis_per15ms_lm <- get.models(Soregonensis_per15ms, delta == 0)[[1]]
  summary(Soregonensis_per15ms_lm)
  
  # 50 - INTERCEPT ONLY
  Soregonensis_per50ms <- CumSumBioLMfun("Soregonensis_per50")
  model.sel(Soregonensis_per50ms)
  Soregonensis_per50ms_lm <- get.models(Soregonensis_per50ms, delta == 0)[[1]]
  summary(Soregonensis_per50ms_lm)
  
  # 85 - biomass only, +
  Soregonensis_per85ms <- CumSumBioLMfun("Soregonensis_per85")
  model.sel(Soregonensis_per85ms)
  Soregonensis_per85ms_lm <- get.models(Soregonensis_per85ms, delta == 0)[[1]]
  summary(Soregonensis_per85ms_lm)
  
  ## Model predictions ----
  PerMass_scales2$Veliger_per15 <- predict(Veliger_per15ms_lm, PerMass_scales2)
  PerMass_scales2$Veliger_per50 <- predict(Veliger_per50ms_lm, PerMass_scales2)
  PerMass_scales2$Veliger_per85 <- predict(Veliger_per85ms_lm, PerMass_scales2)
  
  PerMass_scales2$Mesocyclops_per15 <- predict(Mesocyclops_per15ms_lm, PerMass_scales2)
  PerMass_scales2$Mesocyclops_per50 <- predict(Mesocyclops_per50ms_lm, PerMass_scales2)
  PerMass_scales2$Mesocyclops_per85 <- predict(Mesocyclops_per85ms_lm, PerMass_scales2)
  
  PerMass_scales2$Dretrocurva_per15 <- predict(Dretrocurva_per15ms_lm, PerMass_scales2)
  PerMass_scales2$Dretrocurva_per50 <- predict(Dretrocurva_per50ms_lm, PerMass_scales2)
  PerMass_scales2$Dretrocurva_per85 <- predict(Dretrocurva_per85ms_lm, PerMass_scales2)
  
  PerMass_scales2$Soregonensis_per15 <- predict(Soregonensis_per15ms_lm, PerMass_scales2)
  PerMass_scales2$Soregonensis_per50 <- predict(Soregonensis_per50ms_lm, PerMass_scales2)
  PerMass_scales2$Soregonensis_per85 <- predict(Soregonensis_per85ms_lm, PerMass_scales2)
            
  PerMass_scales3 <- PerMass_scales2 %>% 
                      mutate(Biomass.bmS_pred = case_when(Taxa_moment == "Veliger_per15" ~ Veliger_per15,
                                                      Taxa_moment == "Veliger_per50" ~ Veliger_per50,
                                                      Taxa_moment == "Veliger_per85" ~ Veliger_per85,
                                                      Taxa_moment == "Mesocyclops_per15" ~ Mesocyclops_per15,
                                                      Taxa_moment == "Mesocyclops_per50" ~ Mesocyclops_per50,
                                                      Taxa_moment == "Mesocyclops_per85" ~ Mesocyclops_per85,
                                                      Taxa_moment == "Dretrocurva_per15" ~ Dretrocurva_per15,
                                                      Taxa_moment == "Dretrocurva_per50" ~ Dretrocurva_per50,
                                                      Taxa_moment == "Dretrocurva_per85" ~ Dretrocurva_per85,
                                                      Taxa_moment == "Soregonensis_per15" ~ Soregonensis_per15,
                                                      Taxa_moment == "Soregonensis_per50" ~ Soregonensis_per50,
                                                      Taxa_moment == "Soregonensis_per85" ~ Soregonensis_per85),
                             # model
                             Model = case_when(Taxa_moment == "Veliger_per15" ~ "BM ~ DOY + Y",
                                                          Taxa_moment == "Veliger_per50" ~ "BM ~ DOY * Y",
                                                          Taxa_moment == "Veliger_per85" ~ "BM ~ DOY + Y",
                                                          Taxa_moment == "Mesocyclops_per15" ~ "BM ~ DOY + Y",
                                                          Taxa_moment == "Mesocyclops_per50" ~ "BM ~ DOY + Y",
                                                          Taxa_moment == "Mesocyclops_per85" ~ "BM ~ DOY * Y",
                                                          Taxa_moment == "Dretrocurva_per15" ~ "BM ~ DOY + Y",
                                                          Taxa_moment == "Dretrocurva_per50" ~ "BM ~ DOY",
                                                          Taxa_moment == "Dretrocurva_per85" ~ "BM ~ DOY",
                                                          Taxa_moment == "Soregonensis_per15" ~ "BM ~ DOY + Y",
                                                          Taxa_moment == "Soregonensis_per50" ~ "BM ~ DOY + Y",
                                                          Taxa_moment == "Soregonensis_per85" ~ "BM ~ DOY + Y"),
                             # R2
                             R2adj = case_when(Taxa_moment == "Veliger_per15" ~ "0.177",
                                                          Taxa_moment == "Veliger_per50" ~ "0.332",
                                                          Taxa_moment == "Veliger_per85" ~ "0.305",
                                                          Taxa_moment == "Mesocyclops_per15" ~ "0.049",
                                                          Taxa_moment == "Mesocyclops_per50" ~ "0.176",
                                                          Taxa_moment == "Mesocyclops_per85" ~ "0.105",
                                                          Taxa_moment == "Dretrocurva_per15" ~ "0.102",
                                                          Taxa_moment == "Dretrocurva_per50" ~ "0.031",
                                                          Taxa_moment == "Dretrocurva_per85" ~ "0.229",
                                                          Taxa_moment == "Soregonensis_per15" ~ "0.182",
                                                          Taxa_moment == "Soregonensis_per50" ~ "0.284",
                                                          Taxa_moment == "Soregonensis_per85" ~ "0.244")) %>%
                      ungroup() %>% 
                      mutate(Yearf2 = ifelse(Taxa_moment %in% c("Dretrocurva_per50", "Dretrocurva_per85"),
                                             "1995_2022", as.character(Yearf)),
                             Yearf2 = fct_relevel(as.factor(Yearf2), 
                                                  "1995_2022", "1995_2003", "2000_2014", "2015_2022"),
                             Taxa = fct_relevel(as.factor(Taxa),
                                                "Dretrocurva", "Veliger", "Soregonensis", "Mesocyclops"),
                             Moment = as.factor(Moment)) %>% 
                      select(Taxa:Biomass.bmS, Yearf2, Model, R2adj, Biomass.bmS_pred)

  ## Fig. S6 ----
  Taxa.labels <- c("D. retrocurva", "Veligers", "S. oregonensis", "Mesocyclops")
  names(Taxa.labels) <- c("Dretrocurva", "Veliger", "Soregonensis", "Mesocyclops")
  

  levels(PerMass_scales3$Moment) <- c("DOY[start]",
                                      "DOY[middle]",
                                      "DOY[end]")

  
  png(file.path(here::here("05_Figures"), "05_FigS6_BiomassVdoy.png"), 
      units = "in", height = 9, width = 10, res = 300)
  ggplot() +
    geom_point(data = PerMass_scales3, aes(y = Biomass.bmS, x = DOY.csS, color = Yearf2), alpha = 50/100) +
    geom_line(data = PerMass_scales3, aes(y = Biomass.bmS_pred, x = DOY.csS, color = Yearf2), size = 1.25) +
    geom_text(data = PerMass_scales3, aes( x = -4, y = -3.9, label = Model), hjust = 0) +
    geom_text(data = PerMass_scales3, aes( x = 3, y = -3.9, label = R2adj), hjust = 1) +
    annotate("text", x = 1.6, y = -3.75, hjust = 1, label = "R^2 ", parse = TRUE) +
    annotate("text", x = 1.8, y = -3.85, hjust = 1, label = "=") +
    ylab(expression(paste("Zooplankton biomass (µg ", L^-1,")"))) +
    xlab("Centered and scaled day of year") +
    scale_color_manual(values = c("black","#FC8D59", "#1B9E77", "#E7298A"), name = "Year group") +
    facet_grid(Taxa ~ Moment, 
               labeller = labeller(Taxa = Taxa.labels, Moment = label_parsed)) +
    theme_bw()+
    theme(axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.position = "top",
          strip.text = element_text(face = "bold.italic", size = 16),
          strip.background = element_rect(color = "transparent"),
          panel.background = element_rect(fill = "white", color = "black")) 
dev.off()  

  
  

# Fig S5 ----
# spatial patterns
dFigS5 <- Seas_pred %>% 
        filter(DOY == 200 & Year == 2010) %>%
        # Exclude bimodal estimates for unimodal sites and vice versa
        filter(!(Taxa == "Veliger" & Region == "Bimodal" & Sample_site %in% c("36-873","37-890","29-905","27-918", "14-972"))) %>% 
        filter(!(Taxa == "Veliger" & Region == "Unimodal" & Sample_site %in% c("16-970", "3-996", "8-994"))) %>% 
        mutate(Sample_site = fct_relevel(Sample_site, "36-873","37-890","29-905","27-918", "14-972", "16-970", "3-996", "8-994"),
               Sample_site = fct_recode(Sample_site, 
                                         "36" = "36-873", "37" = "37-890",
                                          "29" = "29-905","27" = "27-918", 
                                         "14" = "14-972","16" = "16-970",
                                         "3" ="3-996", "8" = "8-994"),
               Taxa = as.factor(Taxa),
               Taxa = fct_relevel(Taxa, "Dretrocurva","Veliger","Soregonensis","Mesocyclops"))

TaxaLables2 <- c("Dreissenid veligers", "Mesocyclops", "D. retrocuva", "S. oregonensis")
names(TaxaLables2) <- c("Veliger", "Mesocyclops","Dretrocurva", "Soregonensis")

png(file.path(here::here("05_Figures"), "05_FigS5_Spatial patterns.png"), 
    units = "in", height = 7, width = 10, res = 300)
ggplot(dFigS5, aes(y = ppmedianBM, x = Sample_site)) +
        geom_errorbar(aes(ymin = .lower, ymax = .upper)) +
        geom_point(shape = 21, fill = "grey", size = 5) +
        facet_wrap(~Taxa, scales = "free_y",
                   labeller = labeller(Taxa = TaxaLables2)) +
        ylab(expression(paste("Biomass (µg ",L^-1,")"))) +
        xlab("Site") +
        theme_bw() +
        theme(axis.title = element_text(size = 22),
              axis.text.x = element_text(size = 16),
              axis.text.y = element_text(size = 16),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.text = element_text(size = 14),
              legend.title = element_text(size = 16),
              legend.position = "right",
              strip.text = element_text(face = "bold.italic", size = 16),
              strip.background = element_rect(color = "transparent"),
              panel.background = element_rect(fill = "white", color = "black"))
dev.off()


## Fig 3 - veligers
Fig3a <- ggplot(Seas_pred %>% 
           filter(Taxa == "Veliger" & Year == 2010 & Sample_site %in% c("36-873", "16-970")),
       aes(y = ppmedianBM, x = Date, color = Region)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = Region), alpha = 0.25, color = 'transparent') +
    geom_line(size = 1.5) +
    annotate("text", x = as.Date("2010-05-01", format = "%Y-%m-%d"), y = 95, label = "A)", size = 8) +
    scale_color_manual(values = c("firebrick", "steelblue"), name = "Modality") +
    scale_fill_manual(values = c("firebrick", "steelblue"), name = "Modality") +
    ylab(NULL) +
    xlab("Month") +
    theme_bw() +
    theme(axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.position = c(0.12, 0.75),
          strip.text = element_text(face = "bold.italic", size = 16),
          strip.background = element_rect(color = "transparent"),
          panel.background = element_rect(fill = "white", color = "black"))

Fig3b <- ggplot(Seas_pred %>% 
           filter(Taxa == "Veliger" & DOY == 200 & Sample_site %in% c("36-873", "16-970")) %>% 
           mutate(Year = ifelse(Region == "Bimodal", Year, Year + 0.2)),
       aes(y = ppmedianBM, x = Year, color = Region, fill = Region)) +
    geom_errorbar(aes(ymin = .lower, ymax = .upper)) +
    geom_point(size = 4, shape = 21, color = "black") +
    annotate("text", x = 1995, y = 620, label = "B)", size = 8) +
    scale_color_manual(values = c("firebrick", "steelblue"), name = "Modality") +
    scale_fill_manual(values = c("firebrick", "steelblue"), name = "Modality") +
    ylab(NULL) +
    xlab("Year") +
    theme_bw() +
    theme(axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          legend.position = c(0.12, 0.75),
          strip.text = element_text(face = "bold.italic", size = 16),
          strip.background = element_rect(color = "transparent"),
          panel.background = element_rect(fill = "white", color = "black"))


Fig3yaxisLab <- textGrob(
    expression(paste(italic("Dreissena spp."), " veliger biomass (µg ", L^-1,")")),
    gp = gpar(fontsize = 22),
    rot = 90
)

png(file.path(here::here("05_Figures"), "05_Fig3_Veligers.png"), 
    units = "in", height = 8, width = 8, res = 300)        
grid.arrange(arrangeGrob(plot_grid(Fig3a, Fig3b, 
                                   nrow = 2,
                                   ncol = 1),
                         left = Fig3yaxisLab))
dev.off()


# export phenology data ----
write.csv(ZP_moments, "03_GeneratedData/05_ZPPhenologyMetrics.csv")

# save/load ----
# save.image("04_SavedRImages/05_Fig2_SeasonalityPhenol.RData")
# load("04_SavedRImages/05_Fig2_SeasonalityPhenol.RData")

