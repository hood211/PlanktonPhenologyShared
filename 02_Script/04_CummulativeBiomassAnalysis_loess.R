# Median cummulative accumulation
# 30 Aug 23 - JMH
# 27 Oct 23

# library ----
library(tidyverse)
library(mgcv)
library(mgcViz)
library(tidygam)
library(boot)
library(MuMIn)

# Data ----
ZP0 <-  read.csv(file.path(here::here("03_GeneratedData/"),"01_LEPASzpData.csv" ), row.names = 1) %>% 
  mutate(Date = as.POSIXct(Sample_date, format = "%Y-%m-%d")) %>% 
  filter(DOY > 125 & DOY < 275) %>% 
  filter(Y <= 2022) %>% 
  #excluding 2020 because COVID delayed sampling
  filter(Y != 2020) %>% 
  select(Sample_site, Date, DOY, Y, Veliger, Dretrocurva, Mesocyclops, Soregonensis) %>% 
  pivot_longer(cols = Veliger:Soregonensis, names_to = "Taxa", values_to = "Biomass") %>% 
  arrange(Taxa, Sample_site, Y, DOY) %>% 
  group_by(Taxa, Sample_site, Y, DOY) %>% 
  summarise(Biomass = mean(Biomass, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(Taxa, Sample_site, Y) %>% 
  mutate(CumSum = cumsum(Biomass)) %>% 
  group_by(Taxa, Sample_site, Y) %>% 
  mutate(MaxYearCumSum = max(CumSum),
         MinYearCumSum = min(CumSum),
         PorCumSum = CumSum/MaxYearCumSum,
         Sample_site = as.factor(Sample_site),
         SiteGroup = case_when(Sample_site %in% c("27-918", "29-905") ~ "s2729",
                               Sample_site %in% c("36-873", "37-890") ~ "s3637",
                               Sample_site %in% c("14-972", "16-970") ~ "s1416",
                               Sample_site %in% c("3-996", "8-994") ~ "s38"),
         SiteGroup = as.factor(SiteGroup),
         TaxaSiteY = paste0(Taxa,"_", Sample_site,"_", Y),
         TaxaSiteGroupY = paste0(Taxa,"_", SiteGroup,"_", Y)) %>% 
  # removes 8 oreg Site*Y combos
  filter(!is.na(PorCumSum))


# model cumsum ----
## clean data ----
ZP <- ZP0 %>%
  # delete anything with <4 pts
  group_by(TaxaSiteY) %>% 
  filter(n() >3) %>% 
  # delete anything that doesn't have data before 160 or after 220
  mutate(DOYmin = min(DOY),
         DOYmax = max(DOY),
         DOYviolation = ifelse(DOYmin > 160 | DOYmax < 220, "bad", "good")) %>% 
  filter(DOYviolation != "bad")
  # oreg: common for cum sum to increase to end of smp period

# Taxa*Site*Y lost
# lost 70 of 852
length(unique(ZP0$TaxaSiteY)) - length(unique(ZP$TaxaSiteY))

  # generate matrix
  
## Loess function ----
  LoessFunJMH <- function(TaxaSiteYlist, LoessTSdat, span4fun){
    for(i in 1:length(TaxaSiteYlist)) { #
      # i = 1
      TaxaSiteY_i <- TaxaSiteYlist[i]
      df_i <- ZP[ZP$TaxaSiteY == TaxaSiteY_i,]
      # LoessMod_i <- loess(PorCumSum ~ DOY, data = df_i, span = 0.3, degree = 2)
      
      skip_to_next1 <- FALSE
      
      LoessMod_i <- tryCatch(loess(PorCumSum ~ DOY, data = df_i, 
                                   span = span4fun, degree = 2),
                             error = function(e) {skip_to_next1 <<- TRUE})
      if(skip_to_next1) {next}
      
      LoessTSi <- seq(min(df_i$DOY), max(df_i$DOY), length.out = 200)
      LoessTSi <- as.data.frame(LoessTSi)
      names(LoessTSi) <-  "DOY"
      
      
      # LoessTSi$predLoess <- predict(LoessMod_i, LoessTSi)
      
      skip_to_next2 <- FALSE
      
      LoessTSi$predLoess <- tryCatch(predict(LoessMod_i, LoessTSi), error = function(e) {skip_to_next2 <<- TRUE})
      
      if(skip_to_next2) {next}
      
      LoessTSi$TaxaSiteY <- TaxaSiteY_i
      
      LoessTSdat <- rbind(LoessTSdat,LoessTSi)
    }
    LoessTSdat
  }
  

## Est loess fits ----  
  ### span 0.2 ----
  TaxaSiteYlist_s2 <- unique(ZP$TaxaSiteY)
  LoessTSdat_s2 <- as.data.frame(matrix(nrow = 1, ncol = 3))
  names(LoessTSdat_s2) <- c("DOY", "predLoess", "TaxaSiteY")
  LoessTSdat_s2 <- LoessFunJMH(TaxaSiteYlist_s2, LoessTSdat_s2, 0.2) %>% filter(!is.na(DOY))
  length(unique(ZP$TaxaSiteY)) - length(unique(LoessTSdat_s2$TaxaSiteY)) #402 remaining
  
  ### span 0.3 ----
  TaxaSiteYlist_s3 <- setdiff(unique(ZP$TaxaSiteY), unique(LoessTSdat_s2$TaxaSiteY))
  LoessTSdat_s3 <- as.data.frame(matrix(nrow = 1, ncol = 3))
  names(LoessTSdat_s3) <- c("DOY", "predLoess", "TaxaSiteY")
  LoessTSdat_s3 <- LoessFunJMH(TaxaSiteYlist_s3, LoessTSdat_s3, 0.3)
  length(TaxaSiteYlist_s3) - length(unique(LoessTSdat_s3$TaxaSiteY)) #31 remaining
  
  ### span 0.4 ----
  TaxaSiteYlist_s4 <- setdiff(TaxaSiteYlist_s3, unique(LoessTSdat_s3$TaxaSiteY))
  LoessTSdat_s4 <- as.data.frame(matrix(nrow = 1, ncol = 3))
  names(LoessTSdat_s4) <- c("DOY", "predLoess", "TaxaSiteY")
  LoessTSdat_s4 <- LoessFunJMH(TaxaSiteYlist_s4, LoessTSdat_s4, 0.4) %>% filter(!is.na(DOY))
  length(TaxaSiteYlist_s4) - length(unique(LoessTSdat_s4$TaxaSiteY)) # 4 remaining
  
  ### span 0.5 ----
  TaxaSiteYlist_s5 <- setdiff(TaxaSiteYlist_s4, unique(LoessTSdat_s4$TaxaSiteY))
  LoessTSdat_s5 <- as.data.frame(matrix(nrow = 1, ncol = 3))
  names(LoessTSdat_s5) <- c("DOY", "predLoess", "TaxaSiteY")
  LoessTSdat_s5 <- LoessFunJMH(TaxaSiteYlist_s5, LoessTSdat_s5, 0.5) %>% filter(!is.na(DOY))
  length(TaxaSiteYlist_s5) - length(unique(LoessTSdat_s5$TaxaSiteY)) # 0 remaining
  
  ### combine ----
  LoessTSdat_sAll <- rbind(LoessTSdat_s2, LoessTSdat_s3, LoessTSdat_s4, LoessTSdat_s5) %>% 
    filter(!is.na(DOY))
  
  
  png(file.path(here::here("05_Figures"), "04_AllSitePorCumSum_AllModelsDat.png"), 
      units = "in", height = 30, width = 20, res = 300)
  
  ggplot(LoessTSdat_sAll %>% 
           separate(TaxaSiteY, sep = "_", into = c("Taxa", "Site", "Y")) %>% 
           filter(!(Site %in% c("36-873", "37-890", "27-918", "29-905") &
                      Taxa == "Soregonensis")), 
         aes(y = predLoess, x = DOY, color = Y)) +
    geom_hline(yintercept = 0.15) +
    geom_hline(yintercept = 0.5) +
    geom_hline(yintercept = 0.85) +
    geom_line() +
    geom_point(data = ZP %>% 
                 rename(Site = Sample_site) %>% 
                 filter(!(Site %in% c("36-873", "37-890", "27-918", "29-905") &
                            Taxa == "Soregonensis")), aes(y = PorCumSum, x = DOY, color = as.factor(Y))) +
    facet_grid(Taxa + Site ~ Y) +
    theme_bw()
  dev.off()
  
  ### clean ----
  # Oreg not looking at sites i'm dropping
  # Oreg, 14, 1995 flat
  # Oreg, 3, 1995 almost flat
  # Oreg, 8, 1996 almost flat
  # Vel handful where we miss 15%
  # Vel, 3, 1995 - miss 15% and 50%
  
  LoessTSdat_sAllc <- LoessTSdat_sAll %>% 
    filter(!(TaxaSiteY %in% c("Soregonensis_14-972_1995",
                              "Soregonensis_3-996_1995",
                              "Soregonensis_8-994_1995",
                              "Veliger_3-996_1995"))) %>% 
    separate(TaxaSiteY, sep = "_", into = c("Taxa", "Site", "Y")) %>% 
    filter(!(Site %in% c("36-873", "37-890", "27-918", "29-905") &
               Taxa == "Soregonensis"))
  
# calculate moments ----
  p_cumsum_15per <- LoessTSdat_sAllc %>% 
    group_by(Taxa, Site, Y) %>% 
    slice_min(abs(predLoess-0.15)) %>% 
    slice_min(DOY) %>% 
    mutate(Moment = "per15")
  
  p_cumsum_50per <- LoessTSdat_sAllc %>% 
    group_by(Taxa, Site, Y) %>% 
    slice_min(abs(predLoess-0.5)) %>% 
    slice_min(DOY) %>% 
    mutate(Moment = "per50")
  
  p_cumsum_85per <- LoessTSdat_sAllc %>% 
    group_by(Taxa, Site, Y) %>% 
    slice_min(abs(predLoess-0.85)) %>% 
    slice_min(DOY) %>% 
    mutate(Moment = "per85")

  p_cumsum_comb <- rbind(p_cumsum_15per, p_cumsum_50per, p_cumsum_85per) %>% 
    # exclude cases where %s are not near target - 25 in total
    filter(!((Moment == "per15" & predLoess > 0.2) |
             (Moment == "per50" & predLoess > 0.55) |
             (Moment == "per85" & predLoess > 0.90)))

  ggplot(p_cumsum_comb, aes(x = predLoess, fill = Moment)) +
    geom_histogram() +
    facet_wrap(~Taxa) +
    xlim(0,1)
  
  ggplot(p_cumsum_comb, aes(y = DOY, x = as.numeric(Y), color = Moment)) +
    geom_point() +
    geom_line() +
    facet_grid(Site ~ Taxa) +
    theme_bw()

# Join w/ BM ----
  
  p_cumsum_comb_LT <- p_cumsum_comb %>% 
                        left_join(ZP %>% 
                                    mutate(Y = as.character(Y)) %>% 
                                    ungroup() %>% 
                                    select(Taxa:PorCumSum), 
                                  by = join_by(Taxa == Taxa,
                                                   Site == Sample_site,
                                                   Y == Y,
                                                   closest(DOY >= DOY)),
                                  suffix = c(".cs", ".bmLT"))
  
  p_cumsum_comb_B <- p_cumsum_comb_LT %>% 
                        left_join(ZP %>% 
                                    mutate(Y = as.character(Y)) %>% 
                                    ungroup() %>% 
                                    select(Taxa:PorCumSum), 
                                  by = join_by(Taxa == Taxa,
                                               Site == Sample_site,
                                               Y == Y,
                                               closest(DOY.cs <= DOY)),
                                  suffix = c(".bmLT", ".bmGT")) %>% 
                        rename(DOY.bmGT = DOY) %>% 
                        mutate(DOYdif_GT = DOY.bmGT - DOY.cs,
                               DOYdif_LT = DOY.cs - DOY.bmLT,
                               DOY.bm = ifelse(is.na(DOYdif_GT) | is.na(DOYdif_LT), 
                                               coalesce(DOY.bmLT, DOY.bmGT),
                                               ifelse(DOYdif_GT < DOYdif_LT, 
                                                      DOY.bmGT, DOY.bmLT)),
                               # biomass
                               Biomass.bm = ifelse(is.na(DOYdif_GT) | is.na(DOYdif_LT), 
                                                   coalesce(Biomass.bmLT, Biomass.bmGT),
                                                   ifelse(DOYdif_GT < DOYdif_LT, 
                                                          Biomass.bmGT, Biomass.bmLT)),
                               # Cum sum
                               CumSum.bm = ifelse(is.na(DOYdif_GT) | is.na(DOYdif_LT), 
                                                  coalesce(CumSum.bmLT, CumSum.bmGT),
                                                  ifelse(DOYdif_GT < DOYdif_LT, 
                                                         CumSum.bmGT, CumSum.bmLT)),
                               # Max cum sum
                               MaxYearCumSum.bm = ifelse(is.na(DOYdif_GT) | is.na(DOYdif_LT), 
                                                         coalesce(MaxYearCumSum.bmLT, MaxYearCumSum.bmGT),
                                                         ifelse(DOYdif_GT < DOYdif_LT, 
                                                                MaxYearCumSum.bmGT, MaxYearCumSum.bmLT))) %>% 
                        rowwise() %>% 
                        mutate(DOYdifmin = min(DOYdif_GT, DOYdif_LT, na.rm = T)) %>% 
                        select(-c(DOY.bmLT:DOYdif_LT))


  hist(p_cumsum_comb_B$DOYdifmin)
  ggplot(p_cumsum_comb_B, aes(x = DOYdifmin, fill = Moment)) +
    geom_density(alpha = 0.5)
  
  # 14 with > 2 wk between bm and moment, max is 17 d
  view(p_cumsum_comb_B[p_cumsum_comb_B$DOYdifmin > 14,])

  ggplot(p_cumsum_comb_B, aes(y = DOYdifmin, x = Y, color = Moment)) +
    geom_point()+
    facet_grid(Site~Taxa)
  
  ggplot(p_cumsum_comb_B, aes(y = MaxYearCumSum.bm, x = DOY.cs, color = Moment)) +
    geom_point() +
    facet_wrap(~Taxa)
  
  write.csv(p_cumsum_comb_B, "03_GeneratedData/04_ZPCumSumPercentiles_LOESSEst.csv")
  
# save/load ----
  # save.image(file.path(here::here("04_SavedRImages"),"04_CummulativeBiomassAnalysis_loess.Rdata"))
  # load(file.path(here::here("04_SavedRImages"),"04_CummulativeBiomassAnalysis_loess.Rdata"))
  
  
 