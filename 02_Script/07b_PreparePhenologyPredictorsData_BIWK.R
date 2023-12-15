# combine phenology & predictor data
# Need three different datasets

# ret PB, meso PB, oreg PB, and Oreg %
  # Basin average
# vel PB
  # bimodal ("16-970", "3-996", "8-994") and unimodal ("36-873","37-890","29-905","27-918", "14-972") sites 
# Ret %, Meso %, Veliger %
  # site groups: ("36-873","37-890"), ("29-905","27-918"), ("14-972", "16-970"), and ("3-996", "8-994")


# JMH Oct 2023, Dec 2023
library(tidyverse)
library(ggcorrplot)
library(cowplot)

# data ----
# phenology metrics
ZPphen <- read.csv("03_GeneratedData/05_ZPPhenologyMetrics.csv", row.names = 1)

## fish ----
fish <- read.csv(file.path(here::here("03_GeneratedData"),"06a_WBplanktivoreDensity_NotSummed2Month.csv"), row.names = 1) 

## plankton ----
Plankton0 <- read.csv("03_GeneratedData/06b_PlanktonTempDat_BiWk.csv", row.names = 1) 


Plankton <- Plankton0 %>% 
              select(-Sample_ID) %>% 
              pivot_wider(id_cols = c(Sample_date, Response), names_from = Sample_site, values_from = Value) %>% 
              # filter(Response != "Byth") # ONLY FOR COR PLOTS, byth has 2 sites
              rowwise() %>% 
              # need to come up with Byth data for sites other than 27 and 37 - calc avg of these sites
              mutate(Average = ifelse(Response == "Byth", mean(c(`27-918`,`37-890`), na.rm = T), as.numeric("NA")),
                     `14-972` = ifelse(Response == "Byth", Average, `14-972`),
                     `16-970` = ifelse(Response == "Byth", Average, `16-970`),
                     `29-905` = ifelse(Response == "Byth", Average, `29-905`),
                     `3-996` = ifelse(Response == "Byth", Average, `3-996`),
                     `36-873` = ifelse(Response == "Byth", Average, `36-873`),
                     `8-994` = ifelse(Response == "Byth", Average, `8-994`),
                     # also fill in NA's for FTG sites
                     `27-918` = ifelse(Response == "Byth" & is.na(`27-918`), Average, `27-918`),
                     `37-890` = ifelse(Response == "Byth" & is.na(`37-890`), Average, `37-890`)) %>% 
              select(-Average) %>% 
              pivot_longer(cols = '3-996':`16-970`, names_to = "Sample_site", values_to = "Values") %>% 
              select(Response, Sample_date, Sample_site, Values)

# 
ggplot(Plankton0 %>% 
         filter(Response == "Byth") %>% 
         mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d")),
                aes(y = Value, x = Sample_date)) +
  geom_point() +
  facet_wrap(~Sample_site)

# need to run ONLY first 4 lines of Plankton to make this work
    pdf("05_Figures/07_CorrAmongSites.pdf")
    for(i in 1:(length(unique(Plankton$Response)))){
      Response_i <- unique(Plankton$Response)[i]

      cor_i <- cor(Plankton[Plankton$Response == Response_i,-c(1,2)],
                   use = "pairwise.complete.obs", method = "spearman")
      pval_i <- cor_pmat(Plankton[Plankton$Response == Response_i,-c(1,2)],
                         use = "pairwise.complete.obs")
      corrPlot_i <-  ggcorrplot(cor_i, p.mat = pval_i, outline.col = "white", type = "lower",
                                show.diag = T, lab = T) + #insig = "blank",
        ggtitle(Response_i)
      print(corrPlot_i)
    }
    dev.off()

## ice ----
# basin wide, ANNUAL
ice0 <- read.csv("03_GeneratedData/06c_LErieMeanIce.csv", row.names = 1) %>% 
              mutate(Response = "MeanIce",
                     Sample_site = "none", 
                     Month = "none") %>% 
              select(Response, Y = Year, Sample_site, Month, Values = MeanIce)

  # need to dup for all sites
  ice0.5 <- rbind(ice0 %>% 
                 mutate(Sample_site = "3-996"),
               ice0 %>% 
                 mutate(Sample_site = "8-994"),
               ice0 %>% 
                 mutate(Sample_site = "14-972"),
               ice0 %>% 
                 mutate(Sample_site = "16-970"),
               ice0 %>% 
                 mutate(Sample_site = "27-918"),
               ice0 %>% 
                 mutate(Sample_site = "29-905"),
               ice0 %>% 
                 mutate(Sample_site = "36-873"),
               ice0 %>% 
                 mutate(Sample_site = "37-890")) %>% 
              select(Response, Sample_site, Year = Y, BiWOY = Month, Values) %>% 
              mutate(BiWOY = as.numeric("NA")) %>%
              group_by(Response, Sample_site, BiWOY) %>% 
              mutate(across(Values, \(x) scale(x, center = T, scale = T)))

  ice <- cbind(ice0.5[,1:4], ice0.5$Values[,1])
  names(ice) <- c("Response", "Sample_site", "Year", "BiWOY", "Values")

## combine fish/plankton ----
  
Pred0 <- rbind(fish %>% 
                select(Response, Sample_date, Sample_site, Values = Value),
              Plankton) %>% 
              mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d"),
                     Year = as.numeric(strftime(Sample_date, format = "%Y")),
                     BiWOY = floor(as.numeric(strftime(Sample_date, format = "%U"))/2)) 
  
  Pred0 %>% split(.$Response) %>% map(summary)
  Byth_HalfMin <- min(Pred0[Pred0$Response == "Byth" & Pred0$Values >0 & !is.na(Pred0$Values),]$Values)/2
  Chla_HalfMin <- min(Pred0[Pred0$Response == "Chla" & Pred0$Values >0 & !is.na(Pred0$Values),]$Values)/2
  Edible_HalfMin <- min(Pred0[Pred0$Response == "EdiblePP" & Pred0$Values >0 & !is.na(Pred0$Values),]$Values)/2
  InEdible_HalfMin <- min(Pred0[Pred0$Response == "InEdiblePP" & Pred0$Values >0 & !is.na(Pred0$Values),]$Values)/2
  Lepto_HalfMin <- min(Pred0[Pred0$Response == "Lepto" & Pred0$Values >0 & !is.na(Pred0$Values),]$Values)/2
  MesoZP_HalfMin <- min(Pred0[Pred0$Response == "MesoZPfood" & Pred0$Values >0 & !is.na(Pred0$Values),]$Values)/2
  PPLthen10_HalfMin <- min(Pred0[Pred0$Response == "PPLthen10" & Pred0$Values >0 & !is.na(Pred0$Values),]$Values)/2
  
  Pred <- Pred0%>% 
            # Replace zeros with 1/2 min
            mutate(Values = ifelse(Response == "Lepto" & Values == 0, Lepto_HalfMin,
                                   ifelse(Response == "InEdiblePP" & Values == 0, InEdible_HalfMin,
                                          ifelse(Response == "Chla" & Values == 0, Chla_HalfMin,
                                                 ifelse(Response == "Byth" & Values == 0, Byth_HalfMin, 
                                                    ifelse(Response == "EdiblePP" & Values == 0, Edible_HalfMin,
                                                      ifelse(Response == "PPLthen10" & Values == 0, PPLthen10_HalfMin,
                                                          ifelse(Response == "MesoZPfood" & Values == 0, MesoZP_HalfMin, 
                                                                 Values)))))))) %>%
            # average to biweek
            group_by(Response, Sample_site, Year, BiWOY) %>% 
            summarise(Values = mean(Values,na.rm = T)) %>% 
            # log all but ice and temp
            mutate(Values = ifelse(!(Response %in% c("AvgTemp", "Ice")), log(Values), Values)) %>% 
            # Scale within each BiWeek across years
            group_by(Response, BiWOY) %>%
            mutate(across(Values, \(x) scale(x, center = T, scale = T)))

# this removes the "[,1]'
Pred2 <- cbind(Pred[,1:4], Pred$Values[,1])
names(Pred2)[5] <- "Values"
    

## combine fish/plank/ice ----
Pred3 <- Pred2 %>% 
           pivot_wider(id_cols = c("Sample_site", "Year", "BiWOY"), names_from = Response, values_from = Values) %>% 
          # join with ice, note this generates dups of ice - intentional
          left_join(ice %>% 
                      ungroup() %>% 
                      select(Sample_site, Year, MeanIce = Values), by = c("Sample_site", "Year")) %>% 
          pivot_longer(cols = AvgTemp:MeanIce, names_to = "Response", values_to = "BiWOYanom") %>% 
          select(Response, Sample_site, Year, BiWOY, BiWOYanom) %>% 
          # missing too much data before biwk 9 and after 20
          filter(BiWOY >= 9 & BiWOY <= 20) %>% 
          droplevels()


# makes three df ----
##  basin wide ----
Pred3_WB <- Pred3 %>% 
              # average over sites
              group_by(Response, Year, BiWOY) %>% 
              summarise(across(BiWOYanom, \(x) mean(x, na.rm = T))) %>% 
              droplevels()

    # checks full Pred3
    ggplot(Pred3_WB, aes(y = BiWOYanom, x = Year)) +
      geom_point() +
      geom_line() +
      facet_grid(Response ~ BiWOY)
    
    ### Cor among months ----
    Pred_WB_wide <- Pred3_WB %>% 
      mutate(BiWOY2 = as.factor(paste0("bw", BiWOY)),
             BiWOY2 = fct_relevel(BiWOY2,
                                  "bw9", "bw10", "bw11", "bw12", "bw13", "bw14",
                                  "bw15", "bw16", "bw17", "bw18", "bw19", "bw20")) %>% 
      select(-BiWOY) %>% 
      pivot_wider(id_cols = Response:Year, names_from = BiWOY2, values_from = "BiWOYanom") %>% 
      select(Response:Year, bw9, bw10:bw15, bw16, bw17:bw20)
    
    # cor heatmap
    # errors associated with ice which is lake-wide and annual
    pdf("05_Figures/07b_PredictorsCorrAmongMonths_Pred_WB_wide_BiWeek.pdf")
    for(i in 1:(length(unique(Pred_WB_wide$Response))-1)){
      Response_i <- unique(Pred_WB_wide$Response)[i]
      
      cor_i <- cor(Pred_WB_wide[Pred_WB_wide$Response == Response_i,-c(1,2)], 
                   use = "pairwise.complete.obs", method = "spearman")
      pval_i <- cor_pmat(Pred_WB_wide[Pred_WB_wide$Response == Response_i,-c(1,2)], 
                         use = "pairwise.complete.obs")
      corrPlot_i <-  ggcorrplot(cor_i, p.mat = pval_i, outline.col = "white", type = "lower",
                                show.diag = T, lab = T) + #insig = "blank", 
        ggtitle(Response_i)
      print(corrPlot_i)
    }
    dev.off()


## veliger Modality ----
Pred_VelGroups <- Pred3 %>% 
              mutate(VelGroups = ifelse(Sample_site %in% c("16-970", "3-996", "8-994"), "Bimodal",
                                  ifelse(Sample_site %in% c("36-873","37-890","29-905","27-918", "14-972"), "Unimodal", "NA")),
                     VelGroups = as.factor(VelGroups)) %>% 
              # average to VelGroups
              group_by(Response, Year, BiWOY, VelGroups) %>% 
              summarise(across(BiWOYanom, \(x) mean(x, na.rm = T))) %>% 
      mutate(BiWOY2 = as.factor(paste0("bw", BiWOY)),
             BiWOY2 = fct_relevel(BiWOY2,
                                  "bw9", "bw10", "bw11", "bw12", "bw13", "bw14",
                                  "bw15", "bw16", "bw17", "bw18", "bw19", "bw20")) %>% 
      ungroup() %>% 
      select(-BiWOY) 
      

      # checks full Pred3
      ggplot(Pred_VelGroups, aes(y = BiWOYanom, x = Year, color = VelGroups)) +
        geom_point() +
        geom_line() +
        facet_grid(Response ~ BiWOY2)
      
      ### Cor among months ----
      Pred_VelGroups_wide <- Pred_VelGroups %>% 
        pivot_wider(id_cols = Response:VelGroups, names_from = BiWOY2, values_from = "BiWOYanom") %>%
        select(Response:VelGroups, bw9, bw10:bw15, bw16, bw17, bw18:bw20) 
        
      Pred_VelGroups_wide %>% split(.$Response) %>% map(summary)
    
      # doesn't do fish, fails due to missing data
      pdf("05_Figures/07b_PredictorsCorrAmongMonths_Pred_VelGroups_wide_BiWeek_NoFish.pdf")
      for(i in 1:(length(unique(Pred_VelGroups_wide$Response)))){
        Response_i <- unique(Pred_VelGroups_wide$Response)[i]
        
        cor_i <- cor(Pred_VelGroups_wide[Pred_VelGroups_wide$Response == Response_i,-c(1:3)], 
                     use = "pairwise.complete.obs", method = "spearman")
        pval_i <- cor_pmat(Pred_VelGroups_wide[Pred_VelGroups_wide$Response == Response_i,-c(1:3)], 
                           use = "pairwise.complete.obs")
        corrPlot_i <-  ggcorrplot(cor_i, p.mat = pval_i, outline.col = "white", type = "lower",
                                  show.diag = T, lab = T) + #insig = "blank", 
          ggtitle(Response_i)
        print(corrPlot_i)
      }
      dev.off()
      
      # Fish - DF just for corr heatmap
      Pred_VelGroups_wide_fish <- Pred_VelGroups_wide %>% 
        filter(Response == "plank_per_ha") %>% 
        select(-c(bw9, bw18, bw20, bw11))
      
      pdf("05_Figures/07b_PredictorsCorrAmongMonths_Pred_VelGroups_wide_BiWeek_FishOnly.pdf")
      for(i in 1:(length(unique(Pred_VelGroups_wide_fish$Response)))){
        Response_i <- unique(Pred_VelGroups_wide_fish$Response)[1]
        
        cor_i <- cor(Pred_VelGroups_wide_fish[Pred_VelGroups_wide_fish$Response == Response_i,-c(1:3)], 
                     use = "pairwise.complete.obs", method = "spearman")
        # pval_i <- cor_pmat(Pred_VelGroups_wide_fish[Pred_VelGroups_wide_fish$Response == Response_i,-c(1:3)], 
        #                    use = "pairwise.complete.obs")
        corrPlot_i <-  ggcorrplot(cor_i,  outline.col = "white", type = "lower",
                                  show.diag = T, lab = T) + #insig = "blank", 
          ggtitle(Response_i)
        print(corrPlot_i)
      }
      dev.off()

# ## Sites ----
Pred_Sites <- Pred3 %>%
          # average to Sitegroup
          group_by(Response, Sample_site, Year, BiWOY) %>%
          summarise(across(BiWOYanom, \(x) mean(x, na.rm = T))) %>%
          droplevels() %>%
        mutate(BiWOY2 = as.factor(paste0("bw", BiWOY)),
               BiWOY2 = fct_relevel(BiWOY2,
                                    "bw9", "bw10", "bw11", "bw12", "bw13", "bw14",
                                    "bw15", "bw16", "bw17", "bw18", "bw19", "bw20")) %>%
        ungroup() %>%
        select(-BiWOY)


# checks full Pred3
ggplot(Pred_Sites, aes(y = BiWOYanom, x = Year, color = Sample_site)) +
  geom_point() +
  geom_line() +
  facet_grid(Response ~ BiWOY2)

### Cor among months ----
Pred_Sites_w <- Pred_Sites %>%
  pivot_wider(id_cols = c(Response, Sample_site, Year), names_from = BiWOY2, values_from = BiWOYanom) %>%
  select(Response:Year, bw9, bw10:bw15, bw16, bw17, bw18:bw20)

# doesn't do ice
pdf("05_Figures/07b_PredictorsCorrAmongMonths_Pred_Sites_w_BiWeek.pdf")
for(i in 1:(length(unique(Pred_Sites_w$Response))-1)){
  Response_i <- unique(Pred_Sites_w$Response)[i]

  cor_i <- cor(Pred_Sites_w[Pred_Sites_w$Response == Response_i,-c(1:3)],
                  use = "pairwise.complete.obs", method = "spearman")
  pval_i <- cor_pmat(Pred_Sites_w[Pred_Sites_w$Response == Response_i,-c(1:3)],
                        use = "pairwise.complete.obs")
 corrPlot_i <-  ggcorrplot(cor_i, p.mat = pval_i, outline.col = "white", type = "lower",
                           show.diag = T, lab = T) + #insig = "blank",
    ggtitle(Response_i)
 print(corrPlot_i)
}
dev.off()


# "Autocorrelation" -----
Pred_Sites_cordf <- as.data.frame(matrix(nrow = 1, ncol = 4))
names(Pred_Sites_cordf) <- c("RowBiWk", "ColBiWk", "SpearmanCor", "Predictor")


for(i in 1:(length(unique(Pred_Sites_w$Response)))){
  # i = 1
  Response_i <- unique(Pred_Sites_w$Response)[i]
  
  cor_i <- cor(Pred_Sites_w[Pred_Sites_w$Response == Response_i,-c(1:3)], 
               use = "pairwise.complete.obs", method = "spearman")
  
  cor_i2 <- as.data.frame(cor_i) %>% 
            mutate(RowBiWk = row.names(cor_i)) %>% 
            pivot_longer(cols = bw9:bw20, names_to = "ColBiWk", values_to = "SpearmanCor") %>% 
            mutate(Predictor = Response_i)
  
  Pred_Sites_cordf <- rbind(Pred_Sites_cordf,cor_i2)
}

BiWk2BiWkPairs <- c("bw9_bw10",
                    "bw10_bw11",
                    "bw11_bw12",
                    "bw12_bw13",
                    "bw13_bw14",
                    "bw14_bw15",
                    "bw15_bw16",
                    "bw17_bw18",
                    "bw18_bw19",
                    "bw19_bw20")

# Duplicates with M_t-1 then M_t
BiWkPairs2Keep <- c("bw9_bw10",  "bw9_bw11",  "bw9_bw12",  "bw9_bw13",  "bw9_bw14",  "bw9_bw15",  "bw9_bw16",  "bw9_bw17" ,
  "bw9_bw18",  "bw9_bw19",  "bw9_bw20",  
  "bw10_bw11", "bw10_bw12", "bw10_bw13", "bw10_bw14","bw10_bw15", "bw10_bw16", "bw10_bw17", "bw10_bw18", 
  "bw10_bw19", "bw10_bw20", 
  "bw11_bw12", "bw11_bw13", "bw11_bw14", "bw11_bw15", "bw11_bw16", "bw11_bw17","bw11_bw18", "bw11_bw19", 
  "bw11_bw20",
  "bw12_bw13", "bw12_bw14", "bw12_bw15", "bw12_bw16", "bw12_bw17","bw12_bw18", "bw12_bw19", "bw12_bw20", 
  "bw13_bw14","bw13_bw15", "bw13_bw16", "bw13_bw17", "bw13_bw18", "bw13_bw19", "bw13_bw20" ,
  "bw14_bw15", "bw14_bw16", "bw14_bw17", "bw14_bw18", "bw14_bw19", "bw14_bw20",
  "bw15_bw16", "bw15_bw17", "bw15_bw18", "bw15_bw19", "bw15_bw20",
  "bw16_bw17", "bw16_bw18", "bw16_bw19", "bw16_bw20", 
  "bw17_bw18", "bw17_bw19", "bw17_bw20",
  "bw18_bw19", "bw18_bw20",
  "bw19_bw20" )

One2One2remove <- c("bw9_bw9",   
                    "bw10_bw10", 
                    "bw11_bw11",
                    "bw12_bw12",
                    "bw13_bw13", 
                    "bw14_bw14", 
                    "bw15_bw15", 
                    "bw16_bw16", 
                    "bw17_bw17", 
                    "bw18_bw18",
                    "bw19_bw19", 
                    "bw20_bw20")


Pred_Sites_cordf2 <- as.data.frame(Pred_Sites_cordf[-1,]) %>% 
                  mutate(BiWkPair = paste0(RowBiWk,"_",ColBiWk)) %>% 
                  # removes r = 1 and duplicates
                  filter(BiWkPair %in% BiWk2BiWkPairs | BiWkPair %in% BiWkPairs2Keep) %>% 
                  mutate(PairGroups = ifelse(BiWkPair %in% BiWk2BiWkPairs,"BiWk2tm1s",
                                        ifelse(BiWkPair  %in% BiWkPairs2Keep, "BiWk2Others", "blah")),
                         Predictor = as.factor(Predictor),
                         Predictor = fct_recode(Predictor, 
                                                 "Temp" = "AvgTemp",
                                                "PPngd" = "EdiblePP",
                                                "PPgd" = "InEdiblePP",
                                                "Ice" = "MeanIce",
                                                "PPL10" = "PPLthen10",
                                                "Fish" = "plank_per_ha"))





PPred_Sites_cordf2_t2tm1 <- Pred_Sites_cordf2 %>% 
                  filter(PairGroups == "BiWk2tm1s") 

PPred_Sites_cordf2_Other <- Pred_Sites_cordf2 %>% 
  filter(PairGroups == "BiWk2Others") 




p1_m2m <- ggplot(PPred_Sites_cordf2_t2tm1 %>% 
         filter(Predictor != "Ice") %>% 
         mutate(Predictor = fct_reorder(as.factor(Predictor), SpearmanCor, .fun = median, .desc = T)), 
       aes(y = SpearmanCor, x = Predictor, label = BiWkPair)) +
  geom_boxplot() +
  geom_text()+
  ggtitle("Cor b/w biwk t and t-1") +
  ylim(-1,1) 

p2_allM <- ggplot(PPred_Sites_cordf2_Other %>% 
         filter(Predictor != "Ice") %>% 
         mutate(Predictor = fct_reorder(as.factor(Predictor), SpearmanCor, .fun = median, .desc = T)), 
       aes(y = SpearmanCor, x = Predictor, label = BiWkPair)) +
  geom_boxplot() +
  geom_text() +
  ggtitle("Cor w/ any biwk but t-1") +
  ylim(-1,1)


png("05_Figures/07b_PredictorsCorrBoxPlots_BiWks.png",  
    units = "in", height = 8, width = 12, res = 300)
plot_grid(p1_m2m, p2_allM)
dev.off()

# Save df's ----
write.csv(Pred_WB_wide, file.path(here::here("03_GeneratedData"),
      "07b_Predictors_LogedAnoms_BasinAvg_BiWk.csv"))
write.csv(Pred_VelGroups_wide, file.path(here::here("03_GeneratedData"),
      "07b_Predictors_LogedAnoms_VelGroup_BiWk.csv"))
write.csv(Pred_Sites_w, file.path(here::here("03_GeneratedData"),
      "07b_Predictors_LogedAnoms_Site_BiWk.csv"))

# Save/load ----
# save.image(file.path(here::here("04_SavedRImages"), "07b_PreparePhenologyPredictorsData_BiWk.RData"))
# load(file.path(here::here("04_SavedRImages"), "07b_PreparePhenologyPredictorsData_BiWk.RData"))
