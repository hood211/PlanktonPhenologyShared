# Interannual patterns in predictors
# JMH, Apr 2023, Dec 2023

# Libaries ----
library(tidyverse)
library(MARSS)

# Data ----
PredSites <- read.csv(file.path(here::here("03_GeneratedData"),
                                  "07b_Predictors_LogedAnoms_Site_BiWk.csv"), row.names = 1) %>% 
              pivot_longer(cols = bw9:bw20, names_to = "BiWk", values_to = "Value") %>% 
              mutate(Month = case_when(BiWk == "bw9" ~ "May",
                                       BiWk == "bw10" ~ "May",
                                       BiWk == "bw11" ~ "May",
                                       BiWk == "bw12" ~ "Jun",
                                       BiWk == "bw13" ~ "May",
                                       BiWk == "bw14" ~ "Jul",
                                       BiWk == "bw15" ~ "Jul",
                                       BiWk == "bw16" ~ "Aug",
                                       BiWk == "bw17" ~ "Aug",
                                       BiWk == "bw18" ~ "Sep",
                                       BiWk == "bw19" ~ "Sep",
                                       BiWk == "bw20" ~ "Oct")) %>% 
                group_by(Response, Sample_site, Year, Month) %>% 
                summarise(Value = mean(Value, na.rm = TRUE)) %>% 
                filter(!(Response %in% c("MeanIce", "Chla"))) %>% 
                # exclude all Byth but "27-918" and "37-890",
                filter(!(Response == "Byth" & 
                           Sample_site %in% c("14-972", "16-970", "29-905", "3-996", "36-873",  "8-994"))) %>% 
                arrange(Response, Month, Sample_site, Year)

ggplot(PredSites, aes(y = Value, x = Year, color = Sample_site)) +
  geom_point()+
  facet_grid(Response ~ Month)

# MARS function ----
# unique sites
MARSfunSites <- function(MarssDf, Monthf, Responsef, nbootstraps){
  # TESTING ONLY
  # MarssDf <- PredSites
  # Monthf = "May"
  # Responsef = "Lepto"
  
    # Prepare data
  Pred4Mars_M <- MarssDf %>%
    filter(Month == {{Monthf}}) %>%
    filter(Response == {{Responsef}}) %>%
    pivot_wider(id_cols = c(Response, Year), names_from = Sample_site, values_from = Value)
    
  
  # transpose data
  Pred4Mars_M_t <- t(Pred4Mars_M[,-c(1:2)])
  
  # Build MARS model
  MarssMod <- list(
    B = "identity",
    U = "unconstrained", # allow for bias
    C = "unconstrained",
    c = "zero", # not sure about his
    Q = "diagonal and unequal",
    Z = "identity",
    A = "zero",
    R = "diagonal and equal")
  
  # Run MARS mod
  Pred4Mars_Mod <- MARSS(Pred4Mars_M_t,
                         model = MarssMod,
                         control = list(maxit = 5000),
                         method = "BFGS")

  Pred4Mars_Mod
}

MARSfunSitesSame <- function(MarssDf, Monthf, Responsef, nbootstraps){
  # TESTING ONLY
  # MarssDf <- PredSites
  # Monthf = "May"
  # Responsef = "Lepto"
  
  # Prepare data
  Pred4Mars_M <- MarssDf %>%
    filter(Month == {{Monthf}}) %>%
    filter(Response == {{Responsef}}) %>%
    pivot_wider(id_cols = c(Response, Year), names_from = Sample_site, values_from = Value)
  
  
  # transpose data
  Pred4Mars_M_t <- t(Pred4Mars_M[,-c(1:2)])
  
  # Build MARS model
  MarssMod <- list(
    B = "identity",
    U = "unconstrained", # allow for bias
    C = "unconstrained",
    c = "zero", # not sure about his
    Q = "diagonal and unequal",
    Z = matrix(1, nrow = nrow(Pred4Mars_M_t), ncol = 1),
    A = "zero",
    R = "diagonal and equal")
  
  # Run MARS mod
  Pred4Mars_Mod <- MARSS(Pred4Mars_M_t,
                         model = MarssMod,
                         control = list(maxit = 5000),
                         method = "BFGS")
  
  Pred4Mars_Mod
}

# MARSS all predictors ----
## Unique sites ----
PredictorsList1 <- unique(PredSites$Response)

ListOfMarsModUnique <- list()

for(j in 1:length(PredictorsList1)){
  # j = 1
  Responsef_j <- PredictorsList1[j]
  
  # select month
  for(i in 1:length(unique(PredSites$Month))){
    # i = 1
    Monthf_i <- unique(PredSites$Month)[i]
    
    ModName_ij <- MARSfunSites(MarssDf = PredSites, Monthf = Monthf_i, Responsef = Responsef_j, 2)
    
    ListOfMarsModUnique[[length(ListOfMarsModUnique)+1]] <- ModName_ij
    
    names(ListOfMarsModUnique)[[length(ListOfMarsModUnique)]] <-  paste0("MARSS_._", Responsef_j,"_._", Monthf_i)
    
  }
  ListOfMarsModUnique
}

ListOfMarsModUniqueAICs <- as.data.frame(sapply(ListOfMarsModUnique, AIC))
ListOfMarsModUniqueAICs$Model <- row.names(ListOfMarsModUniqueAICs)
names(ListOfMarsModUniqueAICs) <- c("AIC", "Model")
ListOfMarsModUniqueAICs2 <- ListOfMarsModUniqueAICs %>% 
  separate(Model, sep = "_._", into = c("Mars", "Predictor", "Month"), remove = FALSE)

## Sites same ----

ListOfMarsModSame <- list()

for(j in 1:length(PredictorsList1)){
  # j = 1
  Responsef_j <- PredictorsList1[j]
  
  # select month
  for(i in 1:length(unique(PredSites$Month))){
    # i = 1
    Monthf_i <- unique(PredSites$Month)[i]
    
    ModName_ij <- MARSfunSitesSame(MarssDf = PredSites, Monthf = Monthf_i, Responsef = Responsef_j, 2)
    
    ListOfMarsModSame[[length(ListOfMarsModSame)+1]] <- ModName_ij
    
    names(ListOfMarsModSame)[[length(ListOfMarsModSame)]] <-  paste0("MARSS_._", Responsef_j,"_._", Monthf_i)
    
  }
  ListOfMarsModSame
}

ListOfMarsModSameAICs <- as.data.frame(sapply(ListOfMarsModSame, AIC))
ListOfMarsModSameAICs$Model <- row.names(ListOfMarsModSameAICs)
names(ListOfMarsModSameAICs) <- c("AIC", "Model")
ListOfMarsModSameAICs2 <- ListOfMarsModSameAICs %>% 
                          separate(Model, sep = "_._", into = c("Mars", "Predictor", "Month"), remove = FALSE)

## Combine and find best model
ListOfMarsModSameAllAICs <- ListOfMarsModSameAICs2 %>% 
                              left_join(ListOfMarsModUniqueAICs2, by = c("Mars", "Predictor","Month"), suffix = c("_Same", "_Unique")) %>% 
                              mutate(BestModel = ifelse(AIC_Same <= AIC_Unique,"Same", "Unique"))

ListOfSameStateModels2Keep <- ListOfMarsModSameAllAICs[ListOfMarsModSameAllAICs$BestModel == "Same",]$Model_Same
ListOfMars_BestSame <- ListOfMarsModSame[ListOfSameStateModels2Keep]

ListOfUniqueStateModels2Keep <- ListOfMarsModSameAllAICs[ListOfMarsModSameAllAICs$BestModel == "Unique",]$Model_Unique
ListOfMars_BestUnique <- ListOfMarsModUnique[ListOfUniqueStateModels2Keep]

ListOfMars_Best <- c(ListOfMars_BestSame, ListOfMars_BestUnique)

### Export AIC table ----
write.csv(ListOfMarsModSameAllAICs %>% 
            select(Predictor, Month, AIC_Same, AIC_Unique, BestModel), "07_Tables/11_MARS_AIC_Table.csv")

## bootstrap MARS models ----
library(doParallel)

num_cores <- detectCores() - 1

startTime <-Sys.time()
# do 500
bias_bootstrap_same <- mclapply(ListOfMars_BestSame, MARSSparamCIs, method = "parametric", nboot = 500, mc.cores = num_cores)
EndTime <- Sys.time()
EndTime - startTime

startTime <-Sys.time()
# do 500
bias_bootstrap_unique <- mclapply(ListOfMars_BestUnique, MARSSparamCIs, method = "parametric", nboot = 500, mc.cores = num_cores)
EndTime <- Sys.time()
EndTime - startTime


bias_bootstrap <- c(bias_bootstrap_same, bias_bootstrap_unique)

## Pull bias est ----
  ### Same ----
  MarsUdf_same <- matrix(nrow = length(bias_bootstrap_same), ncol = 4)
  
  # loop
  for(i in 1:length(bias_bootstrap_same)){
    # i = 1
    ModName <- names(bias_bootstrap_same)[[i]]
    U <- bias_bootstrap_same[[i]]$coef[[2]]
    lowCI <- bias_bootstrap_same[[i]]$par.lowCI$U[2,1]
    upCI <- bias_bootstrap_same[[i]]$par.upCI$U[2,1]
    
    MarsUdf_same[i,] <- c(ModName,U, lowCI, upCI)
    MarsUdf_same
  }
  
  MarsUdf_same.df <- as.data.frame(MarsUdf_same) %>% 
                      mutate(State = "Same",
                             Site = "All")
  names(MarsUdf_same.df) <- c("Model", "U", "U_lowCI", "U_upCI", "State", "Site")
  
  ### Unique ----
  MarsUdf_unique <- matrix(ncol = 4)
  
  # loop
  for(i in 1:length(bias_bootstrap_unique)){
    # i = 1
    ModName <- names(bias_bootstrap_unique)[[i]]
    # U <- bias_bootstrap_unique[[i]]$coef[[2]]
    U <- bias_bootstrap_unique[[i]]$coef[grepl("U", names(bias_bootstrap_unique[[i]]$coef)) == TRUE]
    lowCI <- bias_bootstrap_unique[[i]]$par.lowCI$U[1:length(U),1]
    upCI <- bias_bootstrap_unique[[i]]$par.upCI$U[1:length(U),1]
  
    MarsUdf_unique_i <- cbind(rep(ModName, length = length(U)),U, lowCI, upCI)
    MarsUdf_unique <- rbind(MarsUdf_unique, MarsUdf_unique_i)
  }
  
  
  MarsUdf_unique.df <- as.data.frame(MarsUdf_unique[-1,]) %>% 
                        mutate(State = "Unique",
                               Site = row.names(MarsUdf_unique[-1,]),
                               Site = gsub("U.X.","", Site))
  
  names(MarsUdf_unique.df) <- c("Model", "U", "U_lowCI", "U_upCI", "State", "Site")
  
  
  ### Process dat ----
  MarsUdf_Best <- rbind(MarsUdf_same.df, MarsUdf_unique.df) %>% 
    mutate(across(U:U_upCI, as.numeric)) %>% 
    separate(Model, sep = "_._", into = c("MARSS", "Predictor", "Month")) %>% 
    select(-MARSS) %>% 
    mutate(Sig = ifelse((U_lowCI < 0 & U_upCI < 0 )| (U_lowCI >= 0 & U_upCI >= 0),
                        "S",
                        "NS")) 

## Pull data & states ----
### Sites = same ----
  # make matrix
  MarsStates_same <- matrix(nrow = 1, ncol = 12)
  
  # loop
  for(i in 1:length(bias_bootstrap_same)){
    # i = 7
    ModName <- names(bias_bootstrap_same)[[i]]
    IsItBythFun4Dat <- function(ModName){
      BythTest <- grepl("Byth",ModName)
      DatResult_i <- if(BythTest == FALSE) {bias_bootstrap_same[[i]]$marss$data
      } else {rbind(
                                              rep("NA", length = dim(bias_bootstrap_same[[i]]$marss$data)[2]), #14-972
                                              rep("NA", length = dim(bias_bootstrap_same[[i]]$marss$data)[2]), #16-970
                                              bias_bootstrap_same[[i]]$marss$data[1,], #27-918
                                              rep("NA", length = dim(bias_bootstrap_same[[i]]$marss$data)[2]), # 29-905
                                              rep("NA", length = dim(bias_bootstrap_same[[i]]$marss$data)[2]), #3-996
                                              rep("NA", length = dim(bias_bootstrap_same[[i]]$marss$data)[2]), #36-873
                                              bias_bootstrap_same[[i]]$marss$data[2,], #37-890
                                              rep("NA", length = dim(bias_bootstrap_same[[i]]$marss$data)[2])) #8-994
      }
      
      DatResult_i
    }
    
    Dat_i <- IsItBythFun4Dat(ModName)
    
    States_i <- bias_bootstrap_same[[i]]$states
    States.se_i <- bias_bootstrap_same[[i]]$states.se
    
    MarsStates_i <- cbind(rep(ModName, dim(Dat_i)[2]), seq(1995,2023,by = 1), t(Dat_i), t(States_i), t(States.se_i))
    MarsStates_same <- rbind(MarsStates_same, MarsStates_i)
    MarsStates_same
  }


MarsStates_same.df <- as.data.frame(MarsStates_same[-1,]) 
names(MarsStates_same.df) <- c("Model", "Year", paste0("Dat_",names(MarsStates_same.df)[3:10]), "State_all", "State.se_all")
MarsStates_same.df2 <- MarsStates_same.df %>% 
  pivot_longer(cols = 'Dat_14-972':State.se_all, names_to = "DatSite", values_to = "Values") %>% 
  separate(DatSite, sep = "_", into = c("DatType","Site")) %>% 
  pivot_wider(id_cols = c(Model, Year, Site), names_from = DatType, values_from = Values) %>% 
  mutate(SameUnique = "Same") %>%
  separate(Model, sep = "_._", into = c("MARSS", "Predictor", "Month"), remove = FALSE) %>% 
  select(-MARSS)

### Sites = unique ----
  # make matrix
  MarsStates_unique <- matrix(nrow = 1, ncol = 26)
  
  # loop
  for(i in 1:length(bias_bootstrap_unique)){
    # i = 1
    ModName <- names(bias_bootstrap_unique)[[i]]
    IsItBythFun4_Dat <- function(ModName){
      BythTest <- grepl("Byth",ModName)
      DatResult_i <- if(BythTest == FALSE) {bias_bootstrap_unique[[i]]$marss$data
      } else {rbind(
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #14-972
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #16-970
        bias_bootstrap_unique[[i]]$marss$data[1,], #27-918
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), # 29-905
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #3-996
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #36-873
        bias_bootstrap_unique[[i]]$marss$data[2,], #37-890
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2])) #8-994
      }
      
      DatResult_i
    }
    
    Dat_i <- IsItBythFun4_Dat(ModName)
    
    IsItBythFun4_States <- function(ModName){
      BythTest <- grepl("Byth",ModName)
      DatResult_i <- if(BythTest == FALSE) {bias_bootstrap_unique[[i]]$states
      } else {rbind(
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #14-972
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #16-970
        bias_bootstrap_unique[[i]]$states[1,], #27-918
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), # 29-905
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #3-996
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #36-873
        bias_bootstrap_unique[[i]]$states[1,], #37-890
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2])) #8-994
      }
      
      DatResult_i
    }
    
    States_i <- IsItBythFun4_States(ModName)
    
    IsItBythFun4_States.se <- function(ModName){
      BythTest <- grepl("Byth",ModName)
      DatResult_i <- if(BythTest == FALSE) {bias_bootstrap_unique[[i]]$states
      } else {rbind(
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #14-972
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #16-970
        bias_bootstrap_unique[[i]]$states.se[1,], #27-918
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), # 29-905
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #3-996
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2]), #36-873
        bias_bootstrap_unique[[i]]$states.se[1,], #37-890
        rep("NA", length = dim(bias_bootstrap_unique[[i]]$marss$data)[2])) #8-994
      }
      
      DatResult_i
    }
    
    States.se_i <- IsItBythFun4_States.se(ModName)
    
    MarsStates_i <- cbind(rep(ModName, dim(Dat_i)[2]), seq(1995,2023,by = 1), t(Dat_i), t(States_i), t(States.se_i))
    MarsStates_unique <- rbind(MarsStates_unique, MarsStates_i)
    MarsStates_unique
  }
  
  MarsStates_unique.df <- as.data.frame(MarsStates_unique[-1,]) 
  names(MarsStates_unique.df) <- c("Model", "Year",
                                paste0("Dat_",names(MarsStates_unique.df)[3:10]), 
                                 paste0("State_", names(MarsStates_unique.df)[11:18]),
                                 paste0("State.se_",names(MarsStates_unique.df)[19:26]))
  names(MarsStates_unique.df) <- gsub(pattern = "_X.","_", names(MarsStates_unique.df))
  
  MarsStates_unique.df2 <- MarsStates_unique.df %>% 
    pivot_longer(cols = 'Dat_14-972':'State.se_8-994', names_to = "NameSite", values_to =  "Values") %>% 
    separate(NameSite, sep = "_", into = c("DatType", "Site")) %>% 
    pivot_wider(id_cols = c(Model, Year, Site), names_from = DatType, values_from = Values) %>% 
    mutate(SameUnique = "Unique") %>%
    separate(Model, sep = "_._", into = c("MARSS", "Predictor", "Month"), remove = FALSE) %>% 
    select(-MARSS)
    
# combine

  MarsStates <- rbind(MarsStates_unique.df2, MarsStates_same.df2) %>% 
    mutate(Year = as.numeric(Year),
           Dat = as.numeric(Dat),
           State = as.numeric(State),
           State.se = as.numeric(State.se),
           Site = ifelse(Site == "all", "All", Site)) %>% 
    filter(!(is.na(Dat) & is.na(State))) %>% 
    left_join(MarsUdf_Best %>%
                select(-State), by = c("Predictor", "Month", "Site")) %>%
    filter(!(Month == "Oct" | Year == "2023")) %>% 
    mutate(Predictor = as.factor(Predictor),
           Predictor = fct_relevel(Predictor, "AvgTemp", "EdiblePP", "PPLthen10", "InEdiblePP", 
                                   "MesoZPfood", "Byth", "Lepto", "plank_per_ha"),
           Site = as.factor(Site),
           Site = fct_relevel(Site,
                                     "36-873", "37-890", "27-918", "29-905",
                                     "14-972", "16-970", "3-996", "8-994", "All"),
           Month = as.factor(Month),
           Month = fct_relevel(Month, "May", "Jun", "Jul", "Aug", "Sep"))  %>% 
    mutate(State.Lse = State - State.se,
           State.Use = State + State.se)
  
  levels(MarsStates$Predictor) <- c("Temp",
                                    "PP[NGD]",
                                    "PP[L10]",
                                    "PP[GD]",
                                    "ZP[BM]",
                                    expression(italic("B. longimanus")),
                                    expression(italic("L. kindtii")),
                                    "Fish")
  
  MarsStates_dat <- MarsStates %>% 
      filter(!(is.na(Dat))) 
  
  MarsStates_Sig <- MarsStates %>% 
    filter(!is.na(Sig)) %>% 
    filter(State.Lse > -5)
  
  
  
# Plot ----
  SiteColors <- c("#d53e4f","#f46d43","#fdae61","#fee08b","#4292C6", "#9ECAE1","#762A83", "#C2A5CF", "black")
  
  
  png("05_Figures/11_Fig1_PredictorsInterannualPatternsPlot.png",
      units = "in", height = 11, width = 12, res = 300)
  ggplot() +
    geom_point(data = MarsStates_dat, aes(y = Dat, x = Year, color = Site), 
               size = 0.5, alpha = 50/100, shape = 3) + #
    geom_ribbon(data = MarsStates_Sig, aes(ymin = State.Lse, ymax = State.Use, x = Year, fill = Site),
                alpha = 20/100) + #
    geom_line(data = MarsStates_Sig, aes(y = State, x = Year, color = Site, linewidth = Sig)) + 
    facet_grid(Predictor ~ Month, 
               scales = "free_y", 
               labeller = label_parsed) +
    scale_color_manual(values = SiteColors,
                       labels = c("36", "37", "27", "29",
                                  "14", "16", "3", "8", "All"),
                       name = "Site") +
    scale_fill_manual(values = SiteColors,
                       labels = c("36", "37", "27", "29",
                                  "14", "16", "3", "8", "All"),
                       name = "Site") +
    scale_linewidth_manual(values = c(0.25, 1.1)) +
    ylim(-5,4)+
    ylab("Centered and logged response") +
    theme_bw() +
    theme(axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          # legend.position = "none",
          strip.text.y.right = element_text(face = "bold", size = 13),
          strip.text.x.top =  element_text(face = "bold", size = 13),
          strip.background = element_rect(color = "transparent"),
          panel.background = element_rect(fill = "white", color = "black")) 
  dev.off()
  
  MarsUdf_Best2 <- MarsUdf_Best %>% 
    mutate(Predictor = as.factor(Predictor),
           Predictor = fct_relevel(Predictor, "AvgTemp", "EdiblePP", "PPLthen10", "InEdiblePP", 
                                   "MesoZPfood", "Byth", "Lepto", "plank_per_ha"),
           Site = as.factor(Site),
           Site = fct_relevel(Site,
                              "36-873", "37-890", "27-918", "29-905",
                              "14-972", "16-970", "3-996", "8-994", "All"),
           Site = fct_recode(Site,
                              "36" = "36-873", "37" = "37-890", "27" = "27-918", "29" = "29-905",
                              "14" = "14-972", "16" = "16-970", "3" = "3-996", "8" = "8-994"),
           Month = as.factor(Month),
           Month = fct_relevel(Month, "May", "Jun", "Jul", "Aug", "Sep")) %>% 
    filter(Month != "Oct")
    
  levels(MarsUdf_Best2$Predictor) <- c("Temp",
                                    "PP[NGD]",
                                    "PP[L10]",
                                    "PP[GD]",
                                    "ZP[BM]",
                                    expression(italic("B. longimanus")),
                                    expression(italic("L. kindtii")),
                                    "Fish")
  png("05_Figures/11b_FigS4_MARSS_bias_Plot.png",
      units = "in", height = 11.5, width = 14, res = 300)
  ggplot() +
    geom_hline(color = "black", yintercept = 0) +
    geom_pointrange(data = MarsUdf_Best2, aes(y = U, ymin = U_lowCI, ymax = U_upCI, x = Site), 
                    shape =21 , fill = "grey") +
    geom_text(data = MarsUdf_Best2 %>%
                mutate(Sig2 = ifelse(Sig == "S", "*", "")), aes(y = U_upCI, x = Site, label = Sig2),
              size = 8) +
    facet_grid(Predictor ~ Month, 
                scales = "free_y", 
               labeller = label_parsed)+ 

    theme_bw() +
    ylab(expression(atop(paste(italic("u"), " (", y^-1,")"),"(± 95% condfidence interval)"))) +
    theme(axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = 16),
          # panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          # legend.position = "none",
          strip.text.y.right = element_text(face = "bold", size = 13),
          strip.text.x.top =  element_text(face = "bold", size = 13),
          strip.background = element_rect(color = "transparent"),
          panel.background = element_rect(fill = "white", color = "black")) 
dev.off()




# save.image(here::here("04_SavedRImages", "11_PredictorsInterannualPatterns.Rdata"))
# load(here::here("04_SavedRImages", "11_PredictorsInterannualPatterns.Rdata"))
