# Code for extracting and summarizing plankton and temp data
# heatmaps used to better understand correlations among months to match with phenology data
# JMH; Jan 23, updated Sept 23
# updated with fixed Byth on 2 Dec 23

# Libraries ----
library(RSQLite) # automatically installed with dbplyr
library(tidyverse)
library(here)
library(gridExtra)
library(ggcorrplot)

# Load database  ----
  Database <- file.path("01_OriginalData/00_LEPAS_20231202.db")
  lepas <<- dbConnect(SQLite(), dbname=Database)

# Get ZP and PP data ----
  dbListTables(lepas)   
  PPinfo <- dbReadTable(lepas, "Phyto_sample_info") %>% 
              select(Sample_ID:Project)
  PPlg <- dbReadTable(lepas, "Phyto_lengths") %>% 
              select(-ID)
  PPfin <- dbReadTable(lepas, "Phyto_final")
  PPtaxa <- dbReadTable(lepas, "Phyto_taxa_groups") %>% 
              select(Genus_type, Edible, Genus, Type) 
  ZPinfo <- dbReadTable(lepas, "Zoop_sample_info") %>% 
              select(Sample_ID:Project) %>% 
              distinct()
  ZPfin <- dbReadTable(lepas, "Zoop_final") %>% 
              select(Sample_ID:Zoop_biomass)
  ZPtaxa <- dbReadTable(lepas, "Zoop_taxa_groups") %>% 
              select(Genus_sp_lifestage, Order1, Genus) %>% 
              distinct()

  ## Small PP for Vel ----
  PPlg_sum <- PPlg %>% 
    left_join(PPinfo, by = "Sample_ID") %>% 
    left_join(PPtaxa, by = "Genus_type") %>%
    # Just solitary
    filter(Type == "Solitary") %>% 
    # just lepas sites
    filter(Sample_site %in% c("27-918", "37-890", "36-873", "29-905", "16-970", "14-972", "8-994", "3-996")) %>% 
    group_by(Genus_type) %>% 
    summarise(across(c(Length_um:Diameter_cell_um, Diameter_cell_um), \(x) mean(x, na.rm = T))) %>% 
    filter(Length_um < 10 | Diameter_cell_um < 10)
  
  PPTaxaL10 <- PPlg_sum$Genus_type
  
  
  ## Chl a & temp ----
  # see notes at bottom about chla
  Chla <- dbGetQuery(lepas, 
                     "SELECT  
                    PG.Sample_ID,
                    PG.Sample_site,
                    PG.Sample_date,
                    CAST(PG.Chl_a_ugL_L AS text),
                    CAST(PG.Chl_a_ugL_JH AS text)
                    FROM Pigments PG") %>% 
          # both of these columns have an NA that SQLite turns into a ZERO
          # should just use the JH - Jeffrey & Humphrey 1975
          rename(Chl_a_ugL_L2 = "CAST(PG.Chl_a_ugL_L AS text)",
                 Chl_a_ugL_JH2 = "CAST(PG.Chl_a_ugL_JH AS text)") %>% 
          mutate(Chl_a_ugL_L2 = as.numeric(Chl_a_ugL_L2),
                 Chl_a_ugL_JH2 = as.numeric(Chl_a_ugL_JH2)) 
  
  Temp <- dbReadTable(lepas, "Sample_inventory") %>% 
    select(Sample_ID,Sample_site, Sample_date, Surface_temp:Ave_temp)

# Plankton data
  ## ZP data ----
  ZPfin2 <- ZPfin %>% 
                left_join(ZPtaxa, by = "Genus_sp_lifestage") %>% 
                left_join(ZPinfo, by = "Sample_ID") %>% 
                # just crustacean zp
                filter(Order1 %in% c("Calanoida", "Cyclopoida", "Cladocera")) %>% #, "Harpacticoida", "Poecilostomatoida", "Arguloida"
                # just lepas sites
                filter(Sample_site %in% c("27-918", "37-890", "36-873", "29-905", "16-970", "14-972", "8-994", "3-996")) %>% 
                  # fix stupid naming issue
                mutate(Water_body = ifelse(Water_body %in% c("Central Erie", "Cantral Erie", "Central Basin", "Central Basin", "Cental Basin"), "Central_Basin",
                                             ifelse(Water_body %in% c("Western Erie", "Westerm Erie"), "Western_Basin",Water_body))) %>% 
                mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d"),
                       Y = strftime(Sample_date, format = "%Y"),
                       M = strftime(Sample_date, format = "%b")) %>% 
                # remove "B" samples
                filter(!str_detect(Sample_ID, "B"))
  
  ### Leptodora density ----
  # needed for invert predator data
  Lepto <- ZPfin2 %>% 
              filter(Genus_sp_lifestage %in% c("L_kindti_w_eggs", "L_kindti_wo_eggs")) %>% 
              # remove eggs
              filter(Life_stage != "Egg") %>% 
              # sum to smp ID
              group_by(Sample_ID, Genus, Sample_date, Sample_site, Y, M) %>% 
              summarise(Zoop_biomass = sum(Zoop_biomass, na.rm = T)) 

  
  ### Byth ---- 
  # taking this from DB, original version used spreadsheet
  # note that this has all sample sites, need to restrict to WB FTG sites
  Byth <- ZPfin2 %>% 
      # Quality indicates that the whole smp was counted
      # careful there is also a B_longimanus_wo_eggs and ...w_eggs from subsampled data
      filter(Genus_sp_lifestage %in% c("B_longimanus_w_eggs_Quality", "B_longimanus_wo_eggs_Quality")) %>% 
      # remove eggs
      filter(Life_stage != "Egg") %>% 
      # sum to smp ID
      group_by(Sample_ID, Genus, Sample_date, Sample_site, Y, M) %>% 
      summarise(Zoop_biomass = sum(Zoop_biomass, na.rm = T)) 
    
    

      
    # focus on FTG sites with good time series
    Byth2 <- Byth %>% 
      filter(Sample_site %in% c("27-918", "37-890"))
    
    ggplot(Byth2 %>% 
           filter(Y < 2023) %>% 
            mutate(Yg = case_when(Y >= 1995 & Y < 2000 ~ "1995_1999",
                                   Y >= 2000 & Y < 2005 ~ "2000_2004",
                                   Y >= 2005 & Y < 2010 ~ "2005_2009",
                                   Y >= 2010 & Y < 2015 ~ "2010_2014",
                                   Y >= 2015 & Y < 2020 ~ "2015_2020",
                                   Y >= 2020 ~ "2020_2022")) %>% 
             filter(M %in% c("May", "Jun", "Jul", "Aug", "Sep")) %>% 
             mutate(M = as.factor(M),
                    M = fct_relevel(M, "May", "Jun", "Jul", "Aug", "Sep")), 
           aes(y = log10(Zoop_biomass + 0.00312), x = M)) +
      geom_boxplot() +
      facet_wrap(~Yg)
    
### Meso ZP food ----
    MesoZPFood <- ZPfin %>% 
      left_join(ZPtaxa, by = "Genus_sp_lifestage") %>% 
      left_join(ZPinfo, by = "Sample_ID") %>% 
      # just lepas sites
      filter(Sample_site %in% c("27-918", "37-890", "36-873", "29-905", "16-970", "14-972", "8-994", "3-996")) %>% 
      # fix stupid naming issue
      mutate(Water_body = ifelse(Water_body %in% c("Central Erie", "Cantral Erie", "Central Basin", "Central Basin", "Cental Basin"), "Central_Basin",
                                 ifelse(Water_body %in% c("Western Erie", "Westerm Erie"), "Western_Basin",Water_body))) %>% 
      mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d"),
             Y = strftime(Sample_date, format = "%Y"),
             M = strftime(Sample_date, format = "%b")) %>% 
      # remove "B" samples
      filter(!str_detect(Sample_ID, "B")) %>% 
      # want rotifers, veligers, and copepod immatures,excuding parasitic and benthic taxa
      # Based on Williamson 1986
      mutate(TaxaGroups = case_when(Order1 %in% c("Collothecaceae", "Flosculariaceae", "Rotifera spp.", "Ploima") ~ "Rotifer",
                                    Order1 %in% c("Veneroida") ~ "Veligers",
                                    Order1 %in% c("Cyclopoida", "Calanoida", "Copepoda spp.") ~ "Copepods",
                                    Order1 %in% c("Cladocera") ~ "Cladocera",
                                    Order1 %in% c("Harpacticoida", "Arguloida", "Poecilostomatoida") ~ "DROP"),
             TaxaGroups = as.factor(TaxaGroups),
             # keep younger copepods, many copepidites plus removes year-year bias
             Groups2Keep = case_when(TaxaGroups == "Copepods" & Life_stage == "Younger" ~ "Keep",
                                     # Keep immature Calanoids and cyclopoids
                                     Genus_sp_lifestage %in% c("Calanoid_nauplii", "Diaptomidae_immature", "Calanoida_mixed",
                                                               "Cyclopoid_nauplii","Cyclopidae_immature", "Cyclopoida_mixed", "Temoridae_immature",
                                                               "Dreissena_veligers") ~ "Keep",
                                     # Keep rotifers and veligers
                                     TaxaGroups %in% c("Rotifer", "Veligers") ~ "Keep",
                                     # Keep all cladocerans but Byth and Lepto
                                     # Meso eats Lepto, but keeping it seperate for consistancy
                                     !(Genus %in% c("Bythotrephes", "Cercopagis", "Leptodora")) ~ "Keep",
                                     .default = "DROP"),
             Groups2Keep = as.factor(Groups2Keep)) %>% 
      filter(Groups2Keep == "Keep" & TaxaGroups != "Drop") %>% 
      # sum to smp ID
      group_by(Sample_ID, Sample_date, Sample_site, Y, M) %>% 
      summarise(MesoZPfood = sum(Zoop_biomass, na.rm = T)) 
    
## PP data ----
    PPfin2 <-   PPfin %>% 
      left_join(PPtaxa, by = "Genus_type") %>% 
      left_join(PPinfo, by = "Sample_ID") %>% 
      # just lepas sites
      filter(Sample_site %in% c("27-918", "37-890", "36-873", "29-905", "16-970", "14-972", "8-994", "3-996")) %>% 
      # remove microcystis values I do not want
      filter(!(Genus_type %in% c("Microcystis_sm_celled_by_cell", "Microcystis_lg_celled_by_col"))) %>% 
      # fix stupid naming issue
      mutate(Water_body = ifelse(Water_body == "Central Erie", "Central_Basin",
                                 ifelse(Water_body == "Western Erie", "Western_Basin",Water_body))) %>% 
      mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d"),
             Y = strftime(Sample_date, format = "%Y"),
             M = strftime(Sample_date, format = "%b"))%>% 
            # remove "B" samples
            filter(!str_detect(Sample_ID, "B"))

  ### Edible algae biomass ----
    
    EdibleInedible <- PPfin2 %>% 
      # sum for each sample ID
      group_by(Sample_ID, 
               # just retaining these
               Sample_date, Sample_site, Y, M, Edible) %>% 
      # Sum Edible/Inedible taxa, convert to µg
      summarise(Phyto_biomass = sum(Phyto_biomass, na.rm = T)*1000, .groups = "drop") 
  
  EdiblePP <- EdibleInedible %>% 
              filter(Edible == "TRUE")%>% 
    select(-Edible)
    
  
  
  ### Inedible algae biomass ----
  
  InEdiblePP <- EdibleInedible %>% 
    filter(Edible == "FALSE") %>% 
    select(-Edible)
  
  #### Dominant PPGD taxa ----
  InEdibleTaxaSum <- PPfin2 %>% 
    filter(Edible == "FALSE") %>% 
    filter(Y < 2023) %>% 
    mutate(Yg = case_when(Y >= 1995 & Y < 2000 ~ "1995_1999",
                          Y >= 2000 & Y < 2005 ~ "2000_2004",
                          Y >= 2005 & Y < 2010 ~ "2005_2009",
                          Y >= 2010 & Y < 2015 ~ "2010_2014",
                          Y >= 2015 & Y < 2020 ~ "2015_2020",
                          Y >= 2020 ~ "2020_2022")) %>% 
    group_by(Yg, M, Genus) %>% 
    summarise(Phyto_biomass = mean(Phyto_biomass)) %>% 
    filter(M %in% c("May", "Jun", "Jul", "Aug", "Sep")) %>% 
    mutate(M = as.factor(M),
           M = fct_relevel(M, "May", "Jun", "Jul", "Aug", "Sep"),
           Genus = as.factor(Genus),
           Genus = fct_relevel(Genus,
                               "Aphanizomenon", "Aphanocapsa", "Aphanothece", "Chroococcus", "Cylindrospermopsis", "Dolichospermum", 
                               "Merismopedia", "Micractinium", "Microcystis", "Planktothrix", "Pseudanabaena",
                               "Aulacoseira", "Ceratium","Closterium", "Lyngbya", "Spirulina", "")) %>% 
    group_by(Yg, M) %>% 
    mutate(Phyto_biomassSum = sum(Phyto_biomass),
           Rel_PhytoBM = Phyto_biomass/Phyto_biomassSum)
  
 
  
  CommonPPGDTaxa <-  InEdibleTaxaSum %>% 
    group_by(Genus) %>% 
    summarise(max = max(Rel_PhytoBM)) %>% 
    filter(max > 0.1) %>% 
    distinct(Genus) %>% 
    filter(Genus != "")
  
  
  ggplot(InEdibleTaxaSum %>% 
           filter(Genus %in% CommonPPGDTaxa$Genus), 
         aes(y = Phyto_biomass, x = M, fill = Genus)) +
    geom_bar(position = "fill", stat = "identity") +
    facet_wrap(~Yg) 
  
# Microcystis v. InEdible Algae ----
  # need to come back to this and summarise/trim the same as data in models
  PPfin22 <- PPfin2 %>% 
    mutate(Cyanobacteria = ifelse(Genus %in% c("Aphanizomenon", "Aphanocapsa", "Aphanothece", "Chroococcus", "Cylindrospermopsis", "Dolichospermum", 
                                               "Merismopedia", "Micractinium", "Microcystis", "Planktothrix", "Pseudanabaena"),
                                  "Cyanobacteria", "NotCyanobacteria")) %>% 
    filter(Cyanobacteria == "Cyanobacteria") %>% 
    group_by(Sample_ID, Sample_date, Sample_site, Genus) %>% 
    summarise(Phyto_biomass = sum(Phyto_biomass, na.rm = T)) %>% 
    pivot_wider(id_cols = c(Sample_ID, Sample_date, Sample_site), names_from = Genus, values_from = Phyto_biomass, values_fill = 0) %>% 
    left_join(InEdiblePP %>% 
                select(Sample_ID, InEdiblePP = Phyto_biomass), by = "Sample_ID") %>% 
    mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d"),
           Year = as.numeric(strftime(Sample_date, format = "%Y")),
           M = as.character(strftime(Sample_date, format = "%m")),
           BiWOY = floor(as.numeric(strftime(Sample_date, format = "%U"))/2)) %>% 
    filter(BiWOY >= 9 & BiWOY <= 20)  %>% 
    # rowwise() %>%
    mutate(Cyanobacteria = (Aphanizomenon + Aphanocapsa + Aphanothece + Chroococcus + Cylindrospermopsis + Dolichospermum + 
                              Merismopedia + Micractinium + Microcystis + Planktothrix + Pseudanabaena)*1000,
           Microcystis = Microcystis*1000) %>% 
    select(-c("Aphanizomenon", "Aphanocapsa", "Aphanothece", "Chroococcus", "Cylindrospermopsis", "Dolichospermum", 
              "Merismopedia", "Micractinium", "Planktothrix", "Pseudanabaena")) %>% 
    mutate(PerMicrocystis = Microcystis/InEdiblePP *100,
           PerCyanos = Cyanobacteria/InEdiblePP * 100) 
  
  sem <- function(x, na.rm = FALSE) {
    out <-sd(x, na.rm = na.rm)/sqrt(length(x))
    return(out)
  }
  
  PPfin22_mean <- PPfin22 %>% 
    group_by(Year, M) %>% 
    summarise(across(c(Microcystis, InEdiblePP, Cyanobacteria, PerMicrocystis, PerCyanos), \(x) mean(x, na.rm = T), .names = "{.col}_mean"))
  
  PPfin22_se <- PPfin22 %>% 
    group_by(Year, M) %>% 
    summarise(across(c(Microcystis, InEdiblePP, Cyanobacteria, PerMicrocystis, PerCyanos), \(x) sem(x, na.rm = T), .names = "{.col}_se"))
  
  PPfin22_sum <- PPfin22_mean %>% 
                  left_join(PPfin22_se, by = c("Year", "M")) %>% 
                  mutate(M2 = case_when(M == "05" ~ "May",
                                        M == "06" ~ "Jun",
                                        M == "07" ~ "Jul",
                                        M == "08" ~ "Aug",
                                        M == "09" ~ "Sep",
                                        M == "10" ~ "Oct"),
                         M2 = as.factor(M2),
                         M2 = fct_relevel(M2, "May", "Jun", "Jul", "Aug", "Sep", "Oct")) %>% 
                  filter(M2 != "Oct")
  
  ## Fig. S6 ----
  P_PerMicro <- ggplot(PPfin22_sum, aes(y = PerMicrocystis_mean, x = Year, color = M2)) +
    geom_point(size = 4) +
    geom_line(linewidth = 0.25) +
    scale_color_manual(values = c("#084594", "#9ECAE1", "#FD8D3C", "#FC4E2A","#B10026"), # "#E31A1C", 
                       name = "Month") +
    # scale_color_brewer(palette = "YlOrRd") +
    ylab(expression(paste("% Microcystis in ", PP[GD]))) +
    ylim(0,100) +
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
  
  options(scipen = 999)
  p_PPGDbiomass <- ggplot(PPfin22_sum, aes(y = InEdiblePP_mean, x = Year, color = M2)) +
    geom_point(size = 4) +
    geom_line(linewidth = 0.25) +
    # removes 2016 sep value on high side and 1996 Oct value on low sied
    scale_y_log10(labels = scales::comma, limits = c(50, 1000000)) +
    scale_color_manual(values = c("#084594", "#9ECAE1", "#FD8D3C", "#FC4E2A", "#B10026"), #"#E31A1C", 
                       name = "Month") +
    ylab(expression(paste(PP[GD], " biomass (µg ",L^-1,")"))) +
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
    
  
  P_PerCyano <- ggplot(PPfin22_sum, aes(y = PerCyanos_mean, x = Year, color = M2)) +
    geom_point(size = 4) +
    geom_line(linewidth = 0.25)  +
    scale_color_manual(values = c("#084594", "#9ECAE1", "#FD8D3C", "#FC4E2A",  "#B10026"), #"#E31A1C",
                       name = "Month") +
    # scale_color_brewer(palette = "YlOrRd") +
    ylab(expression(paste("% Cyanobacteria in ", PP[GD]))) +
    ylim(0,100) +
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
  
  library(cowplot)
  png("05_Figures/06b_FigS6_MicroCyanosVInedible.png",
      units = "in", height = 16, width = 12, res = 300)
    plot_grid(p_PPGDbiomass, P_PerMicro, P_PerCyano, labels= c("A", "B", "C"), 
              ncol = 1, label_size = 24)
  dev.off()
  
# PP L 10 µm ----
  PPL10 <- PPfin2 %>% 
    mutate(PPSizeGroup = ifelse(Genus_type %in% PPTaxaL10, "PPL10","PPG10")) %>% 
    # sum for each sample ID
    group_by(Sample_ID, 
             # just retaining these
             Sample_date, Sample_site, Y, M, PPSizeGroup) %>% 
    # Sum Edible/Inedible taxa, convert to µg
    summarise(Phyto_biomass = sum(Phyto_biomass, na.rm = T)*1000, .groups = "drop") %>% 
    filter(PPSizeGroup == "PPL10") %>% 
    select(-PPSizeGroup)
  
  
   ## Chla ----
    # 5 of these are zero
    ChlaA <- Chla %>% 
              filter(Sample_site %in% c("27-918", "37-890", "36-873", "29-905", "16-970", "14-972", "8-994", "3-996")) %>% 
              mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d"),
                     M = strftime(Sample_date, format = "%b"),
                     Y = strftime(Sample_date, format = "%Y")) %>% 
              # remove "B" samples
              filter(!str_detect(Sample_ID, "B"))
      
  
  ## Temperature ----
  TempA <- Temp %>% 
              filter(Sample_site %in% c("27-918", "37-890", "36-873", "29-905", "16-970", "14-972", "8-994", "3-996")) %>% 
              mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d"),
                     M = strftime(Sample_date, format = "%b"),
                     Y = strftime(Sample_date, format = "%Y")) %>% 
              # some of these are zero, which is clearly an error
              filter(Ave_temp != 0) %>% 
              # remove "B" samples
              filter(!str_detect(Sample_ID, "B"))
      
  TempA[duplicated(TempA[,-1]) == TRUE,]
    
# Bring together df ----
  
  PlankDat <- rbind(Lepto %>% 
                      ungroup() %>% 
                      mutate(Response = "Lepto") %>% 
                      select(Sample_ID, Sample_date, Sample_site, Response, Value = Zoop_biomass),
                    EdiblePP %>% 
                      ungroup() %>% 
                      mutate(Response = "EdiblePP") %>% 
                      select(Sample_ID, Sample_date, Sample_site, Response, Value = Phyto_biomass),
                    InEdiblePP %>% 
                      ungroup() %>% 
                      mutate(Response = "InEdiblePP")%>% 
                      select(Sample_ID, Sample_date, Sample_site, Response, Value = Phyto_biomass),
                    PPL10%>% 
                      ungroup() %>% 
                      mutate(Response = "PPLthen10")%>% 
                      select(Sample_ID, Sample_date, Sample_site, Response, Value = Phyto_biomass),
                    ChlaA %>% 
                      mutate(Response = "Chla") %>% 
                      select(Sample_ID, Sample_date, Sample_site, Response, Value = Chl_a_ugL_JH2),
                    TempA %>% 
                      mutate(Response = "AvgTemp") %>% 
                      select(Sample_ID, Sample_date, Sample_site, Response, Value = Ave_temp),
                    Byth2 %>% 
                      ungroup() %>% 
                      mutate(Response = "Byth") %>% 
                      select(Sample_ID, Sample_date, Sample_site, Response, Value = Zoop_biomass),
                    MesoZPFood %>% 
                      ungroup() %>% 
                      mutate(Response = "MesoZPfood") %>% 
                      select(Sample_ID, Sample_date, Sample_site, Response, Value = MesoZPfood))

  
  
# save data
  write.csv(PlankDat, "03_GeneratedData/06b_PlanktonTempDat_BiWk.csv")
      
# save image
    # save.image("04_SavedRImages/06b_PredMunge_Plankton_BiWk.RData")
    # load("04_SavedRImages/06b_PredMunge_Plankton_BiWk.RData")
    
              
    
    
    
    
    
    
    
    
    
   
    