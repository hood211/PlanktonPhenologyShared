# Code for Rel biomass table
# JMH; May 23, Dec 23

# Libraries ----
library(RSQLite) # automatically installed with dbplyr
library(tidyverse)
library(here)


# Load database  ----
Database <- paste0("01_OriginalData/00_LEPAS_20231202.db")
lepas <<- dbConnect(SQLite(), dbname=Database)

## ZP data ---
ZPinfo <- dbReadTable(lepas, "Zoop_sample_info") %>% 
  select(Sample_ID:Project) %>% 
  # one dup row rel to sites
  distinct()

ZPfin <- dbReadTable(lepas, "Zoop_final") %>% 
  select(Sample_ID:Zoop_biomass)

ZPtaxa <- dbReadTable(lepas, "Zoop_taxa_groups") %>% 
  distinct(Genus_sp_lifestage, Order1, Genus, Species) %>% 
  select(Genus_sp_lifestage, Order1, Genus, Species)  

## ZP data combined ----
ZPfin2 <- ZPfin %>% 
  left_join(ZPtaxa, by = "Genus_sp_lifestage") %>% 
  left_join(ZPinfo, by = "Sample_ID") %>% 
  # No rotifers, no parasites
  filter(Order1 %in% c("Calanoida", "Cyclopoida", "Cladocera", "Copepoda spp.", "Veneroida")) %>%
  # just lepas sites
  filter(Sample_site %in% c("27-918", "37-890", "36-873", "29-905", "16-970", "14-972", "8-994", "3-996")) %>% 
  # fix stupid naming issue
  mutate(Water_body = ifelse(Water_body %in% c("Central Erie", "Cantral Erie", "Central Basin", "Central Basin", "Cental Basin"), "Central_Basin",
                             ifelse(Water_body %in% c("Western Erie", "Westerm Erie"), "Western_Basin",Water_body))) %>% 
  mutate(Sample_date = as.POSIXct(Sample_date, format = "%Y-%m-%d"),
         Y = strftime(Sample_date, format = "%Y"),
         M = strftime(Sample_date, format = "%b")) %>% 
  # removing duplicates 
  filter(!duplicated(.)) %>% 
  # remove "B" samples
  filter(!str_detect(Sample_ID, "B"))

## select focal taxa ----
ZPfin2_FocalTax <- ZPfin2 %>% 
            filter(grepl("D_retrocurva", Genus_sp_lifestage) | grepl("S_oregonensis", Genus_sp_lifestage) |
                     grepl("M_edax", Genus_sp_lifestage) | grepl("Dreissena", Genus_sp_lifestage) ) %>% 
            # remove "younger" taxa
            filter(Life_stage != "Younger") %>% 
            # sum across life stages
            group_by(Genus, Species, Sample_site, Y, M, Sample_ID) %>% 
            summarise(Zoop_biomass = sum(Zoop_biomass, na.rm = T)) %>% 
            # average for month
            group_by(Genus, Species, Sample_site, Y, M) %>% 
            summarise(Zoop_biomass = mean(Zoop_biomass, na.rm = T)) %>% 
            mutate(Genus_sp = paste0(Genus,"_", Species),
                   Genus_sp = ifelse(Genus_sp == "Dreissena_Dreissena spp.", "Veligers" ,Genus_sp)) %>% 
            pivot_wider(id_cols = Sample_site:M, names_from = Genus_sp, values_from = Zoop_biomass, values_fill = 0)

ZPfin2_TotalZP <- ZPfin2 %>% 
            group_by(Sample_site, Y, M) %>% 
            summarise(TotZP_biomass = sum(Zoop_biomass, na.rm = T))

ZPfin2_TotCladZP <- ZPfin2 %>% 
  filter(Order1 %in% c("Cladocera", "Cyclopoida", "Calanoida", "Copepoda spp.")) %>% 
  group_by(Sample_site, Y, M) %>% 
  summarise(TotCrus_biomass = sum(Zoop_biomass, na.rm = T))

ZPfinPlot_TOT <- ZPfin2_FocalTax %>% 
              left_join(ZPfin2_TotalZP, by = c("Sample_site", "Y", "M")) %>% 
              left_join(ZPfin2_TotCladZP, by = c("Sample_site", "Y", "M")) %>% 
              mutate(TotZP_BM_NF = TotZP_biomass - Daphnia_retrocurva - Veligers - Mesocyclops_edax - Skistodiaptomus_oregonensis) %>% 
              select(-TotZP_biomass, -TotCrus_biomass) %>% 
              pivot_longer(cols = Daphnia_retrocurva:TotZP_BM_NF, names_to = "Taxa", values_to = "Biomass")%>% 
              filter(M %in% c("May", "Jun", "Jul", "Aug", "Sep")) %>%
              ungroup() %>% 
              mutate(M = fct_relevel(as.factor(M), "May", "Jun", "Jul", "Aug", "Sep"),
                     Sample_site = fct_relevel(as.factor(Sample_site), 
                                               "36-873", "37-890", "14-972", "3-996",
                                               "27-918", "29-905", "16-970", "8-994"),
                     Taxa = fct_relevel(as.factor(Taxa), 
                                        "Daphnia_retrocurva", "Veligers", "Mesocyclops_edax", "Skistodiaptomus_oregonensis", "TotZP_BM_NF"),
                     Yg = ifelse(Y >= 1995 & Y <= 1999, "1995-1999",
                            ifelse(Y >= 2000 & Y <= 2004, "2000-2004",
                                ifelse(Y >= 2005 & Y <= 2009, "2005-2009",
                                      ifelse(Y >= 2010 & Y <= 2014, "2010-2014",
                                            ifelse(Y >= 2015, "2015-2020", "blah")))))) %>% 
              group_by(Taxa, Sample_site, Yg, M) %>% 
              summarise(Biomass = mean(Biomass, na.rm = T))


ZPfin2_b <- ZPfin2 %>% 
              group_by(Genus, Species, Sample_site, Y, M, Sample_ID) %>% 
              summarise(Zoop_biomass = sum(Zoop_biomass, na.rm = T)) %>% 
              group_by(Genus, Species, Sample_site, Y, M) %>% 
              summarise(Zoop_biomass = mean(Zoop_biomass, na.rm = T)) %>% 
              mutate(Genus_sp = paste0(Genus,"_", Species),
                     Genus_sp = ifelse(Genus_sp == "Dreissena_Dreissena spp.", "Veligers" ,Genus_sp)) %>%
              filter(M %in% c("May", "Jun", "Jul", "Aug", "Sep")) %>%
              ungroup() %>% 
              mutate(M = fct_relevel(as.factor(M), "May", "Jun", "Jul", "Aug", "Sep"),
                     Sample_site = fct_relevel(as.factor(Sample_site), 
                                               "36-873", "37-890", "14-972", "3-996",
                                               "27-918", "29-905", "16-970", "8-994"),
                     Yg = ifelse(Y >= 1995 & Y <= 1999, "1995-1999",
                                 ifelse(Y >= 2000 & Y <= 2004, "2000-2004",
                                        ifelse(Y >= 2005 & Y <= 2009, "2005-2009",
                                               ifelse(Y >= 2010 & Y <= 2014, "2010-2014",
                                                      ifelse(Y >= 2015, "2015-2020", "blah"))))),  
                    Yg2 = ifelse(Y >= 1995 & Y <= 2004, "1995-2004",
                            ifelse(Y >= 2005 & Y <= 2014, "2005-2014",
                                ifelse(Y >= 2015 & Y <= 2020, "2015-2020", "blah")))) %>% 
              group_by(Genus_sp, Sample_site, Yg2, M) %>% 
              summarise(Zoop_biomass = mean(Zoop_biomass, na.rm = T)) %>% 
              # pivot_wider(id_cols = Sample_site:M, names_from = "Genus_sp", values_from = "Zoop_biomass") %>% 
              # mutate(ID = paste0(Sample_site, "_b_", Yg, "_b_", M)) %>% 
              # ungroup() %>% 
              # select(-Sample_site, -Yg, -M)
              group_by(Sample_site, Yg2, M) %>% 
              mutate(Rank = rank(desc(Zoop_biomass)))  %>% 
              mutate(Relbiomass = Zoop_biomass/max(Zoop_biomass),
                     Genus_sp2 = ifelse(Genus_sp %in% c("Daphnia_retrocurva", "Veligers", "Mesocyclops_edax", "Skistodiaptomus_oregonensis"),
                                        Genus_sp,
                                        "Other"),
                     SizeCode = ifelse(Genus_sp2 == "Other", "S", "B"),
                     SizeCode = fct_relevel(SizeCode, "B", "S"))

# all years averaged
ZPfin2_b_AllYrAvg <- ZPfin2 %>% 
  group_by(Genus, Species, Sample_site, Y, M, Sample_ID) %>% 
  summarise(Zoop_biomass = sum(Zoop_biomass, na.rm = T)) %>% 
  group_by(Genus, Species, Sample_site, Y, M) %>% 
  summarise(Zoop_biomass = mean(Zoop_biomass, na.rm = T)) %>% 
  mutate(Genus_sp = paste0(Genus,"_", Species),
         Genus_sp = ifelse(Genus_sp == "Dreissena_Dreissena spp.", "Veligers" ,Genus_sp)) %>%
  filter(M %in% c("May", "Jun", "Jul", "Aug", "Sep")) %>%
  ungroup() %>% 
  group_by(Genus_sp, Sample_site, Y, M) %>% 
  summarise(Zoop_biomass = mean(Zoop_biomass, na.rm = T)) %>% 
  group_by(Genus_sp, Sample_site, M) %>% 
  summarise(Zoop_biomass = mean(Zoop_biomass, na.rm = T)) %>% 
  group_by(Sample_site, M) %>% 
  mutate(Rank = rank(desc(Zoop_biomass)))  %>% 
  mutate(Relbiomass = Zoop_biomass/max(Zoop_biomass),
         Genus_sp2 = ifelse(Genus_sp %in% c("Daphnia_retrocurva", "Veligers", "Mesocyclops_edax", "Skistodiaptomus_oregonensis"),
                            Genus_sp,
                            "Other"),
         SizeCode = ifelse(Genus_sp2 == "Other", "S", "B"),
         SizeCode = fct_relevel(SizeCode, "B", "S")) %>% 
  ungroup() %>% 
  mutate(M = fct_relevel(M, "May", "Jun", "Jul", "Aug", "Sep"),
         Sample_site = fct_relevel(as.factor(Sample_site), 
                                            "36-873", "37-890", 
                                   "27-918", "29-905", 
                                   "14-972", "16-970", 
                                   "3-996" , "8-994"),
         Sample_site = fct_recode(Sample_site,"36" = "36-873", "37" = "37-890", 
                                  "14" = "14-972", "3" = "3-996",
                                  "27" = "27-918", "29" = "29-905", 
                                  "16" = "16-970", "8" = "8-994"),
         Genus_sp2 = fct_relevel(Genus_sp2, "Daphnia_retrocurva", "Veligers", "Skistodiaptomus_oregonensis",  "Mesocyclops_edax", "Other"))


levels(ZPfin2_b_AllYrAvg$Genus_sp2) <- c("D. retrocurva",
                                         "Mesocyclops spp.",
                                         "S. oregonensis",
                                         "Dreissena sp.",
                                         "Other")

png("05_Figures/12_RelBiomass_RankPlots_FigS1.png", 
    units = "in", height = 8, width = 8, res = 300)
ggplot(ZPfin2_b_AllYrAvg %>% 
         filter(Zoop_biomass >0), aes(y = log(Zoop_biomass), x = Rank)) +
  geom_line(color = "grey20") +
  geom_point(aes(color = Genus_sp2, size = SizeCode)) +
  scale_size_manual(values = c(3,1)) +
  scale_color_manual(values = c("#D55E00", "#CC97A7", "#0072B2", "#009E73", "#000000")) +
  guides(size = "none",  title = "Taxa", color = guide_legend(title = "Taxa")) +
  theme_bw() +
  facet_grid(Sample_site ~ M, scales = "free_y")+ #, labeller = label_parsed
  ylab(expression(paste(log[e]," Zooplankton biomass (µg ",L^-1,")"))) +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold")
  )
dev.off()
  


save.image("04_SavedRImages/12_relBiomass_FigS1.RData")
