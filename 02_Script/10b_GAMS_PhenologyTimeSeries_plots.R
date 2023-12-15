# Predict phenol plots
# JMH Nov 2023, DEc 2023

# libraries ----
library(tidyverse)
library(tidygam)
library(ggeffects)
library(viridis)
library(RColorBrewer)


# Data files ----
# phenology metrics
ZPphen <- read.csv("03_GeneratedData/05_ZPPhenologyMetrics.csv", row.names = 1) %>% 
  mutate(DOY = as.numeric(as.character(strftime(as.POSIXct(DOY_moment, format = "%Y-%m-%d"), format= "%j"))),
         Sample_site = as.factor(Sample_site))


# how do things vary over time
## ret ----
gm_ret_15_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_ret15.RDS"))
gm_ret_50_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_ret50.RDS"))
gm_ret_85_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_ret85.RDS"))

## meso ----
gm_Meso_15_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Meso15.RDS"))
gm_Meso_50_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Meso50.RDS"))
gm_Meso_85_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Meso85.RDS"))

## Veligers ----
gm_Vel_15_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Vel15.RDS"))
gm_Vel_50_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Vel50.RDS"))
gm_Vel_85_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Vel85.RDS"))

## Orego ----
gm_Oreg_15_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Oreg15.RDS"))
gm_Oreg_50_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Oreg50.RDS"))
gm_Oreg_85_ML <- readRDS(file.path(here::here("06_Rdat"), "10_GAMS_seasMoments_ML_Oreg85.RDS"))

# Generate model pred ----

# Generate model pred ----
gge_RDS_Ret_15 <- ggpredict(gm_ret_15_ML, terms = c("Year", "Sample_site")) %>% 
                    as_tibble() %>% 
                    mutate(Taxa = "Ret",
                           Season = "per15")

gge_RDS_Ret_50 <- ggpredict(gm_ret_50_ML, terms = c("Year", "Sample_site"))  %>% 
                    as_tibble() %>% 
                    mutate(Taxa = "Ret",
                           Season = "per50")

# Need to go through and pair sample-site and sample-site2
gge_RDS_Ret_85 <- ggpredict(gm_ret_85_ML, terms = c("Year", "Sample_site", 
                                                    "Sample_site2 [East, West]"))  %>% 
  as_tibble() %>% 
  mutate(Keep = ifelse(group %in% c("36-873", "37-890", "27-918", "29-905") & facet == "West", "Keep",
                  ifelse(group %in% c("14-972","16-970", "3-996", "8-994") & facet == "East", "Keep", "Drop"))) %>% 
  filter(Keep != "Drop") %>% 
  mutate(Taxa = "Ret",
         Season = "per85") %>% 
  select(-c(facet, Keep)) %>% 
  mutate(predicted = case_when(Taxa == "Ret" & Season == "per85" & group == "36-873" ~ predicted,
                               Taxa == "Ret" & Season == "per85" & group == "37-890" ~ predicted+1.5,
                               Taxa == "Ret" & Season == "per85" & group == "27-918" ~ predicted+3,
                               Taxa == "Ret" & Season == "per85" & group == "29-905" ~ predicted+4.5,
                               Taxa == "Ret" & Season == "per85" & group == "14-972" ~ predicted,
                               Taxa == "Ret" & Season == "per85" & group == "16-970" ~ predicted+1.5,
                               Taxa == "Ret" & Season == "per85" & group == "3-996" ~ predicted+3,
                               Taxa == "Ret" & Season == "per85" & group == "8-994" ~ predicted+4.5))


gge_RDS_Meso_15 <- ggpredict(gm_Meso_15_ML, terms = c("Year", "Sample_site", 
                                                      "Sample_site2 [s14162729, s3, s36378]"))  %>% 
  as_tibble() %>% 
  mutate(Keep = ifelse(group %in% c("14-972","16-970", "27-918", "29-905") & facet == "s14162729", "Keep",
                  ifelse(group %in% c("3-996") & facet == "s3", "Keep", 
                    ifelse(group %in% c("36-873", "37-890", "8-994") & facet == "s36378", "Keep", "Drop")))) %>% 
  filter(Keep == "Keep") %>% 
  mutate(Taxa = "Meso",
         Season = "per15") %>% 
  select(-c(facet, Keep))



gge_RDS_Meso_50 <- ggpredict(gm_Meso_50_ML, terms = c("Year", "Sample_site"))  %>% 
  as_tibble() %>% 
  mutate(Taxa = "Meso",
         Season = "per50")



gge_RDS_Meso_85 <- ggpredict(gm_Meso_85_ML, terms = c("Year", "Sample_site", 
                                                      "Sample_site2 [s14162729, s363738]"))  %>% 
  as_tibble() %>% 
  mutate(Keep = ifelse(group %in% c("14-972","16-970",  "27-918", "29-905") & facet == "s14162729", "Keep",
                       ifelse(group %in% c("36-873", "37-890", "3-996", "8-994") & facet == "s363738", "Keep", "Drop"))) %>% 
  filter(Keep == "Keep") %>% 
    mutate(Taxa = "Meso",
           Season = "per85") %>% 
    select(-c(facet, Keep)) %>% 
  mutate(predicted = case_when(Taxa == "Meso" & Season == "per85" & group == "36-873" ~ predicted,
                               Taxa == "Meso" & Season == "per85" & group == "37-890" ~ predicted+1.5,
                               Taxa == "Meso" & Season == "per85" & group == "3-996" ~ predicted+3,
                               Taxa == "Meso" & Season == "per85" & group == "8-994" ~ predicted+4.5,
                               Taxa == "Meso" & Season == "per85" & group == "14-972" ~ predicted,
                               Taxa == "Meso" & Season == "per85" & group == "16-970" ~ predicted+1.5,
                               Taxa == "Meso" & Season == "per85" & group == "27-918" ~ predicted+3,
                               Taxa == "Meso" & Season == "per85" & group == "29-905" ~ predicted+4.5))
  
  
  
  
  
gge_RDS_Vel_15 <- ggpredict(gm_Vel_15_ML, terms = c("Year", "Sample_site"))  %>% 
                    as_tibble()%>% 
                    mutate(Taxa = "Vel",
                           Season = "per15") 




gge_RDS_Vel_50 <- ggpredict(gm_Vel_50_ML, terms = c("Year", "Sample_site", 
                                                    "Sample_site1 [Bimodal, Unimodal]"))  %>% 
  as_tibble() %>% 
  mutate(Keep = ifelse(group %in% c("16-970", "3-996", "8-994") & facet == "Unimodal", "Keep",
                       ifelse(group %in% c("36-873","37-890","29-905","27-918", "14-972") & facet == "Bimodal", "Keep", "Drop"))) %>% 
  filter(Keep == "Keep")%>% 
  mutate(Taxa = "Vel",
         Season = "per50")  %>% 
  select(-c(facet, Keep))




gge_RDS_Vel_85 <- ggpredict(gm_Vel_85_ML, terms = c("Year", "Sample_site", 
                                                    "Sample_site2 [s1416, s2729, s3637, s38]")) %>% 
  as_tibble() %>% 
  mutate(Keep = "Drop",
        Keep = case_when(facet == "s1416" & group %in% c("14-972","16-970") ~ "Keep",
                         facet == "s2729" & group %in% c("27-918", "29-905") ~ "Keep",
                         facet == "s3637" & group %in% c("36-873", "37-890") ~ "Keep",
                         facet == "s38" & group %in% c("3-996", "8-994") ~ "Keep"))  %>% 
        filter(Keep == "Keep")%>% 
        mutate(Taxa = "Vel",
               Season = "per85")  %>% 
        select(-c(facet, Keep))



gge_RDS_Oreg_15 <- ggpredict(gm_Oreg_15_ML, terms = c("Year", "Sample_site"))  %>% 
  as_tibble() %>% 
  mutate(Taxa = "Oreg",
         Season = "per15")

gge_RDS_Oreg_50 <- ggpredict(gm_Oreg_50_ML, terms = c("Year", "Sample_site"))  %>% 
  as_tibble() %>% 
  mutate(Taxa = "Oreg",
         Season = "per50")

gge_RDS_Oreg_85 <- ggpredict(gm_Oreg_85_ML, terms = c("Year", "Sample_site"))  %>% 
  as_tibble() %>% 
  mutate(Taxa = "Oreg",
         Season = "per85")

gge_RDS_dfs <- rbind(gge_RDS_Ret_15, gge_RDS_Ret_50, gge_RDS_Ret_85,
                      gge_RDS_Meso_15, gge_RDS_Meso_50, gge_RDS_Meso_85,
                      gge_RDS_Vel_15, gge_RDS_Vel_50, gge_RDS_Vel_85,
                      gge_RDS_Oreg_15, gge_RDS_Oreg_50, gge_RDS_Oreg_85) %>% 
                rename(Sample_site = group,
                       Year = x) %>% 
                # FIX: ALL OF THE SITES ARE STILL IN HERE
                mutate(Sample_site = ifelse(Taxa == "Vel" & Season == "per15", "All", 
                                            ifelse(Taxa == "Oreg" & Season == "per15", "All",
                                              ifelse(Taxa == "Oreg" & Season == "per50", "All",
                                                     as.character(Sample_site)))),
                       Sample_site = as.factor(Sample_site),
                       Sample_site = fct_relevel(Sample_site,
                                                 "All","36-873", "37-890", "27-918", "29-905",
                                                 "14-972", "16-970", 
                                                 "3-996", "8-994"),
                       Taxa = as.factor(Taxa),
                       Taxa = fct_relevel(Taxa, "Ret", "Vel", "Oreg", "Meso"),
                       PhenologyDate_pred = as.Date(paste("2010-", predicted), format = "%Y-%j"))


# Plot ----
colPal <- c("black","#d53e4f","#f46d43","#fdae61","#fee08b","#4292C6", "#9ECAE1","#762A83", "#C2A5CF")
Taxa.labels <- c("D. retrocurva", "Veligers", "S. oregonensis", "Mesocyclops")
names(Taxa.labels) <- c("Ret", "Vel", "Oreg", "Meso")

library(ggtext)

png("05_Figures/10b_Fig4_GAMSSeasMomentsPlot.png",
    units = "in", height = 9, width = 6, res = 300)
ggplot(gge_RDS_dfs, aes(y = PhenologyDate_pred, x = Year, 
                        # linewidth = Season,
                        linetype = Season,
                        color = Sample_site)) +
  geom_line(size = 1.1) +
  # scale_linewidth_manual(values = c(1.1, 0.9, 1.1)) +
  scale_linetype_manual(values = c("dashed", "dotted", "solid"),
                        labels = c(expression(DOY[start]),
                                   expression(DOY[middle]),
                                   expression(DOY[end]))) +
  scale_color_manual(values = colPal,
                     labels = c("All",
                                "36", "37", "27", "29",
                                "14", "16", "3", "8"),
                     name = "Site")+
  facet_wrap(Taxa ~ ., ncol = 1, labeller = labeller(Taxa = Taxa.labels)) +
  ylab("Month") +
  theme_bw() +
  theme(axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.position = "right",
        legend.key.width = unit(1.25, "cm"),
        strip.text = element_text(face = "bold.italic", size = 16),
        strip.background = element_rect(color = "transparent"),
        panel.background = element_rect(fill = "white", color = "black")) 
dev.off()


# save.image(here::here("04_SavedRImages", "10b_GAMS_seasonalMomentsR_plots.Rdata"))
# load(here::here("04_SavedRImages", "10b_GAMS_seasonalMomentsR_plots.Rdata"))
