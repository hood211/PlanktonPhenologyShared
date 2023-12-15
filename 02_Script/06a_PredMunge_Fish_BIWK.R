# WB fish dat
# Jan 23, JMH

# libaries
library(tidyverse)
library(readxl)
library(here)
library(readr)
library(ggcorrplot)

# western basin ODNR trawl data from Zak Slagle (21 Oct. 22)
# Catches prior to 2002 standardized to new R/V Explorer using Fishing Power Corrections (Tyson et al. 2006). 
# NOTE: some species and age classes may not have an associated FPC due to low catches during that study, resulting in catch numbers that be not be comparable before and after 2002.
fish <- readxl::read_xlsx(file.path(here::here("01_OriginalData"),
                                    "22_10_18 Hood planktivore catch at LEPAS sites.xlsx"),
                           sheet = "data")

# planktivores ----
# based on conversation with Stu Ludsin prior to this project and some literature searches

  unique(fish$common_name) %>% sort()
  
  plank_AllAges <- c("Alewife", 
                     "Brook silverside",
                     "Emerald shiner",
                     "Goldfish",
                     "Mimic shiner",
                     "Notropis sp.",
                     "Rainbow smelt",
                     "Spottail shiner") # 3/4 zp based on gopalan
  
  plankAge_0_only <- c("Bluegill", # age zero only?
                       "Freshwater drum", # age-0 only, ~60% zp based on Gopalan et al. 1998
                       "Gizzard shad",
                       "Orangespotted sunfish", #only age zero in db
                       "Trout-perch", #age-0 only, 1/2 zp based on gopalan
                       "White bass", #age-0 only, >80% zp based on gopalan
                       "White crappie", #age-0 only
                       "White perch", #age-0 only
                       "Yellow perch") # age-0 only
  
  
  plank <- fish %>% 
            filter(common_name %in% plank_AllAges | (common_name %in% plankAge_0_only & ages == "Age-0")) 

  

# Data exploration ----
  plankS <- plank %>% 
            group_by(year, month, day, siteno) %>% 
            summarise(plank_per_ha = sum(fish_per_ha), na.rm = T) %>% 
            mutate(date = as.POSIXct(paste0(year,"-", month, "-", day), format = "%Y-%m-%d"),
                   DOY = as.numeric(strftime(date, format = "%j")),
                   Y = as.numeric(strftime(date, format = "%Y")),
                   M = as.numeric(strftime(date, format = "%m"))) %>% 
            ungroup() %>% 
            select(date, Y, M, DOY, siteno, plank_per_ha) %>% 
            # Oct data is sparse
            filter(M < 10)

  # planktivores v. DOY
  # pretty strong seasonality at some sites
  ggplot(plankS, aes(y = plank_per_ha, x = M, color = as.factor(Y))) +
    geom_point() +
    facet_wrap(vars(siteno)) +
    # no zeros
    scale_y_log10()
  
  # Planktivores v. Year

  ggplot(plankS, aes(y = plank_per_ha, x = Y, color = as.factor(M))) +
    geom_point() +
    geom_line() +
    facet_grid(as.factor(M) ~ siteno, scales = "free_y") +
    # no zeros
    scale_y_log10()
  
  # no NAs in Aug
  plankS2 <-   plankS %>% 
    mutate(Response = "plank_per_ha",
           siteno = case_when(siteno == "3" ~ "3-996",
                              siteno == "8" ~ "8-994",
                              siteno == "14" ~ "14-972",
                              siteno == "16" ~ "16-970",
                              siteno == "27" ~ "27-918",
                              siteno == "29" ~ "29-905",
                              siteno == "36" ~ "36-873",
                              siteno == "37" ~ "37-890")) %>% 
    select(Sample_date = date, Sample_site = siteno, Response, Value = plank_per_ha)
  
  

  

# export data ----  
  write.csv(plankS2, "03_GeneratedData/06a_WBplanktivoreDensity_NotSummed2Month.csv")
  
# save image
  # save.image("04_SavedRImages/06a_FishMungingWB_NotSummed2Month.RData")
  # load("04_SavedRImages/06a_FishMungingWB_NotSummed2Month.RData")
  