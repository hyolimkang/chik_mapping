## load necessary data and library
# load library
##
##

# load functions
source("Functions/age_strat_subinf_func_final.R")
source("Functions/BurdenFunctions_v2.R")

## load data
## lhs data and hosp, mortality data
load("MainData/le_sample.RData")
load("MainData/nh_fatal_sample.RData")
load("MainData/fatal_sample.RData")
load("MainData/hosp_sample.RData")
load("MainData/lhs_sample_young.RData")
load("MainData/lhs_old.RData")
load("MainData/lhs_young.RData")

## FOI data
load("MainData/allfoi_s1.Data")
load("MainData/all_age_infection.RData")
load("MainData/combined_burden_shrink.RData")
## comorbid data
load("MainData/rr_hosp_model.RData")
load("MainData/")

## processing
all_age_infection <- all_age_infection[!all_age_infection$country %in% c("United States of America", "China" ),]

all_age_infection$tot_infec_med_group1 <- rowSums(all_age_infection[, paste0("mid", 1:3)], na.rm = TRUE)
all_age_infection$tot_infec_med_group2 <- rowSums(all_age_infection[, paste0("mid", 4:5)], na.rm = TRUE)
all_age_infection$tot_infec_med_group3 <- rowSums(all_age_infection[, paste0("mid", 6:7)], na.rm = TRUE)
all_age_infection$tot_infec_med_group4 <- rowSums(all_age_infection[, paste0("mid", 8:9)], na.rm = TRUE)
all_age_infection$tot_infec_med_group5 <- rowSums(all_age_infection[, paste0("mid", 10:11)], na.rm = TRUE)
all_age_infection$tot_infec_med_group6 <- rowSums(all_age_infection[, paste0("mid", 12:13)], na.rm = TRUE)
all_age_infection$tot_infec_med_group7 <- rowSums(all_age_infection[, paste0("mid", 14:15)], na.rm = TRUE)
all_age_infection$tot_infec_med_group8 <- rowSums(all_age_infection[, paste0("mid", 16:17)], na.rm = TRUE)
all_age_infection$tot_infec_med_group9 <- rowSums(all_age_infection[, paste0("mid", 18)], na.rm = TRUE)

infection_focal <- combined_burden %>% group_by(country) %>%
  summarise(tot_infec_med  = sum(med_inf),
            tot_infec_lo   = sum(lo_inf),
            tot_infec_hi   = sum(hi_inf),
            iso3           = first(iso3),
            country        = first(country),
            continent      = first(continent),
            tot_pop        = sum(tot_pop)
  ) %>% 
  as.data.frame()
tot_pop <- infection_focal %>% dplyr::select(country, tot_pop)
all_age_infection <- all_age_infection %>% left_join(tot_pop, by = "country")
all_age_infection <- all_age_infection[,-2]
region_classification <- read_excel("MainData/region_classification.xlsx", 
                                    sheet = "region")
continent <- region_classification$continent
all_age_infection <- cbind(continent, all_age_infection)
colnames(all_age_infection)[1] <- "continent"

