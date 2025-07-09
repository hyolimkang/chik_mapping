# age stratified fatality rates 

age_fatality <- hosp_chikv %>% filter(death_chikv == "TRUE")


# load data
## age stratified burden
load("MainData/burden_shrink_1.RData")
load("MainData/burden_shrink_2.RData")
load("MainData/burden_shrink_3.RData")
load("MainData/burden_shrink_4.RData")
load("MainData/combined_burden_shrink.RData")
load("MainData/all_age_infection.RData")


agelist1 <- lapply(burden_shrink_1, function(x) x$infection_per_band)
agelist2 <- lapply(burden_shrink_2, function(x) x$infection_per_band)
agelist3 <- lapply(burden_shrink_3, function(x) x$infection_per_band)
agelist4 <- lapply(burden_shrink_4, function(x) x$infection_per_band)

agelist1 <- lapply(agelist1, as.data.frame)
agelist2 <- lapply(agelist2, as.data.frame)
agelist3 <- lapply(agelist3, as.data.frame)
agelist4 <- lapply(agelist4, as.data.frame)


col_names <- c(1:18)
agelist1 <- lapply(agelist1, function(df){
  colnames(df) <- col_names
  return(df)
})
agelist2 <- lapply(agelist2, function(df){
  colnames(df) <- col_names
  return(df)
})
agelist3 <- lapply(agelist3, function(df){
  colnames(df) <- col_names
  return(df)
})
agelist4 <- lapply(agelist4, function(df){
  colnames(df) <- col_names
  return(df)
})

comb_df_age1 <- list()
comb_df_age2 <- list()
comb_df_age3 <- list()
comb_df_age4 <- list()

for(i in col_names){
  comb_df_age1[[i]] <- do.call(cbind, lapply(agelist1, function(df) df[, i, drop = F]))
}
for(i in col_names){
  comb_df_age2[[i]] <- do.call(cbind, lapply(agelist2, function(df) df[, i, drop = F]))
}
for(i in col_names){
  comb_df_age3[[i]] <- do.call(cbind, lapply(agelist3, function(df) df[, i, drop = F]))
}
for(i in col_names){
  comb_df_age4[[i]] <- do.call(cbind, lapply(agelist4, function(df) df[, i, drop = F]))
}

num_dfs <- length(comb_df_age1)
all_agelist <- vector("list", length = num_dfs)

for(i in 1:num_dfs){
  all_agelist[[i]] <- do.call(cbind, list(comb_df_age1[[i]],
                                          comb_df_age2[[i]],
                                          comb_df_age3[[i]],
                                          comb_df_age4[[i]]))
}

all_agelist <- lapply(all_agelist, function(df) {
  mid_val <- apply(df[,1:100], 1, function(x) quantile(x, 0.5, na.rm = T))
  #lo_val <- apply(df[,1:100], 1, function(x) quantile(x, 0.025, na.rm = T))
  #hi_val <- apply(df[,1:100], 1, function(x) quantile(x, 0.975, na.rm = T))
  #df$lo_val <- lo_val
  #df$hi_val <- hi_val
  df$mid_val <- mid_val
  
  return(df)
})

# extracting common cols 
common_cols <- combined_burden[, c(1:27, 133)]

all_agelist <-lapply(all_agelist, function(df){
  cbind(common_cols, df)
  
})

extract_vals <- lapply(1: length(all_agelist), function(i) {
  df <- all_agelist[[i]]
  common_cols <- df[,1:9]
  age_group <- df[,10:27]
  tot_pop <- df[,28]
  result_df <- data.frame(age_group = i,
                          mid       = df$mid_val
                          #lo        = df$lo,
                          #hi        = df$hi
  )
  result_df <- cbind(common_cols, age_group, tot_pop, result_df)
  return(result_df)
})

# global level summary
all_agelist_global <- lapply(1: length(extract_vals), function(i){
  df <- extract_vals[[i]]
  data.frame(age_group = i,
             mid = sum(df$mid, na.rm = TRUE)
             #lo  = sum(df$lo),
             #hi  = sum(df$hi)
  )
})

# regional level summary
all_agelist_cont <- lapply(seq_along(extract_vals), function(i){
  df <- extract_vals[[i]]
  
  pop_cols <-  as.character(1:18)
  
  sum_df <- df %>% group_by(continent, age_group) %>%
    summarise(age_group = i,
              mid = sum(mid, na.rm = TRUE),
              #lo  = sum(lo)
              #hi  = sum(hi)
    )
  
  pop_sum_df <- df %>%
    group_by(continent) %>%
    summarise(across(all_of(pop_cols), \(x) sum(x, na.rm = TRUE)),
              .groups = 'drop')
  
  final_df <- left_join(sum_df, pop_sum_df, by = "continent")
  
  return(final_df)
})

# country level summary
all_agelist_nat <- lapply(seq_along(extract_vals), function(i){
  df <- extract_vals[[i]]
  
  pop_cols <-  as.character(1:18)
  
  sum_df <- df %>% group_by(country, continent, age_group) %>%
    summarise(age_group = i,
              mid = sum(mid, na.rm = TRUE),
              tot_pop = sum(tot_pop, na.rm = TRUE)
              #lo  = sum(lo)
              #hi  = sum(hi)
    )%>%
    rename(!!paste0("mid", i) := mid)
  
  pop_sum_df <- df %>%
    group_by(country, continent) %>%
    summarise(across(all_of(pop_cols), \(x) sum(x, na.rm = TRUE)),
              .groups = 'drop')
  
  final_df <- left_join(sum_df, pop_sum_df, by = c("country", "continent"))
  
  return(final_df)
})

mid_cols <- paste0("mid", 1:18)
infection_cols <- lapply(seq_along(all_agelist_nat), function(i) {
  df <- all_agelist_nat[[i]]
  df %>% select(starts_with(paste0("mid", i)))
})
combined_df <- bind_cols(infection_cols) %>% select(starts_with("mid"))
common_cols <- all_agelist_nat[[1]][,c(1:2, 5, 6:23)]

all_age_infection <- cbind(common_cols, combined_df)
save(all_age_infection, file = "MainData/all_age_infection.RData")



################################################################################

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

# load combined_burden_shirnk data
load("MainData/combined_burden_shrink.RData")
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
tot_pop <- infection_focal %>% select(country, tot_pop)
all_age_infection <- all_age_infection %>% left_join(tot_pop, by = "country")
all_age_infection <- all_age_infection[,-2]
region_classification <- read_excel("MainData/region_classification.xlsx", 
                                    sheet = "region")
continent <- region_classification$continent
all_age_infection <- cbind(continent, all_age_infection)
colnames(all_age_infection)[1] <- "continent"
# burden under 40 years old 
under40_burden <- age_specific_subinf_under40(all_age_infection, lhs_sample)
results_under40 <- list()

# Loop through each sublist (each group) in the results
for (group in names(under40_burden)) {
  
  # Get the sublist for the current group
  sub_burden_psa <- under40_burden[[group]]
  
  # Run the postprocess_sub_burden function
  results_under40[[group]] <- postprocess_sub_burden(sub_burden_psa, all_age_infection, sub_inf_global_count)
}

# burden over 40 years old 
over40_burden <- age_specific_subinf_over40(all_age_infection, lhs_sample)
results_over40 <- list()

# Loop through each sublist (each group) in the results
for (group in names(over40_burden)) {
  
  # Get the sublist for the current group
  sub_burden_psa <- over40_burden[[group]]
  
  # Run the postprocess_sub_burden function
  results_over40[[group]] <- postprocess_sub_burden(sub_burden_psa, all_age_infection, sub_inf_global_count)
}

save(results_under40, file = "MainData/results_under40_focal.RData")
save(results_over40, file = "MainData/results_over40_focal.RData")


# burden all 

all_burden <-     results_under40$Group_1$all_results_global[,1:3] +
                  results_under40$Group_2$all_results_global[,1:3] +
                  results_under40$Group_3$all_results_global[,1:3] +
                  results_under40$Group_4$all_results_global[,1:3] +
                  results_over40$Group_5$all_results_global[,1:3] +
                  results_over40$Group_6$all_results_global[,1:3] +
                  results_over40$Group_7$all_results_global[,1:3] +
                  results_over40$Group_8$all_results_global[,1:3] +
                  results_over40$Group_9$all_results_global[,1:3]

## daly by age group
results_under40$Group_1$all_results_global$age_group <- "[0,10)"
results_under40$Group_2$all_results_global$age_group <- "[10,20)"
results_under40$Group_3$all_results_global$age_group <- "[20,30)"
results_under40$Group_4$all_results_global$age_group <- "[30,40)"
results_over40$Group_5$all_results_global$age_group  <- "[40,50)"
results_over40$Group_6$all_results_global$age_group  <- "[50,60)"
results_over40$Group_7$all_results_global$age_group  <- "[60,70)"
results_over40$Group_8$all_results_global$age_group  <- "[70,80)"
results_over40$Group_9$all_results_global$age_group  <- "[80,90)"

global_age1 <- results_under40$Group_1$all_results_global
global_age2 <- results_under40$Group_2$all_results_global
global_age3 <- results_under40$Group_3$all_results_global
global_age4 <- results_under40$Group_4$all_results_global
global_age5 <- results_over40$Group_5$all_results_global
global_age6 <- results_over40$Group_6$all_results_global
global_age7 <- results_over40$Group_7$all_results_global
global_age8 <- results_over40$Group_8$all_results_global
global_age9 <- results_over40$Group_9$all_results_global

global_all <- rbind(global_age1, global_age2, global_age3, global_age4,
                    global_age5, global_age6, global_age7, global_age8,
                    global_age9)

write.csv(all_burden, file = "MainData/global_all_case.csv")
write.csv(all_burden_cont, file = "MainData/global_all_case_cont.csv")

## age specific daly graph
daly <- global_all[global_all$type %in% c("yld_acute", "yld_subac", "yld_chronic", "yll"),]
daly <- daly %>%
  group_by(age_group) %>%
  mutate(percent = tot_med / sum(tot_med) * 100)

daly_by_component <- daly %>% group_by(type) %>% 
                             summarise(tot_med = sum(tot_med),
                                       tot_lo  = sum(tot_lo),
                                       tot_hi  = sum(tot_hi))

tot_daly <- global_all[global_all$type %in% c("daly"),]

write.csv(daly_by_component, file = "MainData/daly_component.csv")

color <- c("#00468B99", "#ED000099", "#42B54099", "#FDAF9199")

A1 <- ggplot(daly, aes(x = age_group, y = percent, fill = type))+
  geom_bar(stat = "identity")+
  theme_bw()+
  labs(x = "Age group",
       y = "Proportion of DALY component",
       fill = "DALY component")+
  scale_y_continuous(labels = scales::percent_format(scale = 1))+
  scale_fill_manual(values = color,
                    labels = c("yld_acute" = "Acute morbidity (YLD)", "yld_chronic" = "Chronic morbidity (YLD)", "yld_subac" = "Sub-acute morbidity (YLD)", "yll" = "Years of life lost (YLL)"))+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, hjust = 0.5, vjust = 1.2),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) 

ggsave(filename = "02_Outputs/2_1_Figures/fig_daly_component.jpg", plot = A1,
       width = 7, height = 4)


dodge_width <- 0.9
A2 <- ggplot(daly, aes(x = age_group, y = tot_med, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = tot_lo, ymax = tot_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "DALY",
       fill = "DALY component") +
  scale_fill_manual(values = color,
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, hjust = 0.5),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

## DALY by continent
# List of all the relevant data frames
groups_under40 <- list(results_under40$Group_1$all_results_continent[, 1:6],
                       results_under40$Group_2$all_results_continent[, 1:6],
                       results_under40$Group_3$all_results_continent[, 1:6],
                       results_under40$Group_4$all_results_continent[, 1:6])

groups_over40 <- list(results_over40$Group_5$all_results_continent[, 1:6],
                      results_over40$Group_6$all_results_continent[, 1:6],
                      results_over40$Group_7$all_results_continent[, 1:6],
                      results_over40$Group_8$all_results_continent[, 1:6],
                      results_over40$Group_9$all_results_continent[, 1:6])

# Sum all the data frames together
all_burden_cont <- Reduce(`+`, lapply(c(groups_under40, groups_over40), function(df) {
  df[, 2:4]
}))
all_burden_cont <- cbind(continent = groups_under40[[1]]$continent, all_burden_cont,
                         type = groups_under40[[1]]$type)

age_groups <- c("[0,10)", "[10,20)", "[20,30)", "[30,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)","[80,90)")

# Combine the lists
all_groups <- c(results_under40, results_over40)
for (i in seq_along(all_groups)) {
  all_groups[[i]]$all_results_continent$age_group <- age_groups[i]
}


all_daly_cont <- do.call(rbind, lapply(all_groups, function(group) group$all_results_continent))
daly_cont <- all_daly_cont[all_daly_cont$type %in% c("yld_acute", "yld_subac", "yld_chronic", "yll"),]

B <- ggplot(daly_cont, aes(x = age_group, y = tot_med, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = tot_lo, ymax = tot_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "DALY",
       fill = "DALY component") +
  scale_fill_manual(values = color,
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  #scale_y_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 9),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9))+
  facet_wrap(~continent)

# chronic 
chronic_cont <- all_daly_cont[all_daly_cont$type %in% c("hospitalisation_chronic6m", "hospitalisation_chronic12m", "hospitalisation_chronic30m", "nonhosp_chronic_6m",
                                                        "nonhosp_chronic_12m", "nonhosp_chronic_30m"),]

hosp_chronic <- chronic_cont %>% filter(grepl("hospitalisation", type))
nonhosp_chronic <- chronic_cont %>% filter(grepl("nonhosp_chronic", type))

hosp_chronic <- hosp_chronic %>% group_by(age_group, continent) %>% 
                      summarise(chronic_med = sum(tot_med),
                                chronic_lo  = sum(tot_lo),
                                chronic_hi  = sum(tot_hi))
  
nonhosp_chronic <- nonhosp_chronic %>% group_by(age_group, continent) %>% 
  summarise(chronic_med = sum(tot_med),
            chronic_lo  = sum(tot_lo),
            chronic_hi  = sum(tot_hi))

hosp_chronic_bind <- rbind(hosp_chronic, nonhosp_chronic)

chronic_grouped <- hosp_chronic_bind %>%
  mutate(group = case_when(
    str_detect(type, "nonhosp.*_12m") ~ "chronic_12m",
    str_detect(type, "nonhosp.*_6m")  ~ "chronic_6m",
    str_detect(type, "hosp.*_12m")    ~ "chronic_12m",
    str_detect(type, "hosp.*_6m")     ~ "chronic_6m",
    str_detect(type, "hosp.*_30m")    ~ "chronic_30m",
    str_detect(type, "nonosp.*_30m")  ~ "chronic_30m",
    str_detect(type, "hospitalisation.*_chronic6m")  ~ "chronic_6m",
    str_detect(type, "hospitalisation.*_chronic12m")  ~ "chronic_12m",
    str_detect(type, "hospitalisation.*_chronic30m")  ~ "chronic_30m"
  ))

chronic_grouped <- chronic_grouped %>% group_by(age_group, continent, group) %>% 
  summarise(chronic_med = sum(tot_med),
            chronic_lo  = sum(tot_lo),
            chronic_hi  = sum(tot_hi))


C <- ggplot(chronic_grouped, aes(x = age_group, y = chronic_med, ymin = chronic_lo, ymax = chronic_hi, color = group))+
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.3)+
  #geom_point(aes(x = age_group, y = chronic_med))+
  facet_wrap(~continent)+
  theme_light()+
  scale_y_continuous(labels = scales::comma)+
  ggtitle("Chronic Chikungunya by month/age/continent")+
  theme(plot.title = element_text(size = 10))+
  scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8))+
  labs(x = "Age group",
       y = "Chronic chikungunya (N)",
       fill = "Type")+
  facet_wrap(~continent)

# daly per 10,000 
age_spec_pop <- all_age_infection %>%
  group_by(continent) %>%
  summarise(
    age1 = sum(`1` + `2` + `3`),
    age2 = sum(`4` + `5`),
    age3 = sum(`6` + `7`),
    age4 = sum(`8` + `9`),
    age5 = sum(`10` + `11`),
    age6 = sum(`12` + `13`),
    age7 = sum(`14` + `15`),
    age8 = sum(`16` + `17`),
    age9 = sum(`18`)
  )

age_spec_pop <- melt(as.data.frame(age_spec_pop), id = "continent")

expanded_age <- age_spec_pop %>%
  group_by(variable) %>%
  slice(rep(1:7, times = 26))
expanded_age <- expanded_age[,3]
colnames(expanded_age) <- "tot_pop_age"

all_daly_cont <- cbind(all_daly_cont, expanded_age)

all_daly_cont <- all_daly_cont %>% mutate(per10k_mid = (tot_med/tot_pop_age)*10000,
                                          per10k_lo  = (tot_lo/tot_pop_age)*10000,
                                          per10k_hi  = (tot_hi/tot_pop_age)*10000)

daly_cont <- all_daly_cont[all_daly_cont$type %in% c("yld_acute", "yld_subac", "yld_chronic", "yll"),]

daly_stack <- daly_cont %>% group_by(continent, type) %>%
                          summarise(per10k_mid = sum(per10k_mid),
                                    per10k_lo  = sum(per10k_lo),
                                    per10k_hi  = sum(per10k_hi))

daly_age <- all_daly_cont %>% group_by(age_group, continent) %>% mutate(per10k_mid = (tot_med/tot_pop_age)*10000,
                                                                        per10k_lo  = (tot_lo/tot_pop_age)*10000,
                                                                        per10k_hi  = (tot_hi/tot_pop_age)*10000)

symp_age <- daly_age[daly_age$type %in% c("symptomatic"),]
fatal_age <- daly_age[daly_age$type %in% c("fatal"),]
total_daly_age <- daly_age[daly_age$type %in% c("daly"),]
hosp_age <- daly_age[daly_age$type %in% c("hospitalisation"),]
chronic_age <- daly_age %>% filter(type %in% c("hospitalisation_chronic6m",
                                               "hospitalisation_chronic12m",
                                               "hospitalisation_chronic30m",
                                               "nonhosp_chronic_6m",
                                               "nonhosp_chronic_12m",
                                               "nonhosp_chronic_30m"))%>%
               group_by(continent, age_group)%>%
               summarise(tot_med = sum(tot_med),
                         tot_lo  = sum(tot_lo),
                         tot_hi  = sum(tot_hi),
                         .groups = "keep") %>% ungroup()

symp_age_all <- symp_age %>% group_by(continent) %>% summarise(tot_med = sum(tot_med),
                                                               tot_lo  = sum(tot_lo),
                                                               tot_hi  = sum(tot_hi),
                                                               tot_pop = sum(tot_pop_age),
                                                               .groups = "keep") %>% mutate(per10k_mid = (tot_med/tot_pop)*10000,
                                                                                            per10k_lo  = (tot_lo/tot_pop)*10000,
                                                                                            per10k_hi  = (tot_hi/tot_pop)*10000)
daly_age_all <- daly_age %>% group_by(continent) %>% summarise(tot_med = sum(tot_med),
                                                               tot_lo  = sum(tot_lo),
                                                               tot_hi  = sum(tot_hi),
                                                               tot_pop = sum(tot_pop_age),
                                                               .groups = "keep") %>% mutate(per10k_mid = (tot_med/tot_pop)*10000,
                                                                                            per10k_lo  = (tot_lo/tot_pop)*10000,
                                                                                            per10k_hi  = (tot_hi/tot_pop)*10000)

fatal_age_all <- fatal_age %>% group_by(continent) %>% summarise(tot_med = sum(tot_med),
                                                               tot_lo  = sum(tot_lo),
                                                               tot_hi  = sum(tot_hi),
                                                               tot_pop = sum(tot_pop_age),
                                                               .groups = "keep") %>% mutate(per10k_mid = (tot_med/tot_pop)*10000,
                                                                                            per10k_lo  = (tot_lo/tot_pop)*10000,
                                                                                            per10k_hi  = (tot_hi/tot_pop)*10000)

hosp_age_all <- hosp_age %>% group_by(continent) %>% summarise(tot_med = sum(tot_med),
                                                                 tot_lo  = sum(tot_lo),
                                                                 tot_hi  = sum(tot_hi),
                                                                 tot_pop = sum(tot_pop_age),
                                                                 .groups = "keep") %>% mutate(per10k_mid = (tot_med/tot_pop)*10000,
                                                                                              per10k_lo  = (tot_lo/tot_pop)*10000,
                                                                                              per10k_hi  = (tot_hi/tot_pop)*10000)

symp_age_all$type <- "symptomatic"
fatal_age_all$type <- "fatal"
hosp_age_all$type <- "hosp"
daly_age_all$type <- "daly"

symp_age_global <- symp_age %>%
  group_by(age_group) %>%
  summarise(
    tot_med = sum(tot_med),
    tot_lo  = sum(tot_lo),
    tot_hi  = sum(tot_hi),
    tot_pop = sum(tot_pop_age),
    .groups = "drop"
  ) %>%
  mutate(
    per100k_mid = (tot_med / tot_pop) * 100000,
    per100k_lo  = (tot_lo / tot_pop) * 100000,
    per100k_hi  = (tot_hi / tot_pop) * 100000
  )

hosp_age_global <- hosp_age %>%
  group_by(age_group) %>%
  summarise(
    tot_med = sum(tot_med),
    tot_lo  = sum(tot_lo),
    tot_hi  = sum(tot_hi),
    tot_pop = sum(tot_pop_age),
    .groups = "drop"
  ) %>%
  mutate(
    per100k_mid = (tot_med / tot_pop) * 100000,
    per100k_lo  = (tot_lo / tot_pop) * 100000,
    per100k_hi  = (tot_hi / tot_pop) * 100000
  )

fatal_age_global <- fatal_age %>%
  group_by(age_group) %>%
  summarise(
    tot_med = sum(tot_med),
    tot_lo  = sum(tot_lo),
    tot_hi  = sum(tot_hi),
    tot_pop = sum(tot_pop_age),
    .groups = "drop"
  ) %>%
  mutate(
    per100k_mid = (tot_med / tot_pop) * 100000,
    per100k_lo  = (tot_lo / tot_pop) * 100000,
    per100k_hi  = (tot_hi / tot_pop) * 100000
  )

daly_age_global <- total_daly_age %>%
  group_by(age_group) %>%
  summarise(
    tot_med = sum(tot_med),
    tot_lo  = sum(tot_lo),
    tot_hi  = sum(tot_hi),
    tot_pop = sum(tot_pop_age),
    .groups = "drop"
  ) %>%
  mutate(
    per100k_mid = (tot_med / tot_pop) * 100000,
    per100k_lo  = (tot_lo / tot_pop) * 100000,
    per100k_hi  = (tot_hi / tot_pop) * 100000
  )


allcont_incidence <- rbind(symp_age_all, fatal_age_all, hosp_age_all, daly_age_all)

write.csv(symp_age, file = "MainData/symp_age.csv")

C <- ggplot(daly_stack, aes(x = continent, y = per10k_mid, fill = type)) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per10k_lo, ymax = per10k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Continent",
       y = "DALY per 10,000 people",
       fill = "DALY component") +
  scale_fill_manual(values = color,
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

# symptomatic, fatal, daly by age
D1 <- ggplot(symp_age, aes(x = age_group, y = per10k_mid, fill = "#E64B35B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per10k_lo, ymax = per10k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "Symptomatic per 10,000 people") +
  scale_fill_manual(values = "#E64B35B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")+
  facet_wrap(~continent)

D1 <- ggplot(symp_age_global, aes(x = age_group, y = per100k_mid, fill = "#E64B35B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per100k_lo, ymax = per100k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "Symptomatic per 100,000 people") +
  scale_fill_manual(values = "#E64B35B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")


D2 <- ggplot(hosp_age, aes(x = age_group, y = per10k_mid, fill = "#4DBBD5B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per10k_lo, ymax = per10k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "Hospitalised per 100,000 people") +
  scale_fill_manual(values = "#4DBBD5B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")+
  facet_wrap(~continent)

D2 <- ggplot(hosp_age_global, aes(x = age_group, y = per100k_mid, fill = "#4DBBD5B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per100k_lo, ymax = per100k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "Hospitalised per 100,000 people") +
  scale_fill_manual(values = "#4DBBD5B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")


D3 <- ggplot(fatal_age, aes(x = age_group, y = per10k_mid, fill = "#00A087B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per10k_lo, ymax = per10k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "Fatal case per 10,000 people") +
  scale_fill_manual(values = "#00A087B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")+
  facet_wrap(~continent)

D3 <- ggplot(fatal_age_global, aes(x = age_group, y = per100k_mid, fill = "#00A087B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per100k_lo, ymax = per100k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "Fatal case per 100,000 people") +
  scale_fill_manual(values = "#00A087B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")

D4 <- ggplot(total_daly_age, aes(x = age_group, y = per10k_mid, fill = "#3C5488B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per10k_lo, ymax = per10k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "DALY per 10,000 people") +
  scale_fill_manual(values = "#3C5488B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")+
  facet_wrap(~continent)

D4 <- ggplot(daly_age_global, aes(x = age_group, y = per100k_mid, fill = "#3C5488B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per100k_lo, ymax = per100k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "DALY per 100,000 people") +
  scale_fill_manual(values = "#3C5488B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")


D5 <- ggplot(allcont_incidence, aes(x = continent, y = per10k_mid, fill = "#3C5488B2")) +
  geom_bar(stat = "identity", position = position_dodge(width = dodge_width)) +  # Stacked bar chart with dodging
  geom_errorbar(aes(ymin = per10k_lo, ymax = per10k_hi), 
                width = 0.2, 
                position = position_dodge(width = dodge_width), 
                color = "black") +  # Dodged error bars
  theme_light() +
  labs(x = "Age group",
       y = "DALY per 10,000 people") +
  scale_fill_manual(values = "#3C5488B2",
                    labels = c("yld_acute" = "YLD Acute", "yld_chronic" = "YLD Chronic", "yld_subac" = "YLD Subacute", "yll" = "YLL")) +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_continuous(labels = scales::comma, trans = 'sqrt')+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 8,  hjust = 0.5, angle = 45, vjust = 0.5),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))+
  theme(legend.position = "none")+
  facet_wrap(~type)

pdf("Results_figs/DALY_fig1_new.pdf", height = 8, width = 12)
#ggarrange(A1, A2, 
#          labels = c("A", "B"),
#          ncol = 2, nrow = 2)
A1 <- A1 + theme(legend.position = "none")
A2 <- A2 + theme(legend.position = "none")
ggarrange(A1, A2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 2,  # Two columns, one row
          widths = c(2,2)
)
dev.off()

pdf("Results_figs/DALY_fig2_new.pdf", height = 10, width = 10)
#ggarrange(B, C, 
#          labels = c("C", "D"),
#          ncol = 2, nrow = 2)
ggarrange(B, 
          labels = c("C"),
          ncol = 1, nrow = 1,  # Two columns, one row
          widths = c(1,1),  # Adjust these values to control the space between A1 and A2
          common.legend = TRUE,
          legend = "bottom"
)
dev.off()

pdf("Results_figs/DALY_fig3_new.pdf", height = 8, width = 12)
#ggarrange(A1, A2, 
#          labels = c("A", "B"),
#          ncol = 2, nrow = 2)

D1 <- D1 + theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank()
)

D2 <- D2 + theme(
  axis.title.x = element_blank(),
  axis.text.x  = element_blank(),
  axis.ticks.x = element_blank()
)

D1 <- D1 + ggtitle("Global symptomatic cases per 100K population by age")
D2 <- D2 + ggtitle("Global hospitalised cases per 100K population by age")
D3 <- D3 + ggtitle("Global fatal cases per 100K population by age")
D4 <- D4 + ggtitle("Global DALY per 100K population by age")

D <- ggarrange(D1, D2, D3, D4,
          #labels = c("D", "E", "F" ,"G"),
          ncol = 2, nrow = 2,  # Two columns, one row
          widths = c(2,2)
)
ggsave(filename = "02_Outputs/2_1_Figures/fig_global_burden_age.jpg", plot = D,
       width = 10, height = 6)

dev.off()


## country level 
age_groups <- c("[0,10)", "[10,20)", "[20,30)", "[30,40)", "[40,50)", "[50,60)", "[60,70)", "[70,80)","[80,90)")

# Combine the lists
all_groups <- c(results_under40, results_over40)
for (i in seq_along(all_groups)) {
  all_groups[[i]]$all_results_country$age_group <- age_groups[i]
}
all_countries <- do.call(rbind, lapply(all_groups, function(group) group$all_results_country))
country_daly <- all_countries[all_countries$type %in% c("yld_acute", "yld_subac", "yld_chronic", "yll"),]
country_symp  <- all_countries[all_countries$type %in% c("symptomatic"),]
country_symp  <- country_symp %>% group_by(country) %>% summarise(
  tot_med = sum(tot_med),
  tot_lo  = sum(tot_lo),
  tot_hi  = sum(tot_hi),
  continent = first(continent),
  tot_pop = first(tot_pop)
)

save(all_countries, file = "00_Data/0_2_Processed/allcountries_symp.RData")
save(country_symp, file = "00_Data/0_2_Processed/country_symp.RData")

tot_vals <-country_daly %>%
  group_by(country, continent) %>%  # Group by both country and continent
  summarise(
    tot_med = sum(tot_med, na.rm = TRUE),
    tot_lo  = sum(tot_lo, na.rm = TRUE),
    tot_hi  = sum(tot_hi, na.rm = TRUE),
    .groups = 'drop'  # Optionally drop grouping after summarising
  )

tot_vals_cont <-country_daly %>%
  group_by(continent) %>%  # Group by both country and continent
  summarise(
    tot_med = sum(tot_med, na.rm = TRUE),
    tot_lo  = sum(tot_lo, na.rm = TRUE),
    tot_hi  = sum(tot_hi, na.rm = TRUE),
    .groups = 'drop'  # Optionally drop grouping after summarising
  )

write.csv(tot_vals, file = "MainData/country_daly.csv")


country_rank <- country_daly %>% group_by(country, type) %>% mutate %>%
                  summarise(tot_med = sum(tot_med),
                            tot_lo  = sum(tot_lo),
                            tot_hi  = sum(tot_hi))

country_rank <- country_rank %>%
  mutate(type = case_when(
    type == "yld_acute" ~ "YLD acute",
    type == "yld_subac" ~ "YLD sub-acute",
    type == "yld_chronic" ~ "YLD chronic",
    type == "yll" ~ "YLL"
  ))

country_rank <- country_rank %>%
  group_by(country) %>%
  mutate(rank = rank(-tot_med, ties.method = "first"),
         label = paste0(type, " ",
                        scales::comma(round(tot_med, 0)))) %>%
  mutate(tot_med_bin = cut(tot_med, breaks = 3, labels = c("Low", "Medium", "High"))) %>%
  ungroup()

breaks <- c(0, 100, 1000, 10000, 50000, 100000, 200000, max(tot_vals$tot_med))
labels <- c("0-100", "100-1,000", "1,000-10,000",
            "10,000-50,000", "50,000-100,000",
            "100,000-200,000", ">200,000")

country_rank <- country_rank %>%
  mutate(tot_med_bin = cut(tot_med, 
                           breaks = breaks, 
                           include.lowest = TRUE,
                           labels = labels))

country_rank <- country_rank %>%
  mutate(country = factor(country, levels = tot_vals$country[order(tot_vals$tot_med, decreasing = F)]))

show_col(pal_lancet("lanonc", alpha = 0.6)(9))
color_palette <- c("#00468B99", "#ED000099", "#42B54099", "#0099B499", "#925E9F99",
                   "#FDAF9199", "#AD002A99") 


pdf("Results_figs/DALY_country.pdf", height = 20, width = 12)
p <- ggplot(country_rank, aes(x = factor(rank), y = country, fill = tot_med_bin))+
  geom_tile(color = "black", linewidth = 0.5)+
  geom_text(aes(label = label), size = 2, color = "black", vjust = 0.5) +
  theme_transparent()+
  scale_fill_manual(values = color_palette, name = "Range")+
  theme(axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10))+
  labs(x = "Rank",
       y = "Country")
dev.off()

## country dalys in a log scale y axis bar graph

cum_df_focal <- tot_vals %>% 
  ungroup()%>%
  arrange(desc(tot_med))%>%
  mutate(cum_sum  = cumsum(tot_med),
         cum_perc = cum_sum / sum(tot_med) * 100)

cum_df_focal$country <- factor(cum_df_focal$country, levels = cum_df_focal$country)

color <- pal_lancet("lanonc")(7)


pseudo_log_trans <- function(sigma = 1) {
  scales::trans_new(
    name = "pseudo_log",
    transform = function(x) sign(x) * log1p(abs(x/sigma)),   # Applies a log-like transformation
    inverse = function(x) sign(x) * sigma * (expm1(abs(x))), # Inverse to recover the original scale
    domain = c(0, Inf)
  )
}

cum_graph_focal <- ggplot(cum_df_focal, aes(x = reorder(country, -tot_med), y = tot_med, fill = continent)) +
  geom_bar(stat = "identity", alpha = 0.6) +
  geom_errorbar(aes(ymin = pmax(0.1, tot_lo), ymax = tot_hi), alpha = 0.6)+
  geom_line(aes(y = cum_perc * max(tot_hi) / 100, group = 1), color = "darkgrey", linewidth = 0.6) +
  geom_point(aes(y = cum_perc * max(tot_hi) / 100), color = "darkgrey", size = 1, alpha = 0.6) +
  scale_y_continuous(
    name = "DALYs",
    trans = pseudo_log_trans(sigma = 1000), 
    breaks = c(1e3, 1e4, 1e5, 1e6, 1e7),
    labels = scales::comma,
    sec.axis = sec_axis(~ . / max(cum_df_focal$tot_hi), name = "Cumulative Percentage", labels = scales::percent_format(),
                        breaks = c(0.10, 0.25, 0.5, 0.75, 1))
  ) +
  scale_fill_manual(values = color)+
  #scale_y_log10(scales::comma)+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7), # Tilt and adjust label position
    plot.margin = margin(10, 10, 10, 10), # Adjust plot margins if needed
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    legend.position = "none"         # Remove legend
  )+
  xlab("Country")+
  labs(fill = "Continent")


ggsave("Results_figs/DALY_bargraph_new.jpg", plot = cum_graph, width = 12, height = 9, dpi = 300)
ggsave("Results_figs/DALY_bargraph_focal.pdf", plot = cum_graph, width = 12, height = 9)


write.csv(cum_df, file = "MainData/cumulative_daly_country.csv")

# save into multiple pages
p1 <- ggplot(country_rank, aes(x = factor(rank), y = country, fill = tot_med_bin)) + 
  geom_tile(color = "black", linewidth = 0.5) + 
  geom_text(aes(label = label), size = 2, color = "black", vjust = 0.5) + 
  theme_transparent()+
  scale_fill_manual(values = color_palette, name = "Range") + 
  theme(axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 10),
        plot.title = element_text(size = 10)) +
  labs(x = "Rank", y = "Country")

ggsave("Results_figs/DALY_country_part1.jpg", plot = p1, width = 12, height = 10, dpi = 300)

## country treemap
tot_vals <- tot_vals %>%
  mutate(tot_med_bin = cut(tot_med, 
                           breaks = breaks, 
                           include.lowest = TRUE,
                           labels = labels))
tot_vals <- tot_vals %>%
  group_by(country) %>%
  mutate(rank = rank(-tot_med, ties.method = "first"),
         label = paste0(country, "\n",
                        scales::comma(round(tot_med, 0))))

p <- ggplot(tot_vals, aes(area = tot_med, fill = continent, label = label)) +
  geom_treemap(color = "black", size = 0.5, alpha = 0.6) +  # Add the treemap tiles
  geom_treemap_text(aes(label = label), size = 12, color = "black", place = "centre", grow = FALSE, min.size = 0) +  # Country names
  scale_fill_manual(values = color_palette, name = "Total DALY (95% median)") +  # Apply custom color palette
  theme_minimal()+
  theme(legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
ggsave("Results_figs/DALY_treemap.jpg", plot = p, width = 21, height = 12, dpi = 300)

pdf("Results_figs/DALY_treepmap.pdf", height = 9, width = 16)
ggplot(tot_vals, aes(area = tot_med, fill = continent, label = label)) +
  geom_treemap(color = "black", size = 0.5, alpha = 0.6) +  # Add the treemap tiles
  geom_treemap_text(aes(label = label), size = 9, color = "black", place = "centre", grow = FALSE, min.size = 0) +  # Country names
  scale_fill_manual(values = color_palette, name = "Total DALY (95% median)") +  # Apply custom color palette
  theme_minimal()+
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 6))
dev.off()

## all symptomatic cases by country

symptomatic <- all_countries %>% filter(type == "symptomatic")
symptomatic <- symptomatic %>% group_by(country) %>%
                    summarise(tot_med = sum(tot_med),
                              tot_lo  = sum(tot_lo),
                              tot_hi  = sum(tot_hi))

