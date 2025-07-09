setwd("C:/Users/Hyolim93/OneDrive - London School of Hygiene and Tropical Medicine/chik_mapping")
# open cookie cut data
cookie_cut <- readRDS("MainData/cookie_cut.RDS")
cookie_cut <- cookie_cut %>%
  mutate(presence = if_else(is.na(presence), 0, presence))
cookie_cut_filter <- cookie_cut %>% filter(presence == 1) #unique countries = 109 
# Remove rows with NaN or binary 0 in chik_binmap and chik_bn
cookie_cut_filter <- cookie_cut_filter %>%
  filter(!is.nan(chik_binmap) & !is.nan(chik_bn)) %>%
  filter(chik_binmap != 0) %>%
  filter(!(country %in% c("United States of America", "China")))

unique_countries <- as.data.frame(unique(cookie_cut_filter$country))

cookie_cut_na <- cookie_cut %>%
  filter(!country %in% unique_countries) # 135 countries

# cookie cut map
allfoi = cookie_cut
allfoi_na = cookie_cut_na
foi_cols <- allfoi[,29:128]

# Vectorized quantile calculation
quantile_vectorized <- function(x, probs) {
  sapply(probs, function(p) quantile(x, probs = p, na.rm = TRUE))
}



allfoi <- allfoi %>% 
  mutate(
    foi_mid = ifelse(presence == 0 | chik_binmap == 0 | chik_bn == 0 | country == "United States of America" | country == "China", 0, apply(foi_cols, 1, median, na.rm = TRUE)),
    foi_lo  = ifelse(presence == 0 | chik_binmap == 0 | chik_bn == 0 | country == "United States of America" | country == "China", 0, apply(foi_cols, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))),
    foi_hi  = ifelse(presence == 0 | chik_binmap == 0 | chik_bn == 0 | country == "United States of America" | country == "China", 0, apply(foi_cols, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))),
    foi_sd  = ifelse(presence == 0 | chik_binmap == 0 | chik_bn == 0 | country == "United States of America" | country == "China", 0, apply(foi_cols, 1, sd, na.rm = TRUE))
  )

allfoi <- allfoi %>% 
  mutate(chik_bn = ifelse(country %in% c("United States of America", "China"), 0, chik_bn))

save(allfoi, file = "MainData/allfoi_s2.RData")

allfoi_sf <- st_as_sf(allfoi, coords = c("x", "y"), crs = crs(tsuit))
allfoi_na_sf <- st_as_sf(allfoi_na, coords = c("x", "y"), crs = crs(tsuit))
foi_sf <- st_as_sf(allfoi_sf)
foi_cmask_sf <- st_as_sf(allfoi_na_sf)
foi_cmask_sf$chik_bn <- 0

# make map
tsuit_terra <- rast(tsuit)
foi_rast <- terra::rasterize(foi_sf, tsuit_terra, field = 'foi_mid', filename = "foi_raster.tif", overwrite = TRUE)
foi_cmask <- terra::rasterize(foi_sf, tsuit_terra, field = 'chik_bn', filename = "foi_cmask.tif", overwrite = TRUE)

# convert the raster to SpatRaster
foi_layer1 <- rast("foi_raster.tif")
foi_layer2 <- rast("foi_cmask.tif")
# rasterize SD
sd_rast <- rasterize(foi_sf, tsuit, field = 'foi_sd')
sd_layer1 <- as(sd_rast, "SpatRaster")

# draw map
make_pixel_foi_map(foi_layer1, foi_layer2)

# Extract quantiles and palette
df1 <- as.data.frame(foi_layer1, xy = TRUE) # FOI data
df2 <- as.data.frame(foi_layer2, xy = TRUE) # c_mask data

quantiles <- quantile(df1$last, probs = c(0, 1), na.rm = TRUE)
palette <- rev(paletteer_c("ggthemes::Red-Blue Diverging", 5))
palette <- rev(viridis::rocket(5))
full_palette <- c("lightgrey", "beige", palette)

foi_presence <- classify(foi_layer1, cbind(0, NA))

log_min <- log10(0.001)
log_max <- log10(0.1)
log_breaks <- seq(log_min, log_max, length.out = 6)
#log_breaks <- c(log_min, seq(log_min, log_max, length.out = 6)[-1], log_max)
breaks <- 10^log_breaks
norm_values <- scales::rescale(log_breaks, to = c(0, 1))

# Adjust the breaks to explicitly handle the 0 value
extended_breaks <- c(0, breaks)  # Include a break for 0

# Adjust the normalized values to handle 0 and map the extended breaks
extended_norm_values <- scales::rescale(log_breaks, to = c(0.01, 1))  # Exclude 0 from log scale
extended_norm_values <- c(0, extended_norm_values)  # Add 0 for mapping lightyellow and green
custom_labels <- label_number(accuracy = 0.01)(extended_breaks)
norm_values <- c(0, rescale(log10(extended_breaks[-1]), to = c(0.01, 1)))
# Load world boundaries (country-level)
world_boundaries <- ne_countries(scale = "medium", returnclass = "sf")

# Create the plot
p <- base_map(world) +
  # First, plot df1 where FOI (last == 0) is beige
  geom_tile(data = df1 %>% filter(last == 0), 
            aes(x = x, y = y), fill = "beige", alpha = 1) +
  
  # Then plot df2 where last == 0 (lightgrey)
  geom_tile(data = df2 %>% filter(last == 0), 
            aes(x = x, y = y), fill = "lightgrey", alpha = 1) +
  
  # If you want to handle df2 where last == 1 explicitly, use this:
  #geom_tile(data = df2 %>% filter(last == 1), 
  #          aes(x = x, y = y), fill = NA, alpha = 0) +  # Set alpha = 0 for transparency
  
  # Plot df1 where FOI (last != 0) and map it to a color gradient
  geom_tile(data = df1 %>% filter(last != 0), 
            aes(x = x, y = y, fill = last), alpha = 1) +
  
  # Plot world boundaries
  geom_sf(data = world_boundaries, fill = NA, color = "black", size = 0.5) +
  
  # Set color scale and apply log10 transformation
  scale_fill_gradientn(name = "FOI",
                       colors = full_palette,
                       trans  = "log10", 
                       breaks = extended_breaks,
                       labels = custom_labels,
                       values = norm_values,
                       na.value = "transparent",
                       oob = scales::squish) +
  
  # Set world coordinates
  coord_sf(xlim = c(-180, 180), ylim = c(-50, 75))

ggsave(filename = paste0("Results_figs/FoI_shrink", gsub("-", "_", Sys.Date()), ".jpg"), 
       plot = p, height = 6, width = 12, dpi = 900)


# burden estimation
# filter cookie cut
foi_shrink_df <- allfoi %>% filter(presence == 1)%>% 
                    filter(chik_binmap == 1) %>%
                      filter(chik_bn == 1)

saveRDS(foi_shrink_df, file = "foi_shrink_df.RDS")

foi_shrink_df <- readRDS("MainData/foi_shrink_df.RDS")

#### burden function process (1: mid)
start_time <- Sys.time()

burden_shrink_1 <- age_strat_burden_psa(df = foi_shrink_df, 29,53)

burden_shrink_2 <- age_strat_burden_psa(df = foi_shrink_df, 54,78)

burden_shrink_3 <- age_strat_burden_psa(df = foi_shrink_df, 79,103)

burden_shrink_4 <- age_strat_burden_psa(df = foi_shrink_df, 104,128)

end_time <- Sys.time()

(duration <- end_time - start_time)


all_results <- c(burden_shrink_1, burden_shrink_2, burden_shrink_3, burden_shrink_4)

numbs <- length(all_results)

tot_inf<-do.call(cbind, lapply(all_results, function(df) df$updated_df$total_infection))
colnames(tot_inf) <- paste("tot_inf", seq_len(numbs))
common_cols <- all_results[[1]]$updated_df[, 1:28]
combined_burden <- cbind(common_cols, tot_inf)
tot_inf_cols <- combined_burden[,29:128]
combined_burden$med_inf <- apply(tot_inf_cols, 1, median, na.rm = T)
combined_burden$lo_inf <- apply(tot_inf_cols, 1, function(x) quantile(x, probs = 0.025, na.rm = T))
combined_burden$hi_inf <- apply(tot_inf_cols, 1, function(x) quantile(x, probs = 0.975, na.rm = T))
combined_burden$sd <- apply(tot_inf_cols, 1, function(x) sd(x, na.rm = TRUE))
age_cols <- all_results[[1]]$updated_df[,10:27]
tot_pop <- rowSums(age_cols[, 1:18])
combined_burden <- cbind(combined_burden, tot_pop)
combined_burden <- combined_burden[!is.na(combined_burden$med_inf) & !is.nan(combined_burden$med_inf), ]
combined_burden <- combined_burden[!combined_burden$country %in% c("United States of America" , "China"),]

save(combined_burden, file = "combined_burden_shrink.RData")
save(burden_shrink_1, file = "MainData/burden_shrink_1.RData")
save(burden_shrink_2, file = "MainData/burden_shrink_2.RData")
save(burden_shrink_3, file = "MainData/burden_shrink_3.RData")
save(burden_shrink_4, file = "MainData/burden_shrink_4.RData")

# burden by country (95% mid)
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

infection_focal$scenario <- "Focal"

region_classification <- read_excel("MainData/region_classification.xlsx", 
                                    sheet = "region")

region <- region_classification$continent

infection_focal <- cbind(region, infection_focal)

global_infection_focal <- infection_focal %>% summarise(tot_infec_med = sum(tot_infec_med),
                                                        tot_infec_lo  = sum(tot_infec_lo),
                                                        tot_infec_hi  = sum(tot_infec_hi))

cont_infection_focal <- infection_focal %>% group_by(continent) %>% 
                                             summarise(tot_infec_med = sum(tot_infec_med),
                                                        tot_infec_lo  = sum(tot_infec_lo),
                                                        tot_infec_hi  = sum(tot_infec_hi))

write.csv(global_infection_focal, file = "MainData/global_infection_focal.csv")
write.csv(cont_infection_focal,   file = "MainData/cont_infection_focal.csv")

## comparison
infection_atrisk <- infection_by_iso3[,c(1:7)]
infection_atrisk$scenario <- "At-risk"

country_order <- infection_all %>%
  group_by(country) %>%
  summarise(max_tot_infec_med = max(tot_infec_med)) %>%
  arrange(desc(max_tot_infec_med))

# Rearrange the countries by the largest `tot_infec_med` to the smallest
infection_all <- infection_all %>%
  mutate(country = factor(country, levels = country_order$country))

position_dodge <- position_dodge(width = 0.3)

ggplot(infection_all, aes(x = country, y = tot_infec_med, color = scenario, group = scenario)) +
  geom_point(position = position_dodge, size = 3) +
  geom_errorbar(aes(ymin = tot_infec_lo, ymax = tot_infec_hi), width = 0.2, position = position_dodge) +
  theme_bw() +
  labs(title = "Total Infections by Country and Scenario",
       x = "Country",
       y = "Total Infections")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_y_log10(labels = comma)


### sub-burdens 
lhs_sample <- as.data.frame(lhs_sample)
sub_burden_psa <- age_subinf_lineage_psa(infection_focal, lhs_sample)
results <- postprocess_sub_burden(sub_burden_psa, infection_focal, sub_inf_global_count)

allresults_cont <- results$all_resulst_continent
allresults_global <- results$all_results_global


write.csv(allresults_global, file = "global_burden_s2.csv")


## foi summary
foi_shrink_df <- readRDS("MainData/foi_shrink_df.RDS")
foi_shrink_df <- foi_shrink_df[,-c(7)]
foi_shrink_df <- merge(foi_shrink_df, region_classification, by = "country")

foi_country <- foi_shrink_df %>% filter(foi_mid != 0)
foi_country <- foi_country %>% arrange(desc(foi_hi))

foi_global_mid <- quantile(foi_country$foi_mid, c(0, 0.5, 1))
foi_global_lo  <- quantile(foi_country$foi_lo, c(0.025, 0.5, 0.975))
foi_global_hi  <- quantile(foi_country$foi_hi, c(0.025, 0.5, 0.975))
foi_sd         <- quantile(foi_country$foi_sd, c(0.025, 0.5, 0.975))

foi_global_quantiles <- quantile(foi_country$foi_mid, c(0.25, 0.75))
foi_sd <- quantile(foi_country$foi_sd, c(0.25, 0.75))
# Estimate the IQR
foi_global_iqr <- foi_global_quantiles[2] - foi_global_quantiles[1]
foi_sd_iqr <- foi_sd[2] - foi_sd[1]

foi_continent_mid <- foi_shrink_df %>%
  group_by(continent.y) %>%
  summarise(
    foi_mid = quantile(foi_mid, 0.5),    # Median (50th percentile) within each continent
    foi_lo  = quantile(foi_mid, 0.025),  # 2.5th percentile within each continent
    foi_hi  = quantile(foi_mid, 0.975)   # 97.5th percentile within each continent
  )

foi_ssa <- foi_country %>% filter(continent.y == "Sub-Saharan Africa")  %>% filter(foi_mid != 0)
foi_eu     <- foi_country %>% filter(continent.y == "Europe & Central Asia")  %>% filter(foi_mid != 0)
foi_north_america <- foi_country %>% filter(continent.y == "North America") %>% filter(foi_mid != 0)  
foi_me <- foi_country %>% filter(continent.y == "Middle East & North Africa")  %>% filter(foi_mid != 0)
foi_southasia <- foi_country %>% filter(continent.y == "South Asia")  %>% filter(foi_mid != 0)
foi_latin <- foi_country %>% filter(continent.y == "Latin America & Caribbean")  %>% filter(foi_mid != 0)
foi_eap <- foi_country %>% filter(continent.y == "East Asia & Pacific")  %>% filter(foi_mid != 0)

foi_ssa <- foi_ssa %>% arrange(desc(foi_mid))
foi_eu <- foi_eu %>% arrange(desc(foi_mid))
foi_north_america <- foi_north_america %>% arrange(desc(foi_mid))
foi_me <- foi_me %>% arrange(desc(foi_mid))
foi_southasia <- foi_southasia %>% arrange(desc(foi_mid))
foi_latin <- foi_latin %>% arrange(desc(foi_mid))
foi_eap <- foi_eap %>% arrange(desc(foi_mid))



foi_africa_sum <- foi_ssa %>% group_by(country) %>%
                       summarise(
                         FOI_mid = quantile(foi_mid, 0.5),
                         FOI_lo  = quantile(foi_mid, 0.025),
                         FOI_hi  = quantile(foi_hi, 0.975)
                       )

foi_latin_sum <- foi_latin %>% group_by(country) %>%
  summarise(
    FOI_mid = quantile(foi_mid, 0.5),
    FOI_lo  = quantile(foi_mid, 0.025),
    FOI_hi  = quantile(foi_hi, 0.975)
  )

foi_southasia_sum <- foi_southasia %>% group_by(country) %>%
  summarise(
    FOI_mid = quantile(foi_mid, 0.5),
    FOI_lo  = quantile(foi_mid, 0.025),
    FOI_hi  = quantile(foi_hi, 0.975)
  )

foi_eap_sum <- foi_south_america %>% group_by(country) %>%
  summarise(
    FOI_mid = quantile(foi_mid, 0.5),
    FOI_lo  = quantile(foi_mid, 0.025),
    FOI_hi  = quantile(foi_hi, 0.975)
  )

foi_na_sum <- foi_north_america %>% group_by(country) %>%
  summarise(
    FOI_mid = quantile(foi_mid, 0.5),
    FOI_lo  = quantile(foi_mid, 0.025),
    FOI_hi  = quantile(foi_hi, 0.975)
  )

foi_me_sum <- foi_me %>% group_by(country) %>%
  summarise(
    FOI_mid = quantile(foi_mid, 0.5),
    FOI_lo  = quantile(foi_mid, 0.025),
    FOI_hi  = quantile(foi_hi, 0.975)
  )

foi_eu_sum <- foi_eu %>% group_by(country) %>%
  summarise(
    FOI_mid = quantile(foi_mid, 0.5),
    FOI_lo  = quantile(foi_mid, 0.025),
    FOI_hi  = quantile(foi_hi, 0.975)
  )

foi_africa_mid <- quantile(foi_africa$foi_hi, c(0.025, 0.5, 0.975))
foi_eu_mid     <- quantile(foi_eu$foi_hi, c(0.025, 0.5, 0.975))
foi_na_mid     <- quantile(foi_north_america$foi_hi, c(0.025, 0.5, 0.975))
foi_sa_mid     <- quantile(foi_south_america$foi_hi, c(0.025, 0.5, 0.975))
foi_oce_mid    <- quantile(foi_oceania$foi_hi, c(0.025, 0.5, 0.975))
foi_asia_mid   <- quantile(foi_asia$foi_hi, c(0.025, 0.5, 0.975))
foi_latin_mid   <- quantile(foi_latin$foi_hi, c(0.025, 0.5, 0.975))

sero_age_10 <- function(foi_quant) {
  
  sero_age_10 <- 1 - exp(-foi_quant[2]*10)
  
  return(sero_age_10)
  
}

(sero_age_afr <- sero_age_10(foi_africa_mid))
(sero_age_oce <- sero_age_10(foi_oce_mid))
(sero_age_eu  <- sero_age_10(foi_eu_mid))
(sero_age_asia <- sero_age_10(foi_asia_mid))
(sero_age_latin <- sero_age_10(foi_latin_mid))


foi_country_quantile <- foi_country %>%
  group_by(country) %>%
  summarise(
    Lower = quantile(foi_mid, 0),
    Upper = quantile(foi_mid, 1),
    Mid   = quantile(foi_mid, 0.5),
    country = first(country),
    continent = first(continent)
  ) %>%
  ungroup()

foi_continent_quantile <- foi_country %>%
  group_by(continent) %>%
  summarise(
    Lower = quantile(foi_mid, 0.025),
    Upper = quantile(foi_mid, 0.975),
    Mid   = quantile(foi_mid, 0.5),
    IQR_lo = quantile(foi_mid, 0.25),
    IQR_hi = quantile(foi_mid, 0.75),
    IQR   = IQR_hi - IQR_lo,
    country = first(country),
    continent = first(continent)
  ) %>%
  ungroup()

# total infection
tot_infec_cont <- infection_focal %>% group_by(continent) %>%
                     summarise(tot_med = sum(tot_infec_med),
                               tot_lo  = sum(tot_infec_lo),
                               tot_hi  = sum(tot_infec_hi))

tot_infec_glob <- infection_focal %>% summarise(tot_med = sum(tot_infec_med),
                                                tot_lo  = sum(tot_infec_lo),
                                                tot_hi  = sum(tot_infec_hi))

write.csv(tot_infec_cont, file = "MainData/tot_infec_cont.csv")
