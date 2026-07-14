## 0. Load libraries, shared functions, and all data objects this script
##    depends on (rr_hosp_model, allfoi, all_age_infection, combined_burden,
##    bg_count_dist_wide, hosp_sample, fatal_sample, nh_fatal_sample,
##    lhs_sample_young, etc.)
source("open_data.R")

## 1. Shared ggplot theme for all figures in this script ------------------
theme_lancet_clean <- function(base_size = 11) {
  
  ggplot2::theme_classic(
    base_size = base_size,
    base_family = "Arial"
  ) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(
        size = base_size,
        colour = "black"
      ),
      axis.text = ggplot2::element_text(
        size = base_size - 1,
        colour = "black"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      axis.line = ggplot2::element_line(
        linewidth = 0.45,
        colour = "black"
      ),
      axis.ticks = ggplot2::element_line(
        linewidth = 0.4,
        colour = "black"
      ),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(
        size = base_size - 1
      ),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        face = "bold",
        size = base_size + 1,
        hjust = 0
      ),
      plot.margin = ggplot2::margin(
        8, 10, 8, 8
      )
    )
}

age_crosswalk_9 <- tibble::tibble(
  group = 1:9,
  burden_age_group = c(
    "[0,10)",
    "[10,20)",
    "[20,30)",
    "[30,40)",
    "[40,50)",
    "[50,60)",
    "[60,70)",
    "[70,80)",
    "[80,90)"
  ),
  age_start = seq(0, 80, by = 10),
  age_end = seq(9, 89, by = 10),
  rr_age_group = c(
    "0–19", "0–19",
    "20–39", "20–39",
    "40–59", "40–59",
    "60–79", "60–79",
    "80+"
  )
)

age_crosswalk_9

## 2. Reshape population into 5-year age bands ----------------------------
pop_5yr <- all_age_infection |>
  dplyr::transmute(
    country,
    
    pop_0_4   = `1` + `2`,
    pop_5_9   = `3`,
    pop_10_14 = `4`,
    pop_15_19 = `5`,
    pop_20_24 = `6`,
    pop_25_29 = `7`,
    pop_30_34 = `8`,
    pop_35_39 = `9`,
    pop_40_44 = `10`,
    pop_45_49 = `11`,
    pop_50_54 = `12`,
    pop_55_59 = `13`,
    pop_60_64 = `14`,
    pop_65_69 = `15`,
    pop_70_74 = `16`,
    pop_75_79 = `17`,
    pop_80plus = `18`
  )

pop_5yr_long <- pop_5yr |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("pop_"),
    names_to = "age_band",
    values_to = "population"
  ) |>
  dplyr::mutate(
    age_start = dplyr::case_when(
      age_band == "pop_0_4"    ~ 0,
      age_band == "pop_5_9"    ~ 5,
      age_band == "pop_10_14"  ~ 10,
      age_band == "pop_15_19"  ~ 15,
      age_band == "pop_20_24"  ~ 20,
      age_band == "pop_25_29"  ~ 25,
      age_band == "pop_30_34"  ~ 30,
      age_band == "pop_35_39"  ~ 35,
      age_band == "pop_40_44"  ~ 40,
      age_band == "pop_45_49"  ~ 45,
      age_band == "pop_50_54"  ~ 50,
      age_band == "pop_55_59"  ~ 55,
      age_band == "pop_60_64"  ~ 60,
      age_band == "pop_65_69"  ~ 65,
      age_band == "pop_70_74"  ~ 70,
      age_band == "pop_75_79"  ~ 75,
      age_band == "pop_80plus" ~ 80
    )
  )

## 3. Attach comorbidity prevalence data to population, harmonise on iso3 -
##    (a country_name-keyed join is tried first, then superseded below by
##    an iso3-keyed join once the iso3 crosswalk is available)
bg_count_dist_0_89 <- bg_count_dist_wide |>
  dplyr::filter(!is.na(burden_age_group))

bg_count_with_pop <- bg_count_dist_0_89 |>
  dplyr::left_join(
    pop_5yr_long |>
      dplyr::select(
        country,
        age_start,
        population
      ),
    by = c(
      "country_name" = "country",
      "age_start" = "age_start"
    )
  )

bg_count_dist_analysis <- bg_count_dist_0_89 |>
  dplyr::semi_join(
    pop_5yr_long |>
      dplyr::distinct(country),
    by = c("country_name" = "country")
  )

missing_from_bg <- all_age_infection |>
  dplyr::distinct(country) |>
  dplyr::anti_join(
    bg_count_dist_0_89 |>
      dplyr::distinct(country_name),
    by = c("country" = "country_name")
  )

missing_from_bg

country_iso3 <- combined_burden |>
  dplyr::distinct(
    country,
    iso3
  )

all_age_infection <- all_age_infection |>
  dplyr::left_join(
    country_iso3,
    by = "country",
    relationship = "many-to-one"
  )

all_age_infection <- all_age_infection |>
  dplyr::mutate(
    iso3 = dplyr::if_else(
      country == "France",
      "FRA",
      iso3
    )
  )

pop_5yr_long <- pop_5yr_long |>
  dplyr::left_join(
    all_age_infection |>
      dplyr::distinct(country, iso3),
    by = "country",
    relationship = "many-to-one"
  )


bg_count_with_pop <- bg_count_dist_0_89 |>
  dplyr::left_join(
    pop_5yr_long |>
      dplyr::select(
        iso3,
        age_start,
        population
      ),
    by = c(
      "iso3",
      "age_start"
    )
  )

bg_regional_prev <- bg_count_dist_0_89 |>
  dplyr::group_by(
    region,
    age_start,
    age_end,
    burden_age_group
  ) |>
  dplyr::summarise(
    prev_comorb_0 = mean(prev_comorb_0, na.rm = TRUE),
    prev_comorb_1 = mean(prev_comorb_1, na.rm = TRUE),
    prev_comorb_2 = mean(prev_comorb_2, na.rm = TRUE),
    prev_comorb_3plus = mean(prev_comorb_3plus, na.rm = TRUE),
    .groups = "drop"
  )

bg_prev_10yr <- bg_count_with_pop |>
  dplyr::filter(
    !is.na(population),
    age_start < 80
  ) |>
  dplyr::group_by(
    iso3,
    country_name,
    burden_age_group
  ) |>
  dplyr::summarise(
    population_10yr = sum(population, na.rm = TRUE),
    
    prev_comorb_0 = weighted.mean(
      prev_comorb_0,
      w = population,
      na.rm = TRUE
    ),
    
    prev_comorb_1 = weighted.mean(
      prev_comorb_1,
      w = population,
      na.rm = TRUE
    ),
    
    prev_comorb_2 = weighted.mean(
      prev_comorb_2,
      w = population,
      na.rm = TRUE
    ),
    
    prev_comorb_3plus = weighted.mean(
      prev_comorb_3plus,
      w = population,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

## 4. Reshape total infections into 10-year age bands and attach
##    comorbidity prevalence ------------------------------------------------
infection_10yr_long <- all_age_infection |>
  dplyr::select(
    iso3,
    country,
    continent,
    dplyr::starts_with("tot_infec_med_group")
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("tot_infec_med_group"),
    names_to = "group",
    values_to = "infections"
  ) |>
  dplyr::mutate(
    group_number = as.integer(
      sub("tot_infec_med_group", "", group)
    )
  ) |>
  dplyr::left_join(
    age_crosswalk_9 |>
      dplyr::select(group_number = group, burden_age_group),
    by = "group_number",
    relationship = "many-to-one"
  )

infection_10yr_long_0_79 <- infection_10yr_long |>
  dplyr::filter(group_number <= 8)

infection_with_prev <- infection_10yr_long_0_79 |>
  dplyr::left_join(
    bg_prev_10yr |>
      dplyr::select(
        iso3,
        burden_age_group,
        prev_comorb_0,
        prev_comorb_1,
        prev_comorb_2,
        prev_comorb_3plus
      ),
    by = c(
      "iso3",
      "burden_age_group"
    )
  )

infection_comorb <- infection_with_prev |>
  dplyr::mutate(
    infections_comorb_0 =
      infections * prev_comorb_0,
    
    infections_comorb_1 =
      infections * prev_comorb_1,
    
    infections_comorb_2 =
      infections * prev_comorb_2,
    
    infections_comorb_3plus =
      infections * prev_comorb_3plus
  )

##
infection_comorb <- infection_comorb |>
  dplyr::left_join(
    age_crosswalk_9 |>
      dplyr::select(group_number = group, rr_age_group),
    by = "group_number",
    relationship = "many-to-one"
  )
infection_comorb_long <- infection_comorb |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("infections_comorb_"),
    names_to = "comorbidity",
    values_to = "infections_comorb"
  ) |>
  dplyr::mutate(
    comorbidity = dplyr::recode(
      comorbidity,
      infections_comorb_0 = "0",
      infections_comorb_1 = "1",
      infections_comorb_2 = "2",
      infections_comorb_3plus = "3+"
    )
  )

rr_hosp_model_unique <- rr_hosp_model |>
  dplyr::distinct(
    rr_age_group,
    comorbidity,
    .keep_all = TRUE
  )

infection_comorb_rr <- infection_comorb_long |>
  dplyr::left_join(
    rr_hosp_model_unique,
    by = c(
      "rr_age_group",
      "comorbidity"
    )
  )

###
sample_rr_lhs <- function(rr_df, estimate_col, lower_col, upper_col,
                          runs = 1000, seed = 123) {
  
  set.seed(seed)
  
  n_param <- nrow(rr_df)
  
  U <- lhs::randomLHS(
    n = runs,
    k = n_param
  )
  
  rr_draws <- matrix(
    NA_real_,
    nrow = runs,
    ncol = n_param
  )
  
  for (j in seq_len(n_param)) {
    
    rr_mid <- rr_df[[estimate_col]][j]
    rr_lo  <- rr_df[[lower_col]][j]
    rr_hi  <- rr_df[[upper_col]][j]
    
    if (rr_mid == 1 && rr_lo == 1 && rr_hi == 1) {
      
      # Reference comorbidity group
      rr_draws[, j] <- 1
      
    } else {
      
      meanlog <- log(rr_mid)
      
      sdlog <- (
        log(rr_hi) - log(rr_lo)
      ) / (2 * 1.96)
      
      rr_draws[, j] <- qlnorm(
        p = U[, j],
        meanlog = meanlog,
        sdlog = sdlog
      )
    }
  }
  
  colnames(rr_draws) <- paste(
    rr_df$rr_age_group,
    rr_df$comorbidity,
    sep = "__"
  )
  
  as.data.frame(rr_draws)
}

rr_hosp_sample <- sample_rr_lhs(
  rr_df = rr_hosp_model_unique,
  estimate_col = "rr_hosp",
  lower_col = "rr_hosp_lo",
  upper_col = "rr_hosp_hi",
  runs = 1000,
  seed = 123
)

rr_hosp_sample_long <- rr_hosp_sample |>
  dplyr::mutate(
    run = dplyr::row_number()
  ) |>
  tidyr::pivot_longer(
    cols = -run,
    names_to = c("rr_age_group", "comorbidity"),
    names_sep = "__",
    values_to = "rr_hosp_draw"
  )

prev_comorb_long <- bg_prev_10yr |>
  tidyr::pivot_longer(
    cols = c(
      prev_comorb_0,
      prev_comorb_1,
      prev_comorb_2,
      prev_comorb_3plus
    ),
    names_to = "comorbidity",
    values_to = "prevalence"
  ) |>
  dplyr::mutate(
    comorbidity = dplyr::recode(
      comorbidity,
      prev_comorb_0 = "0",
      prev_comorb_1 = "1",
      prev_comorb_2 = "2",
      prev_comorb_3plus = "3+"
    )
  )

prev_comorb_long <- prev_comorb_long |>
  dplyr::left_join(
    age_crosswalk_9 |>
      dplyr::select(burden_age_group, group_number = group, rr_age_group),
    by = "burden_age_group",
    relationship = "many-to-one"
  )

prev_rr_long <- prev_comorb_long |>
  dplyr::inner_join(
    rr_hosp_sample_long,
    by = c(
      "rr_age_group",
      "comorbidity"
    ),
    relationship = "many-to-many"
  )

hosp_sample_long <- hosp_sample |>
  dplyr::mutate(
    run = dplyr::row_number()
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("hosp_"),
    names_to = "group_number",
    values_to = "marginal_hosp_rate"
  ) |>
  dplyr::mutate(
    group_number = as.integer(
      sub("hosp_", "", group_number)
    )
  )

prev_rr_hosp_long <- prev_rr_long |>
  dplyr::left_join(
    hosp_sample_long,
    by = c(
      "group_number",
      "run"
    ),
    relationship = "many-to-one"
  )

adjusted_hosp_rate <- prev_rr_hosp_long |>
  dplyr::group_by(
    iso3,
    burden_age_group,
    group_number,
    run
  ) |>
  dplyr::mutate(
    weighted_rr = sum(
      prevalence * rr_hosp_draw,
      na.rm = TRUE
    ),
    
    hosp_rate_reference =
      marginal_hosp_rate / weighted_rr,

    adjusted_hosp_rate =
      dplyr::if_else(
        is.finite(hosp_rate_reference * rr_hosp_draw),
        hosp_rate_reference * rr_hosp_draw,
        NA_real_
      )
  ) |>
  dplyr::ungroup()

analysis_iso3 <- all_age_infection |>
  dplyr::distinct(iso3) |>
  dplyr::filter(!is.na(iso3))

prev_comorb_long_103 <- prev_comorb_long |>
  dplyr::semi_join(
    analysis_iso3,
    by = "iso3"
  )

adjusted_hosp_rate_103 <- adjusted_hosp_rate |>
  dplyr::semi_join(
    analysis_iso3,
    by = "iso3"
  )

comorbid_burden_step2 <-
  calculate_comorbid_burden_step2(
    infection_comorb_long =
      infection_comorb_long,
    
    adjusted_hosp_rate =
      adjusted_hosp_rate_103,
    
    lhs_sample =
      lhs_sample_young
  )


###
age_levels <- c(
  "[0,10)",
  "[10,20)",
  "[20,30)",
  "[30,40)",
  "[40,50)",
  "[50,60)",
  "[60,70)",
  "[70,80)"
)

comorb_levels <- c("0", "1", "2", "3+")

population_lookup <- prev_comorb_long |>
  dplyr::select(
    iso3,
    burden_age_group,
    group_number,
    comorbidity,
    population_10yr,
    prevalence
  ) |>
  dplyr::distinct()

comorbid_burden_plot <- comorbid_burden_step2 |>
  dplyr::select(
    -dplyr::any_of(c("population_10yr", "prevalence"))
  ) |>
  dplyr::left_join(
    population_lookup,
    by = c(
      "iso3",
      "burden_age_group",
      "group_number",
      "comorbidity"
    ),
    relationship = "many-to-one"
  )

plot_by_run <- comorbid_burden_plot |>
  dplyr::filter(
    !is.na(adjusted_hosp_rate),
    !is.na(population_10yr),
    !is.na(prevalence),
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      burden_age_group,
      levels = age_levels
    ),
    comorbidity = factor(
      comorbidity,
      levels = comorb_levels
    ),
    subgroup_population =
      population_10yr * prevalence
  ) |>
  dplyr::group_by(
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    symptomatic = sum(
      symptomatic,
      na.rm = TRUE
    ),
    hospitalised = sum(
      hospitalised,
      na.rm = TRUE
    ),
    subgroup_population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    symptomatic_rate_100k =
      symptomatic /
      subgroup_population *
      100000,
    
    hospitalisation_rate_10k =
      hospitalised /
      subgroup_population *
      10000
  )


plot_summary <- plot_by_run |>
  dplyr::group_by(
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    symptomatic_median =
      median(
        symptomatic_rate_100k,
        na.rm = TRUE
      ),
    
    symptomatic_lower =
      quantile(
        symptomatic_rate_100k,
        0.025,
        na.rm = TRUE
      ),
    
    symptomatic_upper =
      quantile(
        symptomatic_rate_100k,
        0.975,
        na.rm = TRUE
      ),
    
    hospitalisation_median =
      median(
        hospitalisation_rate_10k,
        na.rm = TRUE
      ),
    
    hospitalisation_lower =
      quantile(
        hospitalisation_rate_10k,
        0.025,
        na.rm = TRUE
      ),
    
    hospitalisation_upper =
      quantile(
        hospitalisation_rate_10k,
        0.975,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

comorb_colours <- c(
  "0"  = "#ADB6B6",
  "1"  = "#0099B4",
  "2"  = "#00468B",
  "3+" = "#ED0000"
)

p_hospitalisation <- ggplot(
  plot_summary,
  aes(
    x = burden_age_group,
    y = hospitalisation_median,
    colour = comorbidity,
    group = comorbidity
  )
) +
  geom_errorbar(
    aes(
      ymin = hospitalisation_lower,
      ymax = hospitalisation_upper
    ),
    width = 0.12,
    linewidth = 0.45,
    alpha = 0.8
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.5,
    stroke = 0.3
  ) +
  scale_colour_manual(
    values = comorb_colours,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1,
      big.mark = ","
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    title = "Hospitalisation",
    x = "Age group, years",
    y = "Hospitalised cases per 10,000\nsubgroup population",
    colour = "Number of comorbidities"
  ) +
  theme_lancet_clean()


plot_region_by_run <- comorbid_burden_plot |>
  dplyr::filter(
    !is.na(adjusted_hosp_rate),
    !is.na(population_10yr),
    !is.na(prevalence),
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      burden_age_group,
      levels = age_levels
    ),
    comorbidity = factor(
      comorbidity,
      levels = comorb_levels
    ),
    subgroup_population =
      population_10yr * prevalence
  ) |>
  dplyr::group_by(
    continent,
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    symptomatic =
      sum(symptomatic, na.rm = TRUE),
    
    hospitalised =
      sum(hospitalised, na.rm = TRUE),
    
    subgroup_population =
      sum(subgroup_population, na.rm = TRUE),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    hospitalisation_rate_10k =
      hospitalised /
      subgroup_population *
      10000
  )

plot_region_summary <- plot_region_by_run |>
  dplyr::group_by(
    continent,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    median =
      median(
        hospitalisation_rate_10k,
        na.rm = TRUE
      ),
    
    lower =
      quantile(
        hospitalisation_rate_10k,
        0.025,
        na.rm = TRUE
      ),
    
    upper =
      quantile(
        hospitalisation_rate_10k,
        0.975,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


p_region <- ggplot(
  plot_region_summary,
  aes(
    x = burden_age_group,
    y = median,
    colour = comorbidity,
    group = comorbidity
  )
) +
  geom_errorbar(
    aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.1,
    linewidth = 0.3,
    alpha = 0.65
  ) +
  geom_line(
    linewidth = 0.75
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_wrap(
    ~ continent
  ) +
  scale_colour_manual(
    values = comorb_colours,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  labs(
    x = "Age group, years",
    y = "Hospitalised cases per 10,000\nsubgroup population",
    colour = "Number of comorbidities"
  ) +
  theme_lancet_clean(base_size = 10) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_region

country_hosp_by_run <- comorbid_burden_plot |>
  dplyr::filter(
    !is.na(adjusted_hosp_rate),
    !is.na(population_10yr),
    !is.na(prevalence),
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      burden_age_group,
      levels = age_levels
    ),
    comorbidity = factor(
      comorbidity,
      levels = comorb_levels
    ),
    subgroup_population =
      population_10yr * prevalence
  ) |>
  dplyr::group_by(
    iso3,
    country,
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    hospitalised = sum(
      hospitalised,
      na.rm = TRUE
    ),
    subgroup_population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    hospitalisation_rate_10k =
      hospitalised /
      subgroup_population *
      10000
  )

country_hosp_summary <- country_hosp_by_run |>
  dplyr::group_by(
    iso3,
    country,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    median_rate = median(
      hospitalisation_rate_10k,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

global_hosp_by_run <- comorbid_burden_plot |>
  dplyr::filter(
    !is.na(adjusted_hosp_rate),
    !is.na(population_10yr),
    !is.na(prevalence),
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      burden_age_group,
      levels = age_levels
    ),
    comorbidity = factor(
      comorbidity,
      levels = comorb_levels
    ),
    subgroup_population =
      population_10yr * prevalence
  ) |>
  dplyr::group_by(
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    hospitalised = sum(
      hospitalised,
      na.rm = TRUE
    ),
    subgroup_population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    hospitalisation_rate_10k =
      hospitalised /
      subgroup_population *
      10000
  )

global_hosp_summary <- global_hosp_by_run |>
  dplyr::group_by(
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    median_rate = median(
      hospitalisation_rate_10k,
      na.rm = TRUE
    ),
    lower = quantile(
      hospitalisation_rate_10k,
      0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      hospitalisation_rate_10k,
      0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


p_spaghetti_hosp <- ggplot2::ggplot() +
  
  # Individual countries
  ggplot2::geom_line(
    data = country_hosp_summary,
    mapping = ggplot2::aes(
      x = burden_age_group,
      y = median_rate,
      group = iso3
    ),
    colour = "grey65",
    linewidth = 0.3,
    alpha = 0.35
  ) +
  
  # Global 95% uncertainty interval
  ggplot2::geom_ribbon(
    data = global_hosp_summary,
    mapping = ggplot2::aes(
      x = burden_age_group,
      ymin = lower,
      ymax = upper,
      group = 1
    ),
    fill = "grey45",
    alpha = 0.18
  ) +
  
  # Population-weighted global median
  ggplot2::geom_line(
    data = global_hosp_summary,
    mapping = ggplot2::aes(
      x = burden_age_group,
      y = median_rate,
      group = 1
    ),
    colour = "black",
    linewidth = 1
  ) +
  
  ggplot2::facet_wrap(
    ~ comorbidity,
    ncol = 2,
    scales = "free_y",
    labeller = ggplot2::labeller(
      comorbidity = c(
        "0" = "No comorbidity",
        "1" = "One comorbidity",
        "2" = "Two comorbidities",
        "3+" = "Three or more comorbidities"
      )
    )
  ) +
  
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.06)
    )
  ) +
  
  ggplot2::labs(
    x = "Age group, years",
    y = "Hospitalised cases per 10,000\nsubgroup population"
  ) +
  
  theme_lancet_clean(
    base_size = 10
  ) +
  
  ggplot2::theme(
    legend.position = "none",
    
    strip.background =
      ggplot2::element_blank(),
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        size = 10,
        hjust = 0
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
    
    panel.spacing =
      grid::unit(10, "pt")
  )

p_spaghetti_hosp


region_colours <- c(
  "East Asia & Pacific" = "#00468B",
  "Europe & Central Asia" = "#0099B4",
  "Latin America & Caribbean" = "#925E9F",
  "Middle East & North Africa" = "#FDAF91",
  "North America" = "#42B540",
  "South Asia" = "#ED0000",
  "Sub-Saharan Africa" = "#7E6148"
)


#######
fatal_sample_long <- fatal_sample |>
  dplyr::mutate(
    run = dplyr::row_number()
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("fatal_"),
    names_to = "group_number",
    values_to = "fatal_rate_hosp"
  ) |>
  dplyr::mutate(
    group_number = as.integer(
      sub("fatal_", "", group_number)
    )
  )

nh_fatal_sample_long <- nh_fatal_sample |>
  dplyr::mutate(
    run = dplyr::row_number()
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("nh_fatal_"),
    names_to = "group_number",
    values_to = "fatal_rate_nonhosp"
  ) |>
  dplyr::mutate(
    group_number = as.integer(
      sub("nh_fatal_", "", group_number)
    )
  )

comorbid_burden_fatal <- comorbid_burden_step2 |>
  dplyr::left_join(
    fatal_sample_long,
    by = c(
      "group_number",
      "run"
    ),
    relationship = "many-to-one"
  ) |>
  dplyr::left_join(
    nh_fatal_sample_long,
    by = c(
      "group_number",
      "run"
    ),
    relationship = "many-to-one"
  ) |>
  dplyr::mutate(
    fatal_hospitalised =
      hospitalised * fatal_rate_hosp,
    
    fatal_nonhospitalised =
      nonhospitalised * fatal_rate_nonhosp,
    
    fatal =
      fatal_hospitalised +
      fatal_nonhospitalised
  )

fatal_plot_data <- comorbid_burden_fatal |>
  dplyr::select(
    -dplyr::any_of(
      c("population_10yr", "prevalence")
    )
  ) |>
  dplyr::left_join(
    population_lookup,
    by = c(
      "iso3",
      "burden_age_group",
      "group_number",
      "comorbidity"
    ),
    relationship = "many-to-one"
  ) |>
  dplyr::mutate(
    subgroup_population =
      population_10yr * prevalence
  )

fatal_by_run <- fatal_plot_data |>
  dplyr::filter(
    !is.na(fatal),
    !is.na(subgroup_population),
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      burden_age_group,
      levels = age_levels
    ),
    comorbidity = factor(
      comorbidity,
      levels = comorb_levels
    )
  ) |>
  dplyr::group_by(
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    fatal = sum(
      fatal,
      na.rm = TRUE
    ),
    symptomatic = sum(
      symptomatic,
      na.rm = TRUE
    ),
    subgroup_population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    fatal_rate_100k =
      fatal /
      subgroup_population *
      100000,
    
    case_fatality_percent =
      fatal /
      symptomatic *
      100
  )

fatal_plot_summary <- fatal_by_run |>
  dplyr::group_by(
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    median = median(
      fatal_rate_100k,
      na.rm = TRUE
    ),
    lower = quantile(
      fatal_rate_100k,
      0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      fatal_rate_100k,
      0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

p_fatal_burden <- ggplot2::ggplot(
  fatal_plot_summary,
  ggplot2::aes(
    x = burden_age_group,
    y = median,
    colour = comorbidity,
    fill = comorbidity,
    group = comorbidity
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = lower,
      ymax = upper
    ),
    colour = NA,
    alpha = 0.10
  ) +
  ggplot2::geom_line(
    linewidth = 0.95
  ) +
  ggplot2::geom_point(
    size = 2.3
  ) +
  ggplot2::scale_colour_manual(
    values = comorb_colours,
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = comorb_colours,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.01,
      big.mark = ","
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    title = "Fatal disease burden",
    subtitle = "Population-level deaths incorporating infection incidence and age-specific fatality",
    x = "Age group, years",
    y = "Deaths per 100,000\nsubgroup population",
    colour = "Number of comorbidities",
    fill = "Number of comorbidities"
  ) +
  theme_lancet_clean(
    base_size = 11
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    plot.subtitle = ggplot2::element_text(
      size = 9,
      colour = "grey30"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

p_fatal_burden


###
library(data.table)

foi_dt <- as.data.table(allfoi)

foi_cols <- grep(
  "^foi[0-9]+$",
  names(foi_dt),
  value = TRUE
)

length(foi_cols)

foi_dt_valid <- foi_dt[
  !is.na(iso3) &
    is.finite(tot) &
    tot > 0
]

country_foi_draws <- foi_dt_valid[
  ,
  c(
    list(
      country = first(country),
      continent = first(continent),
      region = first(region),
      population = sum(tot, na.rm = TRUE)
    ),
    lapply(
      .SD,
      function(x) {
        weighted.mean(
          x,
          w = tot,
          na.rm = TRUE
        )
      }
    )
  ),
  by = iso3,
  .SDcols = foi_cols
]

country_foi <- country_foi_draws |>
  dplyr::as_tibble() |>
  dplyr::rowwise() |>
  dplyr::mutate(
    foi = median(
      c_across(dplyr::all_of(foi_cols)),
      na.rm = TRUE
    ),
    
    foi_lower = quantile(
      c_across(dplyr::all_of(foi_cols)),
      0.025,
      na.rm = TRUE
    ),
    
    foi_upper = quantile(
      c_across(dplyr::all_of(foi_cols)),
      0.975,
      na.rm = TRUE
    )
  ) |>
  dplyr::ungroup() |>
  dplyr::select(
    iso3,
    country,
    continent,
    region,
    population,
    foi,
    foi_lower,
    foi_upper
  )


surface_data_by_run <- fatal_plot_data |>
  dplyr::filter(
    !is.na(fatal),
    !is.na(population_10yr),
    !is.na(prevalence),
    group_number <= 8
  ) |>
  dplyr::mutate(
    broad_age_group = dplyr::case_when(
      group_number %in% 1:4 ~ "0–39 years",
      group_number %in% 5:8 ~ "40–79 years",
      TRUE ~ NA_character_
    ),
    
    subgroup_population =
      population_10yr * prevalence,
    
    multimorbid_population =
      dplyr::if_else(
        as.character(comorbidity) %in% c("2", "3+"),
        subgroup_population,
        0
      )
  ) |>
  dplyr::filter(
    !is.na(broad_age_group)
  ) |>
  dplyr::group_by(
    iso3,
    country,
    continent,
    run,
    broad_age_group
  ) |>
  dplyr::summarise(
    fatal = sum(
      fatal,
      na.rm = TRUE
    ),
    
    total_population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    
    multimorbid_population = sum(
      multimorbid_population,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    fatal_rate_100k =
      fatal / total_population * 100000,
    
    multimorbidity_prevalence =
      multimorbid_population /
      total_population *
      100
  )


surface_data_by_run <- surface_data_by_run |>
  dplyr::left_join(
    country_foi |>
      dplyr::select(
        iso3,
        foi,
        foi_lower,
        foi_upper
      ),
    by = "iso3",
    relationship = "many-to-one"
  )

surface_country_summary <- surface_data_by_run |>
  dplyr::filter(
    is.finite(foi),
    is.finite(fatal_rate_100k),
    is.finite(multimorbidity_prevalence)
  ) |>
  dplyr::group_by(
    iso3,
    country,
    continent,
    broad_age_group
  ) |>
  dplyr::summarise(
    foi = dplyr::first(foi),
    
    foi_lower =
      dplyr::first(foi_lower),
    
    foi_upper =
      dplyr::first(foi_upper),
    
    multimorbidity_prevalence =
      median(
        multimorbidity_prevalence,
        na.rm = TRUE
      ),
    
    fatal_rate_100k =
      median(
        fatal_rate_100k,
        na.rm = TRUE
      ),
    
    fatal_rate_lower =
      quantile(
        fatal_rate_100k,
        0.025,
        na.rm = TRUE
      ),
    
    fatal_rate_upper =
      quantile(
        fatal_rate_100k,
        0.975,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    broad_age_group = factor(
      broad_age_group,
      levels = c(
        "0–39 years",
        "40–79 years"
      )
    )
  )