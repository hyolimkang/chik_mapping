## 0. Load libraries, shared functions, and all data objects this script
##    depends on (rr_hosp_model, allfoi, all_age_infection, combined_burden,
##    bg_count_dist_wide, hosp_sample, fatal_sample, nh_fatal_sample,
##    lhs_sample_young, etc.)
source("open_data.R")
options(scipen = 999)


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

## ============================================================
## 1. Create a clean country-to-ISO3 lookup
## ============================================================

country_iso3 <- combined_burden |>
  dplyr::transmute(
    country = as.character(country),
    iso3 = as.character(iso3)
  ) |>
  dplyr::filter(
    !is.na(country),
    !is.na(iso3),
    country != "",
    iso3 != ""
  ) |>
  dplyr::distinct()

country_iso3 <- country_iso3 |>
  dplyr::filter(
    country != "France"
  ) |>
  dplyr::bind_rows(
    tibble::tibble(
      country = "France",
      iso3 = "FRA"
    )
  ) |>
  dplyr::distinct()
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

## 3. Attach comorbidity prevalence data to population,
##    harmonising countries using ISO3 -------------------------------
## 3.1 Keep comorbidity-prevalence rows with valid burden age groups
bg_count_dist_0_89 <- bg_count_dist_wide |>
  dplyr::filter(
    !is.na(burden_age_group)
  )


## 3.2 Create a clean country-to-ISO3 lookup
country_iso3 <- combined_burden |>
  dplyr::transmute(
    country = as.character(country),
    iso3 = as.character(iso3)
  ) |>
  dplyr::filter(
    !is.na(country),
    !is.na(iso3),
    country != "",
    iso3 != ""
  ) |>
  dplyr::filter(
    country != "France"
  ) |>
  dplyr::distinct() |>
  dplyr::bind_rows(
    tibble::tibble(
      country = "France",
      iso3 = "FRA"
    )
  ) |>
  dplyr::distinct()


## Check that each country maps to only one ISO3
duplicate_country_iso3 <- country_iso3 |>
  dplyr::count(
    country,
    name = "n_iso3"
  ) |>
  dplyr::filter(
    n_iso3 > 1
  )

if (nrow(duplicate_country_iso3) > 0) {
  print(duplicate_country_iso3)
  
  stop(
    "Some countries are linked to more than one ISO3 code."
  )
}


## 3.3 Attach ISO3 to infection data
##     Do not overwrite all_age_infection
infection_with_iso3 <- all_age_infection |>
  dplyr::select(
    -dplyr::any_of(
      c(
        "iso3",
        "iso3.x",
        "iso3.y"
      )
    )
  ) |>
  dplyr::left_join(
    country_iso3,
    by = "country",
    relationship = "many-to-one"
  )


## Diagnostic 1:
## countries in infection data that could not be assigned an ISO3
missing_iso3 <- infection_with_iso3 |>
  dplyr::filter(
    is.na(iso3)
  ) |>
  dplyr::distinct(
    country
  )

missing_iso3


## Diagnostic 2:
## countries with infection estimates but no comorbidity prevalence
missing_from_bg <- infection_with_iso3 |>
  dplyr::filter(
    !is.na(iso3)
  ) |>
  dplyr::distinct(
    country,
    iso3
  ) |>
  dplyr::anti_join(
    bg_count_dist_0_89 |>
      dplyr::filter(
        !is.na(iso3)
      ) |>
      dplyr::distinct(
        iso3
      ),
    by = "iso3"
  )

missing_from_bg


## 3.4 Attach ISO3 to the five-year population data
##     Do not overwrite pop_5yr_long
pop_5yr_long_with_iso3 <- pop_5yr_long |>
  dplyr::select(
    -dplyr::any_of(
      c(
        "iso3",
        "iso3.x",
        "iso3.y"
      )
    )
  ) |>
  dplyr::left_join(
    infection_with_iso3 |>
      dplyr::distinct(
        country,
        iso3
      ),
    by = "country",
    relationship = "many-to-one"
  )


## 3.5 Create a unique population lookup table
population_lookup <- pop_5yr_long_with_iso3 |>
  dplyr::filter(
    !is.na(iso3),
    !is.na(age_start)
  ) |>
  dplyr::select(
    iso3,
    age_start,
    population
  ) |>
  dplyr::distinct()


## Check that each ISO3-age combination has only one population value
duplicate_population_keys <- population_lookup |>
  dplyr::count(
    iso3,
    age_start,
    name = "n"
  ) |>
  dplyr::filter(
    n > 1
  )

if (nrow(duplicate_population_keys) > 0) {
  print(duplicate_population_keys)
  
  stop(
    "More than one population value exists for some ISO3-age combinations."
  )
}


## 3.6 Attach population to comorbidity-prevalence data
bg_count_with_pop <- bg_count_dist_0_89 |>
  dplyr::left_join(
    population_lookup,
    by = c(
      "iso3",
      "age_start"
    ),
    relationship = "many-to-one"
  )


## Diagnostic 3:
## comorbidity-prevalence rows for which population was not matched
missing_population <- bg_count_with_pop |>
  dplyr::filter(
    is.na(population)
  ) |>
  dplyr::distinct(
    iso3,
    country_name,
    age_start,
    burden_age_group
  )

missing_population


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
    
    prev_comorb_2plus = weighted.mean(
      prev_comorb_2plus,
      w = population,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

## 4. Reshape total infections into 10-year age bands and attach
##    comorbidity prevalence ------------------------------------------------
infection_10yr_long <- infection_with_iso3 |>
  dplyr::filter(!is.na(iso3)) |>
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
      dplyr::transmute(
        group_number = as.integer(group),
        burden_age_group
      ),
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
        prev_comorb_2plus
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
    
    infections_comorb_2plus =
      infections * prev_comorb_2plus
  )

## 5. Split infections by comorbidity count and attach the coarser
##    rr_age_group banding needed to match RR estimates -------------------
infection_comorb <- infection_comorb |>
  dplyr::left_join(
    age_crosswalk_9 |>
      dplyr::transmute(
        group_number = group,
        rr_age_group = burden_age_group
      ),
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
      infections_comorb_2plus = "2+"
    )
  )

rr_hosp_model_unique <- rr_hosp_model |>
  dplyr::distinct(
    rr_age_group,
    comorbidity,
    .keep_all = TRUE
  )

## 6. Draw 1,000 LHS samples of comorbidity-specific hospitalisation RR ----
##    (log-normal fit to point estimate + 95% CI; reference strata fixed at 1)
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
      prev_comorb_2plus
    ),
    names_to = "comorbidity",
    values_to = "prevalence"
  ) |>
  dplyr::mutate(
    comorbidity = dplyr::recode(
      comorbidity,
      prev_comorb_0 = "0",
      prev_comorb_1 = "1",
      prev_comorb_2plus = "2+"
    )
  )

prev_comorb_long <- prev_comorb_long |>
  dplyr::left_join(
    age_crosswalk_9 |>
      dplyr::select(burden_age_group, group_number = group, burden_age_group),
    by = "burden_age_group",
    relationship = "many-to-one"
  )

prev_rr_long <- prev_comorb_long |>
  dplyr::inner_join(
    rr_hosp_sample_long,
    by = c("burden_age_group" = "rr_age_group", "comorbidity"),
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

## 7. Back-calculate a comorbidity-adjusted hospitalisation rate ----------
##    hosp_sample gives a marginal (population-average) rate per age band;
##    weighted_rr reconstructs that average from the comorbidity-specific
##    RR draws, then marginal/weighted_rr backs out an implied "reference"
##    (comorbidity-free) rate that is multiplied by each stratum's own RR.
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

analysis_iso3 <- infection_with_iso3  |>
  dplyr::distinct(iso3) |>
  dplyr::filter(!is.na(iso3))

adjusted_hosp_rate_103 <- adjusted_hosp_rate |>
  dplyr::semi_join(
    analysis_iso3,
    by = "iso3"
  )

## 8. Construct comorbidity-specific fatality rates ----------------------

rr_death_model_unique <- rr_death_model |>
  dplyr::distinct(
    rr_age_group,
    comorbidity,
    .keep_all = TRUE
  )

rr_death_sample <- sample_rr_lhs(
  rr_df = rr_death_model_unique,
  estimate_col = "rr_death",
  lower_col = "rr_death_lo",
  upper_col = "rr_death_hi",
  runs = 1000,
  seed = 456
)

rr_death_sample_long <- rr_death_sample |>
  dplyr::mutate(
    run = dplyr::row_number()
  ) |>
  tidyr::pivot_longer(
    cols = -run,
    names_to = c(
      "rr_age_group",
      "comorbidity"
    ),
    names_sep = "__",
    values_to = "rr_death_draw"
  )

prev_death_rr_long <- prev_comorb_long |>
  dplyr::inner_join(
    rr_death_sample_long,
    by = c(
      "burden_age_group" = "rr_age_group",
      "comorbidity"
    ),
    relationship = "many-to-many"
  )

fatal_sample_long <- fatal_sample |>
  dplyr::mutate(
    run = dplyr::row_number()
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("fatal_"),
    names_to = "group_number",
    values_to = "marginal_fatal_rate_hosp"
  ) |>
  dplyr::mutate(
    group_number = as.integer(
      sub("^fatal_", "", group_number)
    )
  )

nh_fatal_sample_long <- nh_fatal_sample |>
  dplyr::mutate(
    run = dplyr::row_number()
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("nh_fatal_"),
    names_to = "group_number",
    values_to = "marginal_fatal_rate_nonhosp"
  ) |>
  dplyr::mutate(
    group_number = as.integer(
      sub("^nh_fatal_", "", group_number)
    )
  )

prev_death_rr_long <- prev_death_rr_long |>
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
  )

adjusted_fatal_rate <- prev_death_rr_long |>
  dplyr::group_by(
    iso3,
    burden_age_group,
    group_number,
    run
  ) |>
  dplyr::mutate(
    weighted_death_rr = sum(
      prevalence * rr_death_draw,
      na.rm = TRUE
    ),
    
    fatal_rate_hosp_reference =
      marginal_fatal_rate_hosp /
      weighted_death_rr,
    
    fatal_rate_nonhosp_reference =
      marginal_fatal_rate_nonhosp /
      weighted_death_rr,
    
    adjusted_fatal_rate_hosp =
      fatal_rate_hosp_reference *
      rr_death_draw,
    
    adjusted_fatal_rate_nonhosp =
      fatal_rate_nonhosp_reference *
      rr_death_draw
  ) |>
  dplyr::ungroup()

adjusted_fatal_rate_103 <- adjusted_fatal_rate |>
  dplyr::semi_join(
    analysis_iso3,
    by = "iso3"
  )

## 9. Compute the full comorbidity-stratified burden ---------------------

comorbid_burden_full <-
  calculate_comorbid_burden_step3(
    infection_comorb_long =
      infection_comorb_long,
    
    adjusted_hosp_rate =
      adjusted_hosp_rate_103,
    
    adjusted_fatal_rate =
      adjusted_fatal_rate_103,
    
    lhs_sample_young =
      lhs_sample_young,
    
    lhs_old =
      lhs_old,
    
    le_sample =
      le_sample
  )

## 10. Attach subgroup population denominators --------------------------

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

comorbid_burden_analysis <- comorbid_burden_full |>
  dplyr::select(
    -dplyr::any_of(
      c(
        "population_10yr",
        "prevalence",
        "subgroup_population"
      )
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

saveRDS(
  comorbid_burden_full,
  file = "MainData/comorbid_burden_full.rds"
)

saveRDS(
  comorbid_burden_analysis,
  file = "MainData/comorbid_burden_analysis.rds"
)

## Burden calc end --------------------------------------------------------------