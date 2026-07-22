## 0. Load libraries, shared functions, and all data objects this script
##    depends on (rr_hosp_model, allfoi, all_age_infection, combined_burden,
##    bg_count_dist_wide, hosp_sample, fatal_sample, nh_fatal_sample,
##    lhs_sample_young, etc.)

comorbid_burden_analysis <-
  readRDS(
    "MainData/comorbid_burden_analysis.rds"
  )

## -----------------------------------------------------------------------
## Shared comorbidity levels, labels, and colours
## -----------------------------------------------------------------------

comorb_levels <- c(
  "0",
  "1",
  "2+"
)

comorb_labels <- c(
  "0" = "No comorbidity",
  "1" = "One comorbidity",
  "2+" = "Two or more comorbidities"
)

comorb_colours <- c(
  "0" = "#4C78A8",
  "1" = "#F2A541",
  "2+" = "#C95D63"
)

region_colours <- c(
  "East Asia & Pacific" = "#00468B",
  "Europe & Central Asia" = "#0099B4",
  "Latin America & Caribbean" = "#925E9F",
  "Middle East & North Africa" = "#FDAF91",
  "North America" = "#42B540",
  "South Asia" = "#ED0000",
  "Sub-Saharan Africa" = "#7E6148"
)

outcome_levels <- c(
  "Population",
  "Infections",
  "Symptomatic",
  "Hospitalised",
  "Fatal",
  "DALY"
)

outcome_labels <- c(
  "Population" = "Population",
  "Infections" = "Infections",
  "Symptomatic" = "Symptomatic cases",
  "Hospitalised" = "Hospitalisations",
  "Fatal" = "Deaths",
  "DALY" = "DALYs"
)
## =======================================================================
## Hospitalisation burden by age and comorbidity
## Input: comorbid_burden_analysis
## =======================================================================

plot_by_run <- comorbid_burden_analysis |>
  dplyr::filter(
    !is.na(adjusted_hosp_rate),
    !is.na(subgroup_population),
    subgroup_population > 0,
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
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
    symptomatic_median = median(
      symptomatic_rate_100k,
      na.rm = TRUE
    ),
    
    symptomatic_lower = quantile(
      symptomatic_rate_100k,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    symptomatic_upper = quantile(
      symptomatic_rate_100k,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    hospitalisation_median = median(
      hospitalisation_rate_10k,
      na.rm = TRUE
    ),
    
    hospitalisation_lower = quantile(
      hospitalisation_rate_10k,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    hospitalisation_upper = quantile(
      hospitalisation_rate_10k,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


p_hospitalisation <- ggplot2::ggplot(
  plot_summary,
  ggplot2::aes(
    x = burden_age_group,
    y = hospitalisation_median,
    colour = comorbidity,
    group = comorbidity
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = hospitalisation_lower,
      ymax = hospitalisation_upper
    ),
    width = 0.12,
    linewidth = 0.45,
    alpha = 0.8
  ) +
  ggplot2::geom_line(
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    size = 2.5,
    stroke = 0.3
  ) +
  ggplot2::scale_colour_manual(
    values = comorb_colours,
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1,
      big.mark = ","
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    title = "Hospitalisation",
    x = "Age group, years",
    y = "Hospitalised cases per 10,000\nsubgroup population",
    colour = "Number of comorbidities"
  ) +
  theme_lancet_clean()

p_hospitalisation

## =======================================================================
## Hospitalisation attribution and regional/country plots
## Input: comorbid_burden_analysis
## =======================================================================


## 1. Hospitalisation burden attribution ---------------------------------

# Calculate run-specific hospitalisation attribution
attribution_by_run <- plot_by_run |>
  dplyr::group_by(
    run,
    burden_age_group
  ) |>
  dplyr::mutate(
    total_population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    
    total_hospitalised = sum(
      hospitalised,
      na.rm = TRUE
    ),
    
    overall_rate_10k =
      total_hospitalised /
      total_population *
      10000,
    
    contribution_share = dplyr::if_else(
      total_hospitalised > 0,
      hospitalised / total_hospitalised,
      NA_real_
    )
  ) |>
  dplyr::ungroup()


# Summarise the overall age-specific hospitalisation rate
overall_hosp_summary <- attribution_by_run |>
  dplyr::distinct(
    run,
    burden_age_group,
    overall_rate_10k
  ) |>
  dplyr::group_by(
    burden_age_group
  ) |>
  dplyr::summarise(
    overall_median = median(
      overall_rate_10k,
      na.rm = TRUE
    ),
    
    overall_lower = quantile(
      overall_rate_10k,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    overall_upper = quantile(
      overall_rate_10k,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# Summarise each comorbidity group's share of hospitalisations
contribution_share_summary <- attribution_by_run |>
  dplyr::group_by(
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    share_median = median(
      contribution_share,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  dplyr::group_by(
    burden_age_group
  ) |>
  dplyr::mutate(
    share_median_normalised =
      share_median /
      sum(
        share_median,
        na.rm = TRUE
      )
  ) |>
  dplyr::ungroup()


# Convert contribution shares into attributed rates
attribution_plot_summary <- contribution_share_summary |>
  dplyr::left_join(
    overall_hosp_summary,
    by = "burden_age_group"
  ) |>
  dplyr::mutate(
    attributed_rate_10k =
      overall_median *
      share_median_normalised,
    
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
  )


# Prepare overall estimates for plotting
overall_hosp_plot <- overall_hosp_summary |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    )
  )


# Plot hospitalisation burden attribution
p_hospitalisation_attribution <- ggplot2::ggplot(
  attribution_plot_summary,
  ggplot2::aes(
    x = burden_age_group,
    y = attributed_rate_10k,
    fill = comorbidity
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.25,
    position = ggplot2::position_stack(
      reverse = TRUE
    )
  ) +
  ggplot2::geom_errorbar(
    data = overall_hosp_plot,
    mapping = ggplot2::aes(
      x = burden_age_group,
      ymin = overall_lower,
      ymax = overall_upper
    ),
    inherit.aes = FALSE,
    width = 0.14,
    linewidth = 0.55,
    colour = "black"
  ) +
  ggplot2::geom_point(
    data = overall_hosp_plot,
    mapping = ggplot2::aes(
      x = burden_age_group,
      y = overall_median
    ),
    inherit.aes = FALSE,
    shape = 21,
    size = 2.4,
    stroke = 0.55,
    fill = "white",
    colour = "black"
  ) +
  ggplot2::scale_fill_manual(
    values = comorb_colours,
    breaks = comorb_levels,
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1,
      big.mark = ","
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    title = "Hospitalisation burden by comorbidity status",
    x = "Age group, years",
    y = paste0(
      "Attributed hospitalised cases per 10,000\n",
      "total age-group population"
    ),
    fill = "Number of comorbidities"
  ) +
  theme_lancet_clean()

p_hospitalisation_attribution



## =======================================================================
## 2. Hospitalisation plot: region-faceted
## =======================================================================

plot_region_by_run <- comorbid_burden_analysis |>
  dplyr::filter(
    !is.na(hospitalised),
    !is.na(subgroup_population),
    subgroup_population > 0,
    !is.na(continent),
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
  ) |>
  dplyr::group_by(
    continent,
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
    median = median(
      hospitalisation_rate_10k,
      na.rm = TRUE
    ),
    
    lower = quantile(
      hospitalisation_rate_10k,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    upper = quantile(
      hospitalisation_rate_10k,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


p_region <- ggplot2::ggplot(
  plot_region_summary,
  ggplot2::aes(
    x = burden_age_group,
    y = median,
    colour = comorbidity,
    group = comorbidity
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.1,
    linewidth = 0.3,
    alpha = 0.65
  ) +
  ggplot2::geom_line(
    linewidth = 0.75
  ) +
  ggplot2::geom_point(
    size = 1.8
  ) +
  ggplot2::facet_wrap(
    ~ continent
  ) +
  ggplot2::scale_colour_manual(
    values = comorb_colours,
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    x = "Age group, years",
    y = "Hospitalised cases per 10,000\nsubgroup population",
    colour = "Number of comorbidities"
  ) +
  theme_lancet_clean(
    base_size = 10
  ) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(
      face = "bold"
    ),
    
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

p_region



## =======================================================================
## 3. Hospitalisation plot:
##    per-country spaghetti plus global median and uncertainty
## =======================================================================

country_hosp_by_run <- comorbid_burden_analysis |>
  dplyr::filter(
    !is.na(hospitalised),
    !is.na(subgroup_population),
    subgroup_population > 0,
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
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
  ) |>
  dplyr::filter(
    is.finite(median_rate)
  )


global_hosp_by_run <- comorbid_burden_analysis |>
  dplyr::filter(
    !is.na(hospitalised),
    !is.na(subgroup_population),
    subgroup_population > 0,
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
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
      probs = 0.025,
      na.rm = TRUE
    ),
    
    upper = quantile(
      hospitalisation_rate_10k,
      probs = 0.975,
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
      comorbidity = comorb_labels
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
      grid::unit(
        10,
        "pt"
      )
  )

p_spaghetti_hosp



## =======================================================================
## Fatal burden figures
## Input: MainData/comorbid_burden_analysis.rds
## =======================================================================


##  Check required variables -------------------------------------------

required_fatal_columns <- c(
  "iso3",
  "country",
  "continent",
  "run",
  "burden_age_group",
  "group_number",
  "comorbidity",
  "symptomatic",
  "fatal",
  "subgroup_population"
)

missing_fatal_columns <- setdiff(
  required_fatal_columns,
  names(comorbid_burden_analysis)
)

if (length(missing_fatal_columns) > 0) {
  stop(
    paste0(
      "The following required variables are missing from ",
      "comorbid_burden_analysis: ",
      paste(
        missing_fatal_columns,
        collapse = ", "
      )
    )
  )
}


## 4. Prepare common fatal burden dataset --------------------------------

fatal_graph_data <- comorbid_burden_analysis |>
  dplyr::filter(
    group_number <= 8
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
  )


## =======================================================================
## Figure 1. Global fatal burden by age and comorbidity
## =======================================================================

fatal_by_run <- fatal_graph_data |>
  dplyr::filter(
    !is.na(fatal),
    !is.na(subgroup_population),
    subgroup_population > 0
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
      dplyr::if_else(
        symptomatic > 0,
        fatal / symptomatic * 100,
        NA_real_
      )
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
      probs = 0.025,
      na.rm = TRUE
    ),
    
    upper = quantile(
      fatal_rate_100k,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    median_cfr = median(
      case_fatality_percent,
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
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = comorb_colours,
    labels = comorb_labels,
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
    subtitle = paste0(
      "Population-level deaths incorporating infection incidence ",
      "and age-specific fatality"
    ),
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


## =======================================================================
## Figure 2. Regional distribution of country-specific fatal burden
## =======================================================================

country_fatal_by_run <- fatal_graph_data |>
  dplyr::filter(
    !is.na(fatal),
    !is.na(subgroup_population),
    subgroup_population > 0,
    !is.na(continent)
  ) |>
  dplyr::group_by(
    iso3,
    country,
    continent,
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    fatal = sum(
      fatal,
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
      100000
  )


country_fatal_summary <- country_fatal_by_run |>
  dplyr::group_by(
    iso3,
    country,
    continent,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    median_rate = median(
      fatal_rate_100k,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  dplyr::filter(
    is.finite(median_rate)
  )


regional_country_fatal_distribution <-
  country_fatal_summary |>
  dplyr::group_by(
    continent,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    region_median = median(
      median_rate,
      na.rm = TRUE
    ),
    
    country_p10 = quantile(
      median_rate,
      probs = 0.10,
      na.rm = TRUE
    ),
    
    country_p25 = quantile(
      median_rate,
      probs = 0.25,
      na.rm = TRUE
    ),
    
    country_p75 = quantile(
      median_rate,
      probs = 0.75,
      na.rm = TRUE
    ),
    
    country_p90 = quantile(
      median_rate,
      probs = 0.90,
      na.rm = TRUE
    ),
    
    n_countries = dplyr::n_distinct(
      iso3
    ),
    
    .groups = "drop"
  )


global_fatal_by_run <- fatal_graph_data |>
  dplyr::filter(
    !is.na(fatal),
    !is.na(subgroup_population),
    subgroup_population > 0
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
      100000
  )


global_fatal_summary <- global_fatal_by_run |>
  dplyr::group_by(
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    median_rate = median(
      fatal_rate_100k,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


p_regional_fatal_fan <- ggplot2::ggplot(
  regional_country_fatal_distribution,
  ggplot2::aes(
    x = burden_age_group,
    group = continent
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = country_p10,
      ymax = country_p90,
      fill = continent
    ),
    alpha = 0.08,
    colour = NA
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = country_p25,
      ymax = country_p75,
      fill = continent
    ),
    alpha = 0.18,
    colour = NA
  ) +
  ggplot2::geom_line(
    ggplot2::aes(
      y = region_median,
      colour = continent
    ),
    linewidth = 0.9
  ) +
  ggplot2::geom_line(
    data = global_fatal_summary,
    mapping = ggplot2::aes(
      x = burden_age_group,
      y = median_rate,
      group = 1
    ),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 0.9,
    linetype = "22"
  ) +
  ggplot2::facet_wrap(
    ~ comorbidity,
    ncol = 2,
    scales = "fixed",
    labeller = ggplot2::labeller(
      comorbidity = comorb_labels
    )
  ) +
  ggplot2::scale_colour_manual(
    values = region_colours,
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = region_colours,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.07)
    )
  ) +
  ggplot2::labs(
    x = "Age group, years",
    y = "Deaths per 100,000\nsubgroup population",
    colour = "Region",
    fill = "Region"
  ) +
  theme_lancet_clean(
    base_size = 10
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    
    strip.text = ggplot2::element_text(
      face = "bold",
      hjust = 0
    ),
    
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    
    legend.position = "bottom",
    legend.box = "vertical",
    
    panel.spacing = grid::unit(
      10,
      "pt"
    )
  )

p_regional_fatal_fan


## =======================================================================
## Figure 3. Conditional fatality among symptomatic infections
## =======================================================================

conditional_fatality_by_run <- fatal_graph_data |>
  dplyr::filter(
    !is.na(fatal),
    !is.na(symptomatic),
    symptomatic > 0
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
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    conditional_fatality_percent =
      fatal /
      symptomatic *
      100
  )


conditional_fatality_summary <-
  conditional_fatality_by_run |>
  dplyr::group_by(
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    median = median(
      conditional_fatality_percent,
      na.rm = TRUE
    ),
    
    lower = quantile(
      conditional_fatality_percent,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    upper = quantile(
      conditional_fatality_percent,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


p_conditional_fatality <- ggplot2::ggplot(
  conditional_fatality_summary,
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
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = comorb_colours,
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.001
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    title = "Conditional fatality among symptomatic infections",
    subtitle = "Age- and comorbidity-specific fatality risk",
    x = "Age group, years",
    y = "Deaths per 100 symptomatic cases, %",
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

p_conditional_fatality

## =======================================================================
## DALY burden by age and comorbidity
## Input: comorbid_burden_analysis
## =======================================================================

daly_by_run <- comorbid_burden_analysis |>
  dplyr::filter(
    group_number <= 8,
    !is.na(daly),
    !is.na(subgroup_population),
    subgroup_population > 0
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
  ) |>
  dplyr::group_by(
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    yld = sum(
      yld_total,
      na.rm = TRUE
    ),
    
    yll = sum(
      yll,
      na.rm = TRUE
    ),
    
    daly = sum(
      daly,
      na.rm = TRUE
    ),
    
    subgroup_population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    yld_rate_100k =
      yld /
      subgroup_population *
      100000,
    
    yll_rate_100k =
      yll /
      subgroup_population *
      100000,
    
    daly_rate_100k =
      daly /
      subgroup_population *
      100000,
    
    yll_fraction_percent =
      dplyr::if_else(
        daly > 0,
        yll / daly * 100,
        NA_real_
      )
  )


daly_plot_summary <- daly_by_run |>
  dplyr::group_by(
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    daly_mid = median(
      daly_rate_100k,
      na.rm = TRUE
    ),
    
    daly_lo = quantile(
      daly_rate_100k,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    daly_hi = quantile(
      daly_rate_100k,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    yld_mid = median(
      yld_rate_100k,
      na.rm = TRUE
    ),
    
    yll_mid = median(
      yll_rate_100k,
      na.rm = TRUE
    ),
    
    yll_fraction_mid = median(
      yll_fraction_percent,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


p_daly_burden <- ggplot2::ggplot(
  daly_plot_summary,
  ggplot2::aes(
    x = burden_age_group,
    y = daly_mid,
    colour = comorbidity,
    fill = comorbidity,
    group = comorbidity
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = daly_lo,
      ymax = daly_hi
    ),
    alpha = 0.12,
    colour = NA
  ) +
  ggplot2::geom_line(
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    size = 2.2
  ) +
  ggplot2::scale_colour_manual(
    values = comorb_colours,
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = comorb_colours,
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1,
      big.mark = ","
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    title = "Age- and comorbidity-specific chikungunya DALY burden",
    subtitle = "Median and 95% uncertainty interval across model runs",
    x = "Age group, years",
    y = "DALYs per 100,000 subgroup population",
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
      hjust = 1,
      vjust = 1
    )
  )

p_daly_burden


## =======================================================================
## FOI, multimorbidity, and conditional fatality analysis
##
## Inputs:
##   1. allfoi
##   2. comorbid_burden_analysis
## =======================================================================

## =======================================================================
## 1. Country-level FOI aggregation
## =======================================================================

foi_dt <- data.table::as.data.table(
  allfoi
)


# Identify FOI draw columns
foi_cols <- grep(
  pattern = "^foi[0-9]+$",
  x = names(foi_dt),
  value = TRUE
)

if (length(foi_cols) == 0) {
  stop(
    "No FOI draw columns matching '^foi[0-9]+$' were found in allfoi."
  )
}


# Retain rows with valid country and population information
foi_dt_valid <- foi_dt[
  !is.na(iso3) &
    is.finite(tot) &
    tot > 0
]


# Calculate population-weighted FOI for each draw and country
country_foi_draws <- foi_dt_valid[
  ,
  c(
    list(
      country = data.table::first(country),
      continent = data.table::first(continent),
      region = data.table::first(region),
      population = sum(
        tot,
        na.rm = TRUE
      )
    ),
    lapply(
      .SD,
      function(x) {
        stats::weighted.mean(
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


## Convert FOI draws to a country-level median and 95% interval
country_foi <- country_foi_draws |>
  data.table::melt(
    id.vars = c(
      "iso3",
      "country",
      "continent",
      "region",
      "population"
    ),
    measure.vars = foi_cols,
    variable.name = "foi_draw",
    value.name = "foi_value"
  ) |>
  as.data.frame() |>
  dplyr::filter(
    is.finite(foi_value)
  ) |>
  dplyr::group_by(
    iso3,
    country,
    continent,
    region,
    population
  ) |>
  dplyr::summarise(
    foi = median(
      foi_value,
      na.rm = TRUE
    ),
    
    foi_lower = quantile(
      foi_value,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    foi_upper = quantile(
      foi_value,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


## =======================================================================
## 2. Country-specific conditional fatality and multimorbidity prevalence
## =======================================================================

conditional_fatality_country_by_run <-
  comorbid_burden_analysis |>
  dplyr::filter(
    !is.na(fatal),
    !is.na(symptomatic),
    symptomatic > 0,
    !is.na(subgroup_population),
    subgroup_population > 0,
    group_number <= 8
  ) |>
  dplyr::mutate(
    broad_age_group = dplyr::case_when(
      group_number %in% 1:2 ~ "0–19 years",
      group_number %in% 3:6 ~ "20–59 years",
      group_number %in% 7:8 ~ "60–79 years",
      TRUE ~ NA_character_
    ),
    
    multimorbid_population = dplyr::if_else(
      as.character(comorbidity) == "2+",
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
    
    symptomatic = sum(
      symptomatic,
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
    conditional_fatality_percent =
      fatal /
      symptomatic *
      100,
    
    multimorbidity_prevalence =
      multimorbid_population /
      total_population *
      100
  )


## Summarise uncertainty across model runs
conditional_fatality_country <-
  conditional_fatality_country_by_run |>
  dplyr::filter(
    is.finite(
      conditional_fatality_percent
    ),
    is.finite(
      multimorbidity_prevalence
    ),
    is.finite(
      total_population
    ),
    total_population > 0
  ) |>
  dplyr::group_by(
    iso3,
    country,
    continent,
    broad_age_group
  ) |>
  dplyr::summarise(
    conditional_fatality = median(
      conditional_fatality_percent,
      na.rm = TRUE
    ),
    
    conditional_fatality_lower = quantile(
      conditional_fatality_percent,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    conditional_fatality_upper = quantile(
      conditional_fatality_percent,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    multimorbidity_prevalence = median(
      multimorbidity_prevalence,
      na.rm = TRUE
    ),
    
    total_population = median(
      total_population,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


## =======================================================================
## 3. Join fatality estimates with FOI estimates
## =======================================================================

conditional_fatality_country <-
  conditional_fatality_country |>
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
  ) |>
  dplyr::filter(
    is.finite(foi),
    is.finite(
      conditional_fatality
    ),
    is.finite(
      multimorbidity_prevalence
    ),
    is.finite(
      total_population
    ),
    total_population > 0
  ) |>
  dplyr::mutate(
    broad_age_group = factor(
      broad_age_group,
      levels = c(
        "0–19 years",
        "20–59 years",
        "60–79 years"
      )
    )
  )


## =======================================================================
## 4. Prepare bubble-plot data
## =======================================================================

bubble_plot_data <- conditional_fatality_country |>
  dplyr::filter(
    !is.na(broad_age_group),
    is.finite(foi),
    is.finite(
      multimorbidity_prevalence
    ),
    is.finite(
      conditional_fatality
    ),
    is.finite(
      total_population
    ),
    total_population > 0
  )


## Population-weighted reference values within each broad age group
bubble_reference <- bubble_plot_data |>
  dplyr::group_by(
    broad_age_group
  ) |>
  dplyr::summarise(
    reference_foi = stats::weighted.mean(
      foi,
      w = total_population,
      na.rm = TRUE
    ),
    
    reference_multimorbidity =
      stats::weighted.mean(
        multimorbidity_prevalence,
        w = total_population,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


## =======================================================================
## 5. Plot FOI, multimorbidity, and conditional fatality
## =======================================================================

p_foi_multimorbidity_fatality <- ggplot2::ggplot(
  bubble_plot_data,
  ggplot2::aes(
    x = foi,
    y = multimorbidity_prevalence
  )
) +
  
  # Population-weighted reference FOI
  ggplot2::geom_vline(
    data = bubble_reference,
    mapping = ggplot2::aes(
      xintercept = reference_foi
    ),
    inherit.aes = FALSE,
    linetype = "22",
    linewidth = 0.45,
    colour = "grey45"
  ) +
  
  # Population-weighted reference multimorbidity prevalence
  ggplot2::geom_hline(
    data = bubble_reference,
    mapping = ggplot2::aes(
      yintercept =
        reference_multimorbidity
    ),
    inherit.aes = FALSE,
    linetype = "22",
    linewidth = 0.45,
    colour = "grey45"
  ) +
  
  # Country bubbles
  ggplot2::geom_point(
    ggplot2::aes(
      size = total_population,
      fill = conditional_fatality
    ),
    shape = 21,
    colour = "grey20",
    stroke = 0.35,
    alpha = 0.82
  ) +
  
  ggplot2::facet_wrap(
    ~ broad_age_group,
    nrow = 1,
    scales = "free"
  ) +
  
  # Bubble colour: conditional fatality
  ggplot2::scale_fill_viridis_c(
    name = paste0(
      "Conditional fatality\n",
      "among symptomatic cases (%)"
    ),
    option = "magma",
    direction = -1,
    labels = scales::label_number(
      accuracy = 0.001,
      suffix = "%"
    ),
    guide = ggplot2::guide_colourbar(
      title.position = "top",
      title.hjust = 0,
      barheight = grid::unit(
        45,
        "pt"
      ),
      barwidth = grid::unit(
        9,
        "pt"
      ),
      order = 2
    )
  ) +
  
  # Bubble size: broad age-group population
  ggplot2::scale_size_continuous(
    name = "Age-group population",
    range = c(
      2,
      12
    ),
    labels = scales::label_number(
      scale = 1e-6,
      suffix = "M",
      accuracy = 1
    ),
    guide = ggplot2::guide_legend(
      title.position = "top",
      title.hjust = 0,
      order = 1
    )
  ) +
  
  ggplot2::scale_x_continuous(
    labels = scales::label_number(
      accuracy = 0.001
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.04,
        0.08
      )
    )
  ) +
  
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 1,
      suffix = "%"
    ),
    expand = ggplot2::expansion(
      mult = c(
        0.04,
        0.08
      )
    )
  ) +
  
  ggplot2::labs(
    title = paste0(
      "Transmission, multimorbidity, ",
      "and conditional fatality"
    ),
    subtitle = paste(
      "Each bubble represents one country;",
      "dashed lines show age-group",
      "population-weighted averages"
    ),
    x = "Country force of infection",
    y = "Population with ≥2 comorbidities, %"
  ) +
  
  theme_lancet_clean(
    base_size = 10
  ) +
  
  ggplot2::theme(
    legend.position = "right",
    legend.box = "vertical",
    
    legend.title = ggplot2::element_text(
      size = 8.5,
      face = "bold",
      hjust = 0
    ),
    
    legend.text = ggplot2::element_text(
      size = 8
    ),
    
    strip.background =
      ggplot2::element_blank(),
    
    strip.text =
      ggplot2::element_text(
        face = "bold",
        hjust = 0
      ),
    
    axis.text.x =
      ggplot2::element_text(
        angle = 0,
        hjust = 0.5
      ),
    
    panel.spacing =
      grid::unit(
        12,
        "pt"
      ),
    
    plot.margin =
      ggplot2::margin(
        t = 8,
        r = 30,
        b = 8,
        l = 8
      )
  )

p_foi_multimorbidity_fatality


################################################################################
## Comorbidity-stratified burden figures ----------------------------------
## Builds on objects already in memory after running the Step 3
## (comorbid_burden_analysis) pipeline in ASSUMPTIONS.md.
## Requires: infection_comorb_long, adjusted_hosp_rate_103,
##           adjusted_fatal_rate_103, theme_lancet_clean()
## Set to NULL for a global (103-country) aggregate, or e.g. "BRA" for a
## single-country figure.
country_filter <- NULL

filter_country <- function(df) {
  if (is.null(country_filter)) return(df)
  dplyr::filter(df, iso3 == country_filter)
}

## Consistent age-band ordering / labels for the x-axis on both figures ---
age_labels <- infection_comorb_long |>
  dplyr::distinct(group_number, burden_age_group) |>
  dplyr::arrange(group_number)

comorb_colors <- c(
  "0"  = "#B0B0B0",
  "1"  = "#F1A340",
  "2+" = "#B2182B"
)

## ==========================================================================
## FIGURE 2. Comorbidity-stratified hospitalization burden by age
## ==========================================================================
## expected_hosp (per run) = infections_comorb * adjusted_hosp_rate
## infections_comorb is deterministic; adjusted_hosp_rate carries the 1000
## LHS draws, so this propagates RR uncertainty through to expected counts.

expected_hosp_by_stratum <- infection_comorb_long |>
  filter_country() |>
  dplyr::inner_join(
    adjusted_hosp_rate_103 |>
      dplyr::select(iso3, group_number, comorbidity, run, adjusted_hosp_rate),
    by = c("iso3", "group_number", "comorbidity"),
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(
    expected_hosp = infections_comorb * adjusted_hosp_rate
  )

hosp_stack_data <- expected_hosp_by_stratum |>
  dplyr::group_by(group_number, comorbidity, run) |>
  dplyr::summarise(
    total_hosp = sum(expected_hosp, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(group_number, comorbidity) |>
  dplyr::summarise(
    median_hosp = median(total_hosp, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    comorbidity = factor(comorbidity, levels = c("0", "1", "2+"))
  )

fig2_hosp_stack <- ggplot(
  hosp_stack_data,
  aes(x = group_number, y = median_hosp, fill = comorbidity)
) +
  geom_area(position = "stack", colour = "white", linewidth = 0.2) +
  scale_x_continuous(
    breaks = age_labels$group_number,
    labels = age_labels$burden_age_group
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = comorb_colors) +
  labs(
    x = "Age group",
    y = "Expected hospitalizations (median across 1,000 draws)",
    fill = "Comorbidity count",
    title = "Comorbidity-stratified chikungunya hospitalization burden by age"
  ) +
  theme_lancet_clean()

fig2_hosp_stack

ggsave(
  "Figures/fig2_hosp_comorbidity_stack.pdf",
  fig2_hosp_stack,
  width = 7, height = 5, units = "in"
)

## ==========================================================================
## FIGURE 3. Comorbidity-stratified death burden by age
## (hospitalized-pathway vs non-hospitalized-pathway, two panels)
## ==========================================================================
## ASSUMPTION FLAGGED: this treats adjusted_fatal_rate_hosp and
## adjusted_fatal_rate_nonhosp as rates PER INFECTION (i.e. the same
## convention as adjusted_hosp_rate), so both are multiplied by
## infections_comorb directly. If fatal_sample / rr_death_model were
## actually built as a case-fatality ratio CONDITIONAL ON hospitalization
## (i.e. deaths per hospitalized case, not per infection), the hosp-pathway
## term below should instead multiply expected_hosp (from Figure 2), not
## infections_comorb. Please confirm which convention fatal_sample /
## nh_fatal_sample use before trusting this panel.

expected_fatal_by_stratum <- infection_comorb_long |>
  filter_country() |>
  dplyr::inner_join(
    adjusted_fatal_rate_103 |>
      dplyr::select(
        iso3, group_number, comorbidity, run,
        adjusted_fatal_rate_hosp, adjusted_fatal_rate_nonhosp
      ),
    by = c("iso3", "group_number", "comorbidity"),
    relationship = "many-to-many"
  ) |>
  dplyr::mutate(
    expected_death_hosp    = infections_comorb * adjusted_fatal_rate_hosp,
    expected_death_nonhosp = infections_comorb * adjusted_fatal_rate_nonhosp
  )

fatal_stack_data <- expected_fatal_by_stratum |>
  tidyr::pivot_longer(
    cols = c(expected_death_hosp, expected_death_nonhosp),
    names_to = "pathway",
    values_to = "expected_deaths"
  ) |>
  dplyr::mutate(
    pathway = dplyr::recode(
      pathway,
      expected_death_hosp    = "Hospitalized pathway",
      expected_death_nonhosp = "Non-hospitalized pathway"
    )
  ) |>
  dplyr::group_by(group_number, comorbidity, pathway, run) |>
  dplyr::summarise(
    total_deaths = sum(expected_deaths, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::group_by(group_number, comorbidity, pathway) |>
  dplyr::summarise(
    median_deaths = median(total_deaths, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    comorbidity = factor(comorbidity, levels = c("0", "1", "2+"))
  )

fig3_death_stack <- ggplot(
  fatal_stack_data,
  aes(x = group_number, y = median_deaths, fill = comorbidity)
) +
  geom_area(position = "stack", colour = "white", linewidth = 0.2) +
  facet_wrap(~pathway, scales = "free_y", ncol = 2) +
  scale_x_continuous(
    breaks = age_labels$group_number,
    labels = age_labels$burden_age_group
  ) +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = comorb_colors) +
  labs(
    x = "Age group",
    y = "Expected deaths (median across 1,000 draws)",
    fill = "Comorbidity count",
    title = "Comorbidity-stratified chikungunya death burden by age and pathway"
  ) +
  theme_lancet_clean() +
  theme(strip.background = element_blank())

fig3_death_stack

ggsave(
  "Figures/fig3_death_comorbidity_stack.pdf",
  fig3_death_stack,
  width = 9, height = 5, units = "in"
)

################################################################################
## =======================================================================
## Absolute global burden across comorbidity groups
## =======================================================================

absolute_burden_by_run <- comorbid_burden_analysis |>
  dplyr::filter(
    group_number <= 8,
    !is.na(subgroup_population),
    subgroup_population >= 0
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
  ) |>
  dplyr::group_by(
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    Population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    
    Infections = sum(
      infections_comorb,
      na.rm = TRUE
    ),
    
    Symptomatic = sum(
      symptomatic,
      na.rm = TRUE
    ),
    
    Hospitalised = sum(
      hospitalised,
      na.rm = TRUE
    ),
    
    Fatal = sum(
      fatal,
      na.rm = TRUE
    ),
    
    DALY = sum(
      daly,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(
    cols = c(
      Population,
      Infections,
      Symptomatic,
      Hospitalised,
      Fatal,
      DALY
    ),
    names_to = "outcome",
    values_to = "burden"
  )


absolute_burden_summary <- absolute_burden_by_run |>
  dplyr::group_by(
    burden_age_group,
    outcome,
    comorbidity
  ) |>
  dplyr::summarise(
    burden_median = median(
      burden,
      na.rm = TRUE
    ),
    
    burden_lower = quantile(
      burden,
      probs = 0.025,
      na.rm = TRUE
    ),
    
    burden_upper = quantile(
      burden,
      probs = 0.975,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    
    outcome = factor(
      outcome,
      levels = outcome_levels
    ),
    
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
  )


p_absolute_burden <- ggplot2::ggplot(
  absolute_burden_summary,
  ggplot2::aes(
    x = burden_age_group,
    y = burden_median,
    fill = comorbidity
  )
) +
  ggplot2::geom_col(
    width = 0.76,
    colour = "white",
    linewidth = 0.25,
    position = ggplot2::position_stack(
      reverse = TRUE
    )
  ) +
  ggplot2::facet_wrap(
    ~ outcome,
    ncol = 1,
    scales = "free_y",
    labeller = ggplot2::labeller(
      outcome = outcome_labels
    )
  ) +
  ggplot2::scale_fill_manual(
    values = comorb_colours,
    breaks = comorb_levels,
    labels = comorb_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      scale_cut = scales::cut_short_scale()
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.05)
    )
  ) +
  ggplot2::labs(
    title = paste0(
      "Global chikungunya burden ",
      "across comorbidity groups"
    ),
    subtitle = paste0(
      "Median absolute burden across model runs"
    ),
    x = "Age group, years",
    y = "Absolute burden",
    fill = "Comorbidity count"
  ) +
  theme_lancet_clean(
    base_size = 10
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    
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
      grid::unit(
        7,
        "pt"
      )
  )

p_absolute_burden

## =======================================================================
## Age-specific burden attribution plot with 95% uncertainty intervals
##
## Coloured segments:
##   burden occurring within each comorbidity group
##
## Error bars:
##   95% uncertainty interval for the total stacked burden
## =======================================================================

make_burden_attribution_plot <- function(
    data,
    outcome_variable,
    multiplier,
    y_axis_title,
    plot_title,
    plot_subtitle = NULL
) {
  
  outcome_variable <- rlang::ensym(
    outcome_variable
  )
  
  ## ---------------------------------------------------------------------
  ## 1. Calculate comorbidity-specific contributions within each run
  ## ---------------------------------------------------------------------
  
  attribution_by_run <- data |>
    dplyr::filter(
      group_number <= 8,
      !is.na(subgroup_population),
      subgroup_population >= 0,
      !is.na(!!outcome_variable),
      !!outcome_variable >= 0
    ) |>
    dplyr::mutate(
      burden_age_group = factor(
        as.character(burden_age_group),
        levels = age_levels
      ),
      
      comorbidity = factor(
        as.character(comorbidity),
        levels = comorb_levels
      )
    ) |>
    dplyr::group_by(
      run,
      burden_age_group,
      comorbidity
    ) |>
    dplyr::summarise(
      outcome_burden = sum(
        !!outcome_variable,
        na.rm = TRUE
      ),
      
      subgroup_population = sum(
        subgroup_population,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    dplyr::group_by(
      run,
      burden_age_group
    ) |>
    dplyr::mutate(
      total_age_population = sum(
        subgroup_population,
        na.rm = TRUE
      ),
      
      burden_contribution =
        outcome_burden /
        total_age_population *
        multiplier
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(
      is.finite(burden_contribution),
      total_age_population > 0
    )
  
  
  ## ---------------------------------------------------------------------
  ## 2. Median contribution from each comorbidity group
  ## ---------------------------------------------------------------------
  
  attribution_summary <- attribution_by_run |>
    dplyr::group_by(
      burden_age_group,
      comorbidity
    ) |>
    dplyr::summarise(
      contribution_median = median(
        burden_contribution,
        na.rm = TRUE
      ),
      
      contribution_lower = quantile(
        burden_contribution,
        probs = 0.025,
        na.rm = TRUE
      ),
      
      contribution_upper = quantile(
        burden_contribution,
        probs = 0.975,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      burden_age_group = factor(
        as.character(burden_age_group),
        levels = age_levels
      ),
      
      comorbidity = factor(
        as.character(comorbidity),
        levels = comorb_levels
      )
    )
  
  
  ## ---------------------------------------------------------------------
  ## 3. Total stacked burden within each run
  ##
  ## The uncertainty interval must be calculated from the run-specific
  ## total, rather than by summing subgroup-specific quantiles.
  ## ---------------------------------------------------------------------
  
  total_attribution_by_run <- attribution_by_run |>
    dplyr::group_by(
      run,
      burden_age_group
    ) |>
    dplyr::summarise(
      total_burden_rate = sum(
        burden_contribution,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    )
  
  
  total_attribution_summary <- total_attribution_by_run |>
    dplyr::group_by(
      burden_age_group
    ) |>
    dplyr::summarise(
      total_median = median(
        total_burden_rate,
        na.rm = TRUE
      ),
      
      total_lower = quantile(
        total_burden_rate,
        probs = 0.025,
        na.rm = TRUE
      ),
      
      total_upper = quantile(
        total_burden_rate,
        probs = 0.975,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      burden_age_group = factor(
        as.character(burden_age_group),
        levels = age_levels
      )
    )
  
  
  ## ---------------------------------------------------------------------
  ## 4. Plot stacked median contributions with total 95% UI
  ## ---------------------------------------------------------------------
  
  plot_object <- ggplot2::ggplot(
    attribution_summary,
    ggplot2::aes(
      x = burden_age_group,
      y = contribution_median,
      fill = comorbidity
    )
  ) +
    
    ggplot2::geom_col(
      width = 0.76,
      colour = "white",
      linewidth = 0.25,
      position = ggplot2::position_stack(
        reverse = TRUE
      )
    ) +
    
    ## Total burden 95% uncertainty interval
    ggplot2::geom_errorbar(
      data = total_attribution_summary,
      mapping = ggplot2::aes(
        x = burden_age_group,
        ymin = total_lower,
        ymax = total_upper
      ),
      inherit.aes = FALSE,
      width = 0.18,
      linewidth = 0.55,
      colour = "grey20"
    ) +
    
    ## Median total burden marker
    ggplot2::geom_point(
      data = total_attribution_summary,
      mapping = ggplot2::aes(
        x = burden_age_group,
        y = total_median
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 1.7,
      stroke = 0.45,
      colour = "grey20",
      fill = "white"
    ) +
    
    ggplot2::scale_fill_manual(
      values = comorb_colours,
      breaks = comorb_levels,
      labels = comorb_labels,
      drop = FALSE
    ) +
    
    ggplot2::scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.1,
        scale_cut = scales::cut_short_scale()
      ),
      expand = ggplot2::expansion(
        mult = c(0, 0.10)
      )
    ) +
    
    ggplot2::labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Age group, years",
      y = y_axis_title,
      fill = "Comorbidity count"
    ) +
    
    theme_lancet_clean(
      base_size = 10
    ) +
    
    ggplot2::theme(
      legend.position = "bottom",
      
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),
      
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )
  
  
  list(
    by_run = attribution_by_run,
    subgroup_summary = attribution_summary,
    total_by_run = total_attribution_by_run,
    total_summary = total_attribution_summary,
    plot = plot_object
  )
}

fatal_attribution_results <-
  make_burden_attribution_plot(
    data = comorbid_burden_analysis,
    outcome_variable = fatal,
    multiplier = 100000,
    y_axis_title =
      "Deaths per 100,000 age-group population",
    plot_title =
      paste0(
        "Age-specific chikungunya mortality burden ",
        "by comorbidity group"
      ),
    plot_subtitle =
      paste0(
        "Bars show median contributions by comorbidity group; ",
        "error bars show the 95% uncertainty interval ",
        "for total mortality burden"
      )
  )

p_fatal_attribution <-
  fatal_attribution_results$plot

p_fatal_attribution


daly_attribution_results <-
  make_burden_attribution_plot(
    data = comorbid_burden_analysis,
    outcome_variable = daly,
    multiplier = 100000,
    y_axis_title =
      "DALYs per 100,000 age-group population",
    plot_title =
      paste0(
        "Age-specific chikungunya DALY burden ",
        "by comorbidity group"
      ),
    plot_subtitle =
      paste0(
        "Bars show median contributions by comorbidity group; ",
        "error bars show the 95% uncertainty interval ",
        "for total DALY burden"
      )
  )

p_daly_attribution <-
  daly_attribution_results$plot

p_daly_attribution

hospitalisation_attribution_results <-
  make_burden_attribution_plot(
    data = comorbid_burden_analysis,
    outcome_variable = hospitalised,
    multiplier = 10000,
    y_axis_title =
      "Hospitalisations per 10,000 age-group population",
    plot_title =
      "Hospitalisation burden",
    plot_subtitle =
      paste0(
        "Stacked bars show median contributions by comorbidity group; ",
        "error bars show 95% uncertainty intervals for total burden"
      )
  )

p_hospitalisation_attribution <-
  hospitalisation_attribution_results$plot

p_hospitalisation_attribution

library(patchwork)

p_attribution_combined <-
  (
    p_hospitalisation_attribution /
      p_fatal_attribution /
      p_daly_attribution
  ) +
  patchwork::plot_layout(
    guides = "collect",
    heights = c(1, 1, 1)
  ) +
  patchwork::plot_annotation(
    title =
      paste0(
        "Age-specific chikungunya burden ",
        "by comorbidity group"
      ),
    subtitle =
      paste0(
        "Coloured segments represent burden occurring within each ",
        "comorbidity group; error bars show 95% uncertainty intervals ",
        "for the total age-specific burden"
      ),
    tag_levels = "A"
  ) &
  ggplot2::theme(
    legend.position = "bottom"
  )

p_attribution_combined



## =======================================================================
## Subgroup-specific age pattern
##
## Denominator:
##   age- and comorbidity-specific subgroup population
##
## Error bars:
##   95% uncertainty intervals across PSA runs
## =======================================================================

make_subgroup_age_plot <- function(
    data,
    outcome_variable,
    multiplier,
    y_axis_title,
    plot_title
) {
  
  outcome_variable <- rlang::ensym(
    outcome_variable
  )
  
  subgroup_by_run <- data |>
    dplyr::filter(
      group_number <= 8,
      !is.na(subgroup_population),
      subgroup_population > 0,
      !is.na(!!outcome_variable),
      !!outcome_variable >= 0
    ) |>
    dplyr::mutate(
      burden_age_group = factor(
        as.character(burden_age_group),
        levels = age_levels
      ),
      
      comorbidity = factor(
        as.character(comorbidity),
        levels = comorb_levels
      )
    ) |>
    dplyr::group_by(
      run,
      burden_age_group,
      comorbidity
    ) |>
    dplyr::summarise(
      outcome_burden = sum(
        !!outcome_variable,
        na.rm = TRUE
      ),
      
      subgroup_population = sum(
        subgroup_population,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      subgroup_rate =
        outcome_burden /
        subgroup_population *
        multiplier
    ) |>
    dplyr::filter(
      is.finite(subgroup_rate)
    )
  
  
  subgroup_summary <- subgroup_by_run |>
    dplyr::group_by(
      burden_age_group,
      comorbidity
    ) |>
    dplyr::summarise(
      rate_median = median(
        subgroup_rate,
        na.rm = TRUE
      ),
      
      rate_lower = quantile(
        subgroup_rate,
        probs = 0.025,
        na.rm = TRUE
      ),
      
      rate_upper = quantile(
        subgroup_rate,
        probs = 0.975,
        na.rm = TRUE
      ),
      
      .groups = "drop"
    ) |>
    dplyr::mutate(
      burden_age_group = factor(
        as.character(burden_age_group),
        levels = age_levels
      ),
      
      comorbidity = factor(
        as.character(comorbidity),
        levels = comorb_levels
      )
    )
  
  
  plot_object <- ggplot2::ggplot(
    subgroup_summary,
    ggplot2::aes(
      x = burden_age_group,
      y = rate_median,
      fill = comorbidity
    )
  ) +
    
    ggplot2::geom_col(
      width = 0.72,
      colour = "white",
      linewidth = 0.3
    ) +
    
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = rate_lower,
        ymax = rate_upper
      ),
      width = 0.16,
      linewidth = 0.55,
      colour = "grey20"
    ) +
    
    ggplot2::facet_wrap(
      ~ comorbidity,
      nrow = 1,
      labeller = ggplot2::labeller(
        comorbidity = comorb_labels
      )
    ) +
    
    ggplot2::scale_fill_manual(
      values = comorb_colours,
      breaks = comorb_levels,
      labels = comorb_labels,
      drop = FALSE,
      guide = "none"
    ) +
    
    ggplot2::scale_y_continuous(
      labels = scales::label_number(
        accuracy = 0.1,
        scale_cut = scales::cut_short_scale()
      ),
      expand = ggplot2::expansion(
        mult = c(0, 0.08)
      )
    ) +
    
    ggplot2::labs(
      title = plot_title,
      subtitle = paste0(
        "Bars show median subgroup-specific rates; ",
        "error bars show 95% uncertainty intervals"
      ),
      x = "Age group, years",
      y = y_axis_title
    ) +
    
    theme_lancet_clean(
      base_size = 10
    ) +
    
    ggplot2::theme(
      strip.background =
        ggplot2::element_blank(),
      
      strip.text =
        ggplot2::element_text(
          face = "bold",
          size = 9
        ),
      
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
      
      panel.spacing =
        grid::unit(
          10,
          "pt"
        ),
      
      plot.title =
        ggplot2::element_text(
          face = "bold"
        )
    )
  
  
  list(
    by_run = subgroup_by_run,
    summary = subgroup_summary,
    plot = plot_object
  )
}


fatal_subgroup_results <-
  make_subgroup_age_plot(
    data = comorbid_burden_analysis,
    outcome_variable = fatal,
    multiplier = 100000,
    y_axis_title =
      "Deaths per 100,000 subgroup population",
    plot_title =
      paste0(
        "Age-specific chikungunya mortality rates ",
        "within comorbidity groups"
      )
  )

p_fatal_subgroup <-
  fatal_subgroup_results$plot

p_fatal_subgroup


hospitalisation_subgroup_results <-
  make_subgroup_age_plot(
    data = comorbid_burden_analysis,
    outcome_variable = hospitalised,
    multiplier = 10000,
    y_axis_title =
      "Hospitalisations per 10,000 subgroup population",
    plot_title =
      paste0(
        "Age-specific chikungunya hospitalisation rates ",
        "within comorbidity groups"
      )
  )

p_hospitalisation_subgroup <-
  hospitalisation_subgroup_results$plot

p_hospitalisation_subgroup


daly_subgroup_results <-
  make_subgroup_age_plot(
    data = comorbid_burden_analysis,
    outcome_variable = daly,
    multiplier = 100000,
    y_axis_title =
      "DALYs per 100,000 subgroup population",
    plot_title =
      paste0(
        "Age-specific chikungunya DALY rates ",
        "within comorbidity groups"
      )
  )

p_daly_subgroup <-
  daly_subgroup_results$plot

p_daly_subgroup

## ============================================================
## Step 1. Age-specific infection rates by comorbidity group
## ============================================================

infection_rate_by_run <-
  comorbid_burden_analysis |>
  dplyr::filter(
    group_number <= 8,
    !is.na(subgroup_population),
    subgroup_population > 0,
    !is.na(infections_comorb),
    infections_comorb >= 0
  ) |>
  dplyr::mutate(
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    ),
    comorbidity = factor(
      as.character(comorbidity),
      levels = comorb_levels
    )
  ) |>
  dplyr::group_by(
    run,
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    infections = sum(
      infections_comorb,
      na.rm = TRUE
    ),
    subgroup_population = sum(
      subgroup_population,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    infection_rate =
      infections /
      subgroup_population *
      100000
  ) |>
  dplyr::filter(
    is.finite(infection_rate)
  )


infection_rate_summary <-
  infection_rate_by_run |>
  dplyr::group_by(
    burden_age_group,
    comorbidity
  ) |>
  dplyr::summarise(
    rate_median = median(
      infection_rate,
      na.rm = TRUE
    ),
    rate_lower = quantile(
      infection_rate,
      probs = 0.025,
      na.rm = TRUE
    ),
    rate_upper = quantile(
      infection_rate,
      probs = 0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

p_infection_rate_by_subgroup <-
  ggplot2::ggplot(
    infection_rate_summary,
    ggplot2::aes(
      x = burden_age_group,
      y = rate_median,
      fill = comorbidity
    )
  ) +
  ggplot2::geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.3
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = rate_lower,
      ymax = rate_upper
    ),
    width = 0.16,
    linewidth = 0.55,
    colour = "grey20"
  ) +
  ggplot2::facet_wrap(
    ~ comorbidity,
    nrow = 1,
    labeller = ggplot2::labeller(
      comorbidity = comorb_labels
    )
  ) +
  ggplot2::scale_fill_manual(
    values = comorb_colours,
    guide = "none",
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 1,
      scale_cut = scales::cut_short_scale()
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    title =
      "Age-specific chikungunya infection rates within comorbidity groups",
    subtitle =
      "Bars show median subgroup-specific rates; error bars show 95% uncertainty intervals",
    x = "Age group, years",
    y = "Infections per 100,000 subgroup population"
  ) +
  theme_lancet_clean(
    base_size = 10
  ) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )

p_infection_rate_by_subgroup