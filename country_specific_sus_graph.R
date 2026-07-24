## ======================================================================
## 1. Prepare batch 1 plotting data
## ======================================================================
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

archetype_colours <- c(
  "A" = "#E64B35",
  "B" = "#00A087",
  "C" = "#3C5488",
  "D" = "#7A5195"
)

archetype_labels <- c(
  A = "A. Finite circulation since introduction",
  B = "B. Long-term endemic equilibrium",
  C = "C. Episodic transmission with inter-epidemic periods",
  D = "D. Conditional epidemic potential"
)

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(scales)
library(sf)
library(rnaturalearth)
library(ggplot2)

# Identify countries that were removed during post-processing
excluded_countries <- burden_batch_1 |>
  dplyr::filter(
    is.na(model_code) |
      !model_code %in% c("A", "B", "C", "D")
  ) |>
  dplyr::distinct(country) |>
  dplyr::pull(country)

excluded_countries
length(excluded_countries)

burden_plot_data <- burden_ABCD_all |>
  dplyr::mutate(
    country = as.character(country),
    model_code = as.character(model_code),
    age_band = as.character(age_band),
    draw_id = as.integer(draw_id)
  ) |>
  dplyr::filter(
    !is.na(country),
    country != "",
    model_code %in% c("A", "B", "C", "D"),
    !is.na(draw_id),
    !is.na(age_band)
  ) |>
  dplyr::mutate(
    model_code = factor(
      model_code,
      levels = c("A", "B", "C", "D")
    ),
    age_band = factor(
      age_band,
      levels = age_band_labels
    )
  )

country_archetype <- df_for_burden |>
  dplyr::mutate(
    country = as.character(country),
    iso3 = as.character(iso3),
    model_code = as.character(model_code)
  ) |>
  dplyr::filter(
    !is.na(country),
    country != "",
    !country %in% excluded_countries,
    !is.na(iso3),
    iso3 != "",
    model_code %in% c("A", "B", "C", "D")
  ) |>
  dplyr::distinct(
    country,
    iso3,
    model_code
  )

world_sf <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
) |>
  dplyr::mutate(
    map_iso3 = dplyr::if_else(
      iso_a3 == "-99",
      adm0_a3,
      iso_a3
    )
  ) |>
  dplyr::filter(
    continent != "Antarctica"
  )

world_archetype <- world_sf |>
  dplyr::left_join(
    country_archetype,
    by = c("map_iso3" = "iso3")
  ) |>
  dplyr::mutate(
    model_code = factor(
      model_code,
      levels = c("A", "B", "C", "D")
    )
  )

p_archetype_map <- ggplot2::ggplot(
  world_archetype
) +
  ggplot2::geom_sf(
    ggplot2::aes(fill = model_code),
    colour = "black",
    linewidth = 0.12
  ) +
  ggplot2::scale_fill_manual(
    values = archetype_colours,
    labels = archetype_labels,
    na.value = "grey90",
    drop = FALSE
  ) +
  ggplot2::coord_sf(
    xlim = c(-180, 180),
    ylim = c(-58, 84),   # crop Antarctica + high Arctic
    expand = FALSE,
    datum = NA
  ) +
  ggplot2::labs(
    title = "Global chikungunya transmission archetypes",
    caption = paste0(
      "Countries included in the analysis are classified into four transmission archetypes. ",
      "Countries not included are shown in grey."
    )
  ) +
  theme_lancet_clean() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 9),
    legend.key.width = grid::unit(1.8, "cm"),
    plot.caption = ggplot2::element_text(
      size = 9,
      colour = "grey30",
      hjust = 0
    )
  )

p_archetype_map
################################################################################
## ======================================================================
## 1. Libraries and plotting settings
## ======================================================================


age_band_labels <- c(
  "[0,1)",
  "[1,5)",
  "[5,10)",
  "[10,15)",
  "[15,20)",
  "[20,25)",
  "[25,30)",
  "[30,35)",
  "[35,40)",
  "[40,45)",
  "[45,50)",
  "[50,55)",
  "[55,60)",
  "[60,65)",
  "[65,70)",
  "[70,75)",
  "[75,80)",
  "[80,100]"
)

archetype_colours <- c(
  A = "#E64B35",
  B = "#00A087",
  C = "#3C5488",
  D = "#7A5195"
)

archetype_labels <- c(
  A = "A. Finite circulation since introduction",
  B = "B. Long-term endemic equilibrium",
  C = "C. Episodic transmission with inter-epidemic periods",
  D = "D. Conditional epidemic potential"
)

## ======================================================================
## 2. Read all 100 infection draws
## ======================================================================

burden_ABCD_all <- readRDS(
  "Outputs/burden_ABCD_all.rds"
)

burden_plot_data <- burden_ABCD_all |>
  dplyr::mutate(
    country = as.character(country),
    model_code = as.character(model_code),
    age_band = as.character(age_band),
    draw_id = as.integer(draw_id),
    annual_infections = as.numeric(annual_infections)
  ) |>
  dplyr::filter(
    !is.na(country),
    country != "",
    model_code %in% c("A", "B", "C", "D"),
    !is.na(draw_id),
    !is.na(age_band)
  )


## ======================================================================
## 3. Construct country-by-age population denominators
## ======================================================================

age_pop_lookup <- tibble::tibble(
  age_col = as.character(1:18),
  age_band = age_band_labels
)

population_long <- df_for_burden |>
  dplyr::transmute(
    country = as.character(country),
    model_code = as.character(model_code),
    dplyr::across(
      dplyr::all_of(as.character(1:18)),
      as.numeric
    )
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(as.character(1:18)),
    names_to = "age_col",
    values_to = "population"
  ) |>
  dplyr::left_join(
    age_pop_lookup,
    by = "age_col",
    relationship = "many-to-one"
  ) |>
  dplyr::filter(
    model_code %in% c("A", "B", "C", "D"),
    !is.na(country),
    country != "",
    !is.na(age_band)
  ) |>
  dplyr::group_by(
    country,
    model_code,
    age_band
  ) |>
  dplyr::summarise(
    population = sum(
      population,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

## ======================================================================
## 4. Calculate country-specific incidence for each infection draw
## ======================================================================

incidence_country_draw <- burden_plot_data |>
  dplyr::left_join(
    population_long,
    by = c(
      "country",
      "model_code",
      "age_band"
    ),
    relationship = "many-to-one"
  ) |>
  dplyr::filter(
    !is.na(population),
    population > 0,
    !is.na(annual_infections)
  ) |>
  dplyr::mutate(
    annual_incidence =
      annual_infections /
      population
  )

## ======================================================================
## 5. Summarise the typical country within each archetype and draw
## ======================================================================

incidence_archetype_draw <- incidence_country_draw |>
  dplyr::group_by(
    draw_id,
    model_code,
    age_band
  ) |>
  dplyr::summarise(
    annual_incidence =
      stats::median(
        annual_incidence,
        na.rm = TRUE
      ),
    
    n_countries =
      dplyr::n_distinct(country),
    
    .groups = "drop"
  )
## ======================================================================
## 6. Summarise uncertainty across 100 infection draws
## ======================================================================

incidence_arch_summary <- incidence_archetype_draw |>
  dplyr::group_by(
    model_code,
    age_band
  ) |>
  dplyr::summarise(
    estimate =
      stats::median(
        annual_incidence,
        na.rm = TRUE
      ),
    
    lower =
      stats::quantile(
        annual_incidence,
        probs = 0.025,
        na.rm = TRUE
      ),
    
    upper =
      stats::quantile(
        annual_incidence,
        probs = 0.975,
        na.rm = TRUE
      ),
    
    n_draws =
      dplyr::n_distinct(draw_id),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    model_code = factor(
      model_code,
      levels = c("A", "B", "C", "D")
    ),
    
    age_band = factor(
      age_band,
      levels = age_band_labels
    )
  )
## ======================================================================
## 7. Plot country-standardised archetype incidence
## ======================================================================

p_incidence_by_archetype <- ggplot2::ggplot(
  incidence_arch_summary,
  ggplot2::aes(
    x = age_band,
    y = estimate * 100,
    colour = model_code,
    fill = model_code,
    group = model_code
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = lower * 100,
      ymax = upper * 100
    ),
    alpha = 0.14,
    colour = NA
  ) +
  ggplot2::geom_line(
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    size = 2
  ) +
  ggplot2::scale_colour_manual(
    values = archetype_colours,
    labels = archetype_labels,
    drop = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = archetype_colours,
    labels = archetype_labels,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1,
      suffix = "%"
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.06)
    )
  ) +
  ggplot2::labs(
    title = "Age-specific infection incidence by transmission archetype",
    subtitle = "Median country-specific incidence within each archetype",
    x = "Age band",
    y = "Incidence (%)",
    caption = paste0(
      "A–C represent expected annual incidence; ",
      "D represents incidence conditional on epidemic transmission."
    )
  ) +
  theme_lancet_clean() +
  ggplot2::theme(
    legend.position = "right",
    plot.subtitle = ggplot2::element_text(
      size = 10,
      colour = "grey30"
    ),
    plot.caption = ggplot2::element_text(
      size = 9,
      colour = "grey30",
      hjust = 0
    )
  )

p_incidence_by_archetype

## ======================================================================
## 8. Country-specific diagnostic for Archetype D
## ======================================================================

d_country_summary <- incidence_country_draw |>
  dplyr::filter(
    model_code == "D"
  ) |>
  dplyr::group_by(
    country,
    age_band
  ) |>
  dplyr::summarise(
    median_incidence =
      stats::median(
        annual_incidence,
        na.rm = TRUE
      ),
    
    lower =
      stats::quantile(
        annual_incidence,
        0.025,
        na.rm = TRUE
      ),
    
    upper =
      stats::quantile(
        annual_incidence,
        0.975,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) |>
  dplyr::mutate(
    age_band = factor(
      age_band,
      levels = age_band_labels
    )
  )

p_d_country <- ggplot2::ggplot(
  d_country_summary,
  ggplot2::aes(
    x = age_band,
    y = median_incidence * 100,
    group = country
  )
) +
  ggplot2::geom_line(
    alpha = 0.6,
    linewidth = 0.6
  ) +
  ggplot2::geom_point(
    alpha = 0.7,
    size = 1.2
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1,
      suffix = "%"
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.06)
    )
  ) +
  ggplot2::labs(
    title = "Country-specific conditional incidence in Archetype D",
    subtitle = "Median across 100 infection draws",
    x = "Age band",
    y = "Conditional incidence (%)"
  ) +
  theme_lancet_clean()

p_d_country

