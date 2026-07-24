age_map <- tibble::tibble(
  group_number = 1:9,
  burden_age_group = factor(
    c(
      "[0,10)", "[10,20)", "[20,30)",
      "[30,40)", "[40,50)", "[50,60)",
      "[60,70)", "[70,80)", "80+"
    ),
    levels = c(
      "[0,10)", "[10,20)", "[20,30)",
      "[30,40)", "[40,50)", "[50,60)",
      "[60,70)", "[70,80)", "80+"
    )
  )
)

external_hosp_age_summary <-
  tibble::tibble(
    group_number = hosp_band$band_num,
    estimate = hosp_band$hosp_rate,
    alpha = vapply(
      hosp_beta_params,
      function(x) x$alpha,
      numeric(1)
    ),
    beta = vapply(
      hosp_beta_params,
      function(x) x$beta,
      numeric(1)
    )
  ) |>
  dplyr::filter(group_number <= 8) |>
  dplyr::mutate(
    lower = stats::qbeta(
      0.025,
      shape1 = alpha,
      shape2 = beta
    ),
    upper = stats::qbeta(
      0.975,
      shape1 = alpha,
      shape2 = beta
    ),
    source = "Marginal hosp rate input"
  ) |>
  dplyr::left_join(
    age_map,
    by = "group_number",
    relationship = "many-to-one"
  )


brazil_observed_hosp_by_age <-
  readRDS(
    "MainData/brazil_observed_hosp_by_age.rds"
  )

brazil_observed_hosp_summary <-
  brazil_observed_hosp_by_age |>
  dplyr::transmute(
    burden_age_group,
    estimate = observed_hosp_rate,
    lower,
    upper,
    source = "Brazil SINAN observed"
  )

external_hosp_comparison_summary <-
  external_hosp_age_summary |>
  dplyr::transmute(
    burden_age_group,
    estimate,
    lower,
    upper,
    source = "hosp rate (thiago)"
  )

hosp_age_source_comparison <-
  dplyr::bind_rows(
    brazil_observed_hosp_summary,
    external_hosp_comparison_summary
  ) |>
  dplyr::mutate(
    source = factor(
      source,
      levels = c(
        "Brazil SINAN observed",
        "hosp rate (thiago)"
      )
    ),
    burden_age_group = factor(
      as.character(burden_age_group),
      levels = age_levels
    )
  )

p_hosp_age_source_comparison <-
  ggplot2::ggplot(
    hosp_age_source_comparison,
    ggplot2::aes(
      x = burden_age_group,
      y = estimate * 100,
      fill = source
    )
  ) +
  ggplot2::geom_col(
    position = ggplot2::position_dodge(
      width = 0.78
    ),
    width = 0.68,
    colour = "white",
    linewidth = 0.2
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = lower * 100,
      ymax = upper * 100
    ),
    position = ggplot2::position_dodge(
      width = 0.78
    ),
    width = 0.14,
    linewidth = 0.5
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      accuracy = 0.1,
      suffix = "%"
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    title =
      "Age-specific hospitalisation probability",
    x = "Age group, years",
    y = "Hospitalisation probability",
    fill = NULL
  ) +
  theme_lancet_clean(
    base_size = 10
  ) +
  ggplot2::theme(
    legend.position = "bottom",
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )

p_hosp_age_source_comparison

