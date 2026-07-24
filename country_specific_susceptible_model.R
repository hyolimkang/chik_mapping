library(readxl)
library(dplyr)
library(stringr)

lookup_file <-
  "Data/chikungunya_intro_years_model_ready_updated.xlsx"

country_model_lookup_raw <-
  read_excel(
    path = lookup_file,
    sheet = "R model lookup"
  )

load("MainData/foi_comb_all_0707.RData")

country_model_lookup <-
  country_model_lookup_raw |>
  transmute(
    country = str_squish(country),
    
    transmission_type =
      str_squish(transmission_type),
    
    model_class =
      str_squish(model_class),
    
    model_start_year =
      as.integer(model_start_year),
    
    serology_year =
      as.integer(serology_year),
    
    single_start_suitability =
      str_squish(single_start_suitability),
    
    confidence =
      str_squish(confidence)
  )

country_model_lookup <-
  country_model_lookup |>
  mutate(
    transmission_type = case_when(
      country %in% c("France", "Italy") ~
        "Recurrent episodic autochthonous outbreaks",
      
      TRUE ~ transmission_type
    ),
    
    model_class = case_when(
      country %in% c("France", "Italy") ~
        "C — Event-specific immunity",
      
      TRUE ~ model_class
    ),
    
    model_start_year = case_when(
      country %in% c("France", "Italy") ~
        NA_integer_,
      
      TRUE ~ model_start_year
    )
  )

country_model_lookup <-
  country_model_lookup |>
  mutate(
    outbreak_years = case_when(
      country == "France" ~
        "2010;2014;2017;2024",
      
      country == "Italy" ~
        "2007;2017",
      
      TRUE ~ NA_character_
    )
  )


foi_comb_all_clean <- foi_comb_all |>
  mutate(
    country = str_squish(country)
  )

foi_comb_all_clean <- foi_comb_all |>
  sf::st_drop_geometry() |>
  mutate(
    country = stringr::str_squish(country)
  )


country_model_lookup_clean <- country_model_lookup |>
  mutate(
    country = str_squish(country)
  )


n_before_join <- nrow(foi_comb_all_clean)

foi_comb_model <- foi_comb_all_clean |>
  left_join(
    country_model_lookup_clean |>
      select(
        country,
        transmission_type,
        model_class,
        model_start_year,
        serology_year,
        single_start_suitability,
        confidence,
        outbreak_years
      ),
    by = "country",
    relationship = "many-to-one"
  )

foi_comb_model_nogeo <- as.data.frame(foi_comb_model)

sfc_columns <- vapply(
  foi_comb_model_nogeo,
  function(x) inherits(x, "sfc"),
  logical(1)
)

names(foi_comb_model_nogeo)[sfc_columns]

foi_comb_model_nogeo <-
  foi_comb_model_nogeo[
    ,
    !sfc_columns,
    drop = FALSE
  ]

df_for_burden <- foi_comb_model_nogeo

reference_year <- 2025L

df_for_burden <- df_for_burden |>
  dplyr::mutate(
    model_start_year = dplyr::case_when(
      country == "France" & model_code == "C" ~ 2010L,
      country == "Italy"  & model_code == "C" ~ 2007L,
      TRUE ~ as.integer(model_start_year)
    ),
    
    exposure_duration = dplyr::case_when(
      model_code %in% c("A", "C") &
        !is.na(model_start_year) ~
        reference_year - model_start_year,
      
      model_code == "B" ~ Inf,
      
      TRUE ~ NA_real_
    )
  )

### End of create DF fro burden -------------------------------------------------------



## Functions-------------------------------------------------------
################################################################################
calc_incidence_equilibrium <- function(
    FOI,
    l_lim,
    u_lim
) {
  
  if (is.na(FOI)) {
    return(NA_real_)
  }
  
  if (FOI < 0) {
    stop("FOI must be non-negative.")
  }
  
  if (u_lim <= l_lim) {
    stop("u_lim must be greater than l_lim.")
  }
  
  (
    exp(-FOI * l_lim) -
      exp(-FOI * u_lim)
  ) / (u_lim - l_lim)
}


calc_incidence_finite <- function(
    FOI,
    l_lim,
    u_lim,
    exposure_duration
) {
  
  if (
    is.na(FOI) ||
    is.na(exposure_duration)
  ) {
    return(NA_real_)
  }
  
  if (FOI < 0) {
    stop("FOI must be non-negative.")
  }
  
  if (u_lim <= l_lim) {
    stop("u_lim must be greater than l_lim.")
  }
  
  if (exposure_duration < 0) {
    stop("exposure_duration must be non-negative.")
  }
  
  D <- exposure_duration
  
  # Entire age band is younger than the circulation duration.
  if (u_lim <= D) {
    
    return(
      (
        exp(-FOI * l_lim) -
          exp(-FOI * u_lim)
      ) / (u_lim - l_lim)
    )
  }
  
  # Entire age band is older than the circulation duration.
  if (l_lim >= D) {
    
    return(
      FOI * exp(-FOI * D)
    )
  }
  
  # Age band crosses the circulation-duration boundary.
  part_before_D <-
    exp(-FOI * l_lim) -
    exp(-FOI * D)
  
  part_after_D <-
    (u_lim - D) *
    FOI *
    exp(-FOI * D)
  
  (
    part_before_D +
      part_after_D
  ) / (u_lim - l_lim)
}

calc_incidence_by_model <- function(
    FOI,
    l_lim,
    u_lim,
    model_code,
    exposure_duration = NA_real_
) {
  
  if (is.na(model_code)) {
    return(NA_real_)
  }
  
  if (model_code == "A") {
    
    return(
      calc_incidence_finite(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim,
        exposure_duration = exposure_duration
      )
    )
  }
  
  if (model_code == "B") {
    
    return(
      calc_incidence_equilibrium(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim
      )
    )
  }
  
  # C: event-specific model을 나중에 별도로 구현
  # D: primary endemic burden 계산에서 제외
  if (model_code %in% c("C", "D")) {
    return(NA_real_)
  }
  
  stop(
    paste0(
      "Unknown model_code: ",
      model_code
    )
  )
}

###############################################################################
l_lim <- c(
  0, 1, 5, 10, 15, 20, 25, 30, 35,
  40, 45, 50, 55, 60, 65, 70, 75, 80
)

u_lim <- c(
  1, 5, 10, 15, 20, 25, 30, 35, 40,
  45, 50, 55, 60, 65, 70, 75, 80, 100
)

age_pop_cols <- names(df_for_burden)[7:24]

age_pop_cols
length(age_pop_cols)

process_one_foi_draw <- function(
    df,
    foi_col,
    age_pop_cols,
    l_lim,
    u_lim
) {
  
  n_rows <- nrow(df)
  n_age  <- length(l_lim)
  
  infections_per_age_band <- matrix(
    NA_real_,
    nrow = n_rows,
    ncol = n_age
  )
  
  incidence_per_age_band <- matrix(
    NA_real_,
    nrow = n_rows,
    ncol = n_age
  )
  
  total_infection <- rep(
    NA_real_,
    n_rows
  )
  
  for (i in seq_len(n_rows)) {
    
    foi_i   <- as.numeric(df[[foi_col]][i])
    model_i <- as.character(df$model_code[i])
    D_i     <- as.numeric(df$exposure_duration[i])
    
    # FOI 또는 model class가 없는 경우
    if (
      is.na(foi_i) ||
      is.na(model_i)
    ) {
      next
    }
    
    # C와 D는 아직 catalytic calculation을 적용하지 않음
    if (model_i %in% c("C", "D")) {
      next
    }
    
    incidence_rates <- mapply(
      FUN = function(lower_age, upper_age) {
        
        calc_incidence_by_model(
          FOI = foi_i,
          l_lim = lower_age,
          u_lim = upper_age,
          model_code = model_i,
          exposure_duration = D_i
        )
      },
      lower_age = l_lim,
      upper_age = u_lim
    )
    
    age_group_pop <- as.numeric(
      df[i, age_pop_cols, drop = TRUE]
    )
    
    age_group_pop[is.na(age_group_pop)] <- 0
    
    infection_numbers <-
      incidence_rates * age_group_pop
    
    incidence_per_age_band[i, ] <-
      incidence_rates
    
    infections_per_age_band[i, ] <-
      infection_numbers
    
    total_infection[i] <-
      sum(infection_numbers)
  }
  
  updated_df <- df
  updated_df$total_infection <- total_infection
  
  list(
    updated_df = updated_df,
    infection_per_band = infections_per_age_band,
    incidence_per_band = incidence_per_age_band
  )
}

test_result <- process_one_foi_draw(
  df = df_for_burden,
  foi_col = names(df_for_burden)[26],
  age_pop_cols = age_pop_cols,
  l_lim = l_lim,
  u_lim = u_lim
)

age_band_labels <- paste0(
  "[", l_lim, ",", u_lim,
  ifelse(u_lim == 100, "]", ")")
)

incidence_long <- as.data.frame(
  test_result$incidence_per_band
)

colnames(incidence_long) <- age_band_labels

incidence_long <- incidence_long |>
  tibble::as_tibble() |>
  dplyr::mutate(
    row_id = dplyr::row_number(),
    country = df_for_burden$country,
    model_code = df_for_burden$model_code
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(age_band_labels),
    names_to = "age_band",
    values_to = "incidence"
  ) |>
  dplyr::mutate(
    age_band = factor(
      age_band,
      levels = age_band_labels
    )
  )

incidence_ab_summary <- incidence_long |>
  dplyr::filter(
    model_code %in% c("A", "B")
  ) |>
  dplyr::group_by(
    model_code,
    age_band
  ) |>
  dplyr::summarise(
    median_incidence = median(
      incidence,
      na.rm = TRUE
    ),
    lower = quantile(
      incidence,
      0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      incidence,
      0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


p_incidence_ab <- ggplot2::ggplot(
  incidence_ab_summary,
  ggplot2::aes(
    x = age_band,
    y = median_incidence,
    colour = model_code,
    group = model_code
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.12,
    linewidth = 0.4
  ) +
  ggplot2::geom_line(
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    size = 2.2
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
    title = "Age-specific incidence curve by model type",
    x = "Age band",
    y = "Annual incidence",
    colour = "Model"
  ) +
  theme_lancet_clean(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

p_incidence_ab

infection_long <- as.data.frame(
  test_result$infection_per_band
)

colnames(infection_long) <- age_band_labels

infection_long <- infection_long |>
  tibble::as_tibble() |>
  dplyr::mutate(
    row_id = dplyr::row_number(),
    country = df_for_burden$country,
    model_code = df_for_burden$model_code
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(age_band_labels),
    names_to = "age_band",
    values_to = "infections"
  ) |>
  dplyr::mutate(
    age_band = factor(
      age_band,
      levels = age_band_labels
    )
  )

infection_ab_summary <- infection_long |>
  dplyr::filter(
    model_code %in% c("A", "B")
  ) |>
  dplyr::group_by(
    model_code,
    age_band
  ) |>
  dplyr::summarise(
    median_infections = median(
      infections,
      na.rm = TRUE
    ),
    lower = quantile(
      infections,
      0.025,
      na.rm = TRUE
    ),
    upper = quantile(
      infections,
      0.975,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

p_infections_ab <- ggplot2::ggplot(
  infection_ab_summary,
  ggplot2::aes(
    x = age_band,
    y = median_infections,
    colour = model_code,
    group = model_code
  )
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(
      ymin = lower,
      ymax = upper
    ),
    width = 0.12,
    linewidth = 0.4
  ) +
  ggplot2::geom_line(
    linewidth = 0.9
  ) +
  ggplot2::geom_point(
    size = 2.2
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(
      scale_cut = scales::cut_short_scale()
    ),
    expand = ggplot2::expansion(
      mult = c(0, 0.08)
    )
  ) +
  ggplot2::labs(
    title = "Age-specific infection numbers by model type",
    x = "Age band",
    y = "Infections",
    colour = "Model"
  ) +
  theme_lancet_clean(
    base_size = 10
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    )
  )

p_infections_ab

# Functions for episodic (Typd C) countries
calc_incidence_episodic <- function(
    FOI,
    l_lim,
    u_lim,
    exposure_duration,
    tau
) {
  
  # -----------------------------
  # 1. Input checks
  # -----------------------------
  if (
    is.na(FOI) ||
    is.na(exposure_duration) ||
    is.na(tau)
  ) {
    return(NA_real_)
  }
  
  if (FOI < 0) {
    stop("FOI must be non-negative.")
  }
  
  if (u_lim <= l_lim) {
    stop("u_lim must be greater than l_lim.")
  }
  
  if (exposure_duration < 0) {
    stop("exposure_duration must be non-negative.")
  }
  
  if (tau <= 0) {
    stop("tau must be greater than zero.")
  }
  
  D <- exposure_duration
  
  # -----------------------------
  # 2. Number of effective episodes
  # -----------------------------
  # The introduction year is treated as the first episode.
  K_total <- 1L + floor(D / tau)
  
  # Episode times measured as years since introduction:
  # 0, tau, 2*tau, ...
  episode_times_since_intro <- seq(
    from = 0,
    by = tau,
    length.out = K_total
  )
  
  # Age in the current year required to have been alive
  # during each episode.
  #
  # Example:
  # D = 11, tau = 7
  # episodes occurred 11 and 4 years before the present
  # age thresholds = 11 and 4.
  episode_age_thresholds <-
    D - episode_times_since_intro
  
  # -----------------------------
  # 3. Split the age band at episode thresholds
  # -----------------------------
  thresholds_inside_band <-
    episode_age_thresholds[
      episode_age_thresholds > l_lim &
        episode_age_thresholds < u_lim
    ]
  
  break_points <- sort(
    unique(
      c(
        l_lim,
        thresholds_inside_band,
        u_lim
      )
    )
  )
  
  interval_widths <- diff(break_points)
  
  interval_midpoints <-
    head(break_points, -1) +
    interval_widths / 2
  
  # -----------------------------
  # 4. Episodes experienced in each age segment
  # -----------------------------
  K_experienced <- vapply(
    interval_midpoints,
    function(age) {
      sum(age >= episode_age_thresholds)
    },
    integer(1)
  )
  
  # -----------------------------
  # 5. Susceptibility and incidence
  # -----------------------------
  susceptibility <- exp(
    -FOI * D *
      (K_experienced / K_total)
  )
  
  incidence_by_interval <-
    FOI * susceptibility
  
  # Width-weighted mean incidence over the full age band
  incidence_age_band <-
    sum(
      incidence_by_interval *
        interval_widths
    ) /
    (u_lim - l_lim)
  
  return(incidence_age_band)
}

# expand func
calc_incidence_by_model <- function(
    FOI,
    l_lim,
    u_lim,
    model_code,
    exposure_duration = NA_real_,
    tau = NA_real_
) {
  
  if (
    length(model_code) != 1L ||
    is.na(model_code)
  ) {
    return(NA_real_)
  }
  
  model_code <- as.character(model_code)
  
  if (model_code == "A") {
    
    return(
      calc_incidence_finite(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim,
        exposure_duration = exposure_duration
      )
    )
  }
  
  if (model_code == "B") {
    
    return(
      calc_incidence_equilibrium(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim
      )
    )
  }
  
  if (model_code == "C") {
    
    return(
      calc_incidence_episodic(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim,
        exposure_duration = exposure_duration,
        tau = tau
      )
    )
  }
  
  # D remains unimplemented at this stage.
  if (model_code == "D") {
    return(NA_real_)
  }
  
  stop(
    paste0(
      "Unknown model_code: ",
      model_code
    )
  )
}

## LHS sample for inter-epidemic periods (to estimate number of episodes after introduction)
set.seed(1234)

n_draws <- 100

tau_lhs_unit <- lhs::randomLHS(
  n = n_draws,
  k = 1
)[, 1]

tau_draws <- 7 + tau_lhs_unit * (20 - 7)

tau_psa <- tibble::tibble(
  draw_id = seq_len(n_draws),
  tau = tau_draws
)

foi_cols <- names(foi_comb_all)[26:125]

stopifnot(length(foi_cols) == nrow(tau_psa))

tau_psa <- tau_psa |>
  dplyr::mutate(
    foi_col = foi_cols
  )

process_one_foi_draw <- function(
    df,
    foi_col,
    tau_i,
    age_pop_cols,
    l_lim,
    u_lim
) {
  
  # -----------------------------
  # 1. Input checks
  # -----------------------------
  if (length(tau_i) != 1L) {
    stop("tau_i must be a single value for one PSA draw.")
  }
  
  if (
    any(df$model_code == "C", na.rm = TRUE) &&
    (is.na(tau_i) || tau_i <= 0)
  ) {
    stop("A valid tau_i is required when model C is present.")
  }
  
  n_rows <- nrow(df)
  n_age  <- length(l_lim)
  
  infections_per_age_band <- matrix(
    NA_real_,
    nrow = n_rows,
    ncol = n_age
  )
  
  incidence_per_age_band <- matrix(
    NA_real_,
    nrow = n_rows,
    ncol = n_age
  )
  
  total_infection <- rep(
    NA_real_,
    n_rows
  )
  
  # -----------------------------
  # 2. Loop through cells
  # -----------------------------
  for (i in seq_len(n_rows)) {
    
    foi_i <- as.numeric(
      df[[foi_col]][i]
    )
    
    model_i <- as.character(
      df$model_code[i]
    )
    
    D_i <- as.numeric(
      df$exposure_duration[i]
    )
    
    if (
      is.na(foi_i) ||
      is.na(model_i)
    ) {
      next
    }
    
    # D is not yet included in the baseline burden calculation
    if (model_i == "D") {
      next
    }
    
    incidence_rates <- mapply(
      FUN = function(lower_age, upper_age) {
        
        calc_incidence_by_model(
          FOI = foi_i,
          l_lim = lower_age,
          u_lim = upper_age,
          model_code = model_i,
          exposure_duration = D_i,
          tau = tau_i
        )
      },
      lower_age = l_lim,
      upper_age = u_lim
    )
    
    age_group_pop <- as.numeric(
      unlist(
        df[
          i,
          age_pop_cols,
          drop = FALSE
        ],
        use.names = FALSE
      )
    )
    
    age_group_pop[is.na(age_group_pop)] <- 0
    
    infection_numbers <-
      incidence_rates * age_group_pop
    
    incidence_per_age_band[i, ] <-
      incidence_rates
    
    infections_per_age_band[i, ] <-
      infection_numbers
    
    total_infection[i] <-
      sum(
        infection_numbers,
        na.rm = FALSE
      )
  }
  
  updated_df <- df
  updated_df$total_infection <- total_infection
  
  list(
    updated_df = updated_df,
    infection_per_band = infections_per_age_band,
    incidence_per_band = incidence_per_age_band,
    foi_col = foi_col,
    tau = tau_i
  )
}

test_C_df <- df_for_burden |>
  dplyr::filter(
    country == "France",
    model_code == "C"
  ) |>
  dplyr::slice_head(n = 100)

### Function for D: all susceptible (no transmission evidence before 2025) -- transmission potential if introduced in year T
calc_incidence_potential <- function(
    FOI,
    l_lim,
    u_lim,
    time_horizon = 1
) {
  
  if (is.na(FOI)) {
    return(NA_real_)
  }
  
  if (FOI < 0) {
    stop("FOI must be non-negative.")
  }
  
  if (u_lim <= l_lim) {
    stop("u_lim must be greater than l_lim.")
  }
  
  if (time_horizon <= 0) {
    stop("time_horizon must be greater than zero.")
  }
  
  # All ages are fully susceptible at the start.
  annual_risk <- 1 - exp(
    -FOI * time_horizon
  )
  
  return(annual_risk)
}

## Full function (A-D)
calc_incidence_by_model <- function(
    FOI,
    l_lim,
    u_lim,
    model_code,
    exposure_duration = NA_real_,
    tau = NA_real_,
    d_time_horizon = 1
) {
  
  if (
    length(model_code) != 1L ||
    is.na(model_code)
  ) {
    return(NA_real_)
  }
  
  model_code <- as.character(model_code)
  
  # ---------------------------------
  # A: Finite transmission history
  # ---------------------------------
  if (model_code == "A") {
    
    return(
      calc_incidence_finite(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim,
        exposure_duration = exposure_duration
      )
    )
  }
  
  # ---------------------------------
  # B: Equilibrium catalytic model
  # ---------------------------------
  if (model_code == "B") {
    
    return(
      calc_incidence_equilibrium(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim
      )
    )
  }
  
  # ---------------------------------
  # C: Episodic finite history
  # ---------------------------------
  if (model_code == "C") {
    
    return(
      calc_incidence_episodic(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim,
        exposure_duration = exposure_duration,
        tau = tau
      )
    )
  }
  
  # ---------------------------------
  # D: Conditional transmission potential
  # ---------------------------------
  if (model_code == "D") {
    
    return(
      calc_incidence_potential(
        FOI = FOI,
        l_lim = l_lim,
        u_lim = u_lim,
        time_horizon = d_time_horizon
      )
    )
  }
  
  stop(
    paste0(
      "Unknown model_code: ",
      model_code
    )
  )
}


process_one_foi_draw <- function(
    df,
    foi_col,
    tau_i = NA_real_,
    age_pop_cols,
    l_lim,
    u_lim,
    d_time_horizon = 1
) {
  
  # ---------------------------------
  # 1. Input validation
  # ---------------------------------
  if (
    !foi_col %in% names(df)
  ) {
    stop(
      paste0(
        "FOI column not found: ",
        foi_col
      )
    )
  }
  
  if (
    length(l_lim) != length(u_lim)
  ) {
    stop(
      "l_lim and u_lim must have the same length."
    )
  }
  
  if (
    length(age_pop_cols) != length(l_lim)
  ) {
    stop(
      paste0(
        "Number of age population columns (",
        length(age_pop_cols),
        ") must equal number of age bands (",
        length(l_lim),
        ")."
      )
    )
  }
  
  if (
    any(!age_pop_cols %in% names(df))
  ) {
    stop(
      "Some age population columns are missing from df."
    )
  }
  
  if (
    any(df$model_code == "C", na.rm = TRUE)
  ) {
    
    if (
      length(tau_i) != 1L ||
      is.na(tau_i) ||
      tau_i <= 0
    ) {
      stop(
        "A single positive tau_i is required when model C is present."
      )
    }
  }
  
  if (
    length(d_time_horizon) != 1L ||
    is.na(d_time_horizon) ||
    d_time_horizon <= 0
  ) {
    stop(
      "d_time_horizon must be a single positive value."
    )
  }
  
  n_rows <- nrow(df)
  n_age  <- length(l_lim)
  
  infections_per_age_band <- matrix(
    NA_real_,
    nrow = n_rows,
    ncol = n_age
  )
  
  incidence_per_age_band <- matrix(
    NA_real_,
    nrow = n_rows,
    ncol = n_age
  )
  
  total_infection <- rep(
    NA_real_,
    n_rows
  )
  
  # ---------------------------------
  # 2. Calculate each cell
  # ---------------------------------
  for (i in seq_len(n_rows)) {
    
    foi_i <- as.numeric(
      df[[foi_col]][i]
    )
    
    model_i <- as.character(
      df$model_code[i]
    )
    
    D_i <- as.numeric(
      df$exposure_duration[i]
    )
    
    if (
      is.na(foi_i) ||
      is.na(model_i)
    ) {
      next
    }
    
    incidence_rates <- mapply(
      FUN = function(
    lower_age,
    upper_age
      ) {
        
        calc_incidence_by_model(
          FOI = foi_i,
          l_lim = lower_age,
          u_lim = upper_age,
          model_code = model_i,
          exposure_duration = D_i,
          tau = tau_i,
          d_time_horizon = d_time_horizon
        )
      },
    lower_age = l_lim,
    upper_age = u_lim
    )
    
    age_group_pop <- as.numeric(
      unlist(
        df[
          i,
          age_pop_cols,
          drop = FALSE
        ],
        use.names = FALSE
      )
    )
    
    age_group_pop[
      is.na(age_group_pop)
    ] <- 0
    
    infection_numbers <-
      incidence_rates *
      age_group_pop
    
    incidence_per_age_band[i, ] <-
      incidence_rates
    
    infections_per_age_band[i, ] <-
      infection_numbers
    
    if (
      all(
        is.na(infection_numbers)
      )
    ) {
      
      total_infection[i] <- NA_real_
      
    } else {
      
      total_infection[i] <- sum(
        infection_numbers,
        na.rm = TRUE
      )
    }
  }
  
  # ---------------------------------
  # 3. Attach interpretation metadata
  # ---------------------------------
  updated_df <- df
  
  updated_df$total_infection <-
    total_infection
  
  updated_df$estimate_type <-
    ifelse(
      updated_df$model_code == "D",
      "Conditional transmission potential",
      "Estimated annual burden"
    )
  
  updated_df$transmission_condition <- ifelse(
    updated_df$model_code == "D",
    paste0(
      "Local transmission active for ",
      d_time_horizon,
      " year(s)"
    ),
    NA_character_
  )
  
  list(
    updated_df = updated_df,
    infection_per_band =
      infections_per_age_band,
    incidence_per_band =
      incidence_per_age_band,
    foi_col = foi_col,
    tau = tau_i,
    d_time_horizon = d_time_horizon
  )
}

### Full implementation
foi_psa <- tau_psa |>
  dplyr::transmute(
    draw_id = as.integer(draw_id),
    foi_col = as.character(foi_col),
    tau = as.numeric(tau)
  ) |>
  dplyr::arrange(draw_id)

foi_batches <- split(
  foi_psa,
  ceiling(seq_len(nrow(foi_psa)) / 25)
)

sapply(foi_batches, nrow)


## ============================================================
## Run FOI batches 
## ============================================================

dir.create(
  "Outputs",
  showWarnings = FALSE,
  recursive = TRUE
)

group_key <- paste(
  df_for_burden$country,
  df_for_burden$model_code,
  sep = "|||"
)

group_lookup <- tibble::tibble(
  group_key = group_key,
  country = as.character(df_for_burden$country),
  model_code = as.character(df_for_burden$model_code)
) |>
  dplyr::distinct()

for (batch_i in 1:4) {
  
  current_batch <- foi_batches[[batch_i]]
  
  message(
    "\n========================================",
    "\nStarting batch ", batch_i,
    " of 4",
    "\nNumber of draws: ", nrow(current_batch),
    "\n========================================"
  )
  
  current_batch_results <- vector(
    "list",
    nrow(current_batch)
  )
  
  for (j in seq_len(nrow(current_batch))) {
    
    draw_j <- as.integer(
      current_batch$draw_id[j]
    )
    
    foi_j <- as.character(
      current_batch$foi_col[j]
    )
    
    tau_j <- as.numeric(
      current_batch$tau[j]
    )
    
    message(
      "Batch ", batch_i,
      " | Processing draw ", draw_j,
      " (", j, " / ", nrow(current_batch), ")",
      " | tau = ", round(tau_j, 2)
    )
    
    draw_result <- process_one_foi_draw(
      df = df_for_burden,
      foi_col = foi_j,
      tau_i = tau_j,
      age_pop_cols = age_pop_cols,
      l_lim = l_lim,
      u_lim = u_lim,
      d_time_horizon = 1
    )
    
    infection_mat <-
      draw_result$infection_per_band
    
    colnames(infection_mat) <-
      age_band_labels
    
    country_infections <- rowsum(
      infection_mat,
      group = group_key,
      reorder = FALSE,
      na.rm = TRUE
    )
    
    valid_counts <- rowsum(
      1L * !is.na(infection_mat),
      group = group_key,
      reorder = FALSE,
      na.rm = TRUE
    )
    
    country_infections[
      valid_counts == 0
    ] <- NA_real_
    
    current_batch_results[[j]] <-
      tibble::as_tibble(
        country_infections,
        .name_repair = "minimal"
      ) |>
      dplyr::mutate(
        group_key =
          rownames(country_infections),
        .before = 1
      ) |>
      dplyr::left_join(
        group_lookup,
        by = "group_key",
        relationship = "many-to-one"
      ) |>
      dplyr::select(
        -group_key
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(
          age_band_labels
        ),
        names_to = "age_band",
        values_to = "annual_infections"
      ) |>
      dplyr::mutate(
        estimate_type = dplyr::if_else(
          model_code == "D",
          "Conditional transmission potential",
          "Estimated annual burden"
        ),
        
        batch_id = batch_i,
        draw_id = draw_j,
        tau = tau_j,
        
        age_band = factor(
          age_band,
          levels = age_band_labels
        )
      ) |>
      dplyr::select(
        country,
        model_code,
        estimate_type,
        batch_id,
        draw_id,
        tau,
        age_band,
        annual_infections
      )
    
    rm(
      draw_result,
      infection_mat,
      country_infections,
      valid_counts
    )
    
    gc()
  }
  
  burden_batch_i <- dplyr::bind_rows(
    current_batch_results
  )
  
  output_file <- file.path(
    "Outputs",
    sprintf(
      "burden_ABCD_batch_%02d.rds",
      batch_i
    )
  )
  
  saveRDS(
    burden_batch_i,
    file = output_file
  )
  
  message(
    "Completed batch ", batch_i,
    " | Rows saved: ",
    format(
      nrow(burden_batch_i),
      big.mark = ","
    ),
    "\nSaved to: ",
    output_file
  )
  
  rm(
    current_batch,
    current_batch_results,
    burden_batch_i
  )
  
  gc()
}

batch_files <- file.path(
  "Outputs",
  sprintf("burden_ABCD_batch_%02d.rds", 1:4)
)


burden_ABCD_all <- dplyr::bind_rows(
  lapply(batch_files, readRDS)
)

saveRDS(
  burden_ABCD_all,
  "Outputs/burden_ABCD_all.rds"
)