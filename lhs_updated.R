## lhs_updated.R
## Updated version of lhs_samples.R's LHS parameter draws.
##
## `lhs_sample_young`, `lhs_old`, and `le_sample` are reproduced unchanged
## from lhs_samples.R (self-contained here so this script doesn't depend on
## lhs_samples.R's earlier, fragile `overall`/`hosp_chikv`/`non_hosp_chikv`
## objects).
##
## `hosp_sample` and `fatal_sample` are newly derived from raw count data
## (MainData/hosp_chikvna.RDS, MainData/death_chikvna.RDS) using exact
## binomial variance (p(1-p)/n) per 10-year age band — no assumption
## needed, since these come with real case/denominator counts.
##
## `nh_fatal_sample` has no equivalent count data (death_chikvna.RDS is
## in-hospital deaths only, confirmed 2026-07-15), so it still uses the
## KEY METHODOLOGICAL ASSUMPTION from the first pass: chikv_fatal_hosp_rate.RData
## (objects `hosp`, `fatal`, `nh_fatal`, 20 age-group point estimates, no
## CI/variance) supplies the point estimate, and the relative uncertainty
## (CV = sd/mean) is carried over from the old hardcoded beta parameters
## in lhs_samples.R (~line 401-409). This is a judgment call, not a
## statistical fact — inspect the nh_fatal_sample histograms below and
## flag if a real non-hospitalised death count dataset becomes available.

library(lhs)
library(dplyr)

## ---------------------------------------------------------------------
## 1. lhs_sample_young / lhs_old — unchanged from lhs_samples.R
## ---------------------------------------------------------------------
set.seed(123)
runs <- 1000
A <- randomLHS(n = runs, k = 27)

lhs_sample <- matrix(nrow = nrow(A), ncol = ncol(A))

lhs_sample[,1]   <- qbeta(p = A[,1], shape1 = 49.14034, shape2 = 34.14837, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,2]   <- qbeta(p = A[,2], shape1 = 34.21298, shape2 = 31.963, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,3]   <- qbeta(p = A[,3], shape1 = 36.77819, shape2 = 32.7458, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,4]   <- qbeta(p = A[,4], shape1 = 35.84287, shape2 = 32.55955, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,5]   <- qbeta(p = A[,5], shape1 = 539.2823, shape2 = 14152.25, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,6]   <- qbeta(p = A[,6], shape1 = 58.96698, shape2 = 1415.207, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,7]   <- qbeta(p = A[,7], shape1 = 115.4225, shape2 = 111.383, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,8]   <- qlnorm(p = A[,8], meanlog = 2.302585, sdlog = 0.06044082, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,9]   <- qbeta(p = A[,9], shape1 = 4.581639, shape2 = 9.87148, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,10]  <- qlnorm(p = A[,10], meanlog = -0.6301724, sdlog = 0.0852, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,11]  <- qbeta(p = A[,11], shape1 = 22.51835, shape2 = 146.7926, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,12]  <- qlnorm(p = A[,12], meanlog =  -3.734278, sdlog =  0.3877361, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,13]  <- qbeta(p = A[,13], shape1 = 21.45106, shape2 = 399.158, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,14]  <- qlnorm(p = A[,14], meanlog = -4.108138, sdlog =  0.5998406, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,15]  <- qbeta(p = A[,15], shape1 = 393.9252, shape2 = 547.0268, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,16]  <- qbeta(p = A[,16], shape1 = 17875.92, shape2 = 42754.63, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,17]  <- qbeta(p = A[,17], shape1 = 799.2911, shape2 = 3306.976, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,18]  <- qbeta(p = A[,18], shape1 = 77.56885, shape2 = 836.687, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,19]  <- qbeta(p = A[,19], shape1 = 7.605944, shape2 = 1074.944, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,20]  <- qlnorm(p = A[,20], meanlog = -2.145581, sdlog =  0.1815621, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,21]  <- qlnorm(p = A[,21], meanlog = -0.5430045, sdlog =  0.154684, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,22]  <- qlnorm(p = A[,22], meanlog = -2.262311, sdlog =  0.4746817, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,23]  <- qlnorm(p = A[,23], meanlog = -1.148854, sdlog = 0.1815042, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,24]  <- qlnorm(p = A[,24], meanlog =  -0.6931472, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,25]  <- qlnorm(p = A[,25], meanlog = 0, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,26]  <- qlnorm(p = A[,26], meanlog = 0.6931472, sdlog = 0.08380206, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,27]  <- qbeta(p = A[,27], shape1 = 164.6044, shape2 = 618335.4, ncp=0, lower.tail = TRUE, log.p = FALSE)

lhs_cols <- c("symp_asia", "symp_africa", "symp_america", "symp_overall", "fatal_hosp", "hosp",
              "lt", "le_lost", "dw_chronic", "dur_chronic", "dw_hosp", "dur_acute", "dw_nonhosp", "dur_nonhosp",
              "acute", "subac", "chr6m", "chr12m", "chr30m",
              "dw_chronic_mild", "dw_chronic_severe", "dur_subac", "dw_subac",
              "dur_6m", "dur_12m", "dur_30m", "fatal_nonhosp")

colnames(lhs_sample) <- lhs_cols
lhs_sample_young <- as.data.frame(lhs_sample)

## Differential chronic burden for age 40+ (only columns 15-19 differ)
lhs_sample <- matrix(nrow = nrow(A), ncol = ncol(A))
lhs_sample[,1]   <- qbeta(p = A[,1], shape1 = 49.14034, shape2 = 34.14837, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,2]   <- qbeta(p = A[,2], shape1 = 34.21298, shape2 = 31.963, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,3]   <- qbeta(p = A[,3], shape1 = 36.77819, shape2 = 32.7458, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,4]   <- qbeta(p = A[,4], shape1 = 35.84287, shape2 = 32.55955, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,5]   <- qbeta(p = A[,5], shape1 = 539.2823, shape2 = 14152.25, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,6]   <- qbeta(p = A[,6], shape1 = 58.96698, shape2 = 1415.207, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,7]   <- qbeta(p = A[,7], shape1 = 115.4225, shape2 = 111.383, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,8]   <- qlnorm(p = A[,8], meanlog = 2.302585, sdlog = 0.06044082, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,9]   <- qbeta(p = A[,9], shape1 = 4.581639, shape2 = 9.87148, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,10]  <- qlnorm(p = A[,10], meanlog = -0.6301724, sdlog = 0.0852, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,11]  <- qbeta(p = A[,11], shape1 = 22.51835, shape2 = 146.7926, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,12]  <- qlnorm(p = A[,12], meanlog =  -3.734278, sdlog =  0.3877361, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,13]  <- qbeta(p = A[,13], shape1 = 21.45106, shape2 = 399.158, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,14]  <- qlnorm(p = A[,14], meanlog = -4.108138, sdlog =  0.5998406, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,15]  <- qbeta(p = A[,15], shape1 = 388.2874, shape2 = 742.4988, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,16]  <- qbeta(p = A[,16], shape1 = 2487.085, shape2 = 6450.815, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,17]  <- qbeta(p = A[,17], shape1 = 5998.872, shape2 = 21654.86, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,18]  <- qbeta(p = A[,18], shape1 = 184.6674, shape2 = 1216.058, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,19]  <- qbeta(p = A[,19], shape1 = 15.74654, shape2 = 516.3419, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,20]  <- qlnorm(p = A[,20], meanlog = -2.145581, sdlog =  0.1815621, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,21]  <- qlnorm(p = A[,21], meanlog = -0.5430045, sdlog =  0.154684, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,22]  <- qlnorm(p = A[,22], meanlog = -2.262311, sdlog =  0.4746817, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,23]  <- qlnorm(p = A[,23], meanlog = -1.148854, sdlog = 0.1815042, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,24]  <- qlnorm(p = A[,24], meanlog =  -0.6931472, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,25]  <- qlnorm(p = A[,25], meanlog = 0, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,26]  <- qlnorm(p = A[,26], meanlog = 0.6931472, sdlog = 0.08380206, lower.tail = TRUE, log.p = FALSE)
lhs_sample[,27]  <- qbeta(p = A[,27], shape1 = 164.6044, shape2 = 618335.4, ncp=0, lower.tail = TRUE, log.p = FALSE)

colnames(lhs_sample) <- lhs_cols
lhs_old <- as.data.frame(lhs_sample)

## ---------------------------------------------------------------------
## 2. le_sample 
## ---------------------------------------------------------------------
set.seed(123)
D <- randomLHS(n = runs, k = 9)

le_sample <- matrix(nrow = nrow(D), ncol = ncol(D))
le_sample[,1] <- qlnorm(p = D[,1], meanlog = 4.26127, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample[,2] <- qlnorm(p = D[,2], meanlog = 4.119037, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample[,3] <- qlnorm(p = D[,3], meanlog = 3.955082, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample[,4] <- qlnorm(p = D[,4], meanlog = 3.758872, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample[,5] <- qlnorm(p = D[,5], meanlog = 3.523415, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample[,6] <- qlnorm(p = D[,6], meanlog = 3.230804, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample[,7] <- qlnorm(p = D[,7], meanlog = 2.867899, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample[,8] <- qlnorm(p = D[,8], meanlog = 2.397895, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample[,9] <- qlnorm(p = D[,9], meanlog = 1.808289, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)

colnames(le_sample) <- paste0("le_", 1:9)
le_sample <- as.data.frame(le_sample)

## ---------------------------------------------------------------------
## 3. hosp_sample / fatal_sample / nh_fatal_sample
##
##    hosp_sample  <- binomial variance from raw hospitalisation counts
##                    (MainData/hosp_chikvna.RDS)
##    fatal_sample <- binomial variance from raw in-hospital death counts
##                    (MainData/death_chikvna.RDS is hospitalised deaths
##                    only, per user confirmation 2026-07-15)
##    nh_fatal_sample <- UNCHANGED from the earlier CV-carryover approach
##                    (chikv_fatal_hosp_rate.RData's `nh_fatal`, since no
##                    non-hospitalised death count data exists yet)
## ---------------------------------------------------------------------

beta_calc <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu^2
  beta  <- alpha * (1 / mu - 1)
  list(alpha = alpha, beta = beta)
}

beta_mean <- function(shape1, shape2) shape1 / (shape1 + shape2)
beta_var  <- function(shape1, shape2) (shape1 * shape2) / ((shape1 + shape2)^2 * (shape1 + shape2 + 1))

## Binomial variance = p(1-p)/n where p = cases/n. Exact, data-driven —
## not a carried-over CV assumption like nh_fatal_sample below.
binomial_var <- function(numerator, denominator) {
  p <- numerator / denominator
  (p * (1 - p)) / denominator
}

## Single source of truth for the 20-group -> 9-band mapping (band widths
## are irregular: bands 1-2 hold 3 of the 20 groups each, bands 3-9 hold
## 2 each — derived from the age_groups midpoints below, NOT an even
## 20/9 split).
age_groups <- c(mean(0:1), mean(1:4), mean(5:9), mean(10:11), mean(12:17), mean(18:19),
                 mean(20:24), mean(25:29), mean(30:34), mean(35:39), mean(40:44), mean(45:49),
                 mean(50:54), mean(55:59), mean(60:64), mean(65:69), mean(70:74), mean(75:79),
                 mean(80:84), mean(85:89))
band <- pmin(floor(age_groups / 10) + 1, 9L)
age_group_to_band <- data.frame(age_group = 1:20, band_num = band)

## --- nh_fatal_sample: CV carried over from the old hardcoded beta params
## (lhs_samples.R lines ~401-409), applied to chikv_fatal_hosp_rate.RData's
## `nh_fatal` averaged down to the 9 bands. ---
load("MainData/chikv_fatal_hosp_rate.RData")  # -> hosp, fatal, nh_fatal (unused: hosp, fatal)
nh_fatal_mean_9 <- as.numeric(tapply(nh_fatal, band, mean))

old_nh_fatal_shape1 <- c(9.693931, 6.06979, 2.602814, 9.693116, 3.444061, 16.18857, 19.94925, 27.53439, 39.97903)
old_nh_fatal_shape2 <- c(30950.3, 54677.43, 70752.94, 93187.19, 76516.7, 79369.98, 54639.27, 28387.07, 9914.798)
nh_fatal_cv <- sqrt(beta_var(old_nh_fatal_shape1, old_nh_fatal_shape2)) / beta_mean(old_nh_fatal_shape1, old_nh_fatal_shape2)
nh_fatal_var_9 <- (nh_fatal_cv * nh_fatal_mean_9)^2
nh_fatal_beta_params <- mapply(beta_calc, nh_fatal_mean_9, nh_fatal_var_9, SIMPLIFY = FALSE)

## NOTE: hosp_raw/death_raw's `age_group` is single-year age (0, 1, 2, ...
## up to ~100+ — 101 distinct values), NOT the same 20-index scheme as
## chikv_fatal_hosp_rate.RData. Band it directly via decade binning
## (do not reuse `age_group_to_band`/`band` above, which is for the other
## dataset's 20-group index and would silently fail to match here).
raw_age_to_band <- function(age_group) pmin(floor(age_group / 10) + 1, 9L)

## --- hosp_sample: binomial variance from raw hospitalisation counts ---
hosp_raw <- readRDS("MainData/hosp_chikvna.RDS")  # age_group, hosp (logical), n, total_n already present

hosp_band <- hosp_raw %>%
  group_by(age_group) %>%
  summarise(
    hosp_n = sum(n[hosp], na.rm = TRUE),
    denom_n = max(total_n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(band_num = raw_age_to_band(age_group)) %>%
  group_by(band_num) %>%
  summarise(
    hosp_n = sum(hosp_n),
    denom_n = sum(denom_n),
    hosp_rate = hosp_n / denom_n,
    .groups = "drop"
  ) %>%
  arrange(band_num) %>%
  mutate(hosp_var = binomial_var(hosp_n, denom_n))

hosp_beta_params <- mapply(beta_calc, hosp_band$hosp_rate, hosp_band$hosp_var, SIMPLIFY = FALSE)

## --- fatal_sample: binomial variance from raw in-hospital death counts ---
death_raw <- readRDS("MainData/death_chikvna.RDS") %>%
  rename(death = deaht)  # source data's column is misspelled "deaht"

death_band <- death_raw %>%
  group_by(age_group) %>%
  summarise(
    death_n = sum(n[death], na.rm = TRUE),
    denom_n = sum(n, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(band_num = raw_age_to_band(age_group)) %>%
  group_by(band_num) %>%
  summarise(
    death_n = sum(death_n),
    denom_n = sum(denom_n),
    death_rate = death_n / denom_n,
    .groups = "drop"
  ) %>%
  arrange(band_num) %>%
  mutate(death_var = binomial_var(death_n, denom_n))

fatal_beta_params <- mapply(beta_calc, death_band$death_rate, death_band$death_var, SIMPLIFY = FALSE)

cat("Hospitalisation by band:\n"); print(hosp_band)
cat("\nIn-hospital death by band:\n"); print(death_band)

## --- Generate the three 1,000-run LHS samples ---
set.seed(123)
B <- randomLHS(n = runs, k = 9)
hosp_sample <- matrix(nrow = runs, ncol = 9)
for (i in 1:9) {
  hosp_sample[, i] <- qbeta(p = B[, i], shape1 = hosp_beta_params[[i]]$alpha, shape2 = hosp_beta_params[[i]]$beta)
}
colnames(hosp_sample) <- paste0("hosp_", 1:9)
hosp_sample <- as.data.frame(hosp_sample)

set.seed(123)
C <- randomLHS(n = runs, k = 9)
fatal_sample <- matrix(nrow = runs, ncol = 9)
for (i in 1:9) {
  fatal_sample[, i] <- qbeta(p = C[, i], shape1 = fatal_beta_params[[i]]$alpha, shape2 = fatal_beta_params[[i]]$beta)
}
colnames(fatal_sample) <- paste0("fatal_", 1:9)
fatal_sample <- as.data.frame(fatal_sample)

nh_fatal_sample <- matrix(nrow = runs, ncol = 9)
for (i in 1:9) {
  nh_fatal_sample[, i] <- qbeta(p = C[, i], shape1 = nh_fatal_beta_params[[i]]$alpha, shape2 = nh_fatal_beta_params[[i]]$beta)
}
colnames(nh_fatal_sample) <- paste0("nh_fatal_", 1:9)
nh_fatal_sample <- as.data.frame(nh_fatal_sample)

## ---------------------------------------------------------------------
## 4. Diagnostics — inspect these before trusting the new samples.
##    hosp_sample/fatal_sample are now exact (binomial variance from raw
##    counts), so the histogram should center tightly on the target rate.
##    nh_fatal_sample is still the carried-over-CV assumption — sanity
##    check its spread looks reasonable per band.
## ---------------------------------------------------------------------
for (i in 1:9) {
  hist(hosp_sample[[i]], main = paste("hosp_sample, band", i, "- target rate", round(hosp_band$hosp_rate[i], 4)), xlab = "hospitalisation rate")
  print(quantile(hosp_sample[[i]], c(0.025, 0.5, 0.975)))
}
for (i in 1:9) {
  hist(fatal_sample[[i]], main = paste("fatal_sample, band", i, "- target rate", round(death_band$death_rate[i], 4)), xlab = "fatal-if-hospitalised rate")
  print(quantile(fatal_sample[[i]], c(0.025, 0.5, 0.975)))
}
for (i in 1:9) {
  hist(nh_fatal_sample[[i]], main = paste("nh_fatal_sample, band", i, "- target mean", round(nh_fatal_mean_9[i], 4)), xlab = "fatal-if-non-hospitalised rate")
  print(quantile(nh_fatal_sample[[i]], c(0.025, 0.5, 0.975)))
}

## ---------------------------------------------------------------------
## 5. Save — written to new "_v2" filenames so the existing pipeline
## (open_data.R -> age_strat_comorb_burden.R) is NOT affected until you
## decide to promote these. lhs_sample_young/lhs_old/le_sample are
## identical to what's already in MainData/, so not re-saved here.
## ---------------------------------------------------------------------
save(hosp_sample,     file = "MainData/hosp_sample_v2.RData")
save(fatal_sample,    file = "MainData/fatal_sample_v2.RData")
save(nh_fatal_sample, file = "MainData/nh_fatal_sample_v2.RData")
