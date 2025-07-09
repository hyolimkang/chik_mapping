### parameters < 40 years old 

beta_calc <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

var <- function(lo, hi){
  var <- ((hi - lo)/3.92)^2
  return(var = var)
}


expected_arthralgia_prop <- function(rate) {
  
  exp_0.5_mid <- exp(-rate * (14/30))
  exp_3_mid   <- exp(-rate * 3)
  exp_6_mid   <- exp(-rate * 6)
  exp_12_mid  <- exp(-rate * 12)
  exp_30_mid  <- exp(-rate * 30)
  
  sum_mid <- exp_0.5_mid + exp_3_mid + exp_6_mid + exp_12_mid + exp_30_mid

  exp_0.5_norm <- exp_0.5_mid / sum_mid
  exp_3_norm   <- exp_3_mid  / sum_mid
  exp_6_norm   <- exp_6_mid / sum_mid
  exp_12_norm   <- exp_12_mid / sum_mid
  exp_30_norm   <- exp_30_mid / sum_mid
  
  return(list(exp_0.5_norm = exp_0.5_norm,
              exp_3_norm   = exp_3_norm,
              exp_6_norm   = exp_6_norm,
              exp_12_norm  = exp_12_norm,
              exp_30_norm  = exp_30_norm
  ))
}
# for <40 years 
exp_rate_1 <- expected_arthralgia_prop(rate = 0.1384)
exp_rate_2 <- expected_arthralgia_prop(rate = 0.1134)
exp_rate_3 <- expected_arthralgia_prop(rate = 0.1633)

# for > 40 years
exp_rate_1 <- expected_arthralgia_prop(rate = 0.0830)
exp_rate_2 <- expected_arthralgia_prop(rate = 0.0643)
exp_rate_3 <- expected_arthralgia_prop(rate = 0.1017)


var_exp_0.5 <- ((exp_rate_3$exp_0.5_norm - exp_rate_2$exp_0.5_norm)/3.92)^2
med_exp_0.5 <- exp_rate_1$exp_0.5_norm[1]
beta(med_exp_0.5, var_exp_0.5) 

var_exp_3 <- ((exp_rate_3$exp_3_norm - exp_rate_2$exp_3_norm)/3.92)^2
med_exp_3 <- exp_rate_1$exp_3_norm[1]
beta(med_exp_3, var_exp_3) 

var_exp_6 <- ((exp_rate_2$exp_6_norm - exp_rate_3$exp_6_norm)/3.92)^2
med_exp_6 <- exp_rate_1$exp_6_norm
beta(med_exp_6, var_exp_6) 

var_exp_12 <- ((exp_rate_2$exp_12_norm - exp_rate_3$exp_12_norm)/3.92)^2
med_exp_12 <- exp_rate_1$exp_12_norm
beta(med_exp_12, var_exp_12) 

var_exp_30 <- ((exp_rate_2$exp_30_norm - exp_rate_3$exp_30_norm)/3.92)^2
med_exp_30 <- exp_rate_1$exp_30_norm
beta(med_exp_30, var_exp_30) 


med_dw_ch_mild <- 0.117
ui_dw_ch_mild  <- c(0.08, 0.163)

med_dw_ch_sev <- 0.581
ui_dw_ch_sev  <- c(0.403, 0.739)

med_dur_subac <- 38/365
ui_dur_subac  <- c(14/365, 90/365)

med_dw_subac <- 0.317
ui_dw_subac  <- c(0.216, 0.44)

med_dur_6m <- 0.5
ui_dur_6m  <- c(0.5*0.9, 0.5*1.1)

med_dur_12m <- 1
ui_dur_12m  <- c(1*0.9, 1*1.1)

med_dur_30m <- 2
ui_dur_30m  <- c(2*0.9, (30/12))


fatal_hosp <- 0.03670701
ui_fatal_hosp  <- c(0.03378602, 0.03986725)
var_fatal_hosp <- ((ui_fatal_hosp[2] - ui_fatal_hosp[1])/3.92)^2

beta(fatal_hosp, var_fatal_hosp)

fatal_nonhosp <- 0.00026613481
ui_fatal_nonhosp  <- c(0.00022851106, 0.0003098145)
var_fatal_nonhosp <- ((ui_fatal_nonhosp[2] - ui_fatal_nonhosp[1])/3.92)^2

beta(fatal_nonhosp, var_fatal_nonhosp)

# Calculate mu (Log Mean)
calc_log <- function(mid, ui) {
  
  mu <- log(mid)
  
  sigma <- (log(ui[2]) - log(ui[1])) / (2 * 1.96)

  return(list(mu = mu, sigma = sigma))
  
}

calc_log(med_dw_ch_mild, ui_dw_ch_mild)
calc_log(med_dw_ch_sev, ui_dw_ch_sev)
calc_log(med_dur_subac, ui_dur_subac)
calc_log(med_dw_subac, ui_dw_subac)
calc_log(med_dur_6m, ui_dur_6m)
calc_log(med_dur_12m, ui_dur_12m)
calc_log(med_dur_30m, ui_dur_30m)

# hosp rate per age group
total_n <- overall %>% group_by(age_group)%>%
                         mutate(total_n = sum(n)) %>% 
                         mutate(hosp_rate = n / total_n)

hosp_chikv$total_n <- total_n[, 7]
hosp_chikv <- hosp_chikv %>% group_by(age_group)%>%
  mutate(hosp_n = sum(n)) %>% 
  mutate(hosp_rate = hosp_n / total_n)

hosp_chikv <- hosp_chikv %>% mutate(fatal_rate = ifelse(death_chikv == "TRUE", n/hosp_n, 0))

hosp_fatal_all <- hosp_chikv %>% filter(death_chikv == "TRUE") %>%
                    mutate(var =  ((hosp_rate$total_n*1.1 - hosp_rate$total_n*0.9)/3.92)^2 )

hosp_fatal_all <- hosp_fatal_all %>% mutate(fatal_var = ((upper - lower)/3.92)^2 )

values <-  hosp_fatal_all %>% 
  mutate(hosp_rate = hosp_rate$total_n) %>%  
  select(hosp_rate, var)

## nonhosp chikv fatality
non_hosp_chikv <- non_hosp_chikv %>% group_by(age_group)%>%
  mutate(total_n = sum(n))

non_hosp_chikv <- non_hosp_chikv %>% mutate(fatal_rate = ifelse(death_chikv == "TRUE", n/total_n, 0))

non_hosp_fatal <- non_hosp_chikv %>% filter(death_chikv == "TRUE") %>% mutate(fatal_var = ((upper - lower)/3.92)^2 )

## under 40
library(lhs)
set.seed(123)
runs = 1000
A <- randomLHS (n = runs, 
                k = 27) 

lhs_sample <- matrix (nrow = nrow(A), ncol = ncol(A))

lhs_sample [,1]   <- qbeta (p = A[,1], shape1 = 49.14034, shape2 = 34.14837, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,2]   <- qbeta (p = A[,2], shape1 = 34.21298, shape2 = 31.963, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,3]   <- qbeta (p = A[,3], shape1 = 36.77819, shape2 = 32.7458, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,4]   <- qbeta (p = A[,4], shape1 = 35.84287, shape2 = 32.55955, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,5]   <- qbeta (p = A[,5], shape1 = 539.2823, shape2 = 14152.25, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,6]   <- qbeta (p = A[,6], shape1 = 58.96698, shape2 = 1415.207, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,7]   <- qbeta (p = A[,7], shape1 = 115.4225, shape2 = 111.383, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,8]   <- qlnorm (p = A[,8], meanlog = 2.302585, sdlog = 0.06044082, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,9]   <- qbeta (p = A[,9], shape1 = 4.581639, shape2 = 9.87148, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,10]  <- qlnorm (p = A[,10], meanlog = -0.6301724, sdlog = 0.0852, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,11]  <- qbeta (p = A[,11], shape1 = 22.51835, shape2 = 146.7926, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,12]  <- qlnorm (p = A[,12], meanlog =  -3.734278, sdlog =  0.3877361, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,13]  <- qbeta (p = A[,13], shape1 = 21.45106, shape2 = 399.158, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,14]  <- qlnorm (p = A[,14], meanlog = -4.108138, sdlog =  0.5998406, lower.tail = TRUE, log.p = FALSE)

lhs_sample [,15]   <- qbeta (p = A[,15], shape1 = 393.9252, shape2 = 547.0268, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,16]   <- qbeta (p = A[,16], shape1 = 17875.92, shape2 = 42754.63, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,17]   <- qbeta (p = A[,17], shape1 = 799.2911, shape2 = 3306.976, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,18]   <- qbeta (p = A[,18], shape1 = 77.56885, shape2 = 836.687, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,19]   <- qbeta (p = A[,19], shape1 = 7.605944, shape2 = 1074.944, ncp=0, lower.tail = TRUE, log.p = FALSE)

lhs_sample [,20]  <- qlnorm (p = A[,20], meanlog = -2.145581, sdlog =  0.1815621, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,21]  <- qlnorm (p = A[,21], meanlog = -0.5430045, sdlog =  0.154684, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,22]  <- qlnorm (p = A[,22], meanlog = -2.262311, sdlog =  0.4746817, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,23]  <- qlnorm (p = A[,23], meanlog = -1.148854, sdlog = 0.1815042, lower.tail = TRUE, log.p = FALSE)

lhs_sample [,24]  <- qlnorm (p = A[,24], meanlog =  -0.6931472, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,25]  <- qlnorm (p = A[,25], meanlog = 0, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,26]  <- qlnorm (p = A[,26], meanlog = 0.6931472, sdlog = 0.08380206, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,27]   <- qbeta (p = A[,27], shape1 = 164.6044, shape2 = 618335.4, ncp=0, lower.tail = TRUE, log.p = FALSE)


#check 
quantile(lhs_sample[,1], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,2], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,3], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,4], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,5], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,7], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,8], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,9], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,10], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,11], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,12], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,13], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,14], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,15], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,16], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,17], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,18], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,19], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,27], c(0.025, 0.5, 0.975))

cols <- c("symp_asia", "symp_africa", "symp_america", "symp_overall", "fatal_hosp", "hosp", 
          "lt", "le_lost", "dw_chronic", "dur_chronic", "dw_hosp", "dur_acute", "dw_nonhosp", "dur_nonhosp",
          "acute", "subac", "chr6m", "chr12m", "chr30m", 
          "dw_chronic_mild", "dw_chronic_severe", "dur_subac", "dw_subac", 
          "dur_6m", "dur_12m", "dur_30m", "fatal_nonhosp")

colnames (lhs_sample) <- cols
lhs_sample   <- as.data.frame(lhs_sample)

lhs_sample_young <- lhs_sample

saveRDS(lhs_sample_young, "00_Data/0_2_Processed/lhs_sample_young.rds")


### differential chronic burden for age more than 40 years old

lhs_sample <- matrix (nrow = nrow(A), ncol = ncol(A))

lhs_sample [,1]   <- qbeta (p = A[,1], shape1 = 49.14034, shape2 = 34.14837, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,2]   <- qbeta (p = A[,2], shape1 = 34.21298, shape2 = 31.963, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,3]   <- qbeta (p = A[,3], shape1 = 36.77819, shape2 = 32.7458, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,4]   <- qbeta (p = A[,4], shape1 = 35.84287, shape2 = 32.55955, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,5]   <- qbeta (p = A[,5], shape1 = 539.2823, shape2 = 14152.25, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,6]   <- qbeta (p = A[,6], shape1 = 58.96698, shape2 = 1415.207, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,7]   <- qbeta (p = A[,7], shape1 = 115.4225, shape2 = 111.383, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,8]   <- qlnorm (p = A[,8], meanlog = 2.302585, sdlog = 0.06044082, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,9]   <- qbeta (p = A[,9], shape1 = 4.581639, shape2 = 9.87148, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,10]  <- qlnorm (p = A[,10], meanlog = -0.6301724, sdlog = 0.0852, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,11]  <- qbeta (p = A[,11], shape1 = 22.51835, shape2 = 146.7926, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,12]  <- qlnorm (p = A[,12], meanlog =  -3.734278, sdlog =  0.3877361, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,13]  <- qbeta (p = A[,13], shape1 = 21.45106, shape2 = 399.158, ncp = 0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,14]  <- qlnorm (p = A[,14], meanlog = -4.108138, sdlog =  0.5998406, lower.tail = TRUE, log.p = FALSE)

lhs_sample [,15]   <- qbeta (p = A[,15], shape1 = 388.2874, shape2 = 742.4988, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,16]   <- qbeta (p = A[,16], shape1 = 2487.085, shape2 = 6450.815, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,17]   <- qbeta (p = A[,17], shape1 = 5998.872, shape2 = 21654.86, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,18]   <- qbeta (p = A[,18], shape1 = 184.6674, shape2 = 1216.058, ncp=0, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,19]   <- qbeta (p = A[,19], shape1 = 15.74654, shape2 = 516.3419, ncp=0, lower.tail = TRUE, log.p = FALSE)

lhs_sample [,20]  <- qlnorm (p = A[,20], meanlog = -2.145581, sdlog =  0.1815621, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,21]  <- qlnorm (p = A[,21], meanlog = -0.5430045, sdlog =  0.154684, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,22]  <- qlnorm (p = A[,22], meanlog = -2.262311, sdlog =  0.4746817, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,23]  <- qlnorm (p = A[,23], meanlog = -1.148854, sdlog = 0.1815042, lower.tail = TRUE, log.p = FALSE)

lhs_sample [,24]  <- qlnorm (p = A[,24], meanlog =  -0.6931472, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,25]  <- qlnorm (p = A[,25], meanlog = 0, sdlog =  0.0511915, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,26]  <- qlnorm (p = A[,26], meanlog = 0.6931472, sdlog = 0.08380206, lower.tail = TRUE, log.p = FALSE)
lhs_sample [,27]   <- qbeta (p = A[,27], shape1 = 164.6044, shape2 = 618335.4, ncp=0, lower.tail = TRUE, log.p = FALSE)

#check 
quantile(lhs_sample[,1], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,2], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,3], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,4], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,5], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,7], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,8], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,9], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,10], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,11], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,12], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,13], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,14], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,15], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,16], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,17], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,18], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,19], c(0.025, 0.5, 0.975))
quantile(lhs_sample[,27], c(0.025, 0.5, 0.975))

cols <- c("symp_asia", "symp_africa", "symp_america", "symp_overall", "fatal_hosp", "hosp", 
          "lt", "le_lost", "dw_chronic", "dur_chronic", "dw_hosp", "dur_acute", "dw_nonhosp", "dur_nonhosp",
          "acute", "subac", "chr6m", "chr12m", "chr30m", 
          "dw_chronic_mild", "dw_chronic_severe", "dur_subac", "dw_subac", 
          "dur_6m", "dur_12m", "dur_30m", "fatal_nonhosp")

colnames (lhs_sample) <- cols
lhs_sample <- as.data.frame(lhs_sample)

lhs_old <- lhs_sample

saveRDS(lhs_old, "00_Data/0_2_Processed/lhs_old.rds")


## hosp_rate 
library(lhs)
set.seed(123)
runs = 1000
B <- randomLHS (n = runs, 
                k = 9) 

beta_results <- vector("list", length = nrow(values))

# Loop through each row in the dataframe
for (i in 1:nrow(values)) {
  # Store the result, which is a vector of length 2, in each list element
  beta_results[[i]] <- beta(values[i, 2], values[i, 3])
}

hosp_sample <- matrix (nrow = nrow(B), ncol = ncol(B))

hosp_sample [,1]   <- qbeta (p = B[,1], shape1 = 358.111, shape2 = 4936.908, ncp=0, lower.tail = TRUE, log.p = FALSE)
hosp_sample [,2]   <- qbeta (p = B[,2], shape1 = 374.0897, shape2 = 13933.73, ncp=0, lower.tail = TRUE, log.p = FALSE)
hosp_sample [,3]   <- qbeta (p = B[,3], shape1 = 377.5387, shape2 = 21583.82, ncp=0, lower.tail = TRUE, log.p = FALSE)
hosp_sample [,4]   <- qbeta (p = B[,4], shape1 = 377.765, shape2 = 22374.43, ncp=0, lower.tail = TRUE, log.p = FALSE)
hosp_sample [,5]   <- qbeta (p = B[,5], shape1 = 377.9919, shape2 = 23225.18, ncp=0, lower.tail = TRUE, log.p = FALSE)
hosp_sample [,6]   <- qbeta (p = B[,6], shape1 = 377.6536, shape2 = 21978.21, ncp=0, lower.tail = TRUE, log.p = FALSE)
hosp_sample [,7]   <- qbeta (p = B[,7], shape1 = 376.2637, shape2 = 17976.79, ncp=0, lower.tail = TRUE, log.p = FALSE)
hosp_sample [,8]   <- qbeta (p = B[,8], shape1 = 372.6665, shape2 = 12115.83, ncp=0, lower.tail = TRUE, log.p = FALSE)
hosp_sample [,9]   <- qbeta (p = B[,9], shape1 = 363.4367, shape2 = 6391.353, ncp=0, lower.tail = TRUE, log.p = FALSE)

hosp_cols <- paste0("hosp_", 1:9)
colnames (hosp_sample) <- hosp_cols
hosp_sample   <- as.data.frame(hosp_sample)

## fatality rate
beta_results <- vector("list", length = nrow(fatal_vals))

# Loop through each row in the dataframe
for (i in 1:nrow(fatal_vals)) {
  # Store the result, which is a vector of length 2, in each list element
  beta_results[[i]] <- beta(fatal_vals[i, 2], fatal_vals[i, 3])
}

set.seed(123)
runs = 1000
C <- randomLHS (n = runs, 
                k = 9) 

fatal_sample <- matrix (nrow = nrow(C), ncol = ncol(C))

fatal_sample [,1]   <- qbeta (p = C[,5], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,2]   <- qbeta (p = C[,5], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,3]   <- qbeta (p = C[,5], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,4]   <- qbeta (p = C[,5], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,5]   <- qbeta (p = C[,5], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,6]   <- qbeta (p = C[,5], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,7]   <- qbeta (p = C[,7], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,8]   <- qbeta (p = C[,8], shape1 = 593.7136, shape2 = 238748.2, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,9]   <- qbeta (p = C[,9], shape1 = 367.8489, shape2 = 8318.289, ncp=0, lower.tail = TRUE, log.p = FALSE)

fatal_cols <- paste0("fatal_", 1:9)
colnames (fatal_sample) <- fatal_cols
fatal_sample   <- as.data.frame(fatal_sample)


# age-specific fatality for hospitalised
beta_hosp <- list()
for(i in 1:length(hosp_fatal_all)){
  beta_hosp[[i]] <- beta_calc(hosp_fatal_all$fatal_rate[i], hosp_fatal_all$var[i])
}

# age-specific fatality for non-hospitalised
beta_nonhosp <- list()
for(i in 1:length(hosp_fatal_all)){
  beta_nonhosp[[i]] <- beta_calc(non_hosp_fatal$fatal_rate[i], non_hosp_fatal$fatal_var[i])
}

set.seed(123)
runs = 1000
C <- randomLHS (n = runs, 
                k = 9) 

fatal_sample <- matrix (nrow = nrow(C), ncol = ncol(C))

fatal_sample [,1]   <- qbeta (p = C[,1], shape1 = 26.66166, shape2 = 1455.727, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,2]   <- qbeta (p = C[,2], shape1 = 132.9693, shape2 = 8443.553, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,3]   <- qbeta (p = C[,3], shape1 = 407.9357, shape2 = 22412.47, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,4]   <- qbeta (p = C[,4], shape1 = 493.1692, shape2 = 25471.52, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,5]   <- qbeta (p = C[,5], shape1 = 789.5307, shape2 = 33197.89, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,6]   <- qbeta (p = C[,6], shape1 = 944.9872, shape2 = 34239.3, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,7]   <- qbeta (p = C[,7], shape1 = 1805.207, shape2 = 37879.27, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,8]   <- qbeta (p = C[,8], shape1 = 4466.857, shape2 = 36999.06, ncp=0, lower.tail = TRUE, log.p = FALSE)
fatal_sample [,9]   <- qbeta (p = C[,9], shape1 = 3801.148, shape2 = 16439.15, ncp=0, lower.tail = TRUE, log.p = FALSE)

fatal_cols <- paste0("fatal_", 1:9)
colnames (fatal_sample) <- fatal_cols
fatal_sample   <- as.data.frame(fatal_sample)


## lhs params for non-hosp
nh_fatal_sample <- matrix (nrow = nrow(C), ncol = ncol(C))

nh_fatal_sample [,1]   <- qbeta (p = C[,1], shape1 = 9.693931, shape2 = 30950.3, ncp=0, lower.tail = TRUE, log.p = FALSE)
nh_fatal_sample [,2]   <- qbeta (p = C[,2], shape1 = 6.06979, shape2 = 54677.43, ncp=0, lower.tail = TRUE, log.p = FALSE)
nh_fatal_sample [,3]   <- qbeta (p = C[,3], shape1 = 2.602814, shape2 = 70752.94, ncp=0, lower.tail = TRUE, log.p = FALSE)
nh_fatal_sample [,4]   <- qbeta (p = C[,4], shape1 = 9.693116, shape2 = 93187.19, ncp=0, lower.tail = TRUE, log.p = FALSE)
nh_fatal_sample [,5]   <- qbeta (p = C[,5], shape1 = 3.444061, shape2 = 76516.7, ncp=0, lower.tail = TRUE, log.p = FALSE)
nh_fatal_sample [,6]   <- qbeta (p = C[,6], shape1 = 16.18857, shape2 = 79369.98, ncp=0, lower.tail = TRUE, log.p = FALSE)
nh_fatal_sample [,7]   <- qbeta (p = C[,7], shape1 = 19.94925, shape2 = 54639.27, ncp=0, lower.tail = TRUE, log.p = FALSE)
nh_fatal_sample [,8]   <- qbeta (p = C[,8], shape1 = 27.53439, shape2 = 28387.07, ncp=0, lower.tail = TRUE, log.p = FALSE)
nh_fatal_sample [,9]   <- qbeta (p = C[,9], shape1 = 39.97903, shape2 = 9914.798, ncp=0, lower.tail = TRUE, log.p = FALSE)

nh_fatal_cols <- paste0("nh_fatal_", 1:9)
colnames (nh_fatal_sample) <- nh_fatal_cols
nh_fatal_sample   <- as.data.frame(nh_fatal_sample)


# le lost

med_le_1 <- 70.9
ui_le_1  <- c(med_le_1*0.9, med_le_1*1.1)

med_le_2 <- 61.5
ui_le_2  <- c(med_le_2*0.9, med_le_2*1.1)

med_le_3 <- 52.2
ui_le_3  <- c(med_le_3*0.9, med_le_3*1.1)

med_le_4 <- 42.9
ui_le_4 <- c(med_le_4*0.9, med_le_4*1.1)
  
med_le_5 <- 33.9
ui_le_5  <- c(med_le_5*0.9, med_le_5*1.1)

med_le_6 <- 25.3
ui_le_6  <- c(med_le_6*0.9, med_le_6*1.1)

med_le_7 <- 17.6
ui_le_7  <- c(med_le_7*0.9, med_le_7*1.1)

med_le_8 <- 11
ui_le_8  <- c(med_le_8*0.9, med_le_8*1.1)

med_le_9 <- 6.1
ui_le_9  <- c(med_le_9*0.9, med_le_9*1.1)

calc_log(med_le_1, ui_le_1)
calc_log(med_le_2, ui_le_2)
calc_log(med_le_3, ui_le_3)
calc_log(med_le_4, ui_le_4)
calc_log(med_le_5, ui_le_5)
calc_log(med_le_6, ui_le_6)
calc_log(med_le_7, ui_le_7)
calc_log(med_le_8, ui_le_8)
calc_log(med_le_9, ui_le_9)


set.seed(123)
runs = 1000
D <- randomLHS (n = runs, 
                k = 9) 

le_sample <- matrix (nrow = nrow(D), ncol = ncol(D))
le_sample [,1]   <- qlnorm (p = D[,1], meanlog = 4.26127, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample [,2]   <- qlnorm (p = D[,2], meanlog = 4.119037, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample [,3]   <- qlnorm (p = D[,3], meanlog = 3.955082, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample [,4]   <- qlnorm (p = D[,4], meanlog = 3.758872, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample [,5]   <- qlnorm (p = D[,5], meanlog = 3.523415, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample [,6]   <- qlnorm (p = D[,6], meanlog = 3.230804, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample [,7]   <- qlnorm (p = D[,7], meanlog = 2.867899, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample [,8]   <- qlnorm (p = D[,8], meanlog = 2.397895, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)
le_sample [,9]   <- qlnorm (p = D[,9], meanlog = 1.808289, sdlog = 0.0511915, lower.tail = TRUE, log.p = FALSE)

le_cols <- paste0("le_", 1:9)
colnames (le_sample) <- le_cols
le_sample   <- as.data.frame(le_sample)

saveRDS(le_sample, "00_Data/0_2_Processed/le_sample.rds")
