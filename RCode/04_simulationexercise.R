# load libraries
library(tidyverse)
require(readxl)

# source 
source("01_data_adjustment.R")

# load model outputs 
load("../../RData/age_model_fit.rda")
load("../../RData/duration_model_fit.rda")
load("../../RData/model_fits_noslope.rda")

## Load estimates on annual HIV prevalence in FSW from Thembisa model (can be found on the Thembisa webpage)

fsw_prev = read_excel("../../Thembisa_data/fsw_prev_thembisa.xlsx", sheet = 2) %>% 
  select(Prevalence | starts_with("19") | starts_with("20")) %>% 
  pivot_longer(cols = 2:47) %>% 
  pivot_wider(values_from = value, names_from = Prevalence) # first column contains year, second Thembisa's estimated Mean, 3rd 95% LL and 4th 95% UL (LL and UL not used in this exercise)

## Load estimates of annual age-specific HIV prevalence in females from Thembisa model (can be found on the Thembisa webpage) 
f_prev_age_orig <- read_excel("../../Thembisa_data/f_prev_age.xlsx", sheet = 2) %>% 
  as_tibble() %>% 
  pivot_longer(cols = 2:47, values_to = "prevalence", names_to = "year") # First column contains Age-zear, second calendar year, 3rd mean prevalence (from Thembisa) 

f_prev_by_age_year = f_prev_age_orig %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(Age = ifelse(Age == "90+", "90", Age)) %>% # to make ages numeric transform 90+ into 90
  mutate(Age = as.numeric(Age)) %>% 
  mutate(p_agelwr_ylwr = prevalence) %>% # for linear interpolation: p_agelwr_ylwr contains the prevalence in age a (as reported in first column) and year y (as reported in second column)
  group_by(year) %>% 
  arrange(Age, .by_group = TRUE) %>% 
  mutate(p_ageupr_ylwr = lead(prevalence)) %>%  # for linear interpolation: p_ageupr_ylwr = prevalence in age a+1 and year y
  group_by(Age) %>% 
  arrange(year, .by_group = TRUE) %>% 
  mutate(p_agelwr_yupr = lead(prevalence)) # for linear interpolation: p_agelwr_yupr = prevalence in age a and year y+1

data.prev.temp = f_prev_age_orig %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(Age = ifelse(Age == "90+", "90", Age)) %>% 
  mutate(Age = as.numeric(Age)) %>% 
  mutate(Age = Age -1, 
         year = year - 1) %>% 
  rename(p_ageupr_yupr = prevalence) # for linear interpolation: p_ageupr_yupr = prevalence in age a+1 and year y+1

f_prev_by_age_year = f_prev_by_age_year %>% 
  left_join(data.prev.temp) 

#########################################################################################
#### construct rho, allowing for the estimated time trends in SW duration and FSW age ###
#########################################################################################



### posterior predictions of the expected FSW age from the Age model fit:
#########################################################################

sim.fit.age = as.matrix(fit.age) %>% 
  as_tibble()

cutoff = 10

year = c(min(year_age)-0.25, year_age, max(year_age) + 0.25)
year_centred = year - mean(data.age$study_year)
data.new.age = data.frame(intercept = rep(1, length(year)), year_centered = year_centred)

age.long_re = sim.fit.age[, c(36,34,37)] %>% 
  rowwise() %>% 
  mutate(eta = rnorm(1, eta_mu , eta_tau)) %>% select(eta, omega) %>% as.matrix() %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>%
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year, 8000)) %>% 
  rename(age_linpred_re = value)

age.long = as.matrix(sim.fit.age[, c(36,34)]) %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>%
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year, 8000)) %>% 
  rename(age_linpred = value)

age.long = age.long %>% left_join(age.long_re)

age.long = age.long %>% 
  mutate(Age = floor(age_linpred),
         year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year, 
         diffage = age_linpred - Age) 

age.long = age.long %>% 
  left_join(f_prev_by_age_year) # merge posterior predictions of expected age to age-specific HIV prevalence in general female population

# construct p0 based on the predictions excluding random effects (used in the end)

age.long = age.long %>% # make linear interpolation between integer calendar and age years
  mutate(p_agelwr_y = p_agelwr_ylwr + diffyear*(p_agelwr_yupr - p_agelwr_ylwr)) %>% 
  mutate(p_ageupr_y = p_ageupr_ylwr + diffyear*(p_ageupr_yupr - p_ageupr_ylwr)) %>% 
  mutate(p_age_y = p_agelwr_y + diffage*(p_ageupr_y - p_agelwr_y)) 

age.long = age.long %>% 
  select(imp, study_year, age_linpred, age_linpred_re, p0_linpred = p_age_y)

# construct p0 based on the predictions including random effects (not used in the end)

age.long = age.long %>% 
  mutate(study_year = as.numeric(study_year)) %>% 
  mutate(Age = floor(age_linpred_re),
         year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year, 
         diffage = age_linpred_re - Age)

age.long = age.long %>% 
  left_join(f_prev_by_age_year)

age.long = age.long %>% 
  mutate(p_agelwr_y = p_agelwr_ylwr + diffyear*(p_agelwr_yupr - p_agelwr_ylwr)) %>% 
  mutate(p_ageupr_y = p_ageupr_ylwr + diffyear*(p_ageupr_yupr - p_ageupr_ylwr)) %>% 
  mutate(p_age_y = p_agelwr_y + diffage*(p_ageupr_y - p_agelwr_y))


## P0: Prevalence at entry into SW (under assumption that it is the same as the prevalence in the general female population of the same age)
##################################################################################################################################

p0.long = age.long %>% 
  select(imp, study_year, age_linpred, age_linpred_re, p0_linpred, p0_pred = p_age_y) # p0_llinpred is based on age predictions on pop-mean level (excluding random effects variation), while p0_pred includes random effects variation (not used in the end)


### posterior distribution of the expected duration of sW (= 1/rate of leaving sex work)
######################################################################################################

year_centred_dur = year - mean(data.duration$study_year)
data.new.dur = data.frame(intercept = rep(1, length(year)), year_centered = year_centred_dur)

sim.dur = as.matrix(fit.duration) %>% 
  as_tibble()

mu.long = as.matrix(sim.dur[, c(1,3)]) %*% 
  t(as.matrix(data.new.dur)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year, 8000)) %>% 
  rename(mu_linpred = value)


mu.long.re = sim.dur[, c(1,2,3)] %>% 
  rowwise() %>% 
  mutate(theta = rnorm(1, theta_mu , theta_tau)) %>% select(theta, gamma) %>% as.matrix() %*% 
  t(as.matrix(data.new.dur)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year, 8000)) %>% 
  rename(mu_linpred_re = value)


mu.long = mu.long %>% 
  left_join(mu.long.re)

mu.long = mu.long %>% 
  mutate(study_year = as.numeric(study_year))

### FSW prevalence (interpoated between calendar years)

fsw.prev.temp = fsw_prev %>% 
  mutate(year = as.numeric(name), 
         fswprev = Mean) %>% 
  select(year, fswprev) %>% 
  arrange(year) %>% 
  mutate(fswprev_lead = lead(fswprev)) %>% 
  mutate(diff_fswprev = fswprev_lead - fswprev) %>% 
  mutate(diff_fswprevlag = lag(diff_fswprev)) %>% 
  mutate(dpoverdt = (diff_fswprev + diff_fswprevlag)/2) %>% 
  mutate(dpoverdt_lead = lead(dpoverdt)) %>% 
  mutate(diff_dpoverdt = dpoverdt_lead - dpoverdt)

p.long = data.frame(study_year = year) %>% 
  mutate(year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year) %>% 
  left_join(fsw.prev.temp) %>% 
  mutate(fsw_prev_linexp = fswprev + diffyear*diff_fswprev) %>% 
  mutate(dpoverdt_linexp = dpoverdt + diffyear*diff_dpoverdt) %>% 
  select(study_year, p_base = fswprev, p_linexp = fsw_prev_linexp, dpoverdt_linexp) 

## merge all data needed (P0, P, SW duration/rate of leaving SW) 
data.sim = p0.long %>% left_join(mu.long) %>% left_join(p.long)

# calculated the HIV incidence rate in FSW rho(t):
data.sim = data.sim %>% 
  mutate(rho_linpred_re = (dpoverdt_linexp + 1/mu_linpred_re*(p_linexp - p0_pred))/(1-p_linexp), 
         rho_linpred = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred))/(1-p_linexp))

data.sim.r1 = data.sim %>% # plot assuming r = 1 (HIV prevalence at entry into SW equals (r=1) HIV prevalence in general female population of same age)
  filter(study_year >= 1996, study_year <= 2019) %>% 
  group_by(study_year) %>% 
  summarise(rho_med = median(rho_linpred), 
            rho_lwr = quantile(rho_linpred, 0.05), 
            rho_upr = quantile(rho_linpred, 0.95), 
            rho_lwr_re = quantile(rho_linpred_re, 0.05), 
            rho_upr_re = quantile(rho_linpred_re, 0.95))



data.sim.different.r = data.sim %>% # construct data with r = 1.25, 1.5, 1.75, 2 
  mutate(p0_linpred_a = 1.25*p0_linpred, 
         p0_linpred_b = 1.5*p0_linpred, 
         p0_linpred_c = 1.75*p0_linpred, 
         p0_linpred_d = 2*p0_linpred) %>% 
  mutate(rho_linpred_a = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_a))/(1-p_linexp), 
         rho_linpred_b = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_b))/(1-p_linexp),
         rho_linpred_c = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_c))/(1-p_linexp),
         rho_linpred_d = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_d)/(1-p_linexp)))


data.sim.different.r = data.sim.different.r %>% 
  group_by(study_year) %>% 
  summarise(rho_med_0 = median(rho_linpred), 
            rho_lwr_0 = quantile(rho_linpred, 0.05), 
            rho_upr_0 = quantile(rho_linpred, 0.95),
            rho_med_a = median(rho_linpred_a), 
            rho_lwr_a = quantile(rho_linpred_a, 0.05), 
            rho_upr_a = quantile(rho_linpred_a, 0.95), 
            rho_med_b = median(rho_linpred_b), 
            rho_lwr_b = quantile(rho_linpred_b, 0.05), 
            rho_upr_b = quantile(rho_linpred_b, 0.95), 
            rho_med_c = median(rho_linpred_c), 
            rho_lwr_c = quantile(rho_linpred_c, 0.05), 
            rho_upr_c = quantile(rho_linpred_c, 0.95), 
            rho_med_d = median(rho_linpred_d), 
            rho_lwr_d = quantile(rho_linpred_d, 0.05), 
            rho_upr_d = quantile(rho_linpred_d, 0.95)) %>% 
  pivot_longer(cols = contains("rho")) %>% 
  mutate(estimate = str_sub(name, 5, 7)) %>% 
  mutate(f_p0 = str_sub(name, 9,9)) %>% 
  select(-name) %>% 
  pivot_wider(values_from = value, names_from = estimate) %>% 
  mutate(factor_increase = case_when(f_p0 == "0" ~ 1, 
                                     f_p0 == "a" ~ 1.25, 
                                     f_p0 == "b" ~ 1.5, 
                                     f_p0 == "c" ~ 1.75,
                                     f_p0 == "d" ~ 2))


###############################################################################################
#### construct rho, assuming that SW duration and FSW age remained constant over time  ########
###############################################################################################

sim.fit.age.noslope = as.matrix(fit.age.noslope) %>% 
  as_tibble()

age.long_re_noslope = sim.fit.age.noslope[, c(35,36)] %>% 
  rowwise() %>% 
  mutate(eta = rnorm(1, eta_mu , eta_tau)) %>% select(eta) %>% as.matrix() %*% 
  t(as.matrix(data.new.age[,1])) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+10}) %>%
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year, 8000)) %>% 
  rename(age_linpred_re = value)


age.long_noslope = as.matrix(sim.fit.age.noslope[, c(35)]) %*% 
  t(as.matrix(data.new.age[,1])) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+10}) %>%
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year, 8000)) %>% 
  rename(age_linpred = value)


age.long_noslope = age.long_noslope %>% left_join(age.long_re_noslope)

age.long_noslope = age.long_noslope %>% 
  mutate(Age = floor(age_linpred),
         year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year, 
         diffage = age_linpred - Age)

age.long_noslope = age.long_noslope %>% 
  left_join(f_prev_by_age_year)

age.long_noslope = age.long_noslope %>% 
  mutate(p_agelwr_y = p_agelwr_ylwr + diffyear*(p_agelwr_yupr - p_agelwr_ylwr)) %>% 
  mutate(p_ageupr_y = p_ageupr_ylwr + diffyear*(p_ageupr_yupr - p_ageupr_ylwr)) %>% 
  mutate(p_age_y = p_agelwr_y + diffage*(p_ageupr_y - p_agelwr_y))


age.long_noslope = age.long_noslope %>% 
  select(imp, study_year, age_linpred, age_linpred_re, p0_linpred = p_age_y)

age.long_noslope = age.long_noslope %>% 
  mutate(study_year = as.numeric(study_year)) %>% 
  mutate(Age = floor(age_linpred_re),
         year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year, 
         diffage = age_linpred_re - Age)

age.long_noslope = age.long_noslope %>% 
  left_join(f_prev_by_age_year)

age.long_noslope = age.long_noslope %>% 
  mutate(p_agelwr_y = p_agelwr_ylwr + diffyear*(p_agelwr_yupr - p_agelwr_ylwr)) %>% 
  mutate(p_ageupr_y = p_ageupr_ylwr + diffyear*(p_ageupr_yupr - p_ageupr_ylwr)) %>% 
  mutate(p_age_y = p_agelwr_y + diffage*(p_ageupr_y - p_agelwr_y))

p0.long_noslope = age.long_noslope %>% 
  select(imp, study_year, age_linpred, age_linpred_re, p0_linpred, p0_pred = p_age_y)

sim.dur_noslope = as.matrix(fit.duration.noslope) %>% 
  as_tibble()

mu.long_noslope = as.matrix(sim.dur_noslope[, c(1)]) %*% 
  t(as.matrix(data.new.dur[,1])) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year, 8000)) %>% 
  rename(mu_linpred = value)

mu.long.re_noslope = sim.dur_noslope[, c(1,2)] %>% 
  rowwise() %>% 
  mutate(theta = rnorm(1, theta_mu , theta_tau)) %>% select(theta) %>% as.matrix() %*% 
  t(as.matrix(data.new.dur[,1])) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year, 8000)) %>% 
  rename(mu_linpred_re = value)


mu.long_noslope = mu.long_noslope %>% 
  left_join(mu.long.re_noslope)

mu.long_noslope = mu.long_noslope %>% 
  mutate(study_year = as.numeric(study_year))

## combine all data togeterh (P0, P, and duration of SW/rate of leaving SW)
data.sim_noslope = p0.long_noslope %>% left_join(mu.long_noslope) %>% left_join(p.long)

data.sim_noslope = data.sim_noslope %>% 
  mutate(rho_linpred_re = (dpoverdt_linexp + 1/mu_linpred_re*(p_linexp - p0_pred))/(1-p_linexp), 
         rho_linpred = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred))/(1-p_linexp))


data.sim.r1_noslope = data.sim_noslope %>% 
  filter(study_year >= 1996, study_year <= 2019) %>% 
  group_by(study_year) %>% 
  summarise(rho_med = median(rho_linpred), 
            rho_lwr = quantile(rho_linpred, 0.05), 
            rho_upr = quantile(rho_linpred, 0.95), 
            rho_lwr_re = quantile(rho_linpred_re, 0.05), 
            rho_upr_re = quantile(rho_linpred_re, 0.95))



data.sim.different.r_noslope = data.sim_noslope %>% 
  mutate(p0_linpred_a = 1.25*p0_linpred, 
         p0_linpred_b = 1.5*p0_linpred, 
         p0_linpred_c = 1.75*p0_linpred, 
         p0_linpred_d = 2*p0_linpred) %>% 
  mutate(rho_linpred_a = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_a))/(1-p_linexp), 
         rho_linpred_b = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_b))/(1-p_linexp),
         rho_linpred_c = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_c))/(1-p_linexp),
         rho_linpred_d = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_d)/(1-p_linexp)))


data.sim.different.r_noslope = data.sim.different.r_noslope %>% 
  group_by(study_year) %>% 
  summarise(rho_med_0 = median(rho_linpred), 
            rho_lwr_0 = quantile(rho_linpred, 0.05), 
            rho_upr_0 = quantile(rho_linpred, 0.95),
            rho_med_a = median(rho_linpred_a), 
            rho_lwr_a = quantile(rho_linpred_a, 0.05), 
            rho_upr_a = quantile(rho_linpred_a, 0.95), 
            rho_med_b = median(rho_linpred_b), 
            rho_lwr_b = quantile(rho_linpred_b, 0.05), 
            rho_upr_b = quantile(rho_linpred_b, 0.95), 
            rho_med_c = median(rho_linpred_c), 
            rho_lwr_c = quantile(rho_linpred_c, 0.05), 
            rho_upr_c = quantile(rho_linpred_c, 0.95), 
            rho_med_d = median(rho_linpred_d), 
            rho_lwr_d = quantile(rho_linpred_d, 0.05), 
            rho_upr_d = quantile(rho_linpred_d, 0.95)) %>% 
  pivot_longer(cols = contains("rho")) %>% 
  mutate(estimate = str_sub(name, 5, 7)) %>% 
  mutate(f_p0 = str_sub(name, 9,9)) %>% 
  select(-name) %>% 
  pivot_wider(values_from = value, names_from = estimate) %>% 
  mutate(factor_increase = case_when(f_p0 == "0" ~ 1, 
                                     f_p0 == "a" ~ 1.25, 
                                     f_p0 == "b" ~ 1.5, 
                                     f_p0 == "c" ~ 1.75,
                                     f_p0 == "d" ~ 2))


(plot_sim_r1 = data.sim.r1  %>% 
    mutate(type = "time trend") %>% 
    mutate(fac = "r = 1") %>% 
    ggplot(aes(x=study_year, y=rho_med, group = type)) +
    geom_line(data = data.sim.r1_noslope %>% mutate(type = "constant") %>%  mutate(fac = "r = 1"),
              aes(linetype = type, col = type), size = 0.6) + 
    geom_ribbon(data = data.sim.r1_noslope %>% mutate(type = "constant") %>%  mutate(fac = "r = 1"),
                aes(ymin = rho_lwr, ymax = rho_upr, col = type, fill = type), 
                alpha = 0.2, linetype = 2) + 
    geom_line(aes(col = type, linetype = type), size = 0.6) + 
    geom_ribbon(aes(ymin = rho_lwr, ymax = rho_upr, fill = type, col = type, linetype = type), 
                linetype = 1,  alpha = 0.5) + 
    theme_bw() + 
    labs(y=expression(paste("HIV incidence rate ", rho, " in female sex workers")), 
         x = "Calendar year") + 
    scale_color_manual(values =(c("gray40", "seagreen"))) + 
    scale_fill_manual(values =(c("gray40", "seagreen"))) + 
    guides(fill = FALSE, col = FALSE) + 
    facet_wrap(vars(fac)) + 
    scale_linetype_manual(values = c(1,2), name = NULL, breaks = c("time trend", "constant"), 
                          labels = c("time trends in FSW age & SW duration ", "constant FSW age & SW duration")) + 
    scale_x_continuous(breaks = c(1990, 1996, 2000, 2005, 2010, 2015,2019)) + 
    theme(legend.position = "bottom") + 
    # coord_cartesian(ylim = c(0.01,0.35))) 
    coord_cartesian(ylim = c(0.01,0.45))) 

fakedata = data.frame(y = c(0.5, 0.5), ymin = c(0.55, 0.55), 
                       ymax = c(0.65, 0.65), x = c(2000, 2010), fact = "r = 1.25",
                       type = "time trend", factor_increase = 1)

(plot_sim_r1.25_to_1.75 = data.sim.different.r %>% 
    filter(study_year >= 1996, study_year <= 2019) %>% 
    mutate(type = "time trend") %>% 
    mutate(fact = factor(as.character(factor_increase), 
                         labels = c("r = 1", "r = 1.25", "r = 1.5", "r = 1.75", "r = 2"))) %>% 
    filter(!factor_increase %in% c(1,2)) %>% 
    ggplot(aes(x=study_year, y=med, group = type)) +
    geom_line(data = data.sim.different.r_noslope %>% filter(!factor_increase %in% c(1,2)) %>% 
                filter(study_year >= 1996, study_year <= 2019) %>% 
                mutate(type = "constant", 
                       fact =factor(as.character(factor_increase), 
                                    labels = c("r = 1.25", "r = 1.5", "r = 1.75"))), 
              col = "gray40", lty = 2) + 
    geom_ribbon(data = data.sim.different.r_noslope %>% filter(!factor_increase %in% c(1,2)) %>% 
                  filter(study_year >= 1996, study_year <= 2019) %>% 
                  mutate(type = "constant", 
                         fact = factor(as.character(factor_increase), 
                                       labels = c("r = 1.25", "r = 1.5", "r = 1.75"))), 
                aes(ymin = lwr, ymax = upr),fill = "gray40",col = "gray40",alpha = 0.2, lty=2) + 
    geom_line(aes(col = factor_increase)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = factor_increase, col = factor_increase),alpha = 0.5) + 
    geom_line(data = fakedata, aes(y = y, x = x, col = factor_increase)) + 
    geom_ribbon(data = fakedata, aes(ymin = ymin,ymax = ymax, y=y, x = x, fill = factor_increase), alpha = 0.5) + 
    theme_bw() + 
    labs(y=expression(paste("HIV incidence rate ", rho, " in female sex workers")), 
         x = "Calendar year", 
         fill = "Assumed ratio r of HIV prevalence in FSW at entry\ninto sex work compared to general female population", 
         col =  "Assumed ratio r of HIV prevalence in FSW at entry\ninto sex work compared to general female population") + 
    scale_fill_gradient(high = "chartreuse", low = "seagreen", breaks = c(1.0, 1.25, 1.5, 1.75), labels  = c("1", "1.25", "1.5", "1.75")) + 
    scale_color_gradient(high = "chartreuse", low = "seagreen", breaks = c(1.0, 1.25, 1.5, 1.75),
                         labels  = c("1", "1.25", "1.5", "1.75")) +
    facet_wrap(vars(fact)) + 
    scale_x_continuous(breaks = c(1990, 1996, 2000, 2005, 2010, 2015,2019)) + 
    theme(legend.position = "bottom", 
          legend.margin=margin(0,0,0,0), 
          legend.box.margin=margin(-10,0,0,0), 
          axis.text.x = element_text(angle = 45, hjust = 1)) + 
    coord_cartesian(ylim = c(0.01,0.45)) + theme(legend.position = "top")) 


### Final plot simulation exercise:

cowplot::plot_grid(plot_sim_r1, plot_sim_r1.25_to_1.75, nrow = 2, 
                   rel_heights = c(0.6, 0.4), labels = LETTERS)

