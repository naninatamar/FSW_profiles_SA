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
  pivot_wider(values_from = value, names_from = Prevalence)

## Load estimates of annual age-specific HIV prevalence in females from Thembisa model (can be found on the Thembisa webpage) 
f_prev_age_orig <- read_excel("../../Thembisa_data/f_prev_age.xlsx", sheet = 2) %>% 
  as_tibble() %>% 
  pivot_longer(cols = 2:47, values_to = "prevalence", names_to = "year")

f_prev_by_age_year = f_prev_age_orig %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(Age = ifelse(Age == "90+", "90", Age)) %>% 
  mutate(Age = as.numeric(Age)) %>% 
  mutate(p_agelwr_ylwr = prevalence) %>% 
  group_by(year) %>% 
  arrange(Age, .by_group = TRUE) %>% 
  mutate(p_ageupr_ylwr = lead(prevalence)) %>% 
  group_by(Age) %>% 
  arrange(year, .by_group = TRUE) %>% 
  mutate(p_agelwr_yupr = lead(prevalence)) 

data.prev.temp = f_prev_age_orig %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(Age = ifelse(Age == "90+", "90", Age)) %>% 
  mutate(Age = as.numeric(Age)) %>% 
  mutate(Age = Age -1, 
         year = year - 1) %>% 
  rename(p_ageupr_yupr = prevalence)

f_prev_by_age_year = f_prev_by_age_year %>% 
  left_join(data.prev.temp) 

## agelong only mean:
sim.fit.age = as.matrix(fit.age) %>% 
  as_tibble()

year.fin = seq(from = 
                 min(data.age$study_year)-0.25, 
               to = max(data.age$study_year)+0.25, by = 0.25)
year_centred.fin = year.fin - mean(data.age$study_year)
data.new.age = data.frame(intercept = rep(1, length(year.fin)), year_centered = year_centred.fin)

age.long_re = sim.fit.age[, c(36,34,37)] %>% 
  rowwise() %>% 
  mutate(eta = rnorm(1, eta_mu , eta_tau)) %>% select(eta, omega) %>% as.matrix() %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>%
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year.fin, 8000)) %>% 
  rename(age_linpred_re = value)

age.long = as.matrix(sim.fit.age[, c(36,34)]) %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>%
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year.fin, 8000)) %>% 
  rename(age_linpred = value)

age.long = age.long %>% left_join(age.long_re)

age.long.sc2 = age.long %>% 
  mutate(Age = floor(age_linpred),
         year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year, 
         diffage = age_linpred - Age)

age.long.sc2 = age.long.sc2 %>% 
  left_join(f_prev_by_age_year)

age.long.sc2 = age.long.sc2 %>% 
  mutate(p_agelwr_y = p_agelwr_ylwr + diffyear*(p_agelwr_yupr - p_agelwr_ylwr)) %>% 
  mutate(p_ageupr_y = p_ageupr_ylwr + diffyear*(p_ageupr_yupr - p_ageupr_ylwr)) %>% 
  mutate(p_age_y = p_agelwr_y + diffage*(p_ageupr_y - p_agelwr_y)) %>% 
  mutate(p_ylwr_age = p_agelwr_ylwr + diffage*(p_ageupr_ylwr - p_agelwr_ylwr)) %>% 
  mutate(p_yupr_age = p_agelwr_yupr + diffage*(p_ageupr_yupr - p_agelwr_yupr)) %>% 
  mutate(p_age_y2 = p_ylwr_age + diffyear*(p_yupr_age - p_ylwr_age)) ## gives exactly the same! 


age.long.sc2 = age.long.sc2 %>% 
  select(imp, study_year, age_linpred, age_linpred_re, p0_linpred = p_age_y)

age.long.sc2 = age.long.sc2 %>% 
  mutate(study_year = as.numeric(study_year)) %>% 
  mutate(Age = floor(age_linpred_re),
         year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year, 
         diffage = age_linpred_re - Age)

age.long.sc2 = age.long.sc2 %>% 
  left_join(f_prev_by_age_year)

age.long.sc2 = age.long.sc2 %>% 
  mutate(p_agelwr_y = p_agelwr_ylwr + diffyear*(p_agelwr_yupr - p_agelwr_ylwr)) %>% 
  mutate(p_ageupr_y = p_ageupr_ylwr + diffyear*(p_ageupr_yupr - p_ageupr_ylwr)) %>% 
  mutate(p_age_y = p_agelwr_y + diffage*(p_ageupr_y - p_agelwr_y))

p0.long.sc2 = age.long.sc2 %>% 
  select(imp, study_year, age_linpred, age_linpred_re, p0_linpred, p0_pred = p_age_y)

# duration of sW
year_dur = seq(from = min(data.duration$study_year)-0.25, 
               to = max(data.duration$study_year)+0.25, by = 0.25)
year_centred_dur = year_dur - mean(data.duration$study_year)
data.new.dur = data.frame(intercept = rep(1, length(year_dur)), year_centered = year_centred_dur)

sim.dur = as.matrix(fit.duration) %>% 
  as_tibble()

mu.long = as.matrix(sim.dur[, c(1,3)]) %*% 
  t(as.matrix(data.new.dur)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year.fin, 8000)) %>% 
  rename(mu_linpred = value)


mu.long.re = sim.dur[, c(1,2,3)] %>% 
  rowwise() %>% 
  mutate(theta = rnorm(1, theta_mu , theta_tau)) %>% select(theta, gamma) %>% as.matrix() %*% 
  t(as.matrix(data.new.dur)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year.fin, 8000)) %>% 
  rename(mu_linpred_re = value)


mu.long = mu.long %>% 
  left_join(mu.long.re)

mu.long = mu.long %>% 
  mutate(study_year = as.numeric(study_year))

data.new = data.frame(study_year = year.fin)

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

p.long = data.new %>% 
  mutate(year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year) %>% 
  left_join(fsw.prev.temp) %>% 
  mutate(fsw_prev_linexp = fswprev + diffyear*diff_fswprev) %>% 
  mutate(dpoverdt_linexp = dpoverdt + diffyear*diff_dpoverdt) %>% 
  select(study_year, p_base = fswprev, p_linexp = fsw_prev_linexp, dpoverdt_linexp) 



data.sc3 = p0.long.sc2 %>% left_join(mu.long) %>% left_join(p.long)


data.sc3 = data.sc3 %>% 
  mutate(lambda_linpred_re = 1/mu_linpred_re*(p_linexp - p0_pred)/(1-p_linexp), 
         lambda_linpred = 1/mu_linpred*(p_linexp - p0_linpred)/(1-p_linexp))

lambda.sc3 = data.sc3 %>% 
  group_by(study_year) %>% 
  summarise(lambda_med = median(lambda_linpred), 
            lambda_lwr = quantile(lambda_linpred, 0.05), 
            lambda_upr = quantile(lambda_linpred, 0.95), 
            lambda_lwr_re = quantile(lambda_linpred_re, 0.05), 
            lambda_upr_re = quantile(lambda_linpred_re, 0.95))

lambda.sc3_50 = data.sc3 %>% 
  group_by(study_year) %>% 
  summarise(lambda_med = median(lambda_linpred), 
            lambda_lwr = quantile(lambda_linpred, 0.25), 
            lambda_upr = quantile(lambda_linpred, 0.75), 
            lambda_lwr_re = quantile(lambda_linpred_re, 0.25), 
            lambda_upr_re = quantile(lambda_linpred_re, 0.75))

data.sc4 = data.sc3 %>% 
  mutate(p0_linpred_a = 1.25*p0_linpred, 
         p0_linpred_b = 1.5*p0_linpred, 
         p0_linpred_c = 1.75*p0_linpred, 
         p0_linpred_d = 2*p0_linpred) %>% 
  mutate(lambda_linpred_a = 1/mu_linpred*(p_linexp - p0_linpred_a)/(1-p_linexp),
         lambda_linpred_b = 1/mu_linpred*(p_linexp - p0_linpred_b)/(1-p_linexp),
         lambda_linpred_c = 1/mu_linpred*(p_linexp - p0_linpred_c)/(1-p_linexp),
         lambda_linpred_d = 1/mu_linpred*(p_linexp - p0_linpred_d)/(1-p_linexp))


lambda.sc4 = data.sc4 %>% 
  group_by(study_year) %>% 
  summarise(lambda_med_0 = median(lambda_linpred), 
            lambda_lwr_0 = quantile(lambda_linpred, 0.05), 
            lambda_upr_0 = quantile(lambda_linpred, 0.95),
            lambda_med_a = median(lambda_linpred_a), 
            lambda_lwr_a = quantile(lambda_linpred_a, 0.05), 
            lambda_upr_a = quantile(lambda_linpred_a, 0.95), 
            lambda_med_b = median(lambda_linpred_b), 
            lambda_lwr_b = quantile(lambda_linpred_b, 0.05), 
            lambda_upr_b = quantile(lambda_linpred_b, 0.95), 
            lambda_med_c = median(lambda_linpred_c), 
            lambda_lwr_c = quantile(lambda_linpred_c, 0.05), 
            lambda_upr_c = quantile(lambda_linpred_c, 0.95), 
            lambda_med_d = median(lambda_linpred_d), 
            lambda_lwr_d = quantile(lambda_linpred_d, 0.05), 
            lambda_upr_d = quantile(lambda_linpred_d, 0.95)) %>% 
  pivot_longer(cols = contains("lambda")) %>% 
  mutate(estimate = str_sub(name, 8, 10)) %>% 
  mutate(f_p0 = str_sub(name, 12,12)) %>% 
  select(-name) %>% 
  pivot_wider(values_from = value, names_from = estimate) %>% 
  mutate(factor_increase = case_when(f_p0 == "0" ~ 1, 
                                     f_p0 == "a" ~ 1.25, 
                                     f_p0 == "b" ~ 1.5, 
                                     f_p0 == "c" ~ 1.75,
                                     f_p0 == "d" ~ 2))


####
#### the same with the models assuming a constant SW duration and FSW age
####

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
  bind_cols(study_year = rep(year.fin, 8000)) %>% 
  rename(age_linpred_re = value)


age.long_noslope = as.matrix(sim.fit.age.noslope[, c(35)]) %*% 
  t(as.matrix(data.new.age[,1])) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+10}) %>%
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year.fin, 8000)) %>% 
  rename(age_linpred = value)


age.long_noslope = age.long_noslope %>% left_join(age.long_re_noslope)

age.long.sc2_noslope = age.long_noslope %>% 
  mutate(Age = floor(age_linpred),
         year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year, 
         diffage = age_linpred - Age)

age.long.sc2_noslope = age.long.sc2_noslope %>% 
  left_join(f_prev_by_age_year)

age.long.sc2_noslope = age.long.sc2_noslope %>% 
  mutate(p_agelwr_y = p_agelwr_ylwr + diffyear*(p_agelwr_yupr - p_agelwr_ylwr)) %>% 
  mutate(p_ageupr_y = p_ageupr_ylwr + diffyear*(p_ageupr_yupr - p_ageupr_ylwr)) %>% 
  mutate(p_age_y = p_agelwr_y + diffage*(p_ageupr_y - p_agelwr_y))


age.long.sc2_noslope = age.long.sc2_noslope %>% 
  select(imp, study_year, age_linpred, age_linpred_re, p0_linpred = p_age_y)

age.long.sc2_noslope = age.long.sc2_noslope %>% 
  mutate(study_year = as.numeric(study_year)) %>% 
  mutate(Age = floor(age_linpred_re),
         year = floor(study_year)) %>% 
  mutate(diffyear = study_year - year, 
         diffage = age_linpred_re - Age)

age.long.sc2_noslope = age.long.sc2_noslope %>% 
  left_join(f_prev_by_age_year)

age.long.sc2_noslope = age.long.sc2_noslope %>% 
  mutate(p_agelwr_y = p_agelwr_ylwr + diffyear*(p_agelwr_yupr - p_agelwr_ylwr)) %>% 
  mutate(p_ageupr_y = p_ageupr_ylwr + diffyear*(p_ageupr_yupr - p_ageupr_ylwr)) %>% 
  mutate(p_age_y = p_agelwr_y + diffage*(p_ageupr_y - p_agelwr_y))

p0.long.sc2_noslope = age.long.sc2_noslope %>% 
  select(imp, study_year, age_linpred, age_linpred_re, p0_linpred, p0_pred = p_age_y)

sim.dur_noslope = as.matrix(fit.duration.noslope) %>% 
  as_tibble()

mu.long_noslope = as.matrix(sim.dur_noslope[, c(1)]) %*% 
  t(as.matrix(data.new.dur[,1])) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year.fin, 8000)) %>% 
  rename(mu_linpred = value)

mu.long.re_noslope = sim.dur_noslope[, c(1,2)] %>% 
  rowwise() %>% 
  mutate(theta = rnorm(1, theta_mu , theta_tau)) %>% select(theta) %>% as.matrix() %*% 
  t(as.matrix(data.new.dur[,1])) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate(imp  = row_number()) %>% 
  pivot_longer(cols = 1:(93+2)) %>% 
  bind_cols(study_year = rep(year.fin, 8000)) %>% 
  rename(mu_linpred_re = value)


mu.long_noslope = mu.long_noslope %>% 
  left_join(mu.long.re_noslope)

mu.long_noslope = mu.long_noslope %>% 
  mutate(study_year = as.numeric(study_year))

data.sc2_noslope = mu.long_noslope %>% left_join(p0.long.sc2_noslope) %>% select(-name)
data.sc2_noslope = data.sc2_noslope %>% 
  mutate(p1 = 0.2, 
         p2 = 0.25, 
         p3 = 0.3,
         p4 = 0.35, 
         p5 = 0.4, 
         p6 = 0.45, 
         p7 = 0.5) %>% 
  pivot_longer(cols = c(p1,p2,p3,p4,p5,p6,p7), 
               values_to = "p_gen")

data.sc2_noslope = data.sc2_noslope %>% 
  mutate(lambda_pred = 1/mu_linpred_re*(p_gen-p0_pred)/(1-p_gen), 
         lambda_linpred = 1/mu_linpred*(p_gen-p0_pred)/(1-p_gen))

lambda.sc2_noslope = data.sc2_noslope %>% 
  group_by(p_gen, study_year) %>% 
  summarise(tot = n(), 
            lambda_med = median(lambda_linpred), 
            lambda_lwr = quantile(lambda_linpred, 0.05),
            lambda_upr = quantile(lambda_linpred, 0.95), 
            lambda_lwr_pred = quantile(lambda_pred, 0.05), 
            lambda_upr_pred = quantile(lambda_pred, 0.95))

data.sc3_noslope = p0.long.sc2_noslope %>% left_join(mu.long_noslope) %>% left_join(p.long)


data.sc3_noslope = data.sc3_noslope %>% 
  mutate(lambda_linpred_re = 1/mu_linpred_re*(p_linexp - p0_pred)/(1-p_linexp), 
         lambda_linpred = 1/mu_linpred*(p_linexp - p0_linpred)/(1-p_linexp))

lambda.sc3_noslope = data.sc3_noslope %>% 
  group_by(study_year) %>% 
  summarise(lambda_med = median(lambda_linpred), 
            lambda_lwr = quantile(lambda_linpred, 0.05), 
            lambda_upr = quantile(lambda_linpred, 0.95), 
            lambda_lwr_re = quantile(lambda_linpred_re, 0.05), 
            lambda_upr_re = quantile(lambda_linpred_re, 0.95))


lambda.sc3_noslope_50 = data.sc3_noslope %>% 
  group_by(study_year) %>% 
  summarise(lambda_med = median(lambda_linpred), 
            lambda_lwr = quantile(lambda_linpred, 0.25), 
            lambda_upr = quantile(lambda_linpred, 0.75), 
            lambda_lwr_re = quantile(lambda_linpred_re, 0.25), 
            lambda_upr_re = quantile(lambda_linpred_re, 0.75))


data.sc4_noslope = data.sc3_noslope %>% 
  mutate(p0_linpred_a = 1.25*p0_linpred, 
         p0_linpred_b = 1.5*p0_linpred, 
         p0_linpred_c = 1.75*p0_linpred, 
         p0_linpred_d = 2*p0_linpred) %>% 
  mutate(lambda_linpred_a = 1/mu_linpred*(p_linexp - p0_linpred_a)/(1-p_linexp),
         lambda_linpred_b = 1/mu_linpred*(p_linexp - p0_linpred_b)/(1-p_linexp),
         lambda_linpred_c = 1/mu_linpred*(p_linexp - p0_linpred_c)/(1-p_linexp),
         lambda_linpred_d = 1/mu_linpred*(p_linexp - p0_linpred_d)/(1-p_linexp))


lambda.sc4_noslope = data.sc4_noslope %>% 
  group_by(study_year) %>% 
  summarise(lambda_med_0 = median(lambda_linpred), 
            lambda_lwr_0 = quantile(lambda_linpred, 0.05), 
            lambda_upr_0 = quantile(lambda_linpred, 0.95),
            lambda_med_a = median(lambda_linpred_a), 
            lambda_lwr_a = quantile(lambda_linpred_a, 0.05), 
            lambda_upr_a = quantile(lambda_linpred_a, 0.95), 
            lambda_med_b = median(lambda_linpred_b), 
            lambda_lwr_b = quantile(lambda_linpred_b, 0.05), 
            lambda_upr_b = quantile(lambda_linpred_b, 0.95), 
            lambda_med_c = median(lambda_linpred_c), 
            lambda_lwr_c = quantile(lambda_linpred_c, 0.05), 
            lambda_upr_c = quantile(lambda_linpred_c, 0.95), 
            lambda_med_d = median(lambda_linpred_d), 
            lambda_lwr_d = quantile(lambda_linpred_d, 0.05), 
            lambda_upr_d = quantile(lambda_linpred_d, 0.95)) %>% 
  pivot_longer(cols = contains("lambda")) %>% 
  mutate(estimate = str_sub(name, 8, 10)) %>% 
  mutate(f_p0 = str_sub(name, 12,12)) %>% 
  select(-name) %>% 
  pivot_wider(values_from = value, names_from = estimate) %>% 
  mutate(factor_increase = case_when(f_p0 == "0" ~ 1, 
                                     f_p0 == "a" ~ 1.25, 
                                     f_p0 == "b" ~ 1.5, 
                                     f_p0 == "c" ~ 1.75,
                                     f_p0 == "d" ~ 2))



lambda.sc3.tot = lambda.sc3 %>% mutate(type = "time trend") %>% 
  bind_rows(lambda.sc3_noslope %>% mutate(type = "constant"))

lambda.sc4.tot = lambda.sc4 %>% mutate(type = "time trend") %>% 
  bind_rows(lambda.sc4_noslope %>% mutate(type = "constant"))            


fakedata = data.frame(y = c(0.5, 0.5), ymin = c(0.45, 0.45), 
                      ymax = c(0.55, 0.55), x = c(2000, 2010), fact = "d = 1.25",
                      type = "time trend", factor_increase = 1)

## ad the dP/dt term (was not included before)
data.sc3.new = data.sc3 %>% 
  mutate(lambda_linpred_re = (dpoverdt_linexp + 1/mu_linpred_re*(p_linexp - p0_pred))/(1-p_linexp), 
         lambda_linpred = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred))/(1-p_linexp))


lambda.sc3.new = data.sc3.new %>% 
  filter(study_year >= 1996, study_year <= 2019) %>% 
  group_by(study_year) %>% 
  summarise(lambda_med = median(lambda_linpred), 
            lambda_lwr = quantile(lambda_linpred, 0.05), 
            lambda_upr = quantile(lambda_linpred, 0.95), 
            lambda_lwr_re = quantile(lambda_linpred_re, 0.05), 
            lambda_upr_re = quantile(lambda_linpred_re, 0.95))


data.sc3_noslope.new = data.sc3_noslope %>% 
  mutate(lambda_linpred_re = (dpoverdt_linexp + 1/mu_linpred_re*(p_linexp - p0_pred))/(1-p_linexp), 
         lambda_linpred = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred))/(1-p_linexp))


lambda.sc3_noslope.new = data.sc3_noslope.new %>% 
  filter(study_year >= 1996, study_year <= 2019) %>% 
  group_by(study_year) %>% 
  summarise(lambda_med = median(lambda_linpred), 
            lambda_lwr = quantile(lambda_linpred, 0.05), 
            lambda_upr = quantile(lambda_linpred, 0.95), 
            lambda_lwr_re = quantile(lambda_linpred_re, 0.05), 
            lambda_upr_re = quantile(lambda_linpred_re, 0.95))


(plot_sim_d1.new = lambda.sc3.new  %>% 
    mutate(type = "time trend") %>% 
    mutate(fac = "r = 1") %>% 
    ggplot(aes(x=study_year, y=lambda_med, group = type)) +
    geom_line(data = lambda.sc3_noslope.new %>% mutate(type = "constant") %>%  mutate(fac = "r = 1"),
              aes(linetype = type, col = type), size = 0.6) + 
    geom_ribbon(data = lambda.sc3_noslope.new %>% mutate(type = "constant") %>%  mutate(fac = "r = 1"),
                aes(ymin = lambda_lwr, ymax = lambda_upr, col = type, fill = type), 
                alpha = 0.2, linetype = 2) + 
    geom_line(aes(col = type, linetype = type), size = 0.6) + 
    geom_ribbon(aes(ymin = lambda_lwr, ymax = lambda_upr, fill = type, col = type, linetype = type), 
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


data.sc4.new = data.sc3 %>% 
  mutate(p0_linpred_a = 1.25*p0_linpred, 
         p0_linpred_b = 1.5*p0_linpred, 
         p0_linpred_c = 1.75*p0_linpred, 
         p0_linpred_d = 2*p0_linpred) %>% 
  mutate(lambda_linpred_a = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_a))/(1-p_linexp), 
         lambda_linpred_b = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_b))/(1-p_linexp),
         lambda_linpred_c = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_c))/(1-p_linexp),
         lambda_linpred_d = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_d)/(1-p_linexp)))


lambda.sc4.new = data.sc4.new %>% 
  group_by(study_year) %>% 
  summarise(lambda_med_0 = median(lambda_linpred), 
            lambda_lwr_0 = quantile(lambda_linpred, 0.05), 
            lambda_upr_0 = quantile(lambda_linpred, 0.95),
            lambda_med_a = median(lambda_linpred_a), 
            lambda_lwr_a = quantile(lambda_linpred_a, 0.05), 
            lambda_upr_a = quantile(lambda_linpred_a, 0.95), 
            lambda_med_b = median(lambda_linpred_b), 
            lambda_lwr_b = quantile(lambda_linpred_b, 0.05), 
            lambda_upr_b = quantile(lambda_linpred_b, 0.95), 
            lambda_med_c = median(lambda_linpred_c), 
            lambda_lwr_c = quantile(lambda_linpred_c, 0.05), 
            lambda_upr_c = quantile(lambda_linpred_c, 0.95), 
            lambda_med_d = median(lambda_linpred_d), 
            lambda_lwr_d = quantile(lambda_linpred_d, 0.05), 
            lambda_upr_d = quantile(lambda_linpred_d, 0.95)) %>% 
  pivot_longer(cols = contains("lambda")) %>% 
  mutate(estimate = str_sub(name, 8, 10)) %>% 
  mutate(f_p0 = str_sub(name, 12,12)) %>% 
  select(-name) %>% 
  pivot_wider(values_from = value, names_from = estimate) %>% 
  mutate(factor_increase = case_when(f_p0 == "0" ~ 1, 
                                     f_p0 == "a" ~ 1.25, 
                                     f_p0 == "b" ~ 1.5, 
                                     f_p0 == "c" ~ 1.75,
                                     f_p0 == "d" ~ 2))

data.sc4_noslope.new = data.sc3_noslope %>% 
  mutate(p0_linpred_a = 1.25*p0_linpred, 
         p0_linpred_b = 1.5*p0_linpred, 
         p0_linpred_c = 1.75*p0_linpred, 
         p0_linpred_d = 2*p0_linpred) %>% 
  mutate(lambda_linpred_a = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_a))/(1-p_linexp), 
         lambda_linpred_b = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_b))/(1-p_linexp),
         lambda_linpred_c = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_c))/(1-p_linexp),
         lambda_linpred_d = (dpoverdt_linexp + 1/mu_linpred*(p_linexp - p0_linpred_d)/(1-p_linexp)))


lambda.sc4_noslope.new = data.sc4_noslope.new %>% 
  group_by(study_year) %>% 
  summarise(lambda_med_0 = median(lambda_linpred), 
            lambda_lwr_0 = quantile(lambda_linpred, 0.05), 
            lambda_upr_0 = quantile(lambda_linpred, 0.95),
            lambda_med_a = median(lambda_linpred_a), 
            lambda_lwr_a = quantile(lambda_linpred_a, 0.05), 
            lambda_upr_a = quantile(lambda_linpred_a, 0.95), 
            lambda_med_b = median(lambda_linpred_b), 
            lambda_lwr_b = quantile(lambda_linpred_b, 0.05), 
            lambda_upr_b = quantile(lambda_linpred_b, 0.95), 
            lambda_med_c = median(lambda_linpred_c), 
            lambda_lwr_c = quantile(lambda_linpred_c, 0.05), 
            lambda_upr_c = quantile(lambda_linpred_c, 0.95), 
            lambda_med_d = median(lambda_linpred_d), 
            lambda_lwr_d = quantile(lambda_linpred_d, 0.05), 
            lambda_upr_d = quantile(lambda_linpred_d, 0.95)) %>% 
  pivot_longer(cols = contains("lambda")) %>% 
  mutate(estimate = str_sub(name, 8, 10)) %>% 
  mutate(f_p0 = str_sub(name, 12,12)) %>% 
  select(-name) %>% 
  pivot_wider(values_from = value, names_from = estimate) %>% 
  mutate(factor_increase = case_when(f_p0 == "0" ~ 1, 
                                     f_p0 == "a" ~ 1.25, 
                                     f_p0 == "b" ~ 1.5, 
                                     f_p0 == "c" ~ 1.75,
                                     f_p0 == "d" ~ 2))



fakedata2 = data.frame(y = c(0.5, 0.5), ymin = c(0.55, 0.55), 
                       ymax = c(0.65, 0.65), x = c(2000, 2010), fact = "r = 1.25",
                       type = "time trend", factor_increase = 1)

(plot_sim_d1.25_to_1.75.new = lambda.sc4.new %>% 
    filter(study_year >= 1996, study_year <= 2019) %>% 
    mutate(type = "time trend") %>% 
    mutate(fact = factor(as.character(factor_increase), 
                         labels = c("r = 1", "r = 1.25", "r = 1.5", "r = 1.75", "r = 2"))) %>% 
    filter(!factor_increase %in% c(1,2)) %>% 
    ggplot(aes(x=study_year, y=med, group = type)) +
    geom_line(data = lambda.sc4_noslope.new %>% filter(!factor_increase %in% c(1,2)) %>% 
                filter(study_year >= 1996, study_year <= 2019) %>% 
                mutate(type = "constant", 
                       fact =factor(as.character(factor_increase), 
                                    labels = c("r = 1.25", "r = 1.5", "r = 1.75"))), 
              col = "gray40", lty = 2) + 
    geom_ribbon(data = lambda.sc4_noslope.new %>% filter(!factor_increase %in% c(1,2)) %>% 
                  filter(study_year >= 1996, study_year <= 2019) %>% 
                  mutate(type = "constant", 
                         fact = factor(as.character(factor_increase), 
                                       labels = c("r = 1.25", "r = 1.5", "r = 1.75"))), 
                aes(ymin = lwr, ymax = upr),fill = "gray40",col = "gray40",alpha = 0.2, lty=2) + 
    geom_line(aes(col = factor_increase)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, fill = factor_increase, col = factor_increase),alpha = 0.5) + 
    geom_line(data = fakedata2, aes(y = y, x = x, col = factor_increase)) + 
    geom_ribbon(data = fakedata2, aes(ymin = ymin,ymax = ymax, y=y, x = x, fill = factor_increase), alpha = 0.5) + 
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

cowplot::plot_grid(plot_sim_d1.new, plot_sim_d1.25_to_1.75.new, nrow = 2, 
                   rel_heights = c(0.6, 0.4), labels = LETTERS)

