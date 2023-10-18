# load libraries
library(tidyverse)
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("01_data_adjustment.R")


############################################################################
#### Main analysis: SW duration and FSW age during SW by calendar year #####
############################################################################


### Duration model fit: 
####################### 

year_duration = seq(from = min(data.duration$study_year), to = max(data.duration$study_year), by = 0.25)
year_duration_centred = year_duration - mean(data.duration$study_year)
data.new.duration = data.frame(intercept = rep(1, length(year_duration)), year_centered = year_duration_centred)

fit.duration = stan(file = "../Stanmodels/duration_model.stan",
                    iter = 4000, 
                    data = list(N = nrow(data.duration),
                                L = nrow(data.duration), 
                                y = data.duration$mean_duration_adjusted,
                                ll = 1:nrow(data.duration), 
                                x = data.duration$year_centered, 
                                s = data.duration$studysize_duration, 
                                N_new = nrow(data.new.duration), 
                                x_new = data.new.duration[,2]))


save(data.new.duration, year_duration, fit.duration, 
     file = "../../RData/duration_model_fit.rda")

### Age model fit: 
################### 

year_age = seq(from = min(data.age$study_year), to = max(data.age$study_year), by = 0.25)
year_age_centred = year_age - mean(data.age$study_year)
data.new.age = data.frame(intercept = rep(1, length(year_age)), year_centered = year_age_centred)


fit.age = stan(file = "../Stanmodels/age_model.stan",
                    iter = 4000, 
                    data = list(N = nrow(data.age),
                                L = nrow(data.age), 
                                y = data.age$mean_age_adjusted,
                                ll = 1:nrow(data.age), 
                                x = data.age$year_centered, 
                                s = data.age$studysize_age, 
                                c = 10, 
                                N_sd = sum(!is.na(data.age$SD_age_adjusted)), 
                                sd_y = data.age$SD_age_adjusted[!is.na(data.age$SD_age_adjusted)], 
                                ind_sd = which(!is.na(data.age$SD_age_adjusted)), 
                                N_new = nrow(data.new.age), 
                                x_new = data.new.age[,2]))

save(data.new.age, year_age, fit.age, 
     file = "../../RData/age_model_fit.rda")


############################################################################
#### Sensitivity analysis: Age at entry into SW by year of starting SW #####
############################################################################


data.age.sensitivity = data.age %>% 
  left_join(data.duration %>% select(Population_ID, mean_duration_adjusted)) %>% 
  ungroup() %>% 
  mutate(mean_age_entry = mean_age_adjusted - mean_duration_adjusted,  # constructed age at entry into SW
         year_entry = study_year - mean_duration_adjusted) %>% # constructed year of starting SW
  filter(!is.na(mean_age_entry)) %>% 
  mutate(year_entry_centered = year_entry - mean(year_entry)) 

year_age.sens = seq(from = min(data.age.sensitivity$year_entry), to = max(data.age$study_year), by = 0.25)
year_age_centred.sens = year_age.sens - mean(data.age.sensitivity$year_entry)
data.new.age.sens = data.frame(intercept = rep(1, length(year_age.sens)), year_centered = year_age_centred.sens)

fit.age.sens = stan(file = "../Stanmodels/age_model.stan",
                    iter = 4000,
                    data = list(N = nrow(data.age.sensitivity),
                                L = nrow(data.age.sensitivity), 
                                y = data.age.sensitivity$mean_age_entry,
                                ll = 1:nrow(data.age.sensitivity), 
                                x = data.age.sensitivity$year_entry_centered, 
                                s = data.age.sensitivity$studysize_age, 
                                c = 10, 
                                N_sd = sum(!is.na(data.age.sensitivity$SD_age_adjusted)), 
                                sd_y = data.age.sensitivity$SD_age_adjusted[!is.na(data.age.sensitivity$SD_age_adjusted)], 
                                ind_sd = which(!is.na(data.age.sensitivity$SD_age_adjusted)), 
                                N_new = nrow(data.new.age.sens), 
                                x_new = data.new.age.sens[,2]))
                      
save(data.new.age.sens,data.age.sensitivity, year_age.sens, fit.age.sens, 
     file = "../../RData/age_model_fit_sensitivity.rda")


############################################################################
#### Simulation exercise model fits: Constant FSW age and SW duration  #####
############################################################################

fit.duration.noslope = stan(file = "../Stanmodels/duration_model_noslope.stan",
                            iter = 4000, 
                            data = list(N = nrow(data.duration),
                                        L = nrow(data.duration), 
                                        y = data.duration$mean_duration_adjusted,
                                        ll = 1:nrow(data.duration), 
                                        s = data.duration$studysize_duration, 
                                        N_new = nrow(data.new.duration)))

fit.age.noslope = stan(file = "../Stanmodels/age_model_noslope.stan",
                       iter = 4000, 
                       data = list(N = nrow(data.age),
                                   L = nrow(data.age),
                                   y = data.age$mean_age_adjusted,
                                   ll = 1:nrow(data.age),
                                   s = data.age$studysize_age,
                                   c = 10, 
                                   N_sd = sum(!is.na(data.age$SD_age_adjusted)), 
                                   sd_y = data.age$SD_age_adjusted[!is.na(data.age$SD_age_adjusted)], 
                                   ind_sd = which(!is.na(data.age$SD_age_adjusted)), 
                                   N_new = nrow(data.new.age)))


save(fit.duration.noslope, fit.age.noslope, 
     file = "../../RData/model_fits_noslope.rda")


#######################################################
## Sensitivity analysis: Excluding Milovanovic study ##
#######################################################

# duration model fit

data.duration.exclMilo = data.duration %>% filter(Study != "Milovanovic et al., 2021")
year_duration.exclMilo = seq(from = min(data.duration$study_year), to = max(data.duration$study_year), by = 0.25)
year_duration_centred.exclMilo = year_duration.exclMilo - mean(data.duration.exclMilo$study_year)
data.new.duration.exclMilo = data.frame(intercept = rep(1, length(year_duration.exclMilo)), year_centered = year_duration_centred.exclMilo)

fit.duration.exclMilo = stan(file = "../Stanmodels/duration_model.stan",
                             iter = 4000, 
                             data = list(N = nrow(data.duration.exclMilo),
                                         L = nrow(data.duration.exclMilo), 
                                         y = data.duration.exclMilo$mean_duration_adjusted,
                                         ll = 1:nrow(data.duration.exclMilo), 
                                         x = data.duration.exclMilo$year_centered, 
                                         s = data.duration.exclMilo$studysize_duration, 
                                         N_new = nrow(data.new.duration.exclMilo), 
                                         x_new = data.new.duration.exclMilo[,2]))

                    

save(data.new.duration.exclMilo, data.duration.exclMilo, year_duration.exclMilo, fit.duration.exclMilo, 
     file = "../../RData/duration_model_fit_exclMilo.rda")

# age model fit

data.age.exclMilo = data.age %>% filter(Study != "Milovanovic et al., 2021")
year_age.exclMilo = seq(from = min(data.age$study_year), to = max(data.age$study_year), by = 0.25)
year_age_centred.exclMilo = year_age.exclMilo - mean(data.age.exclMilo$study_year)
data.new.age.exclMilo = data.frame(intercept = rep(1, length(year_age.exclMilo)), year_centered = year_age_centred.exclMilo)

fit.age.exclMilo = stan(file = "../Stanmodels/age_model.stan",
               iter = 4000, 
               data = list(N = nrow(data.age.exclMilo),
                           L = nrow(data.age.exclMilo), 
                           y = data.age.exclMilo$mean_age_adjusted,
                           ll = 1:nrow(data.age.exclMilo), 
                           x = data.age.exclMilo$year_centered, 
                           s = data.age.exclMilo$studysize_age, 
                           c = 10, 
                           N_sd = sum(!is.na(data.age.exclMilo$SD_age_adjusted)), 
                           sd_y = data.age.exclMilo$SD_age_adjusted[!is.na(data.age.exclMilo$SD_age_adjusted)], 
                           ind_sd = which(!is.na(data.age.exclMilo$SD_age_adjusted)), 
                           N_new = nrow(data.new.age.exclMilo), 
                           x_new = data.new.age.exclMilo[,2]))


save(data.new.age.exclMilo, data.age.exclMilo, year_age.exclMilo, fit.age.exclMilo, 
     file = "../../RData/age_model_fit_exclMilo.rda")
