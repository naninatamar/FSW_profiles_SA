# load libraries
library(tidyverse)
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


### check differences in choosing an offset of 10 or of 15 years.

source("01_data_adjustment.R")


################################################
### Data for the age model: with new offset ####
################################################

### change medians and other quantiles to means and correct for potential truncation-induced bias
### (based on distributional assumption)

c_offset = 15 # assumed offset of gamma distribution

data.age15 = data %>% 
  filter(!is.na(Mean_age) | !is.na(age_Quantile_1)) %>%   
  select(-contains("duration")) %>%
  select(- Province, -Location) %>% 
  mutate(exclusion_num_imp = ifelse(is.na(exclusion_num), c_offset, exclusion_num))

## Approach 1: Transformation to mean and truncation-correction based on reported mean and standard devation:  

function_optim = function(mean_age, var_age, trunc_age){
  model <- function(x) {
    F1 <- mean_age - 1/x[2]*expint::gammainc(x[1]+1, trunc_age*x[2])/expint::gammainc(x[1], trunc_age * x[2])
    F2 <- var_age -  1/x[2]^2*expint::gammainc(x[1]+2, trunc_age*x[2])/expint::gammainc(x[1], trunc_age * x[2]) +  (1/x[2]*expint::gammainc(x[1]+1, trunc_age*x[2])/expint::gammainc(x[1], trunc_age * x[2]))^2
    c(F1 = F1, F2 = F2)
  }
  ss = multiroot(f = model, start = c(9, 0.4))
  return(list(root1 = ss$root[1], root2 = ss$root[2]))}

data.age15 = data.age15 %>% 
  rowwise() %>% 
  mutate(root1 = ifelse(is.na(SD_age) | is.na(Mean_age), NA_real_, 
                        function_optim(mean_age = Mean_age-c_offset, var_age = SD_age^2, trunc_age = exclusion_num_imp-c_offset)$root1)) %>% 
  mutate(root2 = ifelse(is.na(SD_age)  | is.na(Mean_age), NA_real_, 
                        function_optim(mean_age = Mean_age-c_offset, var_age = SD_age^2, trunc_age = exclusion_num_imp-c_offset)$root2)) %>% 
  mutate(mean_age_adjusted_ap1 = root1/root2+c_offset, 
         SD_age_adjusted_ap1 = sqrt(root1/root2^2))

## Approach 2: Transformation to mean and truncation-correction based on reported quantiles:  

# define CDF of the truncated gamma distribution

gamma_truncated_cdf = function(x, lb=0, shape, rate){
  prob  = (pgamma(x, shape = shape, rate = rate) - pgamma(lb, shape= shape, rate = rate))/(1-pgamma(lb, shape=shape, rate = rate))
  return(prob)}

# select studies that did report quantiles 

data.age.quant.t115 = data.age15 %>% 
  filter(!is.na(age_Quantile_1)) %>% 
  select(Population_ID, exclusion_num_imp, contains("Quantile")) %>% 
  pivot_longer(cols = c(age_Quantile_1, age_Quantile_2, age_Quantile_3, age_Quantile_4), 
               names_to = "whichquant", values_to = "quantile_age") %>% 
  filter(!is.na(quantile_age)) %>% 
  mutate(whichquant = gsub("age_Quantile_", "", whichquant))


data.age.quant.t215 = data.age15 %>% 
  filter(!is.na(age_Quantile_1)) %>% 
  select(Population_ID, exclusion_num_imp, contains("Probability")) %>%
  select(-contains("rds")) %>% 
  pivot_longer(cols = c(age_Probability_1, age_Probability_2, age_Probability_3, age_Probability_4), 
               names_to = "whichquant", values_to = "probability_age") %>% 
  filter(!is.na(probability_age)) %>% 
  mutate(whichquant = gsub("age_Probability_", "", whichquant))         


data.age.quant15 = data.age.quant.t115 %>% left_join(data.age.quant.t215)

# Loss-function 
function_quantopt = function(quantile, prob, trunc_age){
  model <- function(x){
    shape = x[1]
    rate = x[2]
    prob_tr = gamma_truncated_cdf(x = quantile, lb = trunc_age, shape = x[1], rate = x[2])
    diff = prob_tr - prob
    sum(diff^2)}
  opt = optim(c(9, 0.4), model)
  return(list(shape = opt$par[1], rate = opt$par[2], converged = opt$convergence, value = opt$value))}

data.quant.shape15 = data.age.quant15 %>% group_by(Population_ID ) %>% 
  summarise(root1_quant = function_quantopt(quantile = quantile_age-c_offset, prob = probability_age, trunc_age = exclusion_num_imp-c_offset)$shape)

data.quant.rate15 =  data.age.quant15 %>% group_by(Population_ID ) %>% 
  summarise(root2_quant = function_quantopt(quantile = quantile_age-c_offset, prob = probability_age, trunc_age = exclusion_num_imp-c_offset)$rate)

data.quant.conv15 =  data.age.quant15 %>% group_by(Population_ID) %>% 
  summarise(converged_quant = function_quantopt(quantile = quantile_age-c_offset, prob = probability_age, trunc_age = exclusion_num_imp-c_offset)$converged)

# table(data.quant.conv$converged_quant)

data.quant.value15 =  data.age.quant15 %>% group_by(Population_ID) %>% 
  summarise(value_quant = function_quantopt(quantile = quantile_age-c_offset, prob = probability_age, trunc_age = exclusion_num_imp-c_offset)$value)

data.age.quant15 = data.age.quant15 %>% 
  select(Population_ID) %>% 
  group_by(Population_ID) %>% 
  summarise(tot_quantavail = n()) %>% 
  left_join(data.quant.shape) %>% 
  left_join(data.quant.rate) %>% 
  mutate(mean_age_adjusted_ap2 = root1_quant/root2_quant + c_offset, 
         SD_age_adjusted_ap2 = sqrt(root1_quant/root2_quant^2)) %>% 
  mutate(mean_age_adjusted_ap2 = ifelse(tot_quantavail ==1, NA_real_, mean_age_adjusted_ap2), ## does not work if only 1 quantile available
         SD_age_adjusted_ap2 = ifelse(tot_quantavail ==1, NA_real_, SD_age_adjusted_ap2)) %>% 
  select(-tot_quantavail) %>% 
  filter(!is.na(mean_age_adjusted_ap2))

## Approach 3: Transformation to mean and truncation-correction based on reported mean and maximum:

function_optim_max = function(mean_age, max_age, trunc_age, n_study){
  model <- function(x) {
    F1 <- mean_age - 1/x[2]*expint::gammainc(x[1]+1, trunc_age*x[2])/expint::gammainc(x[1], trunc_age * x[2])
    F2 <- max_age - 1/x[2]*((1+(x[1]-1)/qgamma(1-1/n_study, shape = x[1], rate = x[2]))*(-digamma(1))+qgamma(1-1/n_study, shape = x[1], rate=x[2]))
    c(F1 = F1, F2 = F2)
  }
  ss = multiroot(f = model, start = c(9, 0.4))
  return(list(root1 = ss$root[1], root2 = ss$root[2]))}


function_optim_max_2 = function(mean_age, max_age, trunc_age, n_study){
  model <- function(x) {
    F1 <- mean_age - 1/x[2]*expint::gammainc(x[1]+1, trunc_age*x[2])/expint::gammainc(x[1], trunc_age * x[2])
    F2 <- max_age - ((-digamma(1)*(1/x[2]*(1+(x[1]-1)/(x[2]*qgamma(1-1/n_study, shape = x[1], rate = x[2])))))+qgamma(1-1/n_study, shape = x[1], rate=x[2]))
    c(F1 = F1, F2 = F2)
  }
  ss = multiroot(f = model, start = c(9, 0.4))
  return(list(root1 = ss$root[1], root2 = ss$root[2]))}



data.age.max15 = data.age15 %>% 
  filter(!is.na(Max_age) & !is.na(Mean_age)) %>% 
  select(Population_ID, Mean_age, Max_age, studysize_age, exclusion_num_imp) 

data.age.max15 = data.age.max15 %>% 
  rowwise() %>% 
  mutate(root1_max = function_optim_max_2(mean_age = Mean_age-c_offset, max_age = Max_age-c_offset, trunc_age = exclusion_num_imp-c_offset, n_study = studysize_age)$root1) %>% 
  mutate(root2_max = function_optim_max_2(mean_age = Mean_age-c_offset, max_age = Max_age-c_offset, trunc_age = exclusion_num_imp-c_offset, n_study = studysize_age)$root2) %>% 
  mutate(mean_age_adjusted_ap3 = root1_max/root2_max+c_offset, 
         SD_age_adjusted_ap3 = sqrt(root1_max/root2_max^2)) %>% 
  select(-Mean_age, - Max_age, -exclusion_num_imp)



# bring all together
data.age15 = data.age15 %>% 
  select(Study, Population_ID, exclusion_num_imp, study_year, studysize_age,Mean_age,  SD_age, root1, root2, mean_age_adjusted_ap1, SD_age_adjusted_ap1) %>% 
  left_join(data.age.quant) %>% 
  left_join(data.age.max) %>% 
  mutate(mean_age_adjusted = case_when(exclusion_num_imp == 15 & !is.na(Mean_age) ~ Mean_age, ## take reported mean age if no truncation (no exclusion criteria)
                                       !is.na(Mean_age) & !is.na(SD_age) ~ mean_age_adjusted_ap1, ## take app1 if mean and SD are reported and truncation is apparent
                                       !is.na(mean_age_adjusted_ap2)  ~ mean_age_adjusted_ap2,  ## take app1 if no SD (or no mean) is reported but >=2 quantiles are
                                       TRUE ~ mean_age_adjusted_ap3)) %>%  ## take app3 if only mean but no SD and not quantiles are reported
  mutate(SD_age_adjusted =  case_when(exclusion_num_imp == 15 & !is.na(Mean_age) ~ SD_age, ## take reported mean age if no truncation (no exclusion criteria)
                                      !is.na(Mean_age) & !is.na(SD_age) ~ SD_age_adjusted_ap1, ## take app1 if mean and SD are reported and truncation is apparent
                                      !is.na(SD_age_adjusted_ap2)  ~ SD_age_adjusted_ap2, ## take app1 if no SD (or no mean) is reported but >=2 quantiles are
                                      TRUE ~ SD_age_adjusted_ap3)) %>% ## take app3 if only mean but no SD and not quantiles are reported
  mutate(SD_age_adjusted = case_when(is.na(SD_age) ~ NA_real_, 
                                     TRUE ~ SD_age_adjusted)) %>% 
  mutate(adjustment_method = case_when(exclusion_num_imp == 15 & !is.na(Mean_age) ~ "No adjustment (reported mean)", ## take reported mean age if no truncation (no exclusion criteria)
                                       !is.na(Mean_age) & !is.na(SD_age) ~ "Approach 1",  ## take app1 if mean and SD are reported and truncation is apparent
                                       !is.na(mean_age_adjusted_ap2)  ~ "Approach 2",  ## take app1 if no SD (or no mean) is reported but >=2 quantiles are
                                       TRUE ~ "Approach 3")) %>%  ## take app3 if only mean but no SD and not quantiles are reported
  select(-contains("root"), -contains("_ap"), - Mean_age, - SD_age) %>% 
  ungroup() %>% 
  left_join(data.age.quant_rds) %>% 
  mutate(year_centered = study_year - mean(study_year))   


names(data.age15)
data.age15 = data.age15 %>% select(Population_ID, studysize_age, mean_age_adjusted_15 = mean_age_adjusted, SD_age_adjusted_15 = SD_age_adjusted)

data.age.incl15  = data.age %>% left_join(data.age15)

### comparison of means with offset 15 and offset 10
data.age.incl15 %>% ggplot(aes(x=mean_age_adjusted, y = mean_age_adjusted_15))  + 
  geom_point() + theme_bw() + 
  geom_abline(intercept = 0, slope = 1, col = "green", lty = 2) + 
  labs(x = "Adjusted mean age (offset = 10 years)",y=  "Adjusted mean age (offset = 15 years)")


data.age.incl15 %>% 
  select(study_year, mean_age_adjusted, mean_age_adjusted_15) %>% 
  pivot_longer(cols = c(2:3), names_to = "offset", values_to = "mean_age") %>% 
  mutate(offset = case_when(offset == "mean_age_adjusted" ~ "10", 
                            offset == "mean_age_adjusted_15" ~ "15")) %>% 
  ggplot(aes(x=study_year, y= mean_age, group = offset, col = offset))  + 
  geom_point() + theme_bw() 

## model fit with adjustment = 15 data:

year_age15 = seq(from = min(data.age.incl15$study_year), to = max(data.age.incl15$study_year), by = 0.25)
year_age_centred15 = year_age15 - mean(data.age.incl15$study_year)
data.new.age15 = data.frame(intercept = rep(1, length(year_age15)), year_centered = year_age_centred15)


fit.age15 = stan(file = "../Stanmodels/age_model.stan",
               iter = 4000, 
               data = list(N = nrow(data.age.incl15),
                           L = nrow(data.age.incl15), 
                           y = data.age.incl15$mean_age_adjusted_15,
                           ll = 1:nrow(data.age.incl15), 
                           x = data.age.incl15$year_centered, 
                           s = data.age.incl15$studysize_age, 
                           c = 15, 
                           N_sd = sum(!is.na(data.age.incl15$SD_age_adjusted)), 
                           sd_y = data.age.incl15$SD_age_adjusted_15[!is.na(data.age.incl15$SD_age_adjusted_15)], 
                           ind_sd = which(!is.na(data.age.incl15$SD_age_adjusted_15)), 
                           N_new = nrow(data.new.age15), 
                           x_new = data.new.age15[,2]))
## predictions


## age model

sim.fit.age.15 = as.matrix(fit.age15) %>% 
  as_tibble()

cutoff = 15
## posterior predictions of the expected FSW age - population mean (excluding random effect variability)
pred.age.15 = as.matrix(sim.fit.age.15[, c("eta_mu","omega")]) %*% 
  t(as.matrix(data.new.age15)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_age15)

## posterior predictions of the expected FSW age (including random effect variability)
pred.re.age.15 = sim.fit.age.15[, c("eta_mu","omega","eta_tau")] %>% 
  rowwise() %>% 
  mutate(eta = rnorm(1, eta_mu , eta_tau)) %>% select(eta, omega) %>% as.matrix() %*% 
  t(as.matrix(data.new.age15)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_age15) %>% 
  select(year,Q0025_re = Q0025, Q05_re = Q05, Q0975_re = Q0975)


pred.tot.age.15 = pred.age.15 %>% 
  left_join(pred.re.age.15) 

# posterior predictions of expected FSW age for each individual study: 

random_beta.15 = sim.fit.age.15 %>% select(starts_with("beta["))
common_alpha.15 = sim.fit.age.15 %>% select(alpha)

studylevel_age.15 = apply(random_beta.15, MARGIN = 2, function(x){pull(common_alpha.15)/x + cutoff})

studylevel_age_mean.15 = apply(X =studylevel_age.15, 
                                MARGIN = 2,
                                FUN = mean)

studylevel_age_sd.15 = apply(X = studylevel_age.15, 
                              MARGIN = 2, 
                              FUN = sd)

studylevel_age_quant.15  = apply(X = studylevel_age.15, 
                                  MARGIN = 2, 
                                  FUN = quantile, 
                                  probs = c(0.025, 0.5, 0.975))

studylevel_age_quant.15 = data.frame(t(studylevel_age_quant.15))
names(studylevel_age_quant.15) = c("Q2.5", "Q50", "Q97.5")

age_df.15 = data.frame(studylevel_age_mean.15, studylevel_age_sd.15, studylevel_age_quant.15)
age_df.15$year = data.age.incl15$study_year

age_df.15$studysize_age = data.age.incl15$studysize_age


(plot.age.15 = pred.tot.age.15 %>% 
    ggplot(aes(x=year, y = Q05)) +
    geom_line(col = "tomato", size = 0.8)+
    geom_ribbon(aes(ymin = Q0025, ymax = Q0975),
                alpha = 0.3, fill = "tomato") +
    theme_bw() +
    labs(y = "Mean age (95% CrI) of female sex workers [years]", x = "Study year") +
    geom_point(data = data.age.incl15, 
               aes(y = mean_age_adjusted_15, x = study_year, size = studysize_age),
               col = "black") +
    geom_ribbon(data = pred.tot.age.15, 
                aes(ymin = Q0025_re, ymax = Q0975_re), 
                alpha = 0.2, fill = "tomato") +
    theme(legend.position = "none")  +   
    scale_y_continuous(breaks = c(24,26,28,30, 32, 34, 36, 38)) + 
    scale_x_continuous(breaks = c(1996, 2000,2005, 2010, 2015, 2019)) + 
    geom_pointrange(data = age_df.15, 
                    aes(y = Q50,  ymin = Q2.5, ymax = Q97.5, x = year), 
                    position = position_jitter(width = 0.2, height = 0), 
                    col = "gray40", alpha = 0.6) + 
    scale_size_continuous(range= c(0.1,5))
)  #


load("../../RData/plots_main_analysis.rds")

cowplot::plot_grid(plot.age, plot.age.15, labels = LETTERS)

## combine the two in one plot:

load("../../RData/posterior_predictions.rda")
(plot.comb =  pred.tot.age %>% 
  ggplot(aes(x=year, y = Q05)) +
  geom_line(aes(col = "10"), size = 0.8)+
  geom_ribbon(aes(fill = "10", ymin = Q0025, ymax = Q0975),
              alpha = 0.3) +
  theme_bw() +
  labs(y = "Mean age (95% CrI) of female sex workers [years]", x = "Study year") +
  geom_point(data = data.age.incl15, 
             aes(col = "10", y = mean_age_adjusted, x = study_year, size = studysize_age)) +
  geom_ribbon(data = pred.tot.age, 
              aes(ymin = Q0025_re, ymax = Q0975_re, fill = "10"), 
              alpha = 0.2) + 
  geom_line(data=pred.tot.age.15,aes(col = "15"), size = 0.8, lty = 2)+
  geom_ribbon(data = pred.tot.age.15, aes(ymin = Q0025, ymax = Q0975, fill = "15"),
              alpha = 0.3) +
    geom_point(data = data.age.incl15, 
             aes(y = mean_age_adjusted_15, x = study_year, size = studysize_age, col = "15")) +
  geom_ribbon(data = pred.tot.age.15, 
              aes(ymin = Q0025_re, ymax = Q0975_re, fill =  "15"), 
              alpha = 0.2) + 
  scale_size_continuous(range= c(0.1,5))) + 
  scale_color_manual(values = c("tomato", "seagreen"), name = "offset [years]") + 
  scale_fill_manual(values = c("tomato", "seagreen"), name = "offset [years]") + 
  guides(size = FALSE) + 
  theme(legend.position = c(0.1, 0.9))






