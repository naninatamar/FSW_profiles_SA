
# load libraries
library(tidyverse)
require(readxl)
require(rootSolve)

### 

data = read.csv2("../Data/FSW_summary_data_reported.csv")

####################################
### Data for the duration model ####
####################################

### change medians and other quantiles to means (based on distributional assumption)


data.duration.t1 = data %>% 
  filter(!is.na(Mean_duration) | !is.na(duration_Quantile_1)) %>% 
  select(-contains("age")) %>%
  select(- Province, -Location, - exclusion_num) %>% 
  select(- duration_Quantile_2, -duration_Probability_2) %>% 
  rename(quantile_duration = duration_Quantile_1, probability_duration = duration_Probability_1)


data.duration.t2 = data %>% 
  filter(!is.na(Mean_duration) | !is.na(duration_Quantile_1)) %>% 
  select(-contains("age")) %>%
  select(- Province, -Location, - exclusion_num) %>% 
  select(- duration_Quantile_1, -duration_Probability_1) %>% 
  rename(quantile_duration = duration_Quantile_2, probability_duration = duration_Probability_2) %>% 
  filter(!is.na(quantile_duration))


data.duration = data.duration.t1 %>% 
  bind_rows(data.duration.t2) %>% 
  arrange(Population_ID) %>% 
  rowwise() %>% 
  mutate(lambda_dur = case_when(!is.na(Mean_duration) ~ 1/Mean_duration, 
                                !is.na(quantile_duration)  ~ log(1/(1-probability_duration))/quantile_duration, 
                                TRUE ~ NA_real_)) %>% 
  mutate(reported_summarystatistic = case_when(!is.na(Mean_duration) ~ "mean", 
                                  !is.na(quantile_duration) & probability_duration ==0.5 ~ "median", 
                                  !is.na(quantile_duration) & probability_duration !=0.5 ~ "quantile", 
                                  TRUE ~ NA_character_ )) %>% 
  group_by(Population_ID) %>% 
  mutate(lambda_dur = mean(lambda_dur, na.rm = TRUE)) %>% 
  mutate(mean_duration_adjusted = 1/lambda_dur) %>% 
  select(-quantile_duration, -probability_duration, - lambda_dur, - Mean_duration) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(year_centered = study_year - mean(study_year)) 

###############################
### Data for the age model ####
###############################

### change medians and other quantiles to means and correct for potential truncation-induced bias
### (based on distributional assumption)



c_offset = 10 # assumed offset of gamma distribution

data.age = data %>% 
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

data.age = data.age %>% 
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

data.age.quant.t1 = data.age %>% 
  filter(!is.na(age_Quantile_1)) %>% 
  select(Population_ID, exclusion_num_imp, contains("Quantile")) %>% 
  pivot_longer(cols = c(age_Quantile_1, age_Quantile_2, age_Quantile_3, age_Quantile_4, age_Quantile_5, age_Quantile_6, age_Quantile_7), 
               names_to = "whichquant", values_to = "quantile_age") %>% 
  filter(!is.na(quantile_age)) %>% 
  mutate(whichquant = gsub("age_Quantile_", "", whichquant))


data.age.quant.t2 = data.age %>% 
  filter(!is.na(age_Quantile_1)) %>% 
  select(Population_ID, exclusion_num_imp, contains("Probability")) %>% 
  pivot_longer(cols = c(age_Probability_1, age_Probability_2, age_Probability_3, age_Probability_4, age_Probability_5, age_Probability_6, age_Probability_7), 
               names_to = "whichquant", values_to = "probability_age") %>% 
  filter(!is.na(probability_age)) %>% 
  mutate(whichquant = gsub("age_Probability_", "", whichquant))         
         

data.age.quant = data.age.quant.t1 %>% left_join(data.age.quant.t2)

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

data.quant.shape = data.age.quant %>% group_by(Population_ID ) %>% 
  summarise(root1_quant = function_quantopt(quantile = quantile_age-c_offset, prob = probability_age, trunc_age = exclusion_num_imp-c_offset)$shape)

data.quant.rate =  data.age.quant %>% group_by(Population_ID ) %>% 
  summarise(root2_quant = function_quantopt(quantile = quantile_age-c_offset, prob = probability_age, trunc_age = exclusion_num_imp-c_offset)$rate)

data.quant.conv =  data.age.quant %>% group_by(Population_ID) %>% 
  summarise(converged_quant = function_quantopt(quantile = quantile_age-c_offset, prob = probability_age, trunc_age = exclusion_num_imp-c_offset)$converged)

# table(data.quant.conv$converged_quant)

data.quant.value =  data.age.quant %>% group_by(Population_ID) %>% 
  summarise(value_quant = function_quantopt(quantile = quantile_age-c_offset, prob = probability_age, trunc_age = exclusion_num_imp-c_offset)$value)

data.age.quant = data.age.quant %>% 
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



data.age.max = data.age %>% 
  filter(!is.na(Max_age) & !is.na(Mean_age)) %>% 
  select(Population_ID, Mean_age, Max_age, studysize_age, exclusion_num_imp) 

data.age.max = data.age.max %>% 
  rowwise() %>% 
  mutate(root1_max = function_optim_max_2(mean_age = Mean_age-c_offset, max_age = Max_age-c_offset, trunc_age = exclusion_num_imp-c_offset, n_study = studysize_age)$root1) %>% 
  mutate(root2_max = function_optim_max_2(mean_age = Mean_age-c_offset, max_age = Max_age-c_offset, trunc_age = exclusion_num_imp-c_offset, n_study = studysize_age)$root2) %>% 
  mutate(mean_age_adjusted_ap3 = root1_max/root2_max+c_offset, 
         SD_age_adjusted_ap3 = sqrt(root1_max/root2_max^2)) %>% 
  select(-Mean_age, - Max_age, -exclusion_num_imp)

## data.age

data.age = data.age %>% 
  select(Study, Population_ID, exclusion_num_imp, study_year, studysize_age,Mean_age,  SD_age, root1, root2, mean_age_adjusted_ap1, SD_age_adjusted_ap1) %>% 
  left_join(data.age.quant) %>% 
  left_join(data.age.max) %>% 
  mutate(mean_age_adjusted = case_when(exclusion_num_imp == 10 & !is.na(Mean_age) ~ Mean_age, ## take reported mean age if no truncation (no exclusion criteria)
                                       !is.na(Mean_age) & !is.na(SD_age) ~ mean_age_adjusted_ap1, ## take app1 if mean and SD are reported and truncation is apparent
                                       !is.na(mean_age_adjusted_ap2)  ~ mean_age_adjusted_ap2,  ## take app1 if no SD (or no mean) is reported but >=2 quantiles are
                                       TRUE ~ mean_age_adjusted_ap3)) %>%  ## take app3 if only mean but no SD and not quantiles are reported
  mutate(SD_age_adjusted =  case_when(exclusion_num_imp == 10 & !is.na(Mean_age) ~ SD_age, ## take reported mean age if no truncation (no exclusion criteria)
                                             !is.na(Mean_age) & !is.na(SD_age) ~ SD_age_adjusted_ap1, ## take app1 if mean and SD are reported and truncation is apparent
                                             !is.na(SD_age_adjusted_ap2)  ~ SD_age_adjusted_ap2, ## take app1 if no SD (or no mean) is reported but >=2 quantiles are
                                             TRUE ~ SD_age_adjusted_ap3)) %>% ## take app3 if only mean but no SD and not quantiles are reported
  mutate(SD_age_adjusted = case_when(is.na(SD_age) ~ NA_real_, 
                                     TRUE ~ SD_age_adjusted)) %>% 
  mutate(adjustment_method = case_when(exclusion_num_imp == 10 & !is.na(Mean_age) ~ "No adjustment (reported mean)", ## take reported mean age if no truncation (no exclusion criteria)
                                       !is.na(Mean_age) & !is.na(SD_age) ~ "Approach 1",  ## take app1 if mean and SD are reported and truncation is apparent
                                       !is.na(mean_age_adjusted_ap2)  ~ "Approach 2",  ## take app1 if no SD (or no mean) is reported but >=2 quantiles are
                                       TRUE ~ "Approach 3")) %>%  ## take app3 if only mean but no SD and not quantiles are reported
  select(-contains("root"), -contains("_ap"), - Mean_age, - SD_age) %>% 
  ungroup() %>% 
  mutate(year_centered = study_year - mean(study_year))   

