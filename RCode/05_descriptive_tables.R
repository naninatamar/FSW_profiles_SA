
# load libraries
library(tidyverse)

### 

# load data

source("01_data_adjustment.R") #
## make table:

data.tab = data %>% 
  mutate(Location = case_when(is.na(Location)~ "unspecified", 
                              TRUE ~ Location), 
         Nincluded = case_when(is.na(studysize_duration) | studysize_age == studysize_duration ~ as.character(studysize_age), 
                              studysize_age != studysize_duration ~ paste0(studysize_duration, " / ", studysize_age)), 
         meandur = ifelse(is.na(Mean_duration), "-", sprintf("%.2f", Mean_duration)), 
         quantdur = case_when(is.na(duration_Quantile_1) ~ "-", 
                              is.na(duration_Quantile_2) ~ paste0(round(duration_Quantile_1,1), " (", sprintf("%.2f", duration_Probability_1), ")"), 
                              is.na(duration_Quantile_3) ~ paste0(round(duration_Quantile_1, 1), " (", sprintf("%.2f", duration_Probability_1),  "); ", 
                                                                  round(duration_Quantile_2,1), " (", sprintf("%.2f", duration_Probability_2), ")"), 
                              is.na(duration_Quantile_4) ~  paste0(round(duration_Quantile_1, 1), " (", sprintf("%.2f", duration_Probability_1),  "); ", 
                                                                   round(duration_Quantile_2,1), " (", sprintf("%.2f", duration_Probability_2),"); ", 
                                                                   round(duration_Quantile_3,1), " (", sprintf("%.2f", duration_Probability_3), ")"), 
                              TRUE ~ paste0(round(duration_Quantile_1, 1), " (", sprintf("%.2f", duration_Probability_1),  "); ", 
                                            round(duration_Quantile_2,1), " (", sprintf("%.2f", duration_Probability_2),"); ", 
                                            round(duration_Quantile_3,1), " (", sprintf("%.2f", duration_Probability_3), "); ", 
                                            round(duration_Quantile_4,1), " (", sprintf("%.2f", duration_Probability_4), ")")), 
         mean_sdage = case_when(is.na(Mean_age) ~ "-", 
                                is.na(SD_age) ~ paste0(sprintf("%.2f", Mean_age), " (-)"), 
                                TRUE ~ paste0(sprintf("%.2f", Mean_age), " (", sprintf("%.2f", SD_age), ")")), 
         quantage = case_when(is.na(age_Quantile_1) ~ "-", 
                              is.na(age_Quantile_2) ~ paste0(round(age_Quantile_1,1), " (", sprintf("%.2f", age_Probability_1), ")"), 
                              is.na(age_Quantile_3) ~ paste0(round(age_Quantile_1,1), " (", sprintf("%.2f", age_Probability_1), "); ", 
                                                             round(age_Quantile_2,1), " (", sprintf("%.2f", age_Probability_2),")"), 
                              is.na(age_Quantile_4) ~ paste0(round(age_Quantile_1,1), " (", sprintf("%.2f", age_Probability_1), "); ", 
                                                             round(age_Quantile_2,1), " (", sprintf("%.2f", age_Probability_2),"); ", 
                                                             round(age_Quantile_3,1), " (", sprintf("%.2f", age_Probability_3),")"), 
                              TRUE ~ paste0(round(age_Quantile_1,1), " (", sprintf("%.2f", age_Probability_1), "); ", 
                                            round(age_Quantile_2,1), " (", sprintf("%.2f", age_Probability_2),"); ", 
                                            round(age_Quantile_3,1), " (", sprintf("%.2f", age_Probability_3),"); ", 
                                            round(age_Quantile_4,1), " (", sprintf("%.2f", age_Probability_4),")")), 
         maxage = ifelse(is.na(Max_age), "-", as.character(round(Max_age,1))), 
         exclusion_age = case_when(is.na(exclusion_num) ~ "no", 
                                   TRUE ~ paste0("< ", exclusion_num, " years"))) 

          
years_data = data.tab %>% select(Study, study_year, Province)
impyear = data.tab %>% select(Study, studyyear_calculated)

### 

tab1 = data.tab
dupl = duplicated(tab1$Study)
tab1$Study[dupl] = ""

tab1 = tab1 %>% 
  select(Study, Population_ID, Province, Location, study_year, Nincluded, meandur, quantdur, mean_sdage, quantage, maxage, exclusion_age)

names(tab1) = c("Study", "Population ID", "Province", "Location", "Year", "Study size", 
                    "Mean SW duration", "Quantile (cumulative probability) SW duration", 
                    "Mean (SD) FSW age", "Quantile (cumlative probability) FSW age", 
                    "Max FSW age", "Exclusion")
###

data.tab.duration = data.duration %>% select(Study, Population_ID, study_year, studysize_duration, mean_duration_adjusted)

data.age.raw = data %>% select(Population_ID, Mean_age, SD_age)
data.tab.age = data.age %>% 
  select(Study, Population_ID, study_year, studysize_age, mean_age_adjusted, SD_age_adjusted, adjustment_method) %>% 
  left_join(data.age.raw) %>%
  mutate(truncation_bias = case_when(adjustment_method == "No adjustment (reported mean)" ~ "0 (no exclusion)",
                                     is.na(Mean_age) ~ "-", 
                                     TRUE ~ sprintf("%.2f",Mean_age - mean_age_adjusted)), 
         truncation_bias_sd = case_when(adjustment_method == "No adjustment (reported mean)" & !is.na(SD_age) ~  "0 (no exclusion)",
                                        is.na(SD_age)  ~ "-", 
                                        TRUE ~ sprintf("%.2f", SD_age - SD_age_adjusted))) %>% 
  mutate(truncation_bias_num = case_when(adjustment_method == "No adjustment (reported mean)" | is.na(Mean_age) ~ NA_real_, 
                                         TRUE ~ Mean_age - mean_age_adjusted), 
         truncation_bias_sd_num = case_when(adjustment_method == "No adjustment (reported mean)" |is.na(SD_age) ~  NA_real_, 
                                            TRUE ~ SD_age - SD_age_adjusted)) 
  
  

tab.total = full_join(data.tab.duration, data.tab.age) %>% arrange(Population_ID)

tab.total.print = tab.total %>% 
  select(Study, Population_ID, study_year, studysize_duration, mean_duration_adjusted, 
         studysize_age, mean_age_adjusted, SD_age_adjusted, bias_mean = truncation_bias, bias_sd = truncation_bias_sd) %>% 
  mutate(study_year = as.character(round(study_year, 1)), 
         studysize_duration = ifelse(is.na(studysize_duration), "-", as.character(studysize_duration)), 
         mean_duration_adjusted = ifelse(is.na(mean_duration_adjusted), "-", sprintf("%.2f", mean_duration_adjusted)), 
         studysize_age =  ifelse(is.na(studysize_age), "-", as.character(studysize_age)),
         mean_age_adjusted = ifelse(is.na(mean_age_adjusted), "-", sprintf("%.2f", mean_age_adjusted)),
         SD_age_adjusted = ifelse(is.na(SD_age_adjusted), "-", sprintf("%.2f", SD_age_adjusted)), 
         bias_mean = ifelse(is.na(bias_mean), "-", bias_mean), 
         bias_sd = ifelse(is.na(bias_sd), "-", bias_sd))


names(tab.total.print) = c("Study","Population ID", "Year", "Study size (SW duration)", "Adjusted mean SW duration", "Study size (FSW age)", "Adjusted mean FSW age", "Adjusted Standard devation FSW age", 
                           "Truncation bias mean age", "Truncation bias standard deviation age")

save(tab1, data.tab, years_data, impyear, tab.total, tab.total.print, file = "../../RData/decriptive_tables.rda")

