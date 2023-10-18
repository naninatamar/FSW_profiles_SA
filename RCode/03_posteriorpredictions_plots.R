
# load libraries
library(tidyverse)
require(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# source 

source("01_data_adjustment.R")

# load model outputs 

load("../../RData/duration_model_fit.rda")
load("../../RData/age_model_fit.rda")
load("../../RData/age_model_fit_sensitivity.rda")
load("../../RData/age_model_fit_exclMilo.rda")
load("../../RData/duration_model_fit_exclMilo.rda")

######################
### FSW age:  ########
######################

### posterior predictions 
##########################


sim.fit.age = as.matrix(fit.age) %>% 
  as_tibble()

pred.age = as.matrix(sim.fit.age[, c(36,34)]) %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+10}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year.fin)

cutoff = 10

xpred.fin = data.frame(x1 = rep(1, nrow(data.regr.age.fin)), 
                       x2 = data.regr.age.fin$year_centered)

year.fin = seq(from = min(data.regr.age.fin$studymidyear_prov), to = max(data.regr.age.fin$studymidyear_prov), by = 0.25)
year_centred.fin = year.fin - mean(data.regr.age.fin$studymidyear_prov)
data.new.age = data.frame(intercept = rep(1, length(year.fin)), year_centered = year_centred.fin)

sim.fit.age = as.matrix(fit_untr) %>% 
  as_tibble()

## posterior predictions of the expected FSW age - population mean (excluding random effect variability)
pred.age = as.matrix(sim.fit.age[, c(36,34)]) %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_age)

## posterior predictions of the expected FSW age (including random effect variability)
pred.re.age = sim.fit.age[, c(36,34,37)] %>% 
  rowwise() %>% 
  mutate(mu = rnorm(1, a_mu , a_tau)) %>% select(mu, b) %>% as.matrix() %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_age) %>% 
  select(year,Q0025_re = Q0025, Q05_re = Q05, Q0975_re = Q0975)


pred.tot.age = pred.age %>% 
  left_join(pred.re.age) 

# posterior predictions of expected FSW age for each individual study: 

random_beta = sim.fit.age %>% select(starts_with("beta["))
common_alpha = sim.fit.age %>% select(alpha)

studylevel_age = apply(random_beta, MARGIN = 2, function(x){pull(common_alpha)/x + cutoff})

studylevel_age_mean = apply(X =studylevel_age, 
                            MARGIN = 2,
                            FUN = mean)

studylevel_age_sd = apply(X = studylevel_age, 
                          MARGIN = 2, 
                          FUN = sd)

studylevel_age_quant  = apply(X = studylevel_age, 
                              MARGIN = 2, 
                              FUN = quantile, 
                              probs = c(0.025, 0.5, 0.975))

studylevel_age_quant = data.frame(t(studylevel_age_quant))
names(studylevel_age_quant) = c("Q2.5", "Q50", "Q97.5")

age_df = data.frame(studylevel_age_mean, studylevel_age_sd, studylevel_age_quant)
age_df$year = data.age$study_year

age_df$studysize_age = data.age$studysize_age


### posterior predictions of random intercepts

random_interc = sim.fit.age %>% select(starts_with("a["))
a_mean = apply(X = random_interc, 
               MARGIN = 2, 
               FUN = mean)

a_sd = apply(X = random_interc, 
             MARGIN = 2, 
             FUN = sd)

a_quant  = apply(X = random_interc, 
                 MARGIN = 2, 
                 FUN = quantile, 
                 probs = c(0.025, 0.5, 0.975))

a_quant = data.frame(t(a_quant))
names(a_quant) = c("Q2.5", "Q50", "Q97.5")

a_df = data.frame(a_mean, a_sd, a_quant)

a_df$study = as.character(data.age$Population_ID)
a_df = a_df[order(a_df$a_mean), ]
a_df$a_rank = c(1:dim(a_df)[1])

# posterior distribution of random intercepts (pop-mean):
random_interc_mean = sim.fit.age %>% select(a_mu, a_tau) %>% 
  rowwise() %>% 
  mutate(a_mu_ci = rnorm(1, mean = a_mu, sd = a_tau)) 

a_mu_dat = data.frame(amu_Q05 = rep(median(random_interc_mean$a_mu), nrow(a_df)+2),
                      amu_Q0025 = rep(quantile(random_interc_mean$a_mu, 0.025), nrow(a_df)+2), 
                      amu_Q0975 = rep(quantile(random_interc_mean$a_mu, 0.975), nrow(a_df)+2), 
                      amu_Q0025_re = rep(quantile(random_interc_mean$a_mu_ci, 0.025), nrow(a_df)+2), 
                      amu_Q0975_re = rep(quantile(random_interc_mean$a_mu_ci, 0.975), nrow(a_df)+2), 
                      rank = 0:(nrow(a_df)+1))

### posterior predictions of  parameter estimates of initial Gamma(alpha, beta) distribution

## alpha
alpha_est = sim.fit.age[, "alpha"] %>% 
  summarise(med_alpha = median(alpha), 
            lb_alpha = quantile(alpha, 0.025), 
            ub_alpha = quantile(alpha, 0.975)) # same 

data.alpha = data.new.age %>% 
  mutate(med_alpha = alpha_est$med_alpha, 
         lb_alpha = alpha_est$lb_alpha, 
         ub_alpha = alpha_est$ub_alpha) %>% 
  bind_cols(year = year_age)


## beta

cvar = exp(as.matrix(sim.fit.age[, "log_c_var"]))


pred.beta = as.matrix(sim.fit.age[, c(36,34)]) %*% # pop-level predictions (without random effect variability)
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate_all(., function(x){1/(x*cvar^2)}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_age)


pred.beta.re = sim.fit.age[, c(36,34,37)] %>% # predictions including random effect variability
  rowwise() %>% 
  mutate(mu = rnorm(1, a_mu , a_tau)) %>% 
  select(mu, b) %>% 
  as.matrix() %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){1/(x*cvar^2)}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_age) %>% 
  select(year, re_0975 = Q0975, re_med = Q05, re_0025 = Q0025)


pred.beta.tot = pred.beta %>% 
  left_join(pred.beta.re) 


##  posterior predictions of study-level beta:

random_beta = sim.fit.age %>% select(starts_with("beta["))

random_beta_mean = apply(X = random_beta, 
                           MARGIN = 2,
                           FUN = mean)

random_beta_sd = apply(X = random_beta, 
                         MARGIN = 2, 
                         FUN = sd)

random_beta_quant  = apply(X = random_beta, 
                             MARGIN = 2, 
                             FUN = quantile, 
                             probs = c(0.025, 0.5, 0.975))

random_beta_quant = data.frame(t(random_beta_quant))
names(random_beta_quant) = c("Q2.5", "Q50", "Q97.5")

randombeta_df = data.frame(random_beta_mean, random_beta_sd, random_beta_quant)
randombeta_df$study_year = data.age$study_year

randombeta_df$studysize_age = data.age$studysize_age

dat.alpha = data.alpha %>% 
  select(year = year, 
         Q05=med_alpha, 
         Q0025 = lb_alpha, 
         Q0975 = ub_alpha) %>% 
  mutate(Q0025_re = NA, 
         Q0975_re = NA, 
         type = "alpha")


dat.beta = pred.beta.tot %>% 
  select(year, Q05, Q0025, Q0975, Q0025_re = re_0025, Q0975_re = re_0975) %>% 
  mutate(type = "beta")

dat.gamma.dist = dat.alpha %>% bind_rows(dat.beta)

randombeta_df = randombeta_df %>% 
  mutate(type = "beta") %>% 
  mutate(type = factor(type, levels=c("alpha", "beta"), 
                       labels = c(expression(alpha), expression(beta)))) 

## posterior predictions of standard deviation of FSW age:

pred.sd = as.matrix(sim.fit.age[, c(36,34)]) %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  mutate_all(., function(x){x*cvar}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_age)


pred.sd.re = sim.fit.age[, c(36,34,37)] %>% 
  rowwise() %>% 
  mutate(mu = rnorm(1, a_mu , a_tau)) %>% select(mu, b) %>% as.matrix() %*% 
  t(as.matrix(data.new.age)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x*cvar}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_age) %>% 
  select(year, re_0975 = Q0975, re_med = Q05, re_0025 = Q0025)


pred.sd.tot = pred.sd %>% 
  left_join(pred.sd.re) 

## posterior prediciton of SD FSW age for each study:

random_interc = sim.fit.age %>% select(starts_with("a["))

random_beta = sim.fit.age %>% select(starts_with("beta["))
common_alpha = sim.fit.age %>% select(alpha)
cvar_tot = exp(sim.fit.age %>% select(log_c_var) %>% pull())

studylevel_mu = apply(random_beta, MARGIN = 2, function(x){pull(common_alpha)/x })
studylevel_sd = apply(studylevel_mu, MARGIN = 2, function(x){cvar_tot*x})
studylevel_sd_mean = apply(X =studylevel_sd, 
                           MARGIN = 2,
                           FUN = mean)

studylevel_sd_sd = apply(X = studylevel_sd, 
                         MARGIN = 2, 
                         FUN = sd)

studylevel_sd_quant  = apply(X = studylevel_sd, 
                             MARGIN = 2, 
                             FUN = quantile, 
                             probs = c(0.025, 0.5, 0.975))

studylevel_sd_quant = data.frame(t(studylevel_sd_quant))
names(studylevel_sd_quant) = c("Q2.5", "Q50", "Q97.5")

studylevel_sd_df = data.frame(studylevel_sd_mean, studylevel_sd_sd, studylevel_sd_quant)
studylevel_sd_df$study_year = data.age$study_year


### Plots: Age model fit
#########################


# Plot posterior predictions of expected FSW age
(plot.age = pred.tot.age %>% 
   ggplot(aes(x=year, y = Q05)) +
   geom_line(col = "tomato", size = 0.8)+
   geom_ribbon(aes(ymin = Q0025, ymax = Q0975),
               alpha = 0.3, fill = "tomato") +
   theme_bw() +
   labs(y = "Mean age (95% CrI) of female sex workers [years]", x = "Study year") +
   geom_point(data = data.age, 
              aes(y = mean_age_adjusted, x = study_year, size = studysize_age),
              col = "black") +
   geom_ribbon(data = pred.tot.age, 
               aes(ymin = Q0025_re, ymax = Q0975_re), 
               alpha = 0.2, fill = "tomato") +
   theme(legend.position = "none")  +   
   scale_y_continuous(breaks = c(24,26,28,30, 32, 34, 36, 38)) + 
   scale_x_continuous(breaks = c(1996, 2000,2005, 2010, 2015, 2019)) + 
   geom_pointrange(data = age_df, 
                   aes(y = Q50,  ymin = Q2.5, ymax = Q97.5, x = year), 
                   position = position_jitter(width = 0.2, height = 0), 
                   col = "gray40", alpha = 0.6) + 
   scale_size_continuous(range= c(0.1,5))
)  #


# Plot posterior predictions of alpha, beta parameters of assumed gamma distribution
(plot.prameters.gamma = dat.gamma.dist %>%
    mutate(year = if_else(type == "alpha", 2025+(year - 1995)*0.4, year)) %>%
    mutate(type = factor(type, levels=c("alpha", "beta"), 
                         labels = c(expression(alpha), expression(beta)))) %>%
    ggplot(aes(x=year, y = Q05)) +
    geom_line(col = "seagreen", size = 0.8)+
    geom_ribbon(aes(ymin = Q0025, ymax = Q0975), alpha = 0.3, fill = "seagreen") +
    theme_bw() +
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed) +
    scale_x_continuous(breaks = c(1996, 2000, 2005, 2010, 2015, 2019, 2025, 2030, 2035),
                       labels = c("1996", "2000", "2005", "2010", "2015","2019", "1996", "..........", "2019")) +
    theme(strip.text.x = element_text(size = 12)) +
    labs(y = expression("Estimate of parameters "*alpha*", "*beta*" of "*Gamma(alpha, beta)*" distribution"), x = "Year")+
    geom_ribbon(aes(ymin = Q0025_re, ymax = Q0975_re), alpha = 0.2, fill = "seagreen") +
    theme(legend.position = "none")+ 
    geom_pointrange(data = randombeta_df, 
                    aes(y = Q50,  ymin = Q2.5, ymax = Q97.5, x = study_year), 
                    position = position_jitter(width = 0.2, height = 0), 
                    col = "gray40",
                    alpha = 0.6))

# scale facets so that alpha parameter plot takes less space:
gprandom <- ggplotGrob(plot.prameters.gamma)
facet.columnsrandom <- gprandom$layout$l[grepl("panel", gprandom$layout$name)]
x.varrandom <- sapply(ggplot_build(plot.prameters.gamma)$layout$panel_scales_x,
                      function(l) length(l$range$range))
gprandom$widths[facet.columnsrandom] <- gprandom$widths[facet.columnsrandom] * c(1,3)
grid::grid.draw(gprandom)

#### plot posterior predictions standard deviation:
(plot.sd = pred.sd.tot %>% 
    ggplot(aes(x=year, y = Q05)) +
    geom_line(col = "orange", size = 0.8)+
    geom_ribbon(aes(ymin = Q0025, ymax = Q0975), alpha = 0.3, fill = "orange") +
    theme_bw() +
    labs(y = "Standard deviation (95% CrI) of FSW age [years]", x = "Study year") +
    geom_ribbon(data = pred.sd.tot, 
                aes(ymin = re_0025, ymax = re_0975),
                alpha = 0.2, fill = "orange") +
    # geom_point(data = data.regr.age.fin, 
    #            aes(y = sd_age, x = studymidyear_prov, size = StudyN_age), shape = 0, col = "tomato") +
    geom_point(data = data.age, 
               aes(y = SD_age_adjusted, 
                   x = study_year, 
                   size = studysize_age), 
               col = "black") +
    theme(legend.position = "none")   + 
    scale_x_continuous(breaks = c(1996, 2000,2005, 2010, 2015, 2019)) + 
    geom_pointrange(data = studylevel_sd_df, 
                    aes(y = Q50,  ymin = Q2.5, ymax = Q97.5, 
                        x = study_year), 
                    position = position_jitter(width = 0.2, height = 0), 
                    col = "gray40",
                    alpha = 0.4) + 
    scale_size_continuous(range= c(0.1,5))
  
) 


# Final plot Age model:
cowplot::plot_grid(plot.age,
                   cowplot::plot_grid(plot.sd, gprandom, labels = c("B","C")), 
                   nrow = 2, labels = c("A", NULL))

########################
### Duration of SW: ####
########################

### posterior predictions 
##########################

sim.dur = as.matrix(fit.duration) %>% 
  as_tibble()

pred.rate = as.matrix(sim.dur[, c(1,3)]) %*% 
  t(as.matrix(data.new.dur)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_dur) %>% 
  rename(rate_Q0025 = Q0025, 
         rate_Q05 = Q05, 
         rate_Q0975 = Q0975)


pred.rate.re = sim.dur[, c(1,2,3)] %>% 
  rowwise() %>% 
  mutate(mu = rnorm(1, mu , tau)) %>% select(mu, beta) %>% as.matrix() %*% 
  t(as.matrix(data.new.dur)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_dur) %>% 
  rename(rate_Q0025_re = Q0025, 
         rate_Q05_re = Q05, 
         rate_Q0975_re = Q0975)


pred.rate.tot = pred.rate %>% 
  left_join(pred.rate.re) 

### posterior predictions for each study: 

r_lambda = sim.dur %>% select(starts_with("lambda["))


r_lambda_mean = apply(X =r_lambda, 
                      MARGIN = 2,
                      FUN = mean)

r_lambda_sd = apply(X = r_lambda, 
                    MARGIN = 2, 
                    FUN = sd)

r_lambda_quant  = apply(X = r_lambda, 
                        MARGIN = 2, 
                        FUN = quantile, 
                        probs = c(0.025, 0.5, 0.975))

r_lambda_quant = data.frame(t(r_lambda_quant))
names(r_lambda_quant) = c("Q2.5", "Q50", "Q97.5")

r_lambda_df = data.frame(r_lambda_mean, r_lambda_sd, r_lambda_quant)
r_lambda_df$study_year = data.duration$study_year

r_lambda = sim.dur %>% select(starts_with("lambda["))
r_dur = 1/r_lambda

r_dur_mean = apply(X =r_dur, 
                   MARGIN = 2,
                   FUN = mean)

r_dur_sd = apply(X = r_dur, 
                 MARGIN = 2, 
                 FUN = sd)

r_dur_quant  = apply(X = r_dur, 
                     MARGIN = 2, 
                     FUN = quantile, 
                     probs = c(0.025, 0.5, 0.975))

r_dur_quant = data.frame(t(r_dur_quant))
names(r_dur_quant) = c("Q2.5", "Q50", "Q97.5")

r_dur_df = data.frame(r_dur_mean, r_dur_sd, r_dur_quant)
r_dur_df$study_year = data.duration$study_year

# posterior predictions random intercepts duration model 

random_interc_dur = sim.dur %>% select(starts_with("alpha["))

alpha_mean = apply(X = random_interc_dur, 
                   MARGIN = 2, 
                   FUN = mean)

alpha_sd = apply(X = random_interc_dur, 
                 MARGIN = 2, 
                 FUN = sd)

alpha_quant  = apply(X = random_interc_dur, 
                     MARGIN = 2, 
                     FUN = quantile, 
                     probs = c(0.025, 0.5, 0.975))

alpha_quant = data.frame(t(alpha_quant))
names(alpha_quant) = c("Q2.5", "Q50", "Q97.5")

alpha_df = data.frame(alpha_mean, alpha_sd, alpha_quant)

alpha_df$study = as.character(data.duration$Population_ID)
alpha_df = alpha_df[order(alpha_df$alpha_mean), ]
alpha_df$a_rank = c(1:dim(alpha_df)[1])


random_interc_dur_mean = sim.dur.gamma %>% select(mu, tau) %>% 
  rowwise() %>% 
  mutate(dur_mu_ci = rnorm(1, mean = mu, sd = tau)) 

dur_mu_dat = data.frame(amu_Q05 = rep(median(random_interc_dur_mean$mu), nrow(alpha_df)+2),
                        amu_Q0025 = rep(quantile(random_interc_dur_mean$mu, 0.025), nrow(alpha_df)+2), 
                        amu_Q0975 = rep(quantile(random_interc_dur_mean$mu, 0.975), nrow(alpha_df)+2), 
                        amu_Q0025_re = rep(quantile(random_interc_dur_mean$dur_mu_ci, 0.025), nrow(alpha_df)+2), 
                        amu_Q0975_re = rep(quantile(random_interc_dur_mean$dur_mu_ci, 0.975), nrow(alpha_df)+2), 
                        rank = 0:(nrow(alpha_df)+1))

### Plots: Duration model fit
##############################

# plot posterior prediciton of expected SW duration:

(plot.dur  = pred.rate.tot %>% 
    ggplot(aes(x=year, y = rate_Q05)) + 
    geom_line(col = "dodgerblue", size = 0.8)+ 
    geom_ribbon(aes(ymin = rate_Q0025, ymax = rate_Q0975), alpha = 0.3, fill = "dodgerblue") + 
    theme_bw() + 
    labs(y = "Mean duration (95%-CrI) of sex work [years]", x = "Study year") + 
    scale_x_continuous(breaks = c(1996, 2000, 2005, 2010, 2015, 2019)) +
    geom_ribbon(aes(ymin = rate_Q0025_re, ymax = rate_Q0975_re), alpha = 0.2, fill = "dodgerblue") + 
    geom_point(data = data.duration, aes(y = mean_duration_adjusted, x = study_year, size = studysize_duration), col = "black") + 
    theme(legend.position = "none") + 
    scale_size_continuous(range = c(0.5, 5)) +
    scale_y_continuous(breaks =c(2.5, 5, 7.5, 10, 12.5 ,15, 17.5)) + 
    geom_pointrange(data = r_dur_df, 
                    aes(y = Q50,  ymin = Q2.5, ymax = Q97.5, x = study_year), 
                    position = position_jitter(width = 0.2, height = 0), 
                    col = "gray40",
                    alpha = 0.6)) 

# plot posterior prediction of rate parameter of assumed exponential distribution
(plot.rate.lamba = pred.rate.tot %>%
    mutate(type = factor("lambda", labels = expression(lambda))) %>% 
    ggplot(aes(x=year, y = 1/rate_Q05)) + 
    geom_line(col = "orange3", size = 0.8)+ 
    geom_ribbon(aes(ymin = 1/rate_Q0025, ymax = 1/rate_Q0975), alpha = 0.3, fill = "orange3") + 
    theme_bw() + 
    facet_wrap(. ~ type , scales = "free", labeller = label_parsed) +
    theme(strip.text.x = element_text(size = 12)) +
    labs(y = expression("Estimate of parameter "*lambda*" of "*Exp(lambda)*" distribution [ "*year^-1*" ]"), x = "Year")+
    geom_ribbon(aes(ymin = 1/rate_Q0025_re, ymax = 1/rate_Q0975_re), alpha = 0.2, fill = "orange3") + 
    scale_x_continuous(breaks = c(1996, 2000, 2005, 2010, 2015, 2019)) +
     theme(legend.position = "none") +
    geom_pointrange(data = r_lambda_df, 
                    aes(y = Q50,  ymin = Q2.5, ymax = Q97.5, x = study_year), 
                    position = position_jitter(width = 0.2, height = 0), 
                    col = "gray40",
                    alpha = 0.6))



# final plot SW duration

cowplot::plot_grid(plot.dur, plot.rate.lamba, labels = LETTERS, nrow = 1, rel_widths = c(0.6, 0.4))



###################################################
### Sensitivity analysis : Age at entry into SW ###
###################################################

# predition data (1993 - 2013)
year.sens = c(seq(from = min(data.age.sensitivity$year_entry), to = max(data.age.sensitivity$year_entry), by = 0.25), max(data.age.sensitivity$year_entry))
year_centred.sens = year.sens - mean(data.age.sensitivity$year_entry)
data.new.age.sens = data.frame(intercept = rep(1, length(year.sens)), 
                               year_entry_centered = year_centred.sens)

# extrapolation data (2013-2019)
year.sens.pred = seq(from = ceiling(max(data.age.sensitivity$year_entry)), 
                     to = max(data.age.sensitivity$study_year), by = 0.25)
year_centered.sens.pred =  year.sens.pred - mean(data.age.sensitivity$year_entry)
data.new.age.sens.pred = data.frame(intercept = rep(1, length(year.sens.pred)),
                                    year_entry_centered = year_centered.sens.pred)

sim.sens = as.matrix(fit.age.sens) %>% 
  as_tibble()

pred.sens = as.matrix(sim.sens[, c(26,24)]) %*% 
  t(as.matrix(data.new.age.sens)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year.sens)


pred.sens.re = sim.sens[, c(26,24,27)] %>% 
  rowwise() %>% 
  mutate(mu = rnorm(1, a_mu , a_tau)) %>% select(mu, b) %>% as.matrix() %*% 
  t(as.matrix(data.new.age.sens)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year.sens) %>% 
  select(year, Q0975_re = Q0975, Q05_re = Q05, Q0025_re = Q0025)


pred.sens.tot = pred.sens %>% 
  left_join(pred.sens.re) 

pred.sens.pred = as.matrix(sim.sens[, c(26,24)]) %*% 
  t(as.matrix(data.new.age.sens.pred)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year.sens.pred)


pred.sens.re.pred = sim.sens[, c(26,24,27)] %>% 
  rowwise() %>% 
  mutate(mu = rnorm(1, a_mu , a_tau)) %>% select(mu, b) %>% as.matrix() %*% 
  t(as.matrix(data.new.age.sens.pred)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>%
  mutate_all(., function(x){x+cutoff}) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year.sens.pred) %>% 
  select(year, Q0975_re = Q0975, Q05_re = Q05, Q0025_re = Q0025)

pred.sens.tot.pred = pred.sens.pred %>% 
  left_join(pred.sens.re.pred) 

random_betasens = sim.sens %>% select(starts_with("beta["))
common_alphasens = sim.sens %>% select(alpha)

studylevel_sens = apply(random_betasens, MARGIN = 2, function(x){pull(common_alphasens)/x +10})

studylevel_sens_mean = apply(X =studylevel_sens, 
                             MARGIN = 2,
                             FUN = mean)

studylevel_sens_sd = apply(X = studylevel_sens, 
                           MARGIN = 2, 
                           FUN = sd)

studylevel_sens_quant  = apply(X = studylevel_sens, 
                               MARGIN = 2, 
                               FUN = quantile, 
                               probs = c(0.025, 0.5, 0.975))

studylevel_sens_quant = data.frame(t(studylevel_sens_quant))
names(studylevel_sens_quant) = c("Q2.5", "Q50", "Q97.5")

sens_df = data.frame(studylevel_sens_mean, studylevel_sens_sd, studylevel_sens_quant)
sens_df$year_entry = data.age.sensitivity$year_entry

sens_df$studysize_age = data.age.sensitivity$studysize_age

## combine results with main analysis

pred.sens.tot = pred.sens.tot %>% 
  mutate(type = "at entry into SW (sensitivity analysis)") 

pred.tot.age = pred.tot.age %>% 
  mutate(type = "during SW (main analysis)")

pred.sens.comb = pred.sens.tot %>% bind_rows(pred.tot.age)


data.sens.small = data.age.sensitivity %>% 
  select(age = mean_age_entry, Population_ID, year = year_entry, studysize_age) %>% 
  mutate(type = "at entry into SW (sensitivity analysis)")

data.age.small = data.age %>% 
  select(age = mean_age_adjusted, Population_ID, year = study_year, studysize_age) %>% 
  mutate(type = "during SW (main analysis)")

data.sens.comb = data.sens.small %>% bind_rows(data.age.small)

pred.sens.tot.pred = pred.sens.tot.pred %>% 
  mutate(type = "at entry into SW (sensitivity analysis)")

year.max = max(data.age.sensitivity$year_entry)

data.sens.comb2 = data.sens.comb %>% filter(type =="during SW (main analysis)") %>% 
  mutate(included = as.factor(as.numeric(Population_ID %in% data.sens.small$Population_ID)))


# posterior predictions of random intercepts of sensitivity analysis (age at entry model)

## posterior predictions of random intercept of sensitivity analysis model 

random_interc_entry = sim.sens %>% select(starts_with("a["))

a_entry_mean= apply(X = random_interc_entry, 
                    MARGIN = 2, 
                    FUN = mean)

a_entry_sd = apply(X = random_interc_entry, 
                   MARGIN = 2, 
                   FUN = sd)

a_entry_quant  = apply(X = random_interc_entry, 
                       MARGIN = 2, 
                       FUN = quantile, 
                       probs = c(0.025, 0.5, 0.975))

a_entry_quant = data.frame(t(a_entry_quant))
names(a_entry_quant) = c("Q2.5", "Q50", "Q97.5")

a_entry_df = data.frame(a_entry_mean, a_entry_sd, a_entry_quant)

a_entry_df$study = as.character(data.age.sensitivity$Population_ID)
a_entry_df = a_entry_df[order(a_entry_df$a_entry_mean), ]
a_entry_df$a_rank = c(1:dim(a_entry_df)[1])

random_interc_mean_entry = sim.sens %>% select(a_mu, a_tau) %>% 
  rowwise() %>% 
  mutate(a_mu_ci = rnorm(1, mean = a_mu, sd = a_tau)) 

a_mu_entry_dat = data.frame(amu_Q05 = rep(median(random_interc_mean_entry$a_mu), nrow(a_entry_df)+2),
                            amu_Q0025 = rep(quantile(random_interc_mean_entry$a_mu, 0.025), nrow(a_entry_df)+2), 
                            amu_Q0975 = rep(quantile(random_interc_mean_entry$a_mu, 0.975), nrow(a_entry_df)+2), 
                            amu_Q0025_re = rep(quantile(random_interc_mean_entry$a_mu_ci, 0.025), nrow(a_entry_df)+2), 
                            amu_Q0975_re = rep(quantile(random_interc_mean_entry$a_mu_ci, 0.975), nrow(a_entry_df)+2), 
                            rank = 0:(nrow(a_entry_df)+1))

### Plot: sensitivity analysis
#############################

(plot.sens = pred.sens.comb %>% 
    ggplot(aes(x=year, y = Q05, group = type,  fill = type)) +
    geom_line(aes(col=type), size = 0.8)+
    geom_ribbon(aes(ymin = Q0025, ymax = Q0975), alpha = 0.5) +
    theme_bw() +
    labs(y = "Mean age (95% CrI) of female sex workers [years]", x = "Year") + 
    geom_ribbon(data = pred.sens.comb, aes(ymin = Q0025_re, ymax = Q0975_re), alpha = 0.3) + 
    scale_color_manual(values=c("purple", "tomato"), name = "FSW age") +   
    scale_fill_manual(values=c("purple", "tomato"), name = "FSW age") + 
    geom_point(data = data.sens.comb2, 
               aes(x=year, y=age, size = studysize_age, col = type, shape = included)) + 
    scale_size_continuous(range = c(0.1,5)) +
    geom_line(data = pred.sens.tot.pred %>% filter(year>=year.max), 
             aes(x = year, y = Q05, col = type), lty = 2, size = 0.8) + 
   geom_ribbon(data = pred.sens.tot.pred %>% filter(year>=year.max), 
               aes(x = year, y = Q05, col = type, min = Q0025, ymax = Q0975),
               alpha = 0.2, lty= 1, fill = "gray60", size = 0.3) + 
   geom_ribbon(data = pred.sens.tot.pred %>% filter(year>=year.max),  
               aes(x = year, y = Q05, col = type,  min = Q0025_re, ymax = Q0975_re),  
               fill = "gray60", lty = 1, alpha = 0.2, size = 0.3) + 
   geom_pointrange(data = sens_df %>% mutate(type = "at entry into SW (sensitivity analysis)"), 
                   aes(y = Q50,  ymin = Q2.5, ymax = Q97.5, x = year_entry, col = type), 
                   position = position_jitter(width = 0.2, height = 0),  alpha = 0.6) + 
   theme(legend.position = c(0.2,0.9)) + 
   scale_shape_manual(values = c(3,1)) +
   guides(size ="none", shape = "none") +
   scale_x_continuous(breaks = c(1993,1996, 2000,2005, 2010, 2013,2015, 2019)) 
)

####
#### Plot Slopes and random intercepts of main and sensitivity analysis:
########################################################################

b.age = sim.fit.age %>% 
  select(fsw_age_slope = b) %>%
  as_tibble()

b.age.entry = sim.sens %>% 
  select(fsw_age_entry_slope = b) %>% 
  as_tibble() 

b.dur = sim.dur %>% 
  select(duration_slope=beta) %>% 
  as_tibble() 

b.tot = b.age %>% bind_cols(b.age.entry) %>% bind_cols(b.dur)

b.tot.sum = b.tot %>%  
  mutate_all(., function(x){exp(x*cutoff)}) %>%   
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_Q.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) 

(plot_slopes = b.tot.sum %>% 
    mutate(Q0025 = Q0025 -1, 
           Q05 = Q05-1, 
           Q0975 = Q0975 -1) %>% 
    mutate(type = ifelse(name == "duration_slope", "Duration model", "Age model")) %>% 
    mutate(name = factor(name, levels = rev(c("duration_slope", "fsw_age_slope", "fsw_age_entry_slope")), 
                         labels = rev(c("Main analysis", "Main analysis\n(FSW age during SW)", 
                                        "Sensitivity analysis\n(FSW age at entry into SW)")))) %>% 
    mutate(type = factor(type, levels = c("Duration model", "Age model"))) %>%  
    ggplot(aes(x = name, y = Q05)) + 
    geom_pointrange(aes(ymin = Q0025, ymax = Q0975, col = name),size = 0.8, lwd = 0.8) +
    geom_hline(yintercept=0, lty =2) + 
    theme_bw() + 
    facet_grid(type ~., scales = "free_y", space = "free", switch = "y") + 
    scale_color_manual(values = rev(c("dodgerblue", "tomato", "purple"))) + 
    theme(strip.text.y.left = element_text(angle = 0), 
          strip.placement = "outside",
          strip.switch.pad.grid = unit(0.25, "cm")) + 
    coord_flip(clip = "off") + 
    theme(legend.position = "none")  + 
    scale_y_continuous(labels = scales::percent) +
    labs(x=NULL, y = expression(paste("Percentage increase in outcome"^"&", "(SW duration/FSW age) by 10-year calendar time increase")), 
         caption = expression(paste(" "^"&", "Corresponds to exponentiated slope ", exp(gamma), " (Duration model) /", exp(omega), " (Age model)"))))


# random intercepts
a_df.temp = a_df %>% 
  mutate(type = "Age model\nMain analysis\n(FSW age during SW)")

a_entry_df.temp = a_entry_df %>% 
  mutate(type = "Age model\nSensitivity analysis\n(FSW age at entry into SW)") %>% 
  rename("a_mean"="a_entry_mean", "a_sd" = "a_entry_sd") 


alpha_df.temp = alpha_df %>%
  rename("a_mean"="alpha_mean", "a_sd" = "alpha_sd") %>% 
  mutate(type = "Duration model")


intercepts_tot = a_df.temp %>% bind_rows(a_entry_df.temp) %>% bind_rows(alpha_df.temp) %>% 
  filter(!is.na(study)) %>% 
  mutate(type = factor(type, levels = c("Duration model", 
                                        "Age model\nMain analysis\n(FSW age during SW)", 
                                        "Age model\nSensitivity analysis\n(FSW age at entry into SW)")))

a_mu_dat.temp = a_mu_dat %>% 
  mutate(type = "Age model\nMain analysis\n(FSW age during SW)") %>% 
  mutate(col = "1") %>% 
  rename('a_rank' = 'rank') %>% 
  left_join(a_df.temp %>% select(study, a_rank))

a_mu_entry_dat.temp = a_mu_entry_dat %>% 
  mutate(type = "Age model\nSensitivity analysis\n(FSW age at entry into SW)") %>% 
  mutate(col = "2")   %>% 
  rename('a_rank' = 'rank') %>% 
  left_join(a_entry_df.temp %>% select(study, a_rank))

dur_mu_dat.temp = dur_mu_dat %>% 
  mutate(type = "Duration model") %>% 
  mutate(col ="3") %>%
  rename('a_rank' = 'rank') %>% 
  left_join(alpha_df.temp %>% select(study, a_rank))


tot_mu.dat = a_mu_dat.temp %>% bind_rows(a_mu_entry_dat.temp) %>% bind_rows(dur_mu_dat.temp) %>% 
  filter(!is.na(study)) %>% 
  mutate(type = factor(type, levels = c("Duration model", 
                                        "Age model\nMain analysis\n(FSW age during SW)", 
                                        "Age model\nSensitivity analysis\n(FSW age at entry into SW)")))

(plot_intercepts = intercepts_tot %>%  
    ggplot(aes(x = study, y = Q50, group = type)) +
    geom_ribbon(data = tot_mu.dat, 
                aes(y = amu_Q05, ymin  = amu_Q0025, ymax = amu_Q0975, x = study, fill = col), 
                alpha = 0.4) +
    geom_ribbon(data = tot_mu.dat, aes(y = amu_Q05, ymin  = amu_Q0025_re, ymax = amu_Q0975_re, x = study, 
                                       fill = col), 
                alpha = 0.2) +
    geom_line(data = tot_mu.dat, 
              aes(y = amu_Q05, x = study, col =col), 
              size = 1, lty = 2) +
    geom_pointrange(aes(ymin = Q2.5, ymax = Q97.5)) + 
    # , 
    #                 position = position_jitter(width = 0.1, height = 0)) + 
    facet_grid(type ~ ., scales = "free_y") + 
    # scale_y_continuous() + 
    scale_color_manual(values = c("tomato", "purple", "dodgerblue")) +
    scale_fill_manual(values = c("tomato", "purple", "dodgerblue")) + 
    theme_bw() + 
    scale_y_continuous(breaks = seq(from = 0.4, to =4, by = 0.2), 
                       name = expression(paste("Varying intercept ", eta[j], 
                                               " (Duration model) / ", theta[j], " (Age model)"))) + 
    theme(legend.position = "none") + 
    labs(x = "Population ID"))



############
### Sensitivity analysis: Excluding Milovanovic study
######################################################

# posterior predictions:
#############################

# Duration model 

sim.dur.exclMilo = as.matrix(fit.duration.exclMilo) %>% 
  as_tibble()

pred.rate.exclMilo = as.matrix(sim.dur.exclMilo[, c(1,3)]) %*% 
  t(as.matrix(data.new.duration.exclMilo)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_duration.exclMilo) %>% 
  rename(rate_Q0025 = Q0025, 
         rate_Q05 = Q05, 
         rate_Q0975 = Q0975)

pred.rate.re.exclMilo = sim.dur.exclMilo[, c(1,2,3)] %>% 
  rowwise() %>% 
  mutate(mu = rnorm(1, mu , tau)) %>% select(mu, beta) %>% as.matrix() %*% 
  t(as.matrix(data.new.duration.exclMilo)) %>% 
  as_tibble() %>% 
  mutate_all(exp) %>% 
  summarise_all(list(Q0025 = ~ quantile(., probs = 0.025), Q05 = median, Q0975 = ~quantile(., probs = 0.975))) %>% 
  pivot_longer(cols = everything()) %>%
  mutate(quantile = gsub(".*\\_", "", name)) %>% 
  mutate(name = gsub("\\_.*", "", name)) %>% 
  pivot_wider(names_from = quantile, values_from = value) %>% 
  bind_cols(year = year_duration.exclMilo) %>% 
  rename(rate_Q0025_re = Q0025, 
         rate_Q05_re = Q05, 
         rate_Q0975_re = Q0975)


pred.rate.tot.exclMilo = pred.rate.exclMilo %>% 
  left_join(pred.rate.re.exclMilo) 

### posterior predictions for each study: 

r_lambda.exclMilo = sim.dur.exclMilo %>% select(starts_with("lambda["))


r_lambda_mean.exclMilo = apply(X =r_lambda.exclMilo, 
                      MARGIN = 2,
                      FUN = mean)

r_lambda_sd.exclMilo = apply(X = r_lambda.exclMilo, 
                    MARGIN = 2, 
                    FUN = sd)

r_lambda_quant.exclMilo  = apply(X = r_lambda.exclMilo, 
                        MARGIN = 2, 
                        FUN = quantile, 
                        probs = c(0.025, 0.5, 0.975))

r_lambda_quant.exclMilo = data.frame(t(r_lambda_quant.exclMilo))
names(r_lambda_quant.exclMilo) = c("Q2.5", "Q50", "Q97.5")

r_lambda_df.exclMilo = data.frame(r_lambda_mean.exclMilo, r_lambda_sd.exclMilo, r_lambda_quant.exclMilo)
r_lambda_df.exclMilo$study_year = data.duration.excl.Milovanovic$study_year

r_lambda.exclMilo = sim.dur.exclMilo %>% select(starts_with("lambda["))
r_dur.exclMilo = 1/r_lambda.exclMilo

r_dur_mean.exclMilo = apply(X =r_dur.exclMilo, 
                   MARGIN = 2,
                   FUN = mean)

r_dur_sd.exclMilo = apply(X = r_dur.exclMilo, 
                 MARGIN = 2, 
                 FUN = sd)

r_dur_quant.exclMilo  = apply(X = r_dur.exclMilo, 
                     MARGIN = 2, 
                     FUN = quantile, 
                     probs = c(0.025, 0.5, 0.975))

r_dur_quant.exclMilo = data.frame(t(r_dur_quant.exclMilo))
names(r_dur_quant.exclMilo) = c("Q2.5", "Q50", "Q97.5")

r_dur_df.exclMilo = data.frame(r_dur_mean.exclMilo, r_dur_sd.exclMilo, r_dur_quant.exclMilo)
r_dur_df.exclMilo$study_year = data.duration.exclMilo$study_year


# plot posterior prediciton of expected SW duration:

(plot.dur.exclMilo  = pred.rate.tot.exclMilo %>% 
   ggplot(aes(x=year, y = rate_Q05)) + 
   geom_line(col = "dodgerblue", size = 0.8)+ 
   geom_ribbon(aes(ymin = rate_Q0025, ymax = rate_Q0975), alpha = 0.3, fill = "dodgerblue") + 
   theme_bw() + 
   labs(y = "Mean duration (95%-CrI) of sex work [years]", x = "Study year") + 
   scale_x_continuous(breaks = c(1996, 2000, 2005, 2010, 2015, 2019)) +
   geom_ribbon(aes(ymin = rate_Q0025_re, ymax = rate_Q0975_re), alpha = 0.2, fill = "dodgerblue") + 
   geom_point(data = data.duration.exclMilo, aes(y = mean_duration_adjusted, x = study_year, size = studysize_duration), col = "black") + 
   theme(legend.position = "none") + 
   scale_size_continuous(range = c(0.5, 5)) +
   scale_y_continuous(breaks =c(2.5, 5, 7.5, 10, 12.5 ,15, 17.5)) + 
   geom_pointrange(data = r_dur_df.exclMilo, 
                   aes(y = Q50,  ymin = Q2.5, ymax = Q97.5, x = study_year), 
                   position = position_jitter(width = 0.2, height = 0), 
                   col = "gray40",
                   alpha = 0.6)) ### this looks a bit strange check it!

