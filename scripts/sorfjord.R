## The package "brms" need to be installed.
# install.packages("brms")

# Load package, prepare for parallel computing, and load custom empirical prior function
library(brms)
options(mc.cores = min(parallel::detectCores() - 1L, 5L))
source("scripts/empirical_prior.R")

## Load data
sf_data <- read.csv("data/sf_exposure.csv")

## Model 1
### Preparing the data
sf_mod1_data <- sf_data[sf_data$censored != "left", ]
sf_mod1_data$site <- NULL
sf_mod1_data$date <- NULL
sf_mod1_data$value <- ifelse(sf_mod1_data$censored == "interval", sf_mod1_data$loq / 2, sf_mod1_data$value)
sf_mod1_data$loq <- NULL
sf_mod1_data$lod <- NULL
sf_mod1_data$censored <- NULL

### Preparing the formula
sf_mod1_form <- brmsformula(formula = value ~ 0 + iupac_name,
                            family = lognormal(link_sigma = "identity"))

### Preparing the prior
sf_mod1_prior <- empirical_priors(formula = sf_mod1_form, data = sf_mod1_data)

### Fitting the model
sf_mod1_brm <- brm(formula = sf_mod1_form,
                   data = sf_mod1_data,
                   prior = sf_mod1_prior,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE),
                   chains = 5,
                   iter = 4000,
                   control = list(adapt_delta = 0.99))

### Inspecting the model
summary(sf_mod1_brm)

any(rhat(sf_mod1_brm) > 1.01)

plot(sf_mod1_brm)

fixef(sf_mod1_brm)

pp_check(sf_mod1_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(sf_mod1_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(sf_mod1_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 2
### Preparing the data
sf_mod2_data <- sf_data[sf_data$censored != "left", ]
sf_mod2_data$value <- ifelse(sf_mod2_data$censored == "interval", sf_mod2_data$loq / 2, sf_mod2_data$value)
sf_mod2_data$loq <- NULL
sf_mod2_data$lod <- NULL
sf_mod2_data$censored <- NULL

### Preparing the formula
sf_mod2_form <- brmsformula(formula = value ~ 0 + iupac_name + (0 + iupac_name | site) + (0 + iupac_name | date),
                            family = lognormal(link_sigma = "identity"))

### Preparing the prior
sf_mod2_prior <- empirical_priors(formula = sf_mod2_form, data = sf_mod2_data)

### Fitting the model
sf_mod2_brm <- brm(formula = sf_mod2_form,
                   data = sf_mod2_data,
                   prior = sf_mod2_prior,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE),
                   chains = 5,
                   iter = 4000,
                   control = list(adapt_delta = 0.99))

### Inspecting the model
summary(sf_mod2_brm)

any(rhat(sf_mod2_brm) > 1.01)

plot(sf_mod2_brm)

fixef(sf_mod2_brm)

ranef(sf_mod2_brm)

pp_check(sf_mod2_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(sf_mod2_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(sf_mod2_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 3
### Preparing the data
sf_mod3_data <- sf_data
sf_mod3_data$site <- NULL
sf_mod3_data$date <- NULL
sf_mod3_data$y2 <- ifelse(sf_mod3_data$censored == "none", sf_mod3_data$value, ifelse(sf_mod3_data$censored == "left", sf_mod3_data$lod, sf_mod3_data$loq))
sf_mod3_data$value <- ifelse(sf_mod3_data$censored == "interval", sf_mod3_data$lod, sf_mod3_data$value)
sf_mod3_data$loq <- NULL
sf_mod3_data$lod <- NULL

### Preparing the formula
sf_mod3_form <- brmsformula(formula = value | cens(censored, y2) ~ 0 + iupac_name,
                            family = lognormal(link_sigma = "identity"))

### Preparing the prior
sf_mod3_prior <- empirical_priors(formula = sf_mod3_form, data = sf_mod3_data)

### Fitting the model
sf_mod3_brm <- brm(formula = sf_mod3_form,
                   data = sf_mod3_data,
                   prior = sf_mod3_prior,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE),
                   chains = 5,
                   iter = 4000,
                   control = list(adapt_delta = 0.99))

### Inspecting the model
summary(sf_mod3_brm)

any(rhat(sf_mod3_brm) > 1.01)

plot(sf_mod3_brm)

fixef(sf_mod3_brm)

pp_check(sf_mod3_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(sf_mod3_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(sf_mod3_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 4
### Preparing the data
sf_mod4_data <- sf_data
sf_mod4_data$y2 <- ifelse(sf_mod4_data$censored == "none", sf_mod4_data$value, ifelse(sf_mod4_data$censored == "left", sf_mod4_data$lod, sf_mod4_data$loq))
sf_mod4_data$value <- ifelse(sf_mod4_data$censored == "interval", sf_mod4_data$lod, sf_mod4_data$value)
sf_mod4_data$loq <- NULL
sf_mod4_data$lod <- NULL

### Preparing the formula

sf_mod4_form <- brmsformula(formula = value | cens(censored, y2) ~ 0 + iupac_name + (0 + iupac_name | site) + (0 + iupac_name | date),
                            family = lognormal(link_sigma = "identity"))

### Preparing the prior
sf_mod4_prior <- empirical_priors(formula = sf_mod4_form, data = sf_mod4_data)

### Fitting the model
sf_mod4_brm <- brm(formula = sf_mod4_form,
                   data = sf_mod4_data,
                   prior = sf_mod4_prior,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE),
                   chains = 5,
                   iter = 4000,
                   control = list(adapt_delta = 0.99))

### Inspecting the model
summary(sf_mod4_brm)

any(na.omit(rhat(sf_mod4_brm)) > 1.01)

plot(sf_mod4_brm)

fixef(sf_mod4_brm)

ranef(sf_mod4_brm)

pp_check(sf_mod4_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(sf_mod4_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(sf_mod4_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model comparison

loo(sf_mod1_brm, sf_mod2_brm)

loo(sf_mod3_brm, sf_mod4_brm)
