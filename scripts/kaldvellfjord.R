## The package "brms" need to be installed.
# install.packages("brms")

# Load package, prepare for parallel computing, and load custom empirical prior function
library(brms)
options(mc.cores = min(parallel::detectCores() - 1L, 5L))
source("scripts/empirical_prior.R")

## Load data
kvf_data <- read.csv("data/kvf_exposure.csv")

## Model 1
### Preparing the data
kvf_mod1_data <- kvf_data[kvf_data$censored != "left", ]
kvf_mod1_data$area <- NULL
kvf_mod1_data$site <- NULL
kvf_mod1_data$date <- NULL
kvf_mod1_data$depth <- NULL
kvf_mod1_data$value <- ifelse(kvf_mod1_data$censored == "interval", kvf_mod1_data$loq / 2, kvf_mod1_data$value)
kvf_mod1_data$loq <- NULL
kvf_mod1_data$lod <- NULL
kvf_mod1_data$censored <- NULL

### Preparing the formula
kvf_mod1_form <- brmsformula(formula = value ~ 0 + iupac_name,
                            family = lognormal(link_sigma = "identity"))

### Preparing the prior
kvf_mod1_prior <- empirical_priors(formula = kvf_mod1_form, data = kvf_mod1_data)

### Fitting the model
kvf_mod1_brm <- brm(formula = kvf_mod1_form,
                   data = kvf_mod1_data,
                   prior = kvf_mod1_prior,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE),
                   chains = 5,
                   iter = 4000,
                   control = list(adapt_delta = 0.99))

### Inspecting the model
summary(kvf_mod1_brm)

any(rhat(kvf_mod1_brm) > 1.01)

plot(kvf_mod1_brm)

fixef(kvf_mod1_brm)

pp_check(kvf_mod1_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(kvf_mod1_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(kvf_mod1_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 2
### Preparing the data
kvf_mod2_data <- kvf_data[kvf_data$censored != "left", ]
kvf_mod2_data$depth <- as.character(kvf_mod2_data$depth)
kvf_mod2_data$value <- ifelse(kvf_mod2_data$censored == "interval", kvf_mod2_data$loq / 2, kvf_mod2_data$value)
kvf_mod2_data$loq <- NULL
kvf_mod2_data$lod <- NULL
kvf_mod2_data$censored <- NULL

### Preparing the formula
kvf_mod2_form <- brmsformula(formula = value ~ 0 + iupac_name + (0 + iupac_name | area/site) + (0 + iupac_name | date) + (0 + iupac_name | depth),
                             family = lognormal(link_sigma = "identity"))

### Preparing the prior
kvf_mod2_prior <- empirical_priors(formula = kvf_mod2_form, data = kvf_mod2_data)

### Fitting the model
kvf_mod2_brm <- brm(formula = kvf_mod2_form,
                   data = kvf_mod2_data,
                   prior = kvf_mod2_prior,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE),
                   chains = 5,
                   iter = 4000,
                   control = list(adapt_delta = 0.99))

### Inspecting the model
summary(kvf_mod2_brm)

any(rhat(kvf_mod2_brm) > 1.01)

plot(kvf_mod2_brm)

fixef(kvf_mod2_brm)

ranef(kvf_mod2_brm)

pp_check(kvf_mod2_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(kvf_mod2_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(kvf_mod2_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 3
### Preparing the data
kvf_mod3_data <- kvf_data
kvf_mod3_data$area <- NULL
kvf_mod3_data$site <- NULL
kvf_mod3_data$date <- NULL
kvf_mod3_data$depth <- NULL
kvf_mod3_data$y2 <- ifelse(kvf_mod3_data$censored == "none", kvf_mod3_data$value, ifelse(kvf_mod3_data$censored == "left", kvf_mod3_data$lod, kvf_mod3_data$loq))
kvf_mod3_data$value <- ifelse(kvf_mod3_data$censored == "interval", kvf_mod3_data$lod, kvf_mod3_data$value)
kvf_mod3_data$loq <- NULL
kvf_mod3_data$lod <- NULL

### Preparing the formula
kvf_mod3_form <- brmsformula(formula = value | cens(censored, y2) ~ 0 + iupac_name,
                            family = lognormal(link_sigma = "identity"))

### Preparing the prior
kvf_mod3_prior <- empirical_priors(formula = kvf_mod3_form, data = kvf_mod3_data)

### Fitting the model
kvf_mod3_brm <- brm(formula = kvf_mod3_form,
                   data = kvf_mod3_data,
                   prior = kvf_mod3_prior,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE),
                   chains = 5,
                   iter = 4000,
                   control = list(adapt_delta = 0.99))

### Inspecting the model
summary(kvf_mod3_brm)

any(rhat(kvf_mod3_brm) > 1.01)

plot(kvf_mod3_brm)

fixef(kvf_mod3_brm)

pp_check(kvf_mod3_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(kvf_mod3_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(kvf_mod3_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 4
### Preparing the data
kvf_mod4_data <- kvf_data
kvf_mod4_data$depth <- as.character(kvf_mod4_data$depth)
kvf_mod4_data$y2 <- ifelse(kvf_mod4_data$censored == "none", kvf_mod4_data$value, ifelse(kvf_mod4_data$censored == "left", kvf_mod4_data$lod, kvf_mod4_data$loq))
kvf_mod4_data$value <- ifelse(kvf_mod4_data$censored == "interval", kvf_mod4_data$lod, kvf_mod4_data$value)
kvf_mod4_data$loq <- NULL
kvf_mod4_data$lod <- NULL

### Preparing the formula

kvf_mod4_form <- brmsformula(formula = value | cens(censored, y2) ~ 0 + iupac_name + (0 + iupac_name | area/site) + (0 + iupac_name | date) + (0 + iupac_name | depth),
                             family = lognormal(link_sigma = "identity"))

### Preparing the prior
kvf_mod4_prior <- empirical_priors(formula = kvf_mod4_form, data = kvf_mod4_data)

### Fitting the model
kvf_mod4_brm <- brm(formula = kvf_mod4_form,
                   data = kvf_mod4_data,
                   prior = kvf_mod4_prior,
                   sample_prior = "yes",
                   save_pars = save_pars(all = TRUE),
                   chains = 5,
                   iter = 4000,
                   control = list(adapt_delta = 0.99))

### Inspecting the model
summary(kvf_mod4_brm)

any(na.omit(rhat(kvf_mod4_brm)) > 1.01)

plot(kvf_mod4_brm)

fixef(kvf_mod4_brm)

ranef(kvf_mod4_brm)

pp_check(kvf_mod4_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(kvf_mod4_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(kvf_mod4_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model comparison

loo(kvf_mod1_brm, kvf_mod2_brm)

loo(kvf_mod3_brm, kvf_mod4_brm)
