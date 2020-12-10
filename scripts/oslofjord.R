## The package "brms" need to be installed.
# install.packages("brms")

# Load package, prepare for parallel computing, and load custom empirical prior function
library(brms)
options(mc.cores = min(parallel::detectCores() - 1L, 5L))
source("scripts/empirical_prior.R")

## Load data
of_data <- read.csv("data/of_exposure.csv")

## Model 1
### Preparing the data
of_mod1_data <- of_data[of_data$censored != "left", ]
of_mod1_data$value <- ifelse(of_mod1_data$censored == "interval", of_mod1_data$loq / 2, of_mod1_data$value)
of_mod1_data$site <- NULL
of_mod1_data$date <- NULL
of_mod1_data$pubchem_cid <- NULL
of_mod1_data$loq <- NULL
of_mod1_data$lod <- NULL
of_mod1_data$censored <- NULL

### Preparing the formula
of_mod1_form <- brmsformula(formula = value ~ 0 + iupac_name,
                             family = lognormal(link_sigma = "identity"))

### Preparing the prior
of_mod1_prior <- empirical_priors(formula = of_mod1_form, data = of_mod1_data)

### Fitting the model
of_mod1_brm <- brm(formula = of_mod1_form,
                    data = of_mod1_data,
                    prior = of_mod1_prior,
                    sample_prior = "yes",
                    save_pars = save_pars(all = TRUE),
                    chains = 5,
                    iter = 4000,
                    control = list(adapt_delta = 0.99))

### Inspecting the model
summary(of_mod1_brm)

any(rhat(of_mod1_brm) > 1.01)

plot(of_mod1_brm)

fixef(of_mod1_brm)

pp_check(of_mod1_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(of_mod1_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(of_mod1_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 2
### Preparing the data
of_mod2_data <- of_data[of_data$censored != "left", ]
of_mod2_data$value <- ifelse(of_mod2_data$censored == "interval", of_mod2_data$loq / 2, of_mod2_data$value)
of_mod2_data$date <- NULL
of_mod2_data$pubchem_cid <- NULL
of_mod2_data$loq <- NULL
of_mod2_data$lod <- NULL
of_mod2_data$censored <- NULL

### Preparing the formula
of_mod2_form <- brmsformula(formula = value ~ 0 + iupac_name + (0 + iupac_name | site),
                             family = lognormal(link_sigma = "identity"))

### Preparing the prior
of_mod2_prior <- empirical_priors(formula = of_mod2_form, data = of_mod2_data)

### Fitting the model
of_mod2_brm <- brm(formula = of_mod2_form,
                    data = of_mod2_data,
                    prior = of_mod2_prior,
                    sample_prior = "yes",
                    save_pars = save_pars(all = TRUE),
                    chains = 5,
                    iter = 4000,
                    control = list(adapt_delta = 0.99))

### Inspecting the model
summary(of_mod2_brm)

any(rhat(of_mod2_brm) > 1.01)

plot(of_mod2_brm)

fixef(of_mod2_brm)

ranef(of_mod2_brm)

pp_check(of_mod2_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(of_mod2_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(of_mod2_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 3
### Preparing the data
of_mod3_data <- of_data
of_mod3_data$site <- NULL
of_mod3_data$date <- NULL
of_mod3_data$pubchem_cid <- NULL
of_mod3_data$y2 <- ifelse(of_mod3_data$censored == "none", of_mod3_data$value, ifelse(of_mod3_data == "left", of_mod3_data$lod, of_mod3_data$loq))
of_mod3_data$value <- ifelse(of_mod3_data$censored == "interval", of_mod3_data$lod, of_mod3_data$value)
of_mod3_data$loq <- NULL
of_mod3_data$lod <- NULL

### Preparing the formula
of_mod3_form <- brmsformula(formula = value | cens(censored, y2) ~ 0 + iupac_name,
                             family = lognormal(link_sigma = "identity"))

### Preparing the prior
of_mod3_prior <- empirical_priors(formula = of_mod3_form, data = of_mod3_data)

### Fitting the model
of_mod3_brm <- brm(formula = of_mod3_form,
                    data = of_mod3_data,
                    prior = of_mod3_prior,
                    sample_prior = "yes",
                    save_pars = save_pars(all = TRUE),
                    chains = 5,
                    iter = 4000,
                    control = list(adapt_delta = 0.99))

### Inspecting the model
summary(of_mod3_brm)

any(rhat(of_mod3_brm) > 1.01)

plot(of_mod3_brm)

fixef(of_mod3_brm)

pp_check(of_mod3_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(of_mod3_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(of_mod3_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model 4
### Preparing the data
of_mod4_data <- of_data
of_mod4_data$date <- NULL
of_mod4_data$pubchem_cid <- NULL
of_mod4_data$y2 <- ifelse(of_mod4_data$censored == "none", of_mod4_data$value, ifelse(of_mod4_data$censored == "left", of_mod4_data$lod, of_mod4_data$loq))
of_mod4_data$value <- ifelse(of_mod4_data$censored == "interval", of_mod4_data$lod, of_mod4_data$value)
of_mod4_data$loq <- NULL
of_mod4_data$lod <- NULL

### Preparing the formula

of_mod4_form <- brmsformula(formula = value | cens(censored, y2) ~ 0 + iupac_name + (0 + iupac_name | site),
                             family = lognormal(link_sigma = "identity"))

### Preparing the prior
of_mod4_prior <- empirical_priors(formula = of_mod4_form, data = of_mod4_data)

### Fitting the model
of_mod4_brm <- brm(formula = of_mod4_form,
                    data = of_mod4_data,
                    prior = of_mod4_prior,
                    sample_prior = "yes",
                    save_pars = save_pars(all = TRUE),
                    chains = 5,
                    iter = 4000,
                    control = list(adapt_delta = 0.99))

### Inspecting the model
summary(of_mod4_brm)

any(na.omit(rhat(of_mod4_brm)) > 1.01)

plot(of_mod4_brm)

fixef(of_mod4_brm)

ranef(of_mod4_brm)

pp_check(of_mod4_brm, type = "dens_overlay") +
  ggplot2::scale_x_log10()

pp_check(of_mod4_brm, type = "ecdf_overlay") +
  ggplot2::scale_x_log10()

plot(conditional_effects(of_mod4_brm))[[1]] +
  ggplot2::scale_y_log10()

## Model comparison

loo(of_mod1_brm, of_mod2_brm)

loo(of_mod3_brm, of_mod4_brm)
