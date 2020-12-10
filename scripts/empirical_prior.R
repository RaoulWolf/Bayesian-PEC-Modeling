empirical_priors <- function(formula, data) {
  
  data <- as.data.frame(data)
  
  terms <- brms::brmsterms(formula)
  
  if (terms$family$family != "lognormal") {
    stop("Only the lognormal distribution is supported.", call. = FALSE)
  }
  
  message("This function assumes no intercepts are present and the standard deviation is shared.")
  
  vars <- as.character(terms$dpars$mu$formula)
  yvar <- as.character(terms$dpars$mu$respform)[2]
  xvars <- unlist(strsplit(vars[2], split = "\\+"))
  xvars <- trimws(xvars)
  xvar <- xvars[2]
  
  if (!is.null(terms$adforms$cens) && grepl(",", as.character(terms$adforms$cens)[2])) {
    cens <- as.character(terms$adforms$cens)[2]
    cens <- gsub("resp_cens\\(", "", cens)
    cens <- gsub("\\)", "", cens)
    cens <- unlist(strsplit(cens, split = ", "))
    data[(data[, cens[1]] == "interval" | data[, cens[1]] == 2), yvar] <- data[(data[, cens[1]] == "interval" | data[, cens[1]] == 2), cens[2]]
  }
  
  template <- brms::get_prior(formula = formula, data = data)
  template <- template[(template$class == "cor" & nchar(template$group) > 1) | nchar(template$coef) > 0 | template$class == "sigma", ]
  
  mu_priors <- data.frame(prior = "", class = "b", coef = "", group = "", dpar = "", pre_coef = paste0(xvar, unlist(unique(data[, xvar]))), mean = NA_real_, se = NA_real_, stringsAsFactors = FALSE)
  mu_priors$coef <- gsub("[-]", "M", mu_priors$pre_coef)
  mu_priors$coef <- gsub("[+]", "P", mu_priors$coef)
  mu_priors$coef <- gsub("[^[:alnum:]^;^_^.]", "", mu_priors$coef)
  
  for (i in mu_priors$pre_coef) {
    mu_priors[mu_priors$pre_coef == i, ]$mean <- round(mean(log(unlist(data[data[, xvar] == sub(xvar, "", i), yvar]))))
    mu_priors[mu_priors$pre_coef == i, ]$se <- ceiling(4 * sd(log(unlist(data[data[, xvar] == sub(xvar, "", i), yvar]))) / sqrt(length(unlist(data[data[, xvar] == sub(xvar, "", i), yvar]))))
  }
  
  mu_priors$se <- ifelse(mu_priors$se == 0, max(mu_priors$se, na.rm = TRUE), mu_priors$se)
  mu_priors$se <- ifelse(is.na(mu_priors$se), max(mu_priors$se, na.rm = TRUE), mu_priors$se)
  mu_priors$prior <- paste0("normal(", mu_priors$mean, ", ", mu_priors$se, ")")
  mu_priors$mean <- NULL
  mu_priors$se <- NULL
  
  sigma_priors <- data.frame(prior = "", class = "sigma", coef = "", group = "", dpar = "", mean = NA_real_, se = NA_real_, stringsAsFactors = FALSE)

  sigma_mean <- c()

  for (i in mu_priors$pre_coef) {
    sigma_mean[i] <- sd(log(unlist(data[data[, xvar] == sub(xvar, "", i), yvar])), na.rm = TRUE)
  }

  sigma_priors$mean <- round(mean(sigma_mean[is.finite(sigma_mean)], na.rm = TRUE))
  sigma_priors$se <- ceiling(4 * sd(sigma_mean[is.finite(sigma_mean)], na.rm = TRUE))

  sigma_priors$prior <- paste0("normal(", sigma_priors$mean, ", ", sigma_priors$se, ")")
  sigma_priors$mean <- NULL
  sigma_priors$se <- NULL
  
  if (!is.null(terms$dpars$mu$re)) {
    pre_grid <- data.frame(coef = mu_priors$coef, pre_coef = mu_priors$pre_coef, stringsAsFactors = FALSE)
    group_grid <- data.frame(coef = rep(pre_grid$coef, times = length(terms$dpars$mu$re$group)), pre_coef = rep(pre_grid$pre_coef, times = length(terms$dpars$mu$re$group)), group = sort(rep(terms$dpars$mu$re$group, times = length(pre_grid$coef))), stringsAsFactors = FALSE)
    group_prior <- data.frame(prior = "", class = "sd", coef = group_grid$coef, group = group_grid$group, dpar = "", pre_coef = group_grid$pre_coef, sd = NA_real_, stringsAsFactors = FALSE)
      for (i in unique(group_prior$group)) {
        for (j in unique(group_prior$pre_coef)) {
          group_prior[group_prior$group == i & group_prior$pre_coef == j, ]$sd <- ceiling(4 * sd(log(unlist(data[data[, xvar] == sub(xvar, "", j), yvar]))))
        }
      }
      group_prior$sd <- ifelse(group_prior$sd == 0, max(group_prior$sd, na.rm = TRUE), group_prior$sd)
      group_prior$sd <- ifelse(is.na(group_prior$sd), max(group_prior$sd, na.rm = TRUE), group_prior$sd)
      group_prior$prior <- paste0("normal(0, ", group_prior$sd, ")")
      group_prior$pre_coef <- NULL
      group_prior$sd <- NULL
  }
  
  mu_priors$pre_coef <- NULL
  
  if (!is.null(terms$dpars$mu$re) && any(terms$dpars$mu$re$cor)) {
    cor_prior <- data.frame(prior = "lkj(1)", class = "cor", coef = "", group = terms$dpars$mu$re[terms$dpars$mu$re$cor, ]$group, dpar = "", stringsAsFactors = FALSE)
  }
  
  if (is.null(terms$dpars$mu$re)) {
    prior_list <- rbind(mu_priors, sigma_priors)
  } else if (!any(terms$dpars$mu$re$cor)) {
    prior_list <- rbind(mu_priors, group_prior, sigma_priors)
    } else {
      prior_list <- rbind(mu_priors, cor_prior, group_prior, sigma_priors)
  }
  
  pre_result <- template
  pre_result$prior <- prior_list$prior
  pre_result$class <- prior_list$class
  pre_result$coef <- prior_list$coef
  pre_result$group <- prior_list$group
  pre_result$dpar <- prior_list$dpar

  result <- pre_result
  rownames(result) <- NULL
  class(result) <- c("brmsprior", "data.frame")
  result
}
