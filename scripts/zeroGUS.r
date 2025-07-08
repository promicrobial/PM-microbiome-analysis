#!/usr/bin/env Rscript

#######################################
#                                     #
#      Zero-Inflated Model           #
#        Comparison Tools             #
#                                     #
#######################################

##############################################################################
# Author: [Your Name]                                                         #
# Date: [Current Date]                                                        #
# License: MIT                                                                #
# Version: 1.0.0                                                             #
#                                                                            #
# Description: This script provides functionality to compare different        #
# zero-inflated models across multiple distribution families. It performs     #
# model fitting, diagnostics, and comparison using various statistical        #
# criteria.                                                                   #
#                                                                            #
# Dependencies:                                                              #
#   - glmmTMB: For fitting zero-inflated models                             #
#     https://github.com/glmmTMB/glmmTMB                                    #
#   - DHARMa: For residual diagnostics                                      #
#     https://github.com/florianhartig/DHARMa                              #
#   - knitr: For formatted output                                           #
#     https://github.com/yihui/knitr                                        #
#                                                                            #
##############################################################################

#' Fit Zero-Inflated Negative Binomial Models for Multiple Outcomes
#'
#' @description
#' This function fits zero-inflated negative binomial mixed models for multiple outcomes
#' (e.g., OGUs - Operational Genomic Units) with specified predictors. It handles
#' multiple subjects and timepoints, incorporating random effects and adjusting for
#' multiple comparisons.
#'
#' @param outcome A matrix or data frame where columns represent different outcomes
#'   (e.g., OGUs) and rows represent observations. Column names must be unique
#'   identifiers.
#' @param predictors A data frame containing the following required columns:
#'   \itemize{
#'     \item subject: Factor identifying individual subjects
#'     \item time: Numeric time points of measurements
#'     \item arm: Factor indicating treatment/control groups
#'     \item age_B: Numeric baseline age values
#'   }
#'
#' @return A list containing three elements:
#'   \describe{
#'     \item{adj_results}{A data frame with adjusted p-values for each model parameter}
#'     \item{model_summary}{A data frame with model convergence information and fit statistics}
#'     \item{all_model_results}{A list containing full model objects for each outcome}
#'   }
#'
#' @details
#' The function fits zero-inflated negative binomial mixed models using glmmTMB.
#' It includes random intercepts for subjects and tests for treatment effects,
#' time effects, and their interaction. P-values are adjusted for multiple
#' comparisons using the Benjamini-Hochberg method.
#'
#' @section Input Validation:
#' The function performs several validation checks:
#' \itemize{
#'   \item Ensures outcome is a matrix or data frame
#'   \item Checks for NA values in predictors and outcome
#'   \item Validates required columns in predictors
#' }
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' outcome_data <- matrix(rnorm(100), ncol = 5)
#' predictor_data <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   time = rep(c(0,1), 10),
#'   arm = rep(c("treatment", "control"), each = 10),
#'   age_B = rnorm(20)
#' )
#'
#' # Run analysis
#' results <- genZI(outcome_data, predictor_data)
#' }
#'
#' @export
#' @importFrom glmmTMB glmmTMB
#' @importFrom stats p.adjust
#' @importFrom dplyr arrange
genZI <- function(outcome, predictors, family = nbinom2) {
  if (!is.matrix(outcome) && !is.data.frame(outcome)) {
    stop("Outcome must be a matrix or data frame")
  }

  if (anyNA(predictors)) {
    stop("Predictor variables cannot contain NA values")
  }

  # Ensure proper factor levels
  predictors$subject <- factor(predictors$subject)
  predictors$arm <- factor(predictors$arm)
  predictors$time <- as.numeric(predictors$time)
  predictors$age_B <- as.numeric(predictors$age_B)

  # Verify numeric variables
  if (!is.numeric(predictors$time) || !is.numeric(predictors$age_B)) {
    stop("time and age_B must be numeric variables")
  }

  # Initialize variables
  n_outcomes <- ncol(outcome)
  outcome_ids <- colnames(outcome)
  if (is.null(outcome_ids)) {
    outcome_ids <- paste0("OGU_", 1:n_outcomes)
    colnames(outcome) <- outcome_ids
  }

  # Progress counter
  counter <- make_progress_counter(
    total = ncol(outcome),
    title = "Processing OGUs"
  )

  # Initialize results matrices and data frames
  coef_names <- c("intercept", "time", "arm", "age_B", "interaction")
  results_matrix <- matrix(NA, nrow = n_outcomes, ncol = length(coef_names))
  colnames(results_matrix) <- coef_names

  # Initialize model summary
  model_summary <- data.frame(
    outcome_id = outcome_ids,
    converged = logical(n_outcomes),
    AIC = numeric(n_outcomes),
    BIC = numeric(n_outcomes),
    zero_prop = numeric(n_outcomes),
    error_message = character(n_outcomes),
    stringsAsFactors = FALSE
  )

  # Store all model results
  all_models <- vector("list", n_outcomes)
  names(all_models) <- outcome_ids
  for (i in 1:n_outcomes) {
    counter()
    current_outcome <- outcome[, i]

    # Check for NA values in outcome
    if (anyNA(current_outcome)) {
      model_summary$error_message[i] <- "Outcome contains NA values"
      next
    } else {
      model_summary$error_message[i] <- "No errors"
    }

    # Prepare model data
    model_data <- data.frame(
      y = current_outcome,
      predictors
    ) %>%
      arrange(subject, time)

    # Calculate proportion of zeros
    model_summary$zero_prop[i] <- mean(current_outcome == 0)

    # Fit model
    model_fit <- fit_zinb_model(model_data, family)
    all_models[[i]] <- model_fit$model

    # Update model summary
    model_summary$converged[i] <- model_fit$converged
    if (!is.null(model_fit$error)) {
      model_summary$error_message[i] <- model_fit$error
      next
    }

    # Get model diagnostics
    diagnostics <- zinb_diagnostics(model_fit$model, model_data)

    # Update summary statistics
    model_summary$AIC[i] <- diagnostics$fit_stats$AIC
    model_summary$BIC[i] <- diagnostics$fit_stats$BIC

    if (model_fit$converged) {
      results_matrix[i, ] <- diagnostics$coefficients$conditional[1:5, 4]
    }
  }

  # Create results with adjusted p-values
  adj_results <- data.frame(
    outcome_id = outcome_ids,
    setNames(
      lapply(coef_names, function(x) {
        round(p.adjust(results_matrix[, x], method = "BH"), 4)
      }),
      coef_names
    )
  )

  return(list(
    adj_results = adj_results,
    model_summary = model_summary,
    all_models = all_models
  ))
}

fit_zinb_model <- function(data, family) {
  tryCatch(
    {
      model <- glmmTMB(
        y ~ arm * time + age_B + (1 | subject),
        zi = ~ arm * time,
        data = data,
        family = family
      )
      return(list(
        model = model,
        converged = TRUE,
        error = NULL
      ))
    },
    error = function(e) {
      return(list(
        model = NULL,
        converged = FALSE,
        error = as.character(e)
      ))
    }
  )
}

#' Perform comprehensive ZINB model diagnostics
#' @param model Fitted glmmTMB model
#' @param data Original data used to fit the model
#' @return List containing diagnostic results
zinb_diagnostics <- function(model, data) {
  require(DHARMa)
  require(performance)

  # 1. Basic convergence check
  convergence <- list(
    converged = !is.null(model$fit$convergence) &&
      model$fit$convergence == 0,
    optimizer = model$fit$optimizer,
    message = model$fit$message
  )

  # 2. Model coefficients and their stability
  coef_summary <- list(
    conditional = summary(model)$coefficients$cond,
    zero_inflation = summary(model)$coefficients$zi,
    random_effects = VarCorr(model)
  )

  # 3. Create DHARMa residuals
  dharma_residuals <- tryCatch(
    {
      simulationOutput <- simulateResiduals(model, plot = FALSE)
      list(
        residuals = residuals(simulationOutput),
        scaled_residuals = simulationOutput$scaledResiduals,
        tests = list(
          ks_test = testUniformity(simulationOutput, plot = FALSE),
          dispersion = testDispersion(simulationOutput, plot = FALSE),
          outliers = testOutliers(simulationOutput, plot = FALSE),
          zero_inflation = testZeroInflation(simulationOutput, plot = FALSE)
        )
      )
    },
    error = function(e) NULL
  )

  # 4. Model fit statistics
  fit_stats <- list(
    AIC = AIC(model),
    BIC = BIC(model),
    logLik = logLik(model),
    df.residual = df.residual(model)
  )

  # 5. Zero-inflation assessment
  zero_stats <- list(
    observed_zeros = sum(data$y == 0),
    total_obs = nrow(data),
    zero_proportion = mean(data$y == 0)
  )

  # 6. Random effects assessment
  ranef_stats <- tryCatch(
    {
      re <- ranef(model)
      list(
        variance_components = VarCorr(model),
        icc = performance::icc(model),
        ranef_summary = list(
          mean = mean(unlist(re$cond$subject)),
          sd = sd(unlist(re$cond$subject)),
          quantiles = quantile(
            unlist(re$cond$subject),
            probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
          )
        )
      )
    },
    error = function(e) NULL
  )

  # 7. Model predictions
  predictions <- list(
    fitted_values = fitted(model),
    predicted_probs = predict(model, type = "response")
  )

  # Return all diagnostics
  return(list(
    convergence = convergence,
    coefficients = coef_summary,
    dharma_residuals = dharma_residuals,
    fit_stats = fit_stats,
    zero_stats = zero_stats,
    ranef_stats = ranef_stats,
    predictions = predictions,
    model = model
  ))
}

#' Create diagnostic plots for ZINB model
#' @param diagnostic_results Results from zinb_diagnostics function
#' @param data Original data used to fit the model
#' @return List of ggplot objects
plot_zinb_diagnostics <- function(diagnostic_results) {
  require(ggplot2)
  require(patchwork)

  model <- diagnostic_results$model
  data <- diagnostic_results$model$frame

  plot_list <- list()

  # 1. Residuals vs Fitted
  plot_list$resid_fitted <- ggplot(
    data.frame(
      fitted = diagnostic_results$predictions$fitted_values,
      residuals = diagnostic_results$dharma_residuals$scaled_residuals
    ),
    aes(x = fitted, y = residuals)
  ) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") +
    labs(
      title = "Residuals vs Fitted",
      x = "Fitted values",
      y = "Scaled residuals"
    ) +
    theme_minimal()

  # 2. Q-Q plot of residuals
  plot_list$qq <- ggplot(
    data.frame(
      residuals = diagnostic_results$dharma_residuals$scaled_residuals
    ),
    aes(sample = residuals)
  ) +
    stat_qq() +
    stat_qq_line(color = "red") +
    labs(title = "Normal Q-Q Plot of Residuals") +
    theme_minimal()

  # 3. Random effects plot
  if (!is.null(diagnostic_results$ranef_stats)) {
    re_data <- data.frame(
      subject = names(unlist(ranef(model)$cond$subject)),
      effect = unlist(ranef(model)$cond$subject)
    )
    plot_list$ranef <- ggplot(re_data, aes(y = effect, x = subject)) +
      geom_point() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      labs(
        title = "Random Effects by Subject",
        x = "Subject",
        y = "Random effect"
      ) +
      theme(axis.text.x = element_blank())
  }

  # 4. Zero probability plot
  plot_list$zero_prob <- ggplot(
    data.frame(
      observed = data[, 1] == 0,
      predicted = diagnostic_results$predictions$predicted_probs
    ),
    aes(x = predicted, fill = observed)
  ) +
    geom_density(alpha = 0.5) +
    labs(
      title = "Predicted Probabilities by Observed Zeros",
      x = "Predicted probability",
      y = "Density"
    ) +
    theme_minimal()

  # 5. DHARMa simulated residual plots
  # plot_list$simqq <- plot(plotQQunif(simulateResiduals(model)))
  # plot_list$simresid <- plot(plotResiduals(simulateResiduals(model)))

  # Combine plots
  # BUG cannot combine with DHARMa plots
  # cowplot::plot_grid(plot_list)

  # temp solution - remove DHARMa plots for user plots themselves
  (plot_list$resid_fitted + plot_list$qq) /
    (plot_list$ranef + plot_list$zero_prob)
}

#' Print summary of ZINB diagnostics
#' @param diagnostic_results Results from zinb_diagnostics function
#' @return Printed summary of diagnostics
print_zinb_diagnostics <- function(diagnostic_results) {
  cat("ZINB Model Diagnostics Summary\n")
  cat("============================\n\n")

  # Convergence
  cat("1. Convergence:\n")
  cat("   - Model converged:", diagnostic_results$convergence$converged, "\n")
  cat("   - Optimizer:", diagnostic_results$convergence$optimizer, "\n")
  if (!is.null(diagnostic_results$convergence$message)) {
    cat("   - Message:", diagnostic_results$convergence$message, "\n")
  }
  cat("\n")

  # Model fit statistics
  cat("2. Model Fit Statistics:\n")
  cat("   - AIC:", diagnostic_results$fit_stats$AIC, "\n")
  cat("   - BIC:", diagnostic_results$fit_stats$BIC, "\n")
  cat("   - Log-likelihood:", diagnostic_results$fit_stats$logLik, "\n")
  cat("\n")

  # Zero-inflation
  cat("3. Zero-inflation Summary:\n")
  cat(
    "   - Observed zeros:",
    diagnostic_results$zero_stats$observed_zeros,
    "\n"
  )
  cat(
    "   - Zero proportion:",
    round(diagnostic_results$zero_stats$zero_proportion, 3),
    "\n"
  )
  cat(
    "   - Zero-inflation test p-value:",
    round(diagnostic_results$dharma_residuals$tests$zero_inflation$p.value, 3),
    "\n\n"
  )

  # Random effects
  if (!is.null(diagnostic_results$ranef_stats)) {
    cat("4. Random Effects Summary:\n")
    cat("   - ICC:", round(diagnostic_results$ranef_stats$icc, 3), "\n")
    cat(
      "   - Random effects SD:",
      round(
        sqrt(attr(
          diagnostic_results$ranef_stats$variance_components$cond$subject,
          "stddev"
        )),
        3
      ),
      "\n\n"
    )
  }

  # Coefficient summary
  cat("5. Fixed Effects Summary:\n")
  cat("Conditional model:\n")
  print(diagnostic_results$coefficients$conditional)
  cat("\nZero-inflation model:\n")
  print(diagnostic_results$coefficients$zero_inflation)
}

#' Run complete ZINB analysis pipeline
#' @param phyloseq_obj Input phyloseq object
#' @param output_dir Directory for saving results
#' @return List containing all analysis results
run_zinb_analysis <- function(phyloseq_obj, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Prepare data
  prepared_data <- prepare_ogu_data(
    absAbund,
    filter_prevalence = TRUE,
    prevalence_cutoff = 0.4,
    filter_abundance = FALSE
  )
  predictors <- prepare_predictors(prepared_data$meta_data)

  # Store original data with sample names
  predictor_data <- predictors %>%
    rownames_to_column(var = "sample_names") %>%
    as.data.frame()

  # Run ZINB models
  model_results <- run_zinb_models(
    prepared_data$ogu_data,
    predictors,
    colnames(prepared_data$ogu_data)
  )

  # Identify significant OGUs
  significant_ogus <- identify_significant_ogus(
    model_results$adj_results,
    prepared_data$taxa_data
  )

  # Extract model results for significant OGUs
  significant_ogu_models <- model_results$all_model_results[
    significant_ogus$rowname
  ]

  # Create detailed results for significant OGUs
  significant_ogu_details <- lapply(significant_ogu_models, function(x) {
    if (!is.null(x$best_model)) {
      best_model <- x$best_model
      model_fit <- x$model_fits[[x$best_form]]

      list(
        ogu_id = x$ogu_id,
        model = best_model,
        model_form = x$best_form,
        converged = model_fit$converged,
        AIC = model_fit$AIC,
        BIC = model_fit$BIC,
        coef_summary = model_fit$coef_summary,
        diagnostics = zinb_diagnostics(
          best_model,
          data = data.frame(
            y = best_model@frame$y,
            predictors[rownames(best_model@frame), ]
          )
        )
      )
    }
  })

  # Save results
  write.csv(
    model_results$adj_results,
    file.path(output_dir, "OGU-ZINB-adj-p-values.csv"),
    row.names = FALSE
  )

  write.csv(
    model_results$model_summary,
    file.path(output_dir, "OGU-ZINB-model-summary.csv"),
    row.names = FALSE
  )

  write.csv(
    significant_ogus,
    file.path(output_dir, "significant-OGU-ZINB-adj-p-values.csv"),
    row.names = FALSE
  )

  # Save significant OGU model results
  saveRDS(
    significant_ogu_details,
    file.path(output_dir, "significant-OGU-model-details.rds")
  )

  # Save predictor data
  write.csv(
    predictor_data,
    file.path(output_dir, "predictor-variables.csv"),
    row.names = FALSE
  )

  # Generate and save diagnostic plots for significant OGUs
  for (ogu in names(significant_ogu_details)) {
    if (!is.null(significant_ogu_details[[ogu]])) {
      pdf(file.path(output_dir, paste0("diagnostics-", ogu, ".pdf")))
      print(plot_zinb_diagnostics(
        significant_ogu_details[[ogu]]$diagnostics,
        data.frame(
          y = significant_ogu_details[[ogu]]$model@frame$y,
          predictors[rownames(significant_ogu_details[[ogu]]$model@frame), ]
        )
      ))
      dev.off()
    }
  }

  # Return comprehensive results
  return(list(
    zinb_results = model_results$adj_results,
    model_summary = model_results$model_summary,
    significant_ogus = significant_ogus,
    significant_ogu_models = significant_ogu_details,
    prepared_data = prepared_data,
    predictor_data = predictor_data,
    analysis_info = list(
      date = Sys.Date(),
      n_samples = nrow(predictor_data),
      n_ogus = ncol(prepared_data$ogu_data),
      n_significant = nrow(significant_ogus)
    )
  ))
}

#' Compare Zero-Inflated Models Across Different Distributions
#'
#' @description
#' This function fits and compares zero-inflated models using different distribution
#' families. It provides comprehensive model diagnostics, including AIC/BIC
#' comparisons, residual analysis, and convergence checks.
#'
#' @param formula An object of class "formula" specifying the model structure
#' @param data A data frame containing the variables in the model
#' @param zi A one-sided formula for zero-inflation component
#' @param family Character vector specifying distribution families to test.
#' @param plot Logical; whether to generate diagnostic plots (default: TRUE)
#'
#' @return A list containing:
#'   \describe{
#'     \item{models}{List of fitted models for each distribution}
#'     \item{summary}{Data frame with model comparison statistics}
#'     \item{diagnostics}{List of diagnostic information for each model}
#'     \item{plots}{If plot=TRUE, diagnostic plots for each model}
#'   }
#'
#' @details
#' The function performs multiple checks and validations:
#' \itemize{
#'   \item Model convergence via Hessian positive-definiteness
#'   \item Residual diagnostics using DHARMa
#'   \item Zero-inflation assessment
#'   \item Prediction validity checks
#' }
#'
#' @examples
#' \dontrun{
#' data(example_data)
#' result <- compareZI(
#'   count ~ treatment + time + (1|subject),
#'   data = example_data,
#'   zi = ~treatment,
#'   family = c("poisson", "nbinom2")
#' )
#' }
#'
#' @export
#' @importFrom glmmTMB glmmTMB
#' @importFrom DHARMa simulateResiduals
#' @importFrom knitr kable
#'
#' @note
#' This function requires the following packages:
#' \itemize{
#'   \item glmmTMB (>= 1.1.4)
#'   \item DHARMa (>= 0.4.5)
#'   \item knitr (>= 1.42)
#' }
compareZI <- function(
  formula,
  data,
  zi,
  family = available_families,
  plot = TRUE
) {
  # Validate inputs
  if (!plyr::is.formula(formula)) {
    stop("'formula' must be a valid formula object")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!plyr::is.formula(zi)) {
    stop("'zi' must be a valid formula object")
  }

  available_families <- c(
    "gaussian",
    "poisson",
    "Gamma",
    "nbinom1",
    "nbinom2",
    "nbinom12",
    "compois",
    "truncated_compois",
    "genpois",
    "truncated_genpois",
    "truncated_poisson",
    "truncated_nbinom2",
    "truncated_nbinom1",
    "beta_family",
    "betabinomial",
    "tweedie",
    "skewnormal",
    "nbinom1",
    "nbinom12",
    "beta_family",
    "lognormal",
    "ziGamma",
    "t_family",
    "ordbeta",
    "bell"
  )

  invalid_families <- setdiff(family, available_families)
  if (length(invalid_families) > 0) {
    stop(
      "Invalid family specifications: ",
      paste(invalid_families, collapse = ", ")
    )
  }

  suppressMessages(conflicted::conflicts_prefer(glmmTMB::tweedie))
  # Initialize lists to store results
  models <- list()
  fit_stats <- list()
  plots <- list()

  # Try each family
  for (fam_name in family) {
    tryCatch(
      {
        # Get the family function
        fam_fun <- get(fam_name)

        # Fit model
        model <- glmmTMB::glmmTMB(
          formula = formula,
          ziformula = zi,
          family = fam_fun,
          data = data
        )

        # Store model
        models[[fam_name]] <- model

        # Calculate fit statistics
        sim_res <- tryCatch(
          {
            DHARMa::simulateResiduals(model, n = 250)
          },
          error = function(e) {
            warning(paste(
              "Residual simulation failed for",
              fam_name,
              ":",
              e$message
            ))
            return(NULL)
          }
        )

        # Check for infinite/NA predictions
        pred_check <- tryCatch(
          {
            preds <- predict(model, type = "response")
            has_inf <- any(is.infinite(preds))
            has_na <- anyNA(preds)
            if (has_inf || has_na) {
              warning(paste(
                "Model",
                fam_name,
                "produces infinite or NA predictions"
              ))
            }
            list(infinite = has_inf, na = has_na)
          },
          error = function(e) {
            list(infinite = NA, na = NA)
          }
        )

        fit_stats[[fam_name]] <- list(
          AIC = AIC(model),
          BIC = BIC(model),
          converged = !is.null(model$sdr$pdHess) && model$sdr$pdHess,
          simResiduals = sim_res,
          pred_check = pred_check,
          # Check eigenvalues of Hessian
          hessian_check = tryCatch(
            {
              eigen_vals <- eigen(model$sdr$hessian)$values
              list(
                negative_eigenvals = any(eigen_vals < 0),
                condition_number = max(abs(eigen_vals)) / min(abs(eigen_vals))
              )
            },
            error = function(e) {
              list(negative_eigenvals = NA, condition_number = NA)
            }
          )
        )

        # Create diagnostic plots if requested
        if (plot && !is.null(sim_res)) {
          plots[[fam_name]] <- tryCatch(
            {
              list(
                residuals = DHARMa::plot(sim_res, ask = FALSE),
                qq = qqnorm(residuals(model)),
                fitted = plot(fitted(model), residuals(model))
              )
            },
            error = function(e) NULL
          )
        }
      },
      error = function(e) {
        fit_stats[[fam_name]] <- list(
          AIC = NA,
          BIC = NA,
          converged = FALSE,
          simResiduals = NA,
          pred_check = list(infinite = NA, na = NA),
          hessian_check = list(negative_eigenvals = NA, condition_number = NA),
          error = conditionMessage(e)
        )
        plots[[fam_name]] <- NULL
      }
    )
  }

  # Create summary dataframe
  summary_df <- do.call(
    rbind,
    lapply(names(fit_stats), function(fam) {
      data.frame(
        family = fam,
        AIC = fit_stats[[fam]]$AIC,
        BIC = fit_stats[[fam]]$BIC,
        `positive-definite Hessian (converged)` = fit_stats[[fam]]$converged,
        negative_eigenvals = fit_stats[[fam]]$hessian_check$negative_eigenvals,
        condition_number = fit_stats[[fam]]$hessian_check$condition_number,
        infinite_preds = fit_stats[[fam]]$pred_check$infinite,
        na_preds = fit_stats[[fam]]$pred_check$na,
        error = if (!is.null(fit_stats[[fam]]$error)) {
          fit_stats[[fam]]$error
        } else {
          NA
        }
      )
    })
  )

  # Add zero-inflation statistics
  zero_stats <- lapply(names(models), function(mod_name) {
    mod <- models[[mod_name]]
    if (!is.null(mod)) {
      tryCatch(
        {
          obs_zeros <- sum(mod$frame[, 1] == 0)
          pred_zeros <- sum(predict(mod, type = "response") < 0.0001)
          return(c(observed = obs_zeros, predicted = pred_zeros))
        },
        error = function(e) {
          return(c(observed = NA, predicted = NA))
        }
      )
    } else {
      return(c(observed = NA, predicted = NA))
    }
  })

  # Ensure zero_stats has same length as summary_df
  zero_df <- do.call(rbind, zero_stats)
  if (nrow(zero_df) == nrow(summary_df)) {
    summary_df$observed_zeros <- zero_df[, "observed"]
    summary_df$predicted_zeros <- zero_df[, "predicted"]
  }

  # Return results
  result <- list(
    models = models,
    summary = summary_df,
    diagnostics = fit_stats,
    plots = if (plot) plots else NULL
  )

  # Print formatted summary
  cat("\nModel Comparison Summary:\n")
  print(knitr::kable(summary_df, digits = 2))

  # Only print best model if there are valid AIC/BIC values
  if (any(!is.na(summary_df$AIC))) {
    cat(
      "\nBest model by AIC:",
      summary_df$family[which.min(summary_df$AIC)],
      "\n"
    )
  }
  if (any(!is.na(summary_df$BIC))) {
    cat(
      "Best model by BIC:",
      summary_df$family[which.min(summary_df$BIC)],
      "\n"
    )
  }

  return(invisible(result))
}

#' Create Summary Plots for Zero-Inflated Model Results
#'
#' @description
#' Creates a comprehensive set of diagnostic and summary plots for zero-inflated model results,
#' including p-value distributions, convergence status, zero proportions, and model fit criteria.
#'
#' @param results A list containing the following components:
#'   \itemize{
#'     \item adj_results: Data frame with adjusted p-values for different effects
#'     \item model_summary: Data frame with model convergence and fit information
#'   }
#'
#' @return A list of ggplot objects containing:
#'   \itemize{
#'     \item p_values: Distribution of adjusted p-values by effect
#'     \item convergence: Bar plot of model convergence status
#'     \item zero_prop: Distribution of zero proportions
#'     \item fit_criteria: Distribution of AIC and BIC values
#'     \item error_summary: Summary of error messages (if present)
#'     \item significant: Bar plot of significant results by effect
#'     \item combined: Combined plot of all visualizations
#'   }
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#'
#' @examples
#' \dontrun{
#' results <- genZI(outcome, predictors)
#' plots <- create_summary_plots(results)
#' print(plots$combined)
#' }
#'
#' @export
create_summary_plots <- function(results) {
  require(ggplot2)
  require(dplyr)
  require(patchwork)

  # Initialize plot list
  plots <- list()

  # 1. P-value distributions for different effects
  plots$p_values <- results$adj_results %>%
    tidyr::pivot_longer(
      cols = c(time, arm, age_B, interaction),
      names_to = "effect",
      values_to = "p_value"
    ) %>%
    # Remove non-finite values
    filter(is.finite(p_value)) %>%
    {
      ggplot(., aes(x = log10(p_value), fill = effect)) +
        geom_histogram(bins = 50, position = "dodge") +
        facet_wrap(~effect, scales = "free") +
        labs(
          title = "Distribution of Adjusted P-values by Effect",
          x = "Adjusted P-value (log10)",
          y = "Count",
          caption = "Red line indicates a significance cut off equal to 0.05",
        ) +
        geom_vline(
          xintercept = log10(0.05),
          linetype = "dashed",
          color = "red"
        ) +
        annotate(
          "text",
          x = -Inf,
          y = Inf,
          label = paste(
            "Number of features with q ≤ 0.05:",
            sum(.$p_value <= 0.05)
          ),
          hjust = 1,
          vjust = -1
        )
    }

  # 2. Model convergence summary
  plots$convergence <- results$model_summary %>%
    ggplot(aes(x = converged)) +
    geom_bar(aes(fill = converged)) +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    labs(
      title = "Model Convergence Summary",
      x = "Convergence Status",
      y = "Count"
    )

  # 3. Zero proportion distribution
  plots$zero_prop <- results$model_summary %>%
    # Remove non-finite values
    filter(is.finite(zero_prop)) %>%
    ggplot(aes(x = zero_prop)) +
    geom_histogram(bins = 30, fill = "steelblue") +
    labs(
      title = "Distribution of Zero Proportions",
      x = "Proportion of Zeros",
      y = "Count"
    )

  # 4. Model fit criteria
  plots$fit_criteria <- results$model_summary %>%
    tidyr::pivot_longer(
      cols = c(AIC, BIC),
      names_to = "criterion",
      values_to = "value"
    ) %>%
    # Remove non-finite values
    filter(is.finite(value)) %>%
    ggplot(aes(x = value, fill = criterion)) +
    geom_histogram(bins = 30, position = "dodge") +
    facet_wrap(~criterion, scales = "free") +
    labs(
      title = "Distribution of Model Fit Criteria",
      x = "Value",
      y = "Count"
    )

  # 5. Error message summary if present
  if (any(!is.na(results$model_summary$error_message))) {
    error_summary <- results$model_summary %>%
      filter(!is.na(error_message)) %>%
      count(error_message) %>%
      # Ensure n is finite
      filter(is.finite(n))

    if (nrow(error_summary) > 0) {
      plots$error_summary <- error_summary %>%
        ggplot(aes(x = reorder(error_message, n), y = n)) +
        geom_col() +
        coord_flip() +
        # Adjust text size based on number of error messages
        theme(
          axis.text.y = element_text(
            size = min(max(8, 12 - nrow(error_summary) / 10), 12)
          )
        ) +
        labs(
          title = "Summary of Error Messages",
          x = "Error Message",
          y = "Count"
        )
    }
  }

  # 6. Significant results summary
  if (!is.null(results$adj_results)) {
    sig_threshold <- 0.05
    sig_counts <- results$adj_results %>%
      summarise(across(
        c(time, arm, age_B, interaction),
        ~ sum(. < sig_threshold & is.finite(.)) # Only count finite values
      )) %>%
      tidyr::pivot_longer(
        everything(),
        names_to = "effect",
        values_to = "count"
      )

    plots$significant <- ggplot(sig_counts, aes(x = effect, y = count)) +
      geom_col(fill = "steelblue") +
      # Rotate x-axis labels for better readability
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = "Number of Significant Results by Effect",
        x = "Effect",
        y = "Count of Significant Results"
      )
  }

  # Combine plots using patchwork with explicit layout
  # Create base layout with 2 rows, 2 columns for main plots
  combined_plot <- (plots$p_values +
    plots$convergence +
    plots$zero_prop +
    plots$fit_criteria) +
    plot_layout(ncol = 2, nrow = 2)

  # If we have additional plots, add them in new rows
  if (!is.null(plots$error_summary) && !is.null(plots$significant)) {
    # Add both additional plots in a new row
    combined_plot <- combined_plot /
      (plots$error_summary + plots$significant) +
      plot_layout(heights = c(2, 1))
  } else if (!is.null(plots$error_summary)) {
    # Add just error summary
    combined_plot <- combined_plot /
      plots$error_summary +
      plot_layout(heights = c(2, 1))
  } else if (!is.null(plots$significant)) {
    # Add just significant results
    combined_plot <- combined_plot /
      plots$significant +
      plot_layout(heights = c(2, 1))
  }

  # Add overall theme adjustments
  combined_plot <- combined_plot &
    theme(
      plot.title = element_text(size = 12),
      axis.title = element_text(size = 10)
    )

  plots$combined <- combined_plot

  return(plots)
}

#' Create Zero Proportion Plot
#'
#' @description
#' Creates a scatter plot showing the proportion of zeros for each outcome,
#' colored by model convergence status.
#'
#' @param results A list containing model_summary with convergence and zero proportion information
#'
#' @return A ggplot object showing zero proportions by outcome, colored by convergence status
#'
#' @import ggplot2
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' results <- genZI(outcome, predictors)
#' zero_plot(results)
#' }
#'
#' @export
zero_plot <- function(results) {
  conv_text <- sprintf(
    "%d/%d models converged (%.1f%%)",
    sum(results$model_summary$converged),
    nrow(results$model_summary),
    100 * sum(results$model_summary$converged) / nrow(results$model_summary)
  )

  results$model_summary %>%
    mutate(
      convergence_status = ifelse(converged, "Converged", "Failed")
    ) %>%
    ggplot(aes(
      x = outcome_id,
      y = zero_prop,
      color = convergence_status
    )) +
    geom_point(size = 3) +
    scale_color_manual(
      name = "Convergence Status",
      values = c("Failed" = "#e78284", "Converged" = "#8CAAEE")
    ) +
    xlab(NULL) +
    ylab("Proportion of Zeros") +
    theme(
      axis.text.x = element_blank(),
      legend.position = "top"
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = conv_text,
      hjust = 1,
      vjust = -1
    )
}

#' Process DHARMa Diagnostics for a Single Model
#'
#' @description
#' Performs DHARMa diagnostic tests on a single model, including uniformity,
#' dispersion, outliers, and zero-inflation tests.
#'
#' @param model A fitted model object
#' @param model_id Character string identifying the model
#'
#' @return A data frame containing:
#'   \itemize{
#'     \item model_id: Model identifier
#'     \item test: Name of diagnostic test
#'     \item p.value: P-value from test
#'     \item warning: Any warnings generated during testing
#'     \item n_obs: Number of observations
#'     \item n_zeros: Number of zero values
#'   }
#'
#' @import DHARMa
#'
#' @examples
#' \dontrun{
#' model_diagnostics <- process_dharma_diagnostics(fitted_model, "model1")
#' }
#'
#' @export
process_dharma_diagnostics <- function(model, model_id) {
  if (is.null(model)) {
    return(data.frame(
      model_id = model_id,
      test = NA,
      p.value = NA,
      test_status = "Model NULL", # New status column
      warning = "Null model",
      n_obs = NA,
      n_zeros = NA,
      stringsAsFactors = FALSE
    ))
  }

  tryCatch(
    {
      # Capture warnings during simulation
      sim_warnings <- character()
      withCallingHandlers(
        r <- simulateResiduals(model),
        warning = function(w) {
          sim_warnings <<- c(sim_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )

      # Get basic model information
      n_obs <- length(r$observedResponse)
      n_zeros <- sum(r$observedResponse == 0)

      # List of tests to run with specific parameters
      tests <- list(
        Uniformity = function() testUniformity(r, plot = FALSE),
        Dispersion = function() testDispersion(r, plot = FALSE),
        # Use bootstrap for outlier test
        Outliers = function() testOutliers(r, type = "bootstrap", plot = FALSE),
        ZeroInflation = function() testZeroInflation(r, plot = FALSE)
      )

      # Run each test with error handling
      results <- lapply(names(tests), function(test_name) {
        test_warnings <- character()

        test_result <- withCallingHandlers(
          tryCatch(
            {
              result <- tests[[test_name]]()

              # Determine test status
              status <- if (is.null(result$p.value)) {
                "No test performed"
              } else if (is.na(result$p.value)) {
                "Test not applicable"
              } else {
                "Test completed"
              }

              list(
                p.value = result$p.value,
                status = status,
                error = NULL
              )
            },
            error = function(e) {
              list(
                p.value = NA,
                status = "Error in test",
                error = as.character(e)
              )
            }
          ),
          warning = function(w) {
            test_warnings <<- c(test_warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }
        )

        # Create result dataframe
        data.frame(
          model_id = model_id,
          test = test_name,
          p.value = test_result$p.value,
          test_status = test_result$status, # New column
          warning = if (length(test_warnings) > 0) {
            paste(test_warnings, collapse = "; ")
          } else {
            NA
          },
          n_obs = n_obs,
          n_zeros = n_zeros,
          stringsAsFactors = FALSE
        )
      })

      # Combine all test results
      do.call(rbind, results)
    },
    error = function(e) {
      data.frame(
        model_id = model_id,
        test = NA,
        p.value = NA,
        test_status = "Error in model processing", # New status
        warning = as.character(e),
        n_obs = NA,
        n_zeros = NA,
        stringsAsFactors = FALSE
      )
    }
  )
}

#' Process DHARMa Results for Multiple Models
#'
#' @description
#' Processes DHARMa diagnostic results for multiple models and creates
#' comprehensive summaries and visualizations.
#'
#' @param results List of fitted models to be diagnosed
#'
#' @return A list containing:
#'   \itemize{
#'     \item complete: Complete diagnostic results for all models
#'     \item summary: Summary statistics for each diagnostic test
#'     \item warnings: Summary of warnings by test type
#'     \item wide_format: Wide format version of diagnostic results
#'   }
#'
#' @import dplyr
#' @import tidyr
#'
#' @examples
#' \dontrun{
#' dharma_results <- process_dharma_results(fitted_models)
#' print(dharma_results$summary)
#' }
#'
#' @export
process_dharma_results <- function(results) {
  # Process all models
  all_diagnostics <- lapply(seq_along(results), function(i) {
    process_dharma_diagnostics(
      results[[i]],
      names(results)[i]
    )
  })

  # Combine all results
  complete_diagnostics <- do.call(rbind, all_diagnostics)

  # Create enhanced summary statistics
  test_summary <- complete_diagnostics %>%
    group_by(test) %>%
    summarise(
      n_total = n(),
      n_completed = sum(test_status == "Test completed"),
      n_na = sum(test_status == "Test not applicable"),
      n_error = sum(test_status == "Error in test"),
      n_significant = sum(p.value < 0.05, na.rm = TRUE),
      mean_pval = mean(p.value, na.rm = TRUE),
      median_pval = median(p.value, na.rm = TRUE),
      n_warnings = sum(!is.na(warning)),
      .groups = 'drop'
    )

  # Enhanced warning summary
  warning_summary <- complete_diagnostics %>%
    filter(!is.na(warning)) %>%
    group_by(test, test_status) %>%
    summarise(
      warnings = list(unique(warning)),
      .groups = 'drop'
    )

  # Create wide format version with status
  wide_diagnostics <- complete_diagnostics %>%
    select(-warning) %>%
    tidyr::pivot_wider(
      id_cols = c(model_id, n_obs, n_zeros),
      names_from = test,
      values_from = c(p.value, test_status)
    )

  return(list(
    complete = complete_diagnostics,
    summary = test_summary,
    warnings = warning_summary,
    wide_format = wide_diagnostics
  ))
}

#' Plot DHARMa Diagnostic Summary
#'
#' @description
#' Creates a comprehensive visualization of DHARMa diagnostic results,
#' including p-value distributions, test success rates, and sample size information.
#'
#' @param dharma_results Output from process_dharma_results function
#'
#' @return A patchwork object combining four plots:
#'   \itemize{
#'     \item P-value distributions by test
#'     \item Proportion of significant tests
#'     \item Distribution of observation counts
#'     \item Distribution of zero proportions
#'   }
#'
#' @import ggplot2
#' @import patchwork
#'
#' @examples
#' \dontrun{
#' dharma_results <- process_dharma_results(fitted_models)
#' plot_dharma_summary(dharma_results)
#' }
#'
#' @export
plot_dharma_summary <- function(dharma_results) {
  require(ggplot2)
  require(patchwork)

  # Create status summary plot
  p1 <- ggplot(dharma_results$summary, aes(x = test)) +
    geom_col(aes(y = n_total), fill = "lightgrey") +
    geom_col(aes(y = n_completed), fill = "steelblue") +
    geom_text(aes(y = n_total, label = sprintf("NA: %d", n_na)), vjust = -0.5) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Test Completion Status", y = "Count", x = "Test")

  # Modified p-value distribution plot
  p2 <- ggplot(
    dharma_results$complete %>%
      filter(test_status == "Test completed"),
    aes(x = p.value, fill = test)
  ) +
    geom_histogram(bins = 20, position = "dodge") +
    facet_wrap(~test) +
    theme_minimal() +
    labs(title = "Distribution of p-values (completed tests only)")

  # Enhanced observation count plot
  p3 <- ggplot(dharma_results$wide_format, aes(x = n_obs)) +
    geom_histogram(bins = 20) +
    theme_minimal() +
    labs(
      title = "Distribution of observations",
      subtitle = sprintf("Total models: %d", nrow(dharma_results$wide_format))
    )

  # Modified zero proportion plot with test status
  p4 <- dharma_results$complete %>%
    filter(!is.na(n_zeros)) %>%
    ggplot(aes(x = n_zeros / n_obs, fill = test_status)) +
    geom_histogram(bins = 20, position = "stack") +
    theme_minimal() +
    labs(title = "Zero proportions by test status")

  # Combine plots
  (p1 + p2) / (p3 + p4)
}
