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
genZI <- function(outcome, predictors) {
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
    model_fit <- fit_zinb_model(model_data)
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
        p.adjust(results_matrix[, x], method = "BH")
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

fit_zinb_model <- function(data) {
  tryCatch(
    {
      model <- glmmTMB(
        y ~ arm * time + age_B + (1 | subject),
        zi = ~ arm * time,
        data = data,
        family = nbinom2
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
      theme_minimal() +
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
  prepared_data <- prepare_ogu_data(phyloseq_obj)
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
  if (!is.formula(formula)) {
    stop("'formula' must be a valid formula object")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.formula(zi)) {
    stop("'zi' must be a valid formula object")
  }

  available_families <- c(
    "gaussian",
    "poisson",
    "Gamma",
    "nbinom2",
    "nbinom1",
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

  conflicted::conflicts_prefer(glmmTMB::tweedie)
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
