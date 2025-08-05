#' Validate phyloseq object structure
#' @param phyloseq_obj A phyloseq object to validate
#' @return TRUE if valid, stops with error if invalid
validate_phyloseq_input <- function(phyloseq_obj) {
  if (class(phyloseq_obj) != "phyloseq") {
    stop(
      "Missing required phyloseq components: needs OTU table, taxonomy table, and sample data"
    )
  }
  return(TRUE)
}

#' Prepare OGU data from phyloseq object with optional filtering
#' @param phyloseq_obj A phyloseq object containing OTU table, taxonomy, and sample data
#' @param filter_prevalence Logical, whether to apply prevalence filter (default: FALSE)
#' @param prevalence_cutoff Numeric, minimum proportion of non-zero values (default: 0.4)
#' @param filter_abundance Logical, whether to apply abundance filter (default: FALSE)
#' @param abundance_quantile Numeric, quantile to use for abundance filter (default: 0.9)
#' @param abundance_cutoff Numeric, minimum value at specified quantile (default: 1)
#' @param   rownames_as_samples Indicates if the input count data have sample names at rownames (default: FALSE)
#' @return A list containing prepared OGU data, taxa data, and metadata
prepare_ogu_data <- function(
  phyloseq_obj,
  filter_prevalence = FALSE,
  prevalence_cutoff = 0.4,
  filter_abundance = FALSE,
  abundance_quantile = 0.9,
  abundance_cutoff = 1,
  rownames_as_samples = FALSE
) {
  validate_phyloseq_input(phyloseq_obj)

  tryCatch(
    {
      # Transform OGU data
      if(rownames_as_samples == FALSE){
        ogu_data <- as.data.frame(t(otu_table(phyloseq_obj)))
      } else {
        ogu_data <- as.data.frame(otu_table(phyloseq_obj))
      }

      # Store original dimensions
      original_dims <- dim(ogu_data)

      # Initialize filtering summary
      filter_summary <- list(
        original_features = ncol(ogu_data),
        filtered_features = ncol(ogu_data),
        filters_applied = character()
      )

      # Apply prevalence filter if requested
      if (filter_prevalence) {
        prevalence_index <- apply(ogu_data, 2, function(X) {
          sum(X > 0) > prevalence_cutoff * length(X)
        })

        ogu_data <- ogu_data[, prevalence_index]

        filter_summary$filters_applied <- c(
          filter_summary$filters_applied,
          sprintf("Prevalence (cutoff: %s)", prevalence_cutoff)
        )
        filter_summary$features_after_prevalence <- ncol(ogu_data)
      }

      # Apply abundance filter if requested
      if (filter_abundance) {
        abundance_index <- apply(ogu_data, 2, function(X) {
          quantile(X, abundance_quantile) > abundance_cutoff
        })

        ogu_data <- ogu_data[, abundance_index]

        filter_summary$filters_applied <- c(
          filter_summary$filters_applied,
          sprintf(
            "Abundance (quantile: %s, cutoff: %s)",
            abundance_quantile,
            abundance_cutoff
          )
        )
        filter_summary$features_after_abundance <- ncol(ogu_data)
      }

      # Update final number of filtered features
      filter_summary$filtered_features <- ncol(ogu_data)
      filter_summary$features_removed <- original_dims[2] - ncol(ogu_data)
      filter_summary$proportion_retained <- ncol(ogu_data) / original_dims[2]

      # Extract taxonomy and metadata
      taxa_data <- as.data.frame(tax_table(phyloseq_obj))
      meta_data <- as.data.frame(as.matrix(sample_data(phyloseq_obj)))

      # Subset taxonomy data to match filtered OGUs
      taxa_data <- taxa_data[colnames(ogu_data), ]

      # Verify data alignment
      if (!all(rownames(ogu_data) == rownames(meta_data))) {
        stop("Row names of OGU data and metadata do not match")
      }

      # Print filtering summary
      cat("\nOGU Filtering Summary:\n")
      cat("=====================\n")
      cat(
        "Original number of features:",
        filter_summary$original_features,
        "\n"
      )
      cat(
        "Filters applied:",
        paste(filter_summary$filters_applied, collapse = ", "),
        "\n"
      )
      if (filter_prevalence) {
        cat(
          "Features after prevalence filter:",
          filter_summary$features_after_prevalence,
          "\n"
        )
      }
      if (filter_abundance) {
        cat(
          "Features after abundance filter:",
          filter_summary$features_after_abundance,
          "\n"
        )
      }
      cat("Final number of features:", filter_summary$filtered_features, "\n")
      cat("Features removed:", filter_summary$features_removed, "\n")
      cat(
        "Proportion of features retained:",
        round(filter_summary$proportion_retained * 100, 2),
        "%\n\n"
      )

      return(list(
        ogu_data = ogu_data,
        taxa_data = taxa_data,
        meta_data = meta_data,
        filter_summary = filter_summary
      ))
    },
    error = function(e) {
      stop("Error in prepare_ogu_data: ", conditionMessage(e))
    }
  )
}

#' Prepare predictor variables for GEE analysis
#' @param meta_data Metadata dataframe
#' @return Processed predictor variables dataframe
prepare_predictors <- function(meta_data) {
  tryCatch(
    {
      predictors <- meta_data %>%
        mutate(
          time = scale_values(as.numeric(timesincefeedstart))
        ) %>%
        group_by(displayid) %>%
        mutate(agedays_rounded = as.numeric(agedays_rounded)) %>%
        arrange(agedays_rounded, .by_group = TRUE) %>%
        mutate(
          age_B = case_when(
            timepoint == "Baseline" ~ agedays_rounded,
            TRUE ~ NA # Changed FALSE to TRUE for proper case_when syntax
          )
        ) %>%
        fill(age_B) %>%
        ungroup() %>%
        mutate(
          age_B = scale_values(age_B),
          arm = as.factor(arm),
          subject = str_remove(displayid, "-G|-B"),
          timepoint = timepoint
        ) %>%
        column_to_rownames(var = "sample_names") %>%
        select(
          subject,
          arm,
          timepoint,
          time,
          age_B,
          zlen
        )

      return(predictors)
    },
    error = function(e) {
      stop("Error in prepare_predictors: ", conditionMessage(e))
    }
  )
}

# #' Run GEE model with convergence checking for a single OGU
# #' @param df Data frame containing outcome and predictor variables
# #' @param ogu_id Identifier for the OGU being analyzed
# #' @return List containing model results and convergence information
# run_single_gee <- function(df, ogu_id) {
#   # Try different correlation structures
#   models <- compare_gee_structures(formula = y ~ time * arm + age_B, data = df)

#   return(list(
#     models = models$models,
#     convergence_info = models$summary,
#     best_structure = models$best_structure,
#     best_model = if(!is.null(models$best_structure)) models$models[[models$best_structure]] else NULL,
#     ogu_id = ogu_id
#   ))
# }
#' Run GEE model with convergence checking for a single OGU
#' @param df Data frame containing outcome and predictor variables
#' @param ogu_id Identifier for the OGU being analyzed
#' @param link_family Family object for defining link and variance functions. See stats::family.
#' @return List containing model results and convergence information
run_single_gee <- function(
  df,
  ogu_id,
  link_family = gaussian(link = "identity")
) {
  # Try different correlation structures
  correlation_structures <- c(
    "Exchangeable",
    "Independence",
    "Unstructured",
    "Stationary-M-dependent(1)"
  )

  models <- list()
  convergence_info <- list()

  for (corstr in correlation_structures) {
    tryCatch(
      {
        model <- glmtoolbox::glmgee(
          y ~ time * arm + age_B,
          id = subject,
          data = df,
          corstr = corstr,
          family = link_family,
          verbose = FALSE
        )

        # Store model
        models[[corstr]] <- model

        # Check convergence using different criteria
        convergence_info[[corstr]] <- list(
          converged = !is.null(model$converged) && model$converged,
          QIC = tryCatch(glmtoolbox::QIC(model)$QIC, error = function(e) NA),
          CIC = tryCatch(glmtoolbox::CIC(model)$CIC, error = function(e) NA),
          coef_summary = summary(model)$coefficients,
          warnings = model$messages
        )
      },
      error = function(e) {
        convergence_info[[corstr]] <- list(
          converged = FALSE,
          error_message = conditionMessage(e)
        )
      }
    )
  }

  # Select best model based on QIC
  qic_values <- sapply(convergence_info, function(x) x$QIC)
  best_structure <- names(which.min(qic_values))

  # If no QIC available, use first converged model
  if (all(is.na(qic_values))) {
    converged_models <- names(which(sapply(convergence_info, function(x) {
      x$converged
    })))
    best_structure <- if (length(converged_models) > 0) {
      converged_models[1]
    } else {
      NULL
    }
  }

  return(list(
    models = models,
    convergence_info = convergence_info,
    best_structure = best_structure,
    best_model = if (!is.null(best_structure)) {
      models[[best_structure]]
    } else {
      NULL
    },
    ogu_id = ogu_id
  ))
}

#' Run GEE models for each OGU with convergence diagnostics
#' @param outcome Matrix of OGU abundances
#' @param predictors Predictor variables dataframe
#' @param ogu_ids Vector of OGU identifiers
#' @return List containing model results and convergence information
run_gee_models <- function(
  outcome,
  predictors,
  ogu_ids,
  link_family = gaussian(link = "identity")
) {
  # Validate inputs
  if (!is.matrix(outcome) && !is.data.frame(outcome)) {
    stop("Outcome must be a matrix or data frame")
  }
  if (ncol(outcome) != length(ogu_ids)) {
    stop("Number of OGU IDs must match number of outcome columns")
  }

  n_ogus <- ncol(outcome)
  results <- matrix(NA, nrow = n_ogus, ncol = 5)
  colnames(results) <- c("intercept", "time", "arm", "age_B", "interaction")

  # Store all model results
  all_model_results <- list()

  convergence_summary <- data.frame(
    ogu_id = ogu_ids,
    best_correlation = character(n_ogus),
    converged = logical(n_ogus),
    QIC = numeric(n_ogus),
    CIC = numeric(n_ogus),
    stringsAsFactors = FALSE
  )

  for (j in 1:n_ogus) {
    y <- outcome[, j]
    df <- data.frame(y = y, predictors) %>%
      arrange(subject, time)

    # Run models with convergence checking
    model_results <- run_single_gee(df, ogu_ids[j], link_family)

    # Store all model results
    all_model_results[[ogu_ids[j]]] <- model_results

    # Store convergence information
    convergence_summary$best_correlation[j] <- model_results$best_structure

    if (!is.null(model_results$best_structure)) {
      conv_info <- model_results$convergence_info[[
        model_results$best_structure
      ]]

      convergence_summary$converged[j] <- conv_info$converged
      convergence_summary$QIC[j] <- conv_info$QIC
      convergence_summary$CIC[j] <- conv_info$CIC

      # Store coefficient p-values for converged model
      if (conv_info$converged) {
        results[j, ] <- conv_info$coef_summary[1:5, 4]
      }
    }
  }

  # Create results dataframe with adjusted p-values
  adj_results <- data.frame(
    ogu_id = ogu_ids,
    intercept = p.adjust(results[, "intercept"], method = "BH"),
    time = p.adjust(results[, "time"], method = "BH"),
    arm = p.adjust(results[, "arm"], method = "BH"),
    age_B = p.adjust(results[, "age_B"], method = "BH"),
    interaction = p.adjust(results[, "interaction"], method = "BH")
  )

  return(list(
    adj_results = adj_results,
    convergence_summary = convergence_summary,
    all_model_results = all_model_results
  ))
}

#' Perform comprehensive GEE diagnostics
#' @param model Fitted GEE model
#' @param data Original data used to fit the model
#' @param id_var Name of the ID variable
#' @param time_var Name of the time variable (if applicable)
#' @return List containing diagnostic results
gee_diagnostics <- function(model, data, id_var, time_var = NULL) {
  # 1. Basic convergence check
  convergence <- list(
    converged = !is.null(model$converged) && model$converged,
    iterations = model$geese$iterations
  )

  # 2. Coefficients and their stability
  coef_summary <- summary(model)$coefficients

  # 3. Residual diagnostics
  residuals <- list(
    pearson = residuals(model, type = "pearson"),
    deviance = residuals(model, type = "deviance"),
    mahalanobis = residuals(model, type = "mahalanobis")
  )

  # 4. Correlation structure assessment
  working_correlation <- model$geese$working.correlation

  # 5. Calculate various diagnostic measures
  diagnostics <- list(
    QIC = tryCatch(glmtoolbox::QIC(model)$QIC, error = function(e) NA),
    CIC = tryCatch(glmtoolbox::CIC(model)$CIC, error = function(e) NA)
  )

  # 6. Check for influential observations
  influence <- tryCatch(
    {
      cook_d <- geepack::cook.distance(model)
      list(
        cook_d = cook_d,
        influential_obs = which(cook_d > 4 / nrow(data))
      )
    },
    error = function(e) NULL
  )

  # 7. Temporal correlation check (if time variable provided)
  temporal_correlation <- if (!is.null(time_var)) {
    tryCatch(
      {
        resid_by_time <- split(residuals$pearson, data[[time_var]])
        acf_result <- acf(residuals$pearson, plot = FALSE)
        list(
          acf = acf_result$acf,
          acf_pvalues = acf_result$acf > 0.2 # Flag potentially problematic autocorrelations
        )
      },
      error = function(e) NULL
    )
  } else {
    NULL
  }

  # 8. Create diagnostic plots data
  plot_data <- list(
    fitted_vs_residuals = data.frame(
      fitted = fitted(model),
      residuals = residuals$pearson
    ),
    qq_data = qqnorm(residuals$pearson, plot = FALSE)
  )

  # 9. Summary statistics
  summary_stats <- list(
    sample_size = nrow(data),
    n_subjects = length(unique(data[[id_var]])),
    observations_per_subject = table(data[[id_var]]),
    missing_pattern = table(is.na(fitted(model)))
  )

  # Return all diagnostics
  return(list(
    convergence = convergence,
    coefficients = coef_summary,
    residuals = residuals,
    working_correlation = working_correlation,
    diagnostics = diagnostics,
    influence = influence,
    temporal_correlation = temporal_correlation,
    plot_data = plot_data,
    summary_stats = summary_stats
  ))
}

#' Create diagnostic plots for GEE model
#' @param diagnostic_results Results from gee_diagnostics function
#' @return List of ggplot objects
plot_gee_diagnostics <- function(
  diagnostic_results,
  plot_theme = theme_minimal()
) {
  require(ggplot2)
  if (is.null(diagnostic_results)) {
    stop("Diagnostic results cannot be NULL")
  }

  plot_list <- list()

  tryCatch(
    {
      # 1. Residuals vs Fitted
      plot_list$residuals_vs_fitted <- ggplot(
        diagnostic_results$plot_data$fitted_vs_residuals,
        aes(x = V1, y = pearson)
      ) +
        geom_point(alpha = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_smooth(method = "loess", se = FALSE, color = "blue") +
        labs(
          title = "Residuals vs Fitted Values",
          x = "Fitted values",
          y = "Pearson residuals"
        ) +
        plot_theme

      # 2. Q-Q plot
      plot_list$qq_plot <- data.frame(
        x = diagnostic_results$plot_data$qq_data$x,
        pearson = diagnostic_results$plot_data$qq_data$y
      ) %>%
        ggplot(aes(x = x, y = pearson)) +
        geom_point(alpha = 0.5) +
        geom_abline(
          intercept = 0,
          slope = 1,
          color = "red",
          linetype = "dashed"
        ) +
        labs(
          title = "Normal Q-Q Plot",
          x = "Theoretical Quantiles",
          y = "Sample Quantiles"
        ) +
        plot_theme

      # 3. Residuals histogram
      plot_list$residuals_hist <- ggplot(
        data.frame(pearson = diagnostic_results$residuals$pearson),
        aes(x = pearson)
      ) +
        geom_histogram(bins = 30, fill = "lightblue", color = "black") +
        labs(
          title = "Histogram of Residuals",
          x = "Pearson Residuals",
          y = "Count"
        ) +
        plot_theme

      # 4. Cook's distance plot (if available)
      if (!is.null(diagnostic_results$influence)) {
        plot_list$cooks_distance <- ggplot(
          data.frame(
            index = seq_along(diagnostic_results$influence$cook_d),
            cook_d = diagnostic_results$influence$cook_d
          ),
          aes(x = index, y = cook_d)
        ) +
          geom_point() +
          geom_hline(
            yintercept = diagnostic_results$influence$threshold,
            color = "red",
            linetype = "dashed"
          ) +
          labs(
            title = "Cook's Distance",
            x = "Observation Index",
            y = "Cook's distance"
          ) +
          plot_theme
      }

      # 5. ACF plot (if available)
      if (!is.null(diagnostic_results$temporal_correlation)) {
        plot_list$acf_plot <- ggplot(
          data.frame(
            lag = seq_along(diagnostic_results$temporal_correlation$acf) - 1,
            acf = diagnostic_results$temporal_correlation$acf
          ),
          aes(x = lag, y = acf)
        ) +
          geom_col() +
          geom_hline(
            yintercept = c(0.2, -0.2),
            color = "red",
            linetype = "dashed"
          ) +
          labs(title = "Autocorrelation Function", x = "Lag", y = "ACF") +
          plot_theme
      }

      plot <- ggarrange(plotlist = plot_list, labels = "AUTO")

      return(plot)
    },
    error = function(e) {
      stop("Error in plot_gee_diagnostics: ", conditionMessage(e))
    }
  )
}

#' Print summary of GEE diagnostics
#' @param diagnostic_results Results from gee_diagnostics function
#' @return Printed summary of diagnostics
print_gee_diagnostics <- function(diagnostic_results) {
  cat("GEE Model Diagnostics Summary\n")
  cat("============================\n\n")

  # Convergence
  cat("1. Convergence:\n")
  cat("   - Model converged:", diagnostic_results$convergence$converged, "\n")
  cat(
    "   - Number of iterations:",
    diagnostic_results$convergence$iterations,
    "\n\n"
  )

  # Sample size information
  cat("2. Sample Size Information:\n")
  cat(
    "   - Total observations:",
    diagnostic_results$summary_stats$sample_size,
    "\n"
  )
  cat(
    "   - Number of subjects:",
    diagnostic_results$summary_stats$n_subjects,
    "\n\n"
  )

  # Model fit statistics
  cat("3. Model Fit Statistics:\n")
  cat("   - QIC:", diagnostic_results$diagnostics$QIC, "\n")
  cat("   - CIC:", diagnostic_results$diagnostics$CIC, "\n\n")

  # Influential observations
  if (!is.null(diagnostic_results$influence)) {
    cat("4. Influential Observations:\n")
    cat(
      "   - Number of influential points:",
      length(diagnostic_results$influence$influential_obs),
      "\n\n"
    )
  }

  # Correlation structure
  cat("5. Working Correlation Structure:\n")
  print(diagnostic_results$working_correlation)
  cat("\n")

  # Temporal correlation (if available)
  if (!is.null(diagnostic_results$temporal_correlation)) {
    cat("6. Temporal Correlation:\n")
    cat(
      "   - Significant autocorrelation present:",
      any(diagnostic_results$temporal_correlation$acf_pvalues),
      "\n\n"
    )
  }

  # Coefficient summary
  cat("7. Coefficient Summary:\n")
  print(diagnostic_results$coefficients)
}

#' Identify significant OGUs and merge with taxonomy data
#' @param gee_results GEE model results
#' @param taxa_data Taxonomy data frame
#' @param significance_threshold P-value threshold for significance
#' @return Data frame of significant OGUs with taxonomy information
identify_significant_ogus <- function(
  gee_results,
  taxa_data,
  significance_threshold = 0.05,
  include_intercept = FALSE
) {
  if(include_intercept == FALSE){
    significant_ogus <- left_join(
      rownames_to_column(taxa_data, var = "rowname"),
      gee_results,
      by = c("rowname" = "outcome_id")
    ) %>%
      filter(if_any(
        !`(Intercept)`,
        ~ . <= significance_threshold
      ))
  } else {
        significant_ogus <- left_join(
      rownames_to_column(taxa_data, var = "rowname"),
      gee_results,
      by = c("rowname" = "outcome_id")
    ) %>%
      filter(if_any(
        everything(), ~ . <= significance_threshold
      ))
  }

  return(significant_ogus)
}

#' Run complete GEE analysis pipeline with convergence diagnostics
#' @param phyloseq_obj Input phyloseq object
#' @param output_dir Directory for saving results
#' @return List containing all analysis results including predictor variables
run_gee_analysis <- function(
  phyloseq_obj,
  output_dir,
  link_family = gaussian(link = "identity")
) {
  # Prepare data
  prepared_data <- prepare_ogu_data(phyloseq_obj)
  predictors <- prepare_predictors(prepared_data$meta_data)

  # Store original data with sample names
  predictor_data <- predictors %>%
    rownames_to_column(var = "sample_names") %>%
    as.data.frame()

  # Run GEE models with convergence diagnostics
  model_results <- run_gee_models(
    prepared_data$ogu_data,
    predictors,
    colnames(prepared_data$ogu_data),
    link_family = gaussian(link = "identity")
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
      conv_info <- x$convergence_info[[x$best_structure]]

      list(
        ogu_id = x$ogu_id,
        model = best_model,
        correlation_structure = x$best_structure,
        convergence = conv_info$converged,
        QIC = conv_info$QIC,
        CIC = conv_info$CIC,
        coefficients = conv_info$coef_summary,
        diagnostics = gee_diagnostics(
          best_model,
          data = best_model$data,
          id_var = "subject",
          time_var = "time"
        )
      )
    }
  })

  # Save results
  write.csv(
    model_results$adj_results,
    file.path(output_dir, "OGU-GEE-adj-p-values.csv"),
    row.names = FALSE
  )

  write.csv(
    model_results$convergence_summary,
    file.path(output_dir, "OGU-GEE-convergence-summary.csv"),
    row.names = FALSE
  )

  write.csv(
    significant_ogus,
    file.path(output_dir, "significant-OGU-GEE-adj-p-values.csv"),
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

  # Return comprehensive results including predictor data
  return(list(
    gee_results = model_results$adj_results,
    convergence_summary = model_results$convergence_summary,
    significant_ogus = significant_ogus,
    significant_ogu_models = significant_ogu_details,
    prepared_data = prepared_data,
    predictor_data = predictor_data, # Added predictor data to results
    analysis_info = list(
      date = Sys.Date(),
      n_samples = nrow(predictor_data),
      n_ogus = ncol(prepared_data$ogu_data),
      n_significant = nrow(significant_ogus)
    )
  ))
}

#' Compare GEE models with different correlation structures
#' @param formula Model formula
#' @param data Data frame containing all variables
#' @param id_var Name of the ID variable for repeated measures
#' @param family Distribution family (default: gaussian)
#' @param correlation_structures Vector of correlation structures to try
#' @return List containing model comparison results
compare_gee_structures <- function(
  formula,
  data,
  id_var,
  family = gaussian,
  correlation_structures = c(
    "independence",
    "exchangeable",
    "unstructured",
    "ar1"
  )
) {
  # Initialize lists to store results
  models <- list()
  fit_stats <- list()

  # Try each correlation structure
  for (corstr in correlation_structures) {
    tryCatch(
      {
        # Fit model
        model <- geepack::geeglm(
          formula = formula,
          family = family,
          data = data,
          id = as.formula(paste0("~", id_var)),
          corstr = corstr
        )

        # Store model
        models[[corstr]] <- model

        # Calculate fit statistics
        fit_stats[[corstr]] <- list(
          QIC = tryCatch(glmtoolbox::QIC(model)$QIC, error = function(e) NA),
          CIC = tryCatch(glmtoolbox::CIC(model)$CIC, error = function(e) NA),
          converged = !is.null(model$converged) && model$converged
        )
      },
      error = function(e) {
        fit_stats[[corstr]] <- list(
          QIC = NA,
          CIC = NA,
          converged = FALSE,
          error = conditionMessage(e)
        )
      }
    )
  }

  # Create summary dataframe
  summary_df <- do.call(
    rbind,
    lapply(names(fit_stats), function(corstr) {
      data.frame(
        correlation = corstr,
        QIC = fit_stats[[corstr]]$QIC,
        CIC = fit_stats[[corstr]]$CIC,
        converged = fit_stats[[corstr]]$converged,
        error = if (!is.null(fit_stats[[corstr]]$error)) {
          fit_stats[[corstr]]$error
        } else {
          NA
        }
      )
    })
  )

  return(list(
    models = models,
    summary = summary_df
  ))
}

#' Compare GEE models with different correlation structures
#' @param formula Model formula
#' @param data Data frame containing all variables
#' @param id_var Name of the ID variable for repeated measures
#' @param family Distribution family (default: gaussian)
#' @param correlation_structures Vector of correlation structures to try
#' @return List containing model comparison results
compare_gee_structures.glmtoolbox <- function(
  formula,
  data,
  id_var = "subject",
  family = gaussian,
  correlation_structures = c(
    "Exchangeable",
    "Independence",
    "Unstructured",
    "Stationary-M-dependent",
    "Stationary-M-dependent(1)"
  )
) {
  # Initialize lists to store results
  models <- list()
  fit_stats <- list()

  # Try each correlation structure
  for (corstr in correlation_structures) {
    tryCatch(
      {
        # Fit model
        model <- glmtoolbox::glmgee(
          formula = formula,
          family = family,
          data = data,
          id = id_var,
          corstr = corstr,
          verbose = FALSE
        )

        # Store model
        models[[corstr]] <- model

        # Calculate fit statistics
        fit_stats[[corstr]] <- list(
          converged = !is.null(model$converged) && model$converged,
          QIC = tryCatch(glmtoolbox::QIC(model)$QIC, error = function(e) NA),
          CIC = tryCatch(glmtoolbox::CIC(model)$CIC, error = function(e) NA),
          coef_summary = summary(model)$coefficients,
          warnings = model$messages
        )
      },
      error = function(e) {
        fit_stats[[corstr]] <- list(
          QIC = NA,
          CIC = NA,
          converged = FALSE,
          error = conditionMessage(e)
        )
      }
    )
  }

  # Create summary dataframe
  summary_df <- do.call(
    rbind,
    lapply(names(fit_stats), function(corstr) {
      data.frame(
        correlation = corstr,
        QIC = fit_stats[[corstr]]$QIC,
        CIC = fit_stats[[corstr]]$CIC,
        converged = fit_stats[[corstr]]$converged,
        error = if (!is.null(fit_stats[[corstr]]$error)) {
          fit_stats[[corstr]]$error
        } else {
          NA
        }
      )
    })
  )

  # Select best model based on QIC
  qic_values <- sapply(fit_stats, function(x) x$QIC)
  best_structure <- names(which.min(qic_values))

  # If no QIC available, use first converged model
  if (all(is.na(qic_values))) {
    converged_models <- names(which(sapply(fit_stats, function(x) x$converged)))
    best_structure <- if (length(converged_models) > 0) {
      converged_models[1]
    } else {
      NULL
    }
  }

  return(list(
    models = models,
    summary = summary_df,
    best_structure = best_structure,
    best_model = if (!is.null(best_structure)) {
      models[[best_structure]]
    } else {
      NULL
    }
  ))
}

# start gee-plots
#' Helper function to extract data for a specific OGU
#' @param results Results object from run_gee_analysis
#' @param ogu_id OGU identifier
#' @return List containing OGU-specific data and model results
extract_ogu_data <- function(results, ogu_id) {
  # Get OGU abundances
  ogu_abundance <- results$prepared_data$ogu_data[, ogu_id]

  # Combine with predictor data
  full_data <- results$predictor_data %>%
    mutate(abundance = ogu_abundance)

  # Get model results if significant
  model_results <- results$significant_ogu_models[[ogu_id]]

  # Get taxonomy information
  taxonomy <- results$significant_ogus %>%
    filter(rowname == ogu_id)

  return(list(
    data = full_data,
    model = model_results,
    taxonomy = taxonomy
  ))
}

#' Create summary plot for OGU
#' @param ogu_data Result from extract_ogu_data
#' @return ggplot object
plot_ogu_summary <- function(ogu_data) {
  require(ggplot2)

  # Create abundance plot
  p <- ggplot(
    ogu_data$data,
    aes(x = time, y = abundance, color = arm, group = subject)
  ) +
    geom_line(alpha = 0.3) +
    geom_smooth(aes(group = arm), method = "loess", se = TRUE) +
    theme_minimal() +
    labs(
      title = paste("OGU:", ogu_data$taxonomy$Species),
      subtitle = paste("Phylum:", ogu_data$taxonomy$Phylum),
      x = "Time",
      y = "Abundance"
    )

  return(p)
}
#end

# Analyze a specific significant OGU
analyze_significant_ogu <- function(ogu_details) {
  # Print basic information
  cat("OGU ID:", ogu_details$ogu_id, "\n")
  cat("Correlation structure:", ogu_details$correlation_structure, "\n")
  cat("Converged:", ogu_details$convergence, "\n")
  cat("QIC:", ogu_details$QIC, "\n\n")

  # Print coefficient summary
  cat("Coefficient Summary:\n")
  print(ogu_details$coefficients)
  cat("\n")

  # Create diagnostic plots
  plots <- plot_gee_diagnostics(ogu_details$diagnostics)

  return(plots)
}
