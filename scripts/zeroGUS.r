compareZI <- function(
  formula,
  data,
  zi,
  family = c(
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
  ),
  plot = TRUE
) {
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
        fit_stats[[fam_name]] <- list(
          AIC = AIC(model),
          BIC = BIC(model),
          converged = !is.null(model$sdr$pdHess) && model$sdr$pdHess,
          simResiduals = sim_res,
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
        if (plot) {
          plots[[fam_name]] <- list(
            residuals = DHARMa::plot(sim_res, ask = FALSE),
            qq = qqnorm(residuals(model)),
            fitted = plot(fitted(model), residuals(model))
          )
        }
      },
      error = function(e) {
        fit_stats[[fam_name]] <- list(
          AIC = NA,
          BIC = NA,
          converged = FALSE,
          simResiduals = NA,
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
        converged = fit_stats[[fam]]$converged,
        negative_eigenvals = fit_stats[[
          fam
        ]]$hessian_check$negative_eigenvals,
        condition_number = fit_stats[[fam]]$hessian_check$condition_number,
        error = if (!is.null(fit_stats[[fam]]$error)) {
          fit_stats[[fam]]$error
        } else {
          NA
        }
      )
    })
  )

  # Add zero-inflation statistics
  zero_stats <- lapply(models, function(mod) {
    if (!is.null(mod)) {
      obs_zeros <- sum(mod$frame[, 1] == 0)
      pred_zeros <- sum(predict(mod, type = "response") < 0.0001)
      return(c(observed = obs_zeros, predicted = pred_zeros))
    } else {
      return(c(observed = NA, predicted = NA))
    }
  })

  zero_df <- do.call(rbind, zero_stats)
  summary_df$observed_zeros <- zero_df[, "observed"]
  summary_df$predicted_zeros <- zero_df[, "predicted"]

  # Print formatted summary
  cat("\nModel Comparison Summary:\n")
  print(knitr::kable(summary_df, digits = 2))

  cat(
    "\nBest model by AIC:",
    summary_df$family[which.min(summary_df$AIC)],
    "\n"
  )
  cat(
    "Best model by BIC:",
    summary_df$family[which.min(summary_df$BIC)],
    "\n"
  )

  # Return results
  result <- list(
    models = models,
    summary = summary_df,
    diagnostics = fit_stats,
    plots = if (plot) plots else NULL
  )

  return(invisible(result))
}
