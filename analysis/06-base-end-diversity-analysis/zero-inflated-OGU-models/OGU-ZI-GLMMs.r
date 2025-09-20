# Function to run EMM analysis
basic_emm_analysis <- function(fits, taxa_data, results_dir, emm_adjust = "sidak") {
  
  # Combine all models into single list
  # fits <- list_c(final_models)
  
  # if (length(fits) == 0) {
  #   warning("No models available for EMM analysis")
  #   return(NULL)
  # }

  # Differnece in EMMs between timepoints (no grouping factor)
  emm_time <- lapply(fits, function(x) {
    emm <- emmeans(
      x,
      ~time,
      at = list(time = c(0, 1)),
      adjust = emm_adjust,
      weights = "proportional",
      type = "scale"
    )
    levels(emm)$time <- c("Baseline", "Endpoint")
    return(emm)
  })
  
  # Time contrasts
  time_cont <- lapply(emm_time, function(x){
      contrast(
        x,
        by = c("time"),
        method = "revpairwise",
        type = "scale",
        adjust = emm_adjust,
        weights = "proportional",
        parens = FALSE
      )

  })
  
  time_cont_df <- emm_bind(time_cont, kable = F) %>%
    left_join(rownames_to_column(taxa_data, "Model"), by = "Model") %>%
    mutate(Species = ifelse(Species == "", Genus, Species))
  
  sigTime <- time_cont_df %>%
    filter(p.value <= 0.05)
  
  # Effect sizes
  emmEff <- lapply(1:length(time_cont), function(i) {
    eff_size(
      time_cont[[i]],
      sigma = attr(VarCorr(fits[[i]])$cond$subject, "stddev") * sqrt(2),
      edf = df.residual(fits[[i]]),
      method = "identity",
      parens = FALSE
    )
  })
  
  names(emmEff) <- names(time_cont)
  
  emmEff_df <- emm_bind(emmEff, kable = F)
  
  emmEff_df <- emmEff_df %>%
    filter(Model %in% sigTime$Model) %>%
    mutate(
      magnitude = case_when(
        abs(effect.size) < 0.2 ~ "negligible",
        abs(effect.size) >= 0.2 & abs(effect.size) < 0.5 ~ "small",
        abs(effect.size) >= 0.5 & abs(effect.size) < 0.8 ~ "medium",
        abs(effect.size) >= 0.8 ~ "large",
      )
    ) %>%
    arrange(abs(effect.size)) %>%
    select(Model, effect.size, magnitude)
  
  sigTimeEff <- left_join(sigTime, emmEff_df, by = c("Model"))
  
  # Save results
  write.csv(
    sigTimeEff,
    file.path(results_dir, "significant-ogu-emms-time-ZI.csv"),
    row.names = FALSE
  )
  
  write.csv(
    time_cont_df,
    file.path(results_dir, "ogu-emms-time-ZI.csv"),
    row.names = FALSE
  )
  
  return(list(
    sigTimeEff = sigTimeEff,
    time_cont_df = time_cont_df,
    fits = fits
  ))
}

# Complete analysis wrapper function
run_complete_zi_analysis <- function(phyloseq_obj,
                                   analysis_name = "ZI-analysis",
                                   filter_prevalence = TRUE,
                                   prevalence_cutoff = 0.4,
                                   filter_abundance = FALSE,
                                   formula_str = "time + age_B + age_B:time + (1 | subject)",
                                   zi_formula = "~ age_B * time",
                                   fit_new_models = FALSE) {
  
  cat("Starting ZI analysis for:", analysis_name, "\n")
  
  # Step 1: Prepare data
  cat("Preparing data...\n")
  prepared_data <- run_zi_analysis(
    phyloseq_obj = phyloseq_obj,
    analysis_name = analysis_name,
    filter_prevalence = filter_prevalence,
    prevalence_cutoff = prevalence_cutoff,
    filter_abundance = filter_abundance,
    formula_str = formula_str,
    zi_formula = zi_formula
  )
  
  # Step 2: Fit models
  cat("Fitting ZI models...\n")
  results_list <- fit_zi_distributions(
    outcome = prepared_data$outcome,
    predictors = prepared_data$predictors,
    formula_str = formula_str,
    zi_formula = zi_formula,
    results_dir = prepared_data$results_dir,
    fit_models = fit_new_models
  )
  
  # Step 3: Select best models
  cat("Selecting best fitting models...\n")
  model_selection <- select_best_models(results_list)
  
  # Step 4: Extract final models
  cat("Extracting final models...\n")
  final_models <- extract_final_models(results_list, model_selection$bestDist_filtered)
  
  # Step 5: Run EMM analysis
  cat("Running EMM analysis...\n")
  emm_results <- run_emm_analysis(
    final_models = final_models,
    taxa_data = prepared_data$taxa_data,
    results_dir = prepared_data$results_dir
  )
  
  cat("Analysis complete!\n")
  
  return(list(
    prepared_data = prepared_data,
    results_list = results_list,
    model_selection = model_selection,
    final_models = final_models,
    emm_results = emm_results
  ))
}
