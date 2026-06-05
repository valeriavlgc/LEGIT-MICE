library(LEGIT)
library(ggplot2)
library(writexl)


# ------------------------- #
#  FUNCTION: run_analysis   #
#         IMPUTATED         # 
# ------------------------- #

run_analysis_imputed <- function(gene, e, out, dt_imp, i, analysis_type, multilevel_id = NULL, verbose = TRUE) {

  fits <- vector("list", nimp)
  fit_df_list <- vector("list", nimp)
  crossover_int <- ""
  crossover_int2 <- ""
  tie <- logical()
  
  config <- get_analysis_config(analysis_type, out, multilevel_id)
  params <- list(genes = gene, env = e)
  
  # Write name of the interaction
  #name_combination <- paste(gene, e, out, sep = "_")
  
  name_combination <- paste(
    paste(gene, collapse = "_"),
    paste(e, collapse = "_"),
    out, sep = "_"
  )

  tryCatch({
    for (j in 1: nimp) { 
      data_complete <- complete(dt_imp, action = j)
      
      env <- data_complete[, params$env, drop = FALSE]
      genes <- data_complete[, params$genes, drop = FALSE]
      
      fit <- LEGIT(data = data_complete, genes = genes, env = env, 
                   formula = config$formula_dynamic, lme4 = config$use_lme4, rescale = config$use_lme4) 
      
      fits[[j]] <- fit  # Store the model
    }
    
    
    # Pooling fits
    for (x in 1:nimp) {
      fit_df <- fits[[x]]
      
      fixed_effects <- summary(fit_df$fit_main)$coefficients
      
      # Create a data frame with the necessary columns
      fit_df_imputed <- data.frame(
        term = rownames(fixed_effects),  
        estimate = fixed_effects[, "Estimate"],
        std.error = fixed_effects[, "Std. Error"], 
        residual.df = rep(NA, nrow(fixed_effects)),
        row.names = NULL
    )
      
      fit_df_list[[x]] <- fit_df_imputed
    }
    
    w <- do.call(rbind, fit_df_list)
    rownames(w) <- NULL

    # Now use the pool.table() function to pool the results
    pooled_fits <- pool.table(
      w, 
      type = "all",  # Default to all estimates
      conf.int = TRUE,  # Include confidence intervals
      conf.level = 0.95,  # 95% confidence level
      exponentiate = FALSE,  # Don't exponentiate unless it's logistic regression
      #dfcom = nrow(data_complete) - 6,  
      dfcom = nrow(data_complete) - length(coef(fits[[1]]$fit_main)),
      rule = "rubin1987"  # Default rule for pooling
    )

        # Calculate Chi-squared statistics and p-values for pooled fixed effects
    anova_results <- pooled_fits %>%
      dplyr::mutate(
        chisq = (estimate / std.error)^2,  # Chi-squared statistic
        p.value = pchisq(chisq, df = 1, lower.tail = FALSE)  # P-value from Chi-squared distribution
      ) %>%
      dplyr::select(term, estimate, std.error, chisq, p.value)
    
    rownames(anova_results) <- anova_results[[1]]
    
    
    if (verbose) {
      if (!dir.exists("plots")) dir.create("plots")
      # Plot uses the first imputation as representative
      png(paste0("plots/", i, "_", name_combination, ".png"), width = 800, height = 600)
      fit <- fits[[1]]
      plot(fit, ylab = fit$fit_genes$formula, xlab = fit$fit_env$formula, cex.leg = 1.2)
      dev.off()
      cat("Plot saved as", name_combination)
    }
    
    
    GxE_test_AIC_results <- lapply(1:dt_imp$m, function(k) {
      imp_data <- complete(dt_imp, k) 
      
      env <- imp_data[, params$env, drop = FALSE]
      genes <- imp_data[, params$genes, drop = FALSE]
      
      # Do i want to keep the supress warning?
      suppressWarnings(
        GxE_interaction_test(
          data = imp_data,
          genes = genes,
          env = env,
          formula_noGxE = config$formula_dynamic2,
          crossover = c("min", "max"),
          criterion = "AIC",
          lme4 = config$use_lme4,
          rescale = TRUE
        ))
    })
    
    model1 <- GxE_test_AIC_results[[1]]
    models <- names(model1$fits) 
    
    aic_values <- matrix(NA, nrow = 10, ncol = nimp)
    rownames(aic_values) <- models
    crossover <- matrix(NA, nrow = 2, ncol = nimp)
    
    for (l in 1:nimp) {
      aic <- sapply(models, function(model_name) {
        fit_result <- GxE_test_AIC_results[[l]]$fits[[model_name]]
        
        if(grepl("model", model_name)){
          model_aic <- fit_result$aic
        } else {
          model_aic <- fit_result$true_model_parameters$AIC
        }
        
        return(model_aic)
      })
      
      # crossover DS model weak 
      crossover[1, l] <- GxE_test_AIC_results[[l]]$results["Differential susceptibility WEAK", 3]
      # crossover DS model strong 
      crossover[2, l] <- GxE_test_AIC_results[[l]]$results["Differential susceptibility STRONG", 3]
      
      aic_values[, l] <- aic
    }
    
    pooled_aic <- apply(aic_values, 1, mean, na.rm = TRUE)
    sorted_aic <- sort(pooled_aic)
    
    mod1 <- revert_model(names(sorted_aic)[1], e)
    mod2 <- revert_model(names(sorted_aic)[2], e)
    
    # This could be improved check other script _copia
    
    tie_model <- FALSE
    
    if (abs(sorted_aic[1] - sorted_aic[2]) <= tolerance) {
      tie_model <- TRUE
    }
    
    tie <- c(tie, tie_model)
    
    # Also might be improved, check script.
    if (names(sorted_aic)[1] == "diff_suscept_WEAK" || (names(sorted_aic)[1] == "diff_suscept_STRONG")) {
      
      res <- calculate_crossover(crossover)
      mean_min <- res$mean_min
      mean_max <- res$mean_max
      
      if (names(sorted_aic)[1] == "diff_suscept_WEAK") {
        crossover_int <- paste0(crossover_int, paste0(mean_min[1], "/", mean_max[1]))
      } else if (names(sorted_aic)[1] == "diff_suscept_STRONG") {
        crossover_int <- paste0(crossover_int, paste0(mean_min[2], "/", mean_max[2]))
      }
      
      
    } else {
      crossover_int <- paste0(crossover_int, "")  
    }
    
    if (names(sorted_aic)[2] == "diff_suscept_WEAK" || names(sorted_aic)[2] == "diff_suscept_STRONG") {
      
      res <- calculate_crossover(crossover)
      mean_min <- res$mean_min
      mean_max <- res$mean_max
      
      if (names(sorted_aic)[2] == "diff_suscept_WEAK") {
        crossover_int2 <- paste0(crossover_int2, paste0(mean_min[1], "/", mean_max[1]))
      } else if (names(sorted_aic)[2] == "diff_suscept_STRONG") {
        crossover_int2 <- paste0(crossover_int2, paste0(mean_min[2], "/", mean_max[2]))
      }
      
    } else {
      crossover_int2 <- paste0(crossover_int2, "")  
    }
    
    
    return(list(
      PRSxEnvxOut = name_combination,
      P_value_GxE = anova_results["G:E", "p.value"],
      P_value_G = anova_results["G", "p.value"],
      P_value_C1 = anova_results["C1", "p.value"],
      P_value_C2 = anova_results["C2", "p.value"],
      F_ratio = anova_results["G:E", "chisq"],
      GxE_estimate = anova_results["G:E", "estimate"],  
      Type_of_interaction = mod1,
      AIC = sorted_aic[1],
      crossover_int = crossover_int,
      tie = tie_model,
      Type_of_interaction2 = mod2,
      AIC2 = sorted_aic[2],
      crossover_int2 = crossover_int2
    ))
    
    
  } , warning = function(w) {
    cat("Warning in GxE_interaction_test for", name_combination, ":", w$message, "\n")
    append_na_results()
    
  }, error = function(e) {
    cat("Error in GxE_interaction_test for", name_combination, ":", e$message, "\n")
    append_na_results()
  })
  
}

run_analysis <- function(gene, e, out, data, i, analysis_type, multilevel_id = NULL, verbose = TRUE) {

  config <- get_analysis_config(analysis_type, out, multilevel_id)
  params <- list(genes = gene, env = e)

  genes <- data[, params$genes, drop = FALSE]
  env <- data[, params$env, drop = FALSE] 
  
  # Write name of the interaction
  #name_combination <- paste(gene, e, out, sep = "_")
  
  name_combination <- paste(
    paste(gene, collapse = "_"),
    paste(e, collapse = "_"),
    out, sep = "_"
  )

  # GxE model
  fit = LEGIT(data, genes, env, formula = config$formula_dynamic, 
              lme4 = config$use_lme4, rescale = config$use_lme4)
  
  if (verbose) {
    if (!dir.exists("plots")) dir.create("plots")
    tryCatch({
      png(paste0("plots/", i, "_", name_combination, ".png"), width = 800, height = 600)
      plot(fit, cex.leg = 1.2)
      dev.off()
      cat("Plot saved as", name_combination)
    }, error = function(e) {
      dev.off()
      cat("Plot failed for", name_combination, ":", e$message, "\n")
    })
  }
  
  if(analysis_type == ANALYSIS_TYPE$UNILEVEL) {
    fit_sum <- summary(fit)
    coefficients <- fit_sum$fit_main$coefficients
    # Use Wald t² (partial F) for F_ratio — consistent with the imputed path
    # which also uses (estimate/std.error)^2. The previous anova(fit_glm, test="F")
    # used sequential (Type I) SS which inflates F relative to the partial test.
    pval_col <- "Pr(>|t|)"
    fval_col <- "t value"  # we square this below
  } else if (analysis_type == ANALYSIS_TYPE$MULTILEVEL) {
    coefficients <- car::Anova(fit$fit_main)
    pval_col <- "Pr(>Chisq)"
    fval_col <- "Chisq"
  }

  # Testing the type of interaction
  GxE_test_AIC <- GxE_interaction_test(data = data, genes = genes, env = env, 
                                       formula_noGxE = config$formula_dynamic2, 
                                       crossover = c("min", "max"), 
                                       criterion = "AIC", 
                                       lme4 = config$use_lme4, rescale = TRUE)
  
  
  mod1 <- revert_model(names(GxE_test_AIC$fits)[1], e)

    # could be somewhere else
    tie_model <- FALSE
  
  if (abs(as.numeric(GxE_test_AIC$results[1]) - as.numeric(GxE_test_AIC$results[2])) <= tolerance) {
    tie_model <- TRUE
  }

  mod2 <- revert_model(names(GxE_test_AIC$fits)[2], e)

  if (mod2 == "diff_suscept_WEAK"|| mod2== "diff_suscept_STRONG") {
    crossover2 <- GxE_test_AIC$results[2, 3]
  } else {
    crossover2 <- "NONE"
  }
  
  
  return(list(
    PRSxEnvxOut = name_combination,
    P_value_GxE = coefficients["G:E", pval_col],
    P_value_G = coefficients["G", pval_col],
    P_value_C1 = coefficients["C1", pval_col],
    P_value_C2 = coefficients["C2", pval_col],
    F_ratio = if (analysis_type == ANALYSIS_TYPE$UNILEVEL)
                coefficients["G:E", fval_col]^2   # t² = Wald partial F
              else
                coefficients["G:E", fval_col],     # chi-squared for multilevel
    GxE_estimate = fit_sum$fit_main$coefficients["G:E", "Estimate"],
    Type_of_interaction = mod1,
    AIC = as.numeric(GxE_test_AIC$results[1]),
    crossover_int = GxE_test_AIC$results[1, 3],
    tie = tie_model,
    Type_of_interaction2 = mod2,
    AIC2 = as.numeric(GxE_test_AIC$results[2]),
    crossover_int2 = crossover2
  ))
}

run_analysis_LEGIT <- function(genes_values, env_values, out_values, dt,
                               analysis_type, imputated, multilevel_id = NULL,
                               save_results = TRUE, verbose = TRUE) {
  
  # if(analysis_type == ANALYSIS_TYPE$MULTILEVEL && dt) change ID twin thing???
  i <- 0
  results_list <- list()
  
  cat("Running", 
      if (imputated) "IMPUTED" else "STANDARD", 
      if (analysis_type == ANALYSIS_TYPE$MULTILEVEL) "MULTILEVEL" else "UNILEVEL", 
      "LEGIT analysis...\n")  
  
  
    if(imputated){
      for (out in out_values) {
        i <- i + 1
        result <- run_analysis_imputed(genes_values, env_values, out, dt, i, analysis_type, multilevel_id, verbose = verbose)
        results_list[[i]] <- result
      }
      res <- "/results_MI_"
    } else {
      for (out in out_values) {
        i <- i + 1
        result <- run_analysis(genes_values, env_values, out, dt, i, analysis_type, multilevel_id, verbose = verbose)
        results_list[[i]] <- result
      }
      res <- "/results_"
    }

# TO PASS 1 GENE AT A TIME
  
    # if(imputated){
    #   for (gene in genes_values) {
    #     for (e in env_values) {
    #       for (out in out_values) {
    #         i <- i + 1
    #         result <- run_analysis_imputed(gene, e, out, dt, i, analysis_type, multilevel_id, verbose = verbose)
    #         results_list[[i]] <- result
    #       }
    #     }
    #   }
    #    res <- "/results_MI_"
    # } else {
    #   for (gene in genes_values) {
    #     for (e in env_values) {
    #       for (out in out_values) {
    #         i <- i + 1
    #         result <- run_analysis(gene, e, out, dt, i, analysis_type, multilevel_id, verbose = verbose)
    #         results_list[[i]] <- result
    #       }
    #     }
    #   }
    #   res <- "/results_"
    # }
      
    # Convert to data.frame
    results <- do.call(rbind, lapply(results_list, as.data.frame))
    
    if (save_results) {
      if (!dir.exists("results")) dir.create("results")
        write_xlsx(results, paste0("results", res , analysis_type, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"))
        cat("\nFinished", 
            if (imputated) "IMPUTED" else "STANDARD", 
            if (analysis_type == ANALYSIS_TYPE$MULTILEVEL) "MULTILEVEL" else "UNILEVEL", 
            "LEGIT analysis. Results saved in /results\n")
        
        print(results)  
    }
    
    return(results)
}

