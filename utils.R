append_na_results <- function() {
  return(list(
    PRSxEnvxOut = "",
    P_value_GxE = "",
    P_value_G = "",
    P_value_C1 = "",
    P_value_C2 = "",
    F_ratio = "",
    Type_of_interaction = "",
    AIC = "",
    crossover_int = "",
    tie = "",
    Type_of_interaction2 = "",
    AIC2 = "",
    crossover_int2 = ""
  ))
}

# # Split and extract values from crossover for DS models
calculate_crossover <- function(crossover) {
  for (row in 1:2) {
    for (col in 1:nimp) {
      crossover_value <- as.character(crossover[row, col])
      
      if (!is.na(crossover_value) && grepl("\\(.*\\)", crossover_value)) {
        clean_value <- gsub("[()]", "", crossover_value)
        clean_value <- trimws(clean_value)
        
        values <- as.numeric(strsplit(clean_value, " / ")[[1]])
        min_values[row, col] <- values[1]
        max_values[row, col] <- values[2]
      } else {
        min_values[row, col] <- NA
        max_values[row, col] <- NA
      }
    }
  }
  
  # crossover for DS weak mean_min[1] + max[1]
  mean_min <- rowMeans(min_values, na.rm = TRUE)
  mean_max <- rowMeans(max_values, na.rm = TRUE)
  
  return(list(mean_min = mean_min, mean_max = mean_max))
}


# Revert model diathesis/vantage when environment is negative
revert_model <- function(mod, e) {
  if(e %in% neg_env & mod %in% rev_models) {
    
    if(mod == "diathesis_stress_WEAK") 
      mod_revert <- "vantage_sensitivity_WEAK"
    if(mod == "diathesis_stress_STRONG")
      mod_revert <- "vantage_sensitivity_STRONG"
    if(mod == "vantage_sensitivity_WEAK")
      mod_revert <- "diathesis_stress_WEAK"
    if(mod == "vantage_sensitivity_STRONG")
      mod_revert <- "diathesis_stress_STRONG"
    
    return(mod_revert)
    
  } else {
    return(mod)
  }
}

get_analysis_config <- function(analysis_type, out) {
  
if (analysis_type == ANALYSIS_TYPE$MULTILEVEL) {
  # multilevel variable is in config
  formula_dynamic <- as.formula(paste(out, "~", interaction_term, "+", covariate1, "+", covariate2, "+", idtw))
  formula_dynamic2 <- as.formula(paste(out, "~", 1, "+", idtw))
  use_lme4 <- TRUE
  
} else if (analysis_type == ANALYSIS_TYPE$UNILEVEL) {
  formula_dynamic <- as.formula(paste(out, "~", interaction_term, "+", covariate1, "+", covariate2))
  formula_dynamic2 <- as.formula(paste(out, "~", 1))
  use_lme4 <- FALSE
  
} else {
  stop("Invalid analysis type")
}
  return(list(
    formula_dynamic = formula_dynamic,
    formula_dynamic2 = formula_dynamic2,
    use_lme4 = use_lme4
  ))
  
}

check_vars <- function(col_names, vars) {
  missing_vars <- setdiff(vars, col_names)
  
  if (length(missing_vars) > 0) {
    stop(paste("Error: Missing variables in dataset:", paste(missing_vars, collapse = ", ")))
  }
}

