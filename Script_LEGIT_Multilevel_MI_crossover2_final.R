library(LEGIT)
library(foreign)
library(writexl)
library(dplyr)
library(naniar)
library(mice)


# Add Rubin's rules.
# I tried to imputate categorical variables as such converting them to factors data[vars_to_convert] <- lapply(data[vars_to_convert], as.factor) but then It was problematic taking them back to numeric once imputated.
# Adittionally, pmm imputated according values into this vars. 

data <- read.spss("Master_TwinssCan2.3__GxE_ONLYtwins_0602025_def_woID.sav", to.data.frame = TRUE)

#vis_miss(data)

na_data_df <- sapply(data, function(x) sum(is.na(x)))
print(na_data_df)

# To test data is MCAR, tested in db without PRS mcar_test <- mcar_test(data)

# polr for ordered categorical data and polyreg unordered. Even with categorical vars: Imputation method polyreg is for categorical data.
methods <- make.method(data)
methods[c("PRS_ES_S11", "PRS_ES_S16", "PRScs_ES_auto", "C1", "C2", "idtw")] <- ""

predictor_matrix <- make.predictorMatrix(data)

predictor_matrix[,] <- 0



out_values <- c("cape_pos_freq_1", "si_referen", "si_refrem", "si_suspic","si_magic", "si_illus", "si_psych", "si_dereal", "si_hypersens", 
                "si_socisol", "si_introve",  "si_bluntaf", "si_focus", "si_associa", "si_povspee", "si_excentr", "scl90_par_9f_1", "scl90_psy_9f_1",
                "cape_dep_freq_1",  "scl90_dep_9f_1", "scl90_anx_9f_1", "scl90_pho_9f_1", "scl90_sen_9f_1", "scl90_som_9f_1", "scl90_oc_9f_1",
                "scl90_hos_9f_1", "rse_tot_Rpos", "aws1_total", "aws1_happy", "aws1_env_mast", "aws1_purpose", "aws1_selfaccept", "aws1_autonomy", 
                "aws1_growth", "aws1_relations", "aws1_vitality", "aws1_strength", "ucl_active", "ucl_palliative","ucl_avoid","ucl_social","ucl_passive",
                "ucl_emotions","ucl_thoughts", "gafpsy")


Environment <- c("ctq1_tot", "J1_tot")
#"F1O_Posschz3","si_pos_df","F1_Sipos_2fac","F1_Sipos_3fac",,"si_pos"
Pos_schiz <- c("cape_pos_freq_1", "si_referen", "si_refrem","si_suspic","si_magic", "si_illus", "si_psych", "si_dereal", "si_hypersens", "scl90_par_9f_1", "scl90_psy_9f_1")
#"si_neg_df","F2_SIneg_2fac","F2_SIassoc_3fac", "F3_SIneg_3fac", "si_neg", 
Neg_schiz <- c("si_socisol", "si_introve",  "si_bluntaf", "si_focus", "si_associa", "si_povspee", "si_excentr")
#"F1O_intern",
Internalizing <- c("cape_dep_freq_1",  "scl90_dep_9f_1","scl90_anx_9f_1","scl90_pho_9f_1","scl90_sen_9f_1", 
                   "scl90_som_9f_1","scl90_oc_9f_1", "scl90_hos_9f_1")
#"F1_2fWB", "F2_2fSelfsteem","F1_WB_1f",
Wellbeing <- c("rse_tot_Rpos","aws1_total","aws1_happy","aws1_env_mast","aws1_purpose",
               "aws1_selfaccept","aws1_autonomy","aws1_growth","aws1_relations","aws1_vitality","aws1_strength")
#"F1UCL_confront_2fac","F2UCL_nonconfr_2fac","F1_confront_3fac", "F2_nonconfr_3fac", "F3_socialCoping_3fac", 
Coping <- c("ucl_active", "ucl_palliative","ucl_avoid","ucl_social","ucl_passive","ucl_emotions","ucl_thoughts")


groups <- list(Environment, Pos_schiz, Neg_schiz, Internalizing, Wellbeing, Coping)


for(group in groups) {
  predictor_matrix[group, group] <- 1
  print(predictor_matrix[group, group])
}


gaf_var <- "gafpsy"


predictor_matrix[gaf_var, ] <- 1
#predictor_matrix[gaf_var, c("PRS_ES_S11", "PRS_ES_S16", "PRScs_ES_auto", "C1", "C2", "idtw")] <- 0  

predictor_matrix[, c("PRS_ES_S11", "PRS_ES_S16", "PRScs_ES_auto", "C1", "C2", "idtw")] <- 0

print(predictor_matrix)  






# Number of imputated datasets
nimp <- 5

# Imputate the data
data_imp <- mice(data, m = nimp, maxit = 5, method = methods, predictorMatrix = predictor_matrix, seed = 123)
plot(data_imp)

data_complete_nas <- complete(data_imp)
na_data_comp <- sapply(data_complete_nas, function(x) sum(is.na(x)))
print(na_data_comp)

#genes_values <- c("PRS_ES_S11", "PRS_ES_S16", "PRScs_ES_auto")
genes_values <- c("PRScs_ES_auto")
env_values <- c("ctq1_tot", "J1_tot")
neg_env <- c("ctq1_tot")



rev_models <- c("diathesis_stress_WEAK", "diathesis_stress_STRONG", "vantage_sensitivity_WEAK", "vantage_sensitivity_STRONG")
interaction_term <- "G*E"
covariate1 <- "C1"
covariate2 <- "C2"
idtw <- "(1|idtw)"

coef_list <- list()
fits <- list()  
pooled_coef <- matrix(NA, nrow = 5, ncol = 3)
rownames(pooled_coef) <- c("G", "E", "C1", "C2", "G:E")
colnames(pooled_coef) <- c("Chisq", "Df", "Pr(>Chisq)")


# Minimum difference in the second model AIC so it's not considered tied between the first model and the second.
tolerance <- 1

# All values stored at each iteration of the loop.
name <- character()
pG_E <- numeric()
p_G <- numeric()
p_c1 <- numeric()
p_c2 <- numeric()
FG_E <- numeric()
int <- character()
AIC <- numeric()
crossover_int <- character()
crossover_int2 <- character()
tie <- logical()
int2 <- character()
AIC2 <- numeric()
min_values <- matrix(NA, nrow = 2, ncol = nimp)
max_values <- matrix(NA, nrow = 2, ncol = nimp)
failed_combinations <- character()
fit_df_list <- list()
i = 0

append_na_results <- function() {
  pG_E <<- c(pG_E, NA)
  p_G <<- c(p_G, NA)
  p_c1 <<- c(p_c1, NA)
  p_c2 <<- c(p_c2, NA)
  FG_E <<- c(FG_E, NA)
  int <<- c(int, NA)
  AIC <<- c(AIC, NA)
  crossover_int <<- c(crossover_int, NA)
  crossover_int2 <<- c(crossover_int2, NA)
  tie <<- c(tie, NA)
  int2 <<- c(int2, NA)
  AIC2 <<- c(AIC2, NA)
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




# Loop starts
for (gene in genes_values) {
  for (e in env_values) {
    for(out in out_values){
      
      i <- i + 1
      params <- list(genes = gene, env = e)
      
      formula_dynamic <- as.formula(paste(out, "~", interaction_term, "+", covariate1, "+", covariate2, "+", idtw))
      formula_dynamic2 <- as.formula(paste(out, "~", 1, "+", idtw))
      
      # Write name of the interaction
      name_combination <- paste(i, gene, e, out, sep = "_")
      name <- c(name, name_combination)

    tryCatch({
      for (j in 1: nimp) { 
        data_complete <- complete(data_imp, action = j)

        env <- data_complete[, params$env, drop = FALSE]
        genes <- data_complete[, params$genes, drop = FALSE]

        fit <- LEGIT(data = data_complete, genes = genes, env = env, 
                    formula = formula_dynamic,
                    lme4 = TRUE, rescale = TRUE)
        
        fits[[j]] <- fit  # Store the model
        
        # Extract coefficients and standard errors
        coef_list[[j]] <- car::Anova(fit$fit_main)
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
          residual.df = rep(NA, nrow(fixed_effects)))
        
        fit_df_list[[x]] <- fit_df_imputed
      }

      w <- do.call(rbind, fit_df_list)

      # Now use the pool.table() function to pool the results
      pooled_fits <- pool.table(
        w, 
        type = "all",  # Default to all estimates
        conf.int = TRUE,  # Include confidence intervals
        conf.level = 0.95,  # 95% confidence level
        exponentiate = FALSE,  # Don't exponentiate unless it's logistic regression
        dfcom = 634,  # Use a large sample (infinite degrees of freedom)
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

        # Check to improve the plot, this one is for the first imputation
        name_png <- paste(name_combination, ".png")
        png(name_png, width = 800, height = 600)
        fit <- fits[[1]]
        #interaction_name <- paste(fit$fit_genes$formula, fit$fit_env$formula, sep = "_")
        plot <- plot(fit,ylab = fit$fit_genes$formula , xlab = fit$fit_env$formula, cex.leg=1.2)
        #text(x = 0.5, y = 0.5, labels = interaction_name, cex = 1.5, col = "blue", pos = 4)
        dev.off()
        cat("Plot saved as", name_combination)


        GxE_test_AIC_results <- lapply(1:data_imp$m, function(k) {
          imp_data <- complete(data_imp, k) 
          
          # vars_to_rescale <- c(params$genes, params$env, all.vars(formula_dynamic2))
          # imp_data[vars_to_rescale] <- scale(imp_data[vars_to_rescale])

          env <- imp_data[, params$env, drop = FALSE]
          genes <- imp_data[, params$genes, drop = FALSE]

          suppressWarnings(
          GxE_interaction_test(
            data = imp_data,
            genes = genes,
            env = env,
            formula_noGxE = formula_dynamic2,
            crossover = c("min", "max"),
            criterion = "AIC",
            lme4 = TRUE,
            rescale = TRUE
          ) )
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
          
          cross_weak <- GxE_test_AIC_results[[l]]$results["Differential susceptibility WEAK", 3]
          cross_strong <- GxE_test_AIC_results[[l]]$results["Differential susceptibility STRONG", 3]
          
          # Store the AIC values for this imputation in the matrix
          aic_values[, l] <- aic
          
          crossover[1, l] <-  cross_weak
          crossover[2, l] <- cross_strong
        }

        pooled_aic <- apply(aic_values, 1, mean, na.rm = TRUE)
        sorted_aic <- sort(pooled_aic)

        # p-values extraction
        pG_E <- c(pG_E, anova_results["G:E", "p.value"])
        p_G <- c(p_G, anova_results["G", "p.value"])
        p_c1 <- c(p_c1, anova_results["C1", "p.value"])
        p_c2 <- c(p_c2, anova_results["C2", "p.value"])
        
        FG_E <- c(FG_E, anova_results["G:E", "chisq"])
        

        mod1 <- revert_model(names(sorted_aic)[1], e)
        int <- c(int, mod1)
        AIC_int1 <- sorted_aic[1]
        AIC <- c(AIC, AIC_int1)

        tie_model <- FALSE
        
        AIC_int2 <- sorted_aic[2]
        
        if (abs(AIC_int1 - AIC_int2) <= tolerance) {
          tie_model <- TRUE
        }
        
        tie <- c(tie, tie_model)

        mod2 <- revert_model(names(sorted_aic)[2], e)
        int2 <- c(int2, mod2)
        AIC2 <- c(AIC2, AIC_int2)
        
        
        #cross_int <- GxE_test_AIC$results[1, 3]
        #crossover_int <- c(crossover_int, cross_int)
  
        
        if (names(sorted_aic)[1] == "diff_suscept_WEAK" || (names(sorted_aic)[1] == "diff_suscept_STRONG")) {
          
          res <- calculate_crossover(crossover)
          mean_min <- res$mean_min
          mean_max <- res$mean_max
          
          if (names(sorted_aic)[1] == "diff_suscept_WEAK") {
            crossover_int <- c(crossover_int, paste0(mean_min[1], "/", mean_max[1]))
          } else if (names(sorted_aic)[1] == "diff_suscept_STRONG") {
            crossover_int <- c(crossover_int, paste0(mean_min[2], "/", mean_max[2]))
          }
          

        } else {
          crossover_int <- c(crossover_int, NA)  
        }
        
        if (names(sorted_aic)[2] == "diff_suscept_WEAK" || names(sorted_aic)[2] == "diff_suscept_STRONG") {
          
          res <- calculate_crossover(crossover)
          mean_min <- res$mean_min
          mean_max <- res$mean_max

          if (names(sorted_aic)[2] == "diff_suscept_WEAK") {
            crossover_int2 <- c(crossover_int2, paste0(mean_min[1], "/", mean_max[1]))
          } else if (names(sorted_aic)[2] == "diff_suscept_STRONG") {
            crossover_int2 <- c(crossover_int2, paste0(mean_min[2], "/", mean_max[2]))
          }

        } else {
          crossover_int2 <- c(crossover_int2, NA)  
        }


      }, warning = function(w) {
        # Catch warnings specifically from this test
        cat("Warning in GxE_interaction_test for", name_combination, ":", w$message, "\n")
        append_na_results()
        
      }, error = function(e) {
        # Catch errors specifically from this test
        cat("Error in GxE_interaction_test for", name_combination, ":", e$message, "\n")
        append_na_results()
      })
      
    }
    
  }
}


# variable_list <- list(
#   PRSxEnvxOut = name,
#   P_value_GxE = pG_E,
#   P_value_G = p_G,
#   P_value_C1 = p_c1,
#   P_value_C2 = p_c2,
#   F_ratio = FG_E,
#   Type_of_interaction = int,
#   int_original = int_original,
#   AIC = AIC,
#   crossover_int = crossover_int,
#   Tie = tie,
#   Type_of_interaction2 = int2,
#   int2_original = int2_original,
#   AIC2 = AIC2,
#   crossover_int2 = crossover_int2
# )
# 
# lengths <- sapply(variable_list, length)
# print(lengths)


results <- data.frame (PRSxEnvxOut = name , P_value_GxE = pG_E , P_value_G= p_G,
                       P_value_C1 = p_c1, P_value_C2 = p_c2, F_ratio = FG_E,
                       Type_of_interaction = int, AIC = AIC, 
                       crossover_int = crossover_int, Tie = tie, 
                       Type_of_interaction2 = int2,  
                       AIC2 = AIC2, crossover_int2 = crossover_int2)

write_xlsx(results, "results_interactions_MI_PRS-ESauto.xlsx")


