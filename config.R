DATA_PATH <- "data/Master_TwinssCan2.3__GxE_ONLYtwins_0602025_def_1.sav"
rev_models <- c("diathesis_stress_WEAK", "diathesis_stress_STRONG", "vantage_sensitivity_WEAK", "vantage_sensitivity_STRONG")
interaction_term <- "G*E"
covariate1 <- "C1"
covariate2 <- "C2"
idtw <- "(1|idtw)"
nimp <- 5
coef_list <- list()
fits <- list()  
pooled_coef <- matrix(NA, nrow = 5, ncol = 3)
rownames(pooled_coef) <- c("G", "E", "C1", "C2", "G:E")
colnames(pooled_coef) <- c("Chisq", "Df", "Pr(>Chisq)")
# Minimum difference in the second model AIC so it's not considered tied between the first model and the second.
tolerance <- 1
results_list <- list()
crossover_int <- "" 
crossover_int2 <- "" 
crossover2 <- ""
tie <- logical()
min_values <- matrix(NA, nrow = 2, ncol = nimp)
max_values <- matrix(NA, nrow = 2, ncol = nimp)
failed_combinations <- character()
fit_df_list <- list()
i = 0
ANALYSIS_TYPE <- list(
  MULTILEVEL = "multilevel",
  UNILEVEL = "unilevel"
)