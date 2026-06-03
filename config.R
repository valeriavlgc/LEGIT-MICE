DATA_PATH <- "data/Master_TwinssCan2.3__GxE_ONLYtwins_0602025_def_1.sav"
rev_models <- c("diathesis_stress_WEAK", "diathesis_stress_STRONG", "vantage_sensitivity_WEAK", "vantage_sensitivity_STRONG")
interaction_term <- "G*E"
covariate1 <- "C1"
covariate2 <- "C2"
nimp <- 10
tolerance <- 1
ANALYSIS_TYPE <- list(
  MULTILEVEL = "multilevel",
  UNILEVEL = "unilevel"
)