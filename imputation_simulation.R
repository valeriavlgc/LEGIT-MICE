library(foreign)
library(dplyr)
source("config.R")
source("utils.R")
source("data_simulation.R")
source("data_imputation.R")
source("data_analysis_LEGIT.R")

# -----------------------------------------------------------------------------
# Simulation parameters
# -----------------------------------------------------------------------------
# n_values   <- c(250, 500, 1000)
# pct_values <- c(0.15, 0.25)
# n_reps     <- 50

n_values   <- c(250)
pct_values <- c(0.15, 0.2)
n_reps     <- 50

genes_values  <- c("G1", "G2")
env_values    <- c("E1", "E2")
out_values    <- c("Y1", "Y2")
neg_env       <- c("E2")
analysis_type <- ANALYSIS_TYPE$UNILEVEL

all_results <- list()

# -----------------------------------------------------------------------------
# Simulation loop
# -----------------------------------------------------------------------------
for (n in n_values) {
  for (pct in pct_values) {
    for (rep in 1:n_reps) {
      
      cat("\n--- n =", n, "| pct =", pct, "| rep =", rep, "---\n")
      
      seed <- rep * 1000 + n + round(pct * 100)
      
      dt           <- simulate_data(n = n, seed = seed)
      dt_miss_mcar <- missings_data_MCAR(dt, percentage = pct, seed = seed)
      dt_miss_mar  <- missings_data_MAR( dt, percentage = pct, seed = seed)
      dt_imp_mcar  <- imputate_data(dt_miss_mcar, seed = seed, verbose = FALSE)
      dt_imp_mar   <- imputate_data(dt_miss_mar,  seed = seed, verbose = FALSE)

      # ── Complete data (reference) ───────────────────────────────────────────
      res_complete           <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                   dt, analysis_type, imputated = FALSE,
                                                   save_results = FALSE, verbose = FALSE)
      res_complete$condition <- "complete"
      res_complete$mech      <- "none"
      res_complete$n         <- n
      res_complete$pct       <- pct
      res_complete$rep       <- rep

      # ── MCAR: listwise deletion ─────────────────────────────────────────────
      res_listwise_mcar           <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                        dt_miss_mcar, analysis_type, imputated = FALSE,
                                                        save_results = FALSE, verbose = FALSE)
      res_listwise_mcar$condition <- "listwise"
      res_listwise_mcar$mech      <- "MCAR"
      res_listwise_mcar$n         <- n
      res_listwise_mcar$pct       <- pct
      res_listwise_mcar$rep       <- rep

      # ── MCAR: MICE imputed ──────────────────────────────────────────────────
      res_imputed_mcar           <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                       dt_imp_mcar, analysis_type, imputated = TRUE,
                                                       save_results = FALSE, verbose = FALSE)
      res_imputed_mcar$condition <- "imputed"
      res_imputed_mcar$mech      <- "MCAR"
      res_imputed_mcar$n         <- n
      res_imputed_mcar$pct       <- pct
      res_imputed_mcar$rep       <- rep

      # ── MAR: listwise deletion ──────────────────────────────────────────────
      res_listwise_mar           <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                       dt_miss_mar, analysis_type, imputated = FALSE,
                                                       save_results = FALSE, verbose = FALSE)
      res_listwise_mar$condition <- "listwise"
      res_listwise_mar$mech      <- "MAR"
      res_listwise_mar$n         <- n
      res_listwise_mar$pct       <- pct
      res_listwise_mar$rep       <- rep

      # ── MAR: MICE imputed ───────────────────────────────────────────────────
      res_imputed_mar           <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                      dt_imp_mar, analysis_type, imputated = TRUE,
                                                      save_results = FALSE, verbose = FALSE)
      res_imputed_mar$condition <- "imputed"
      res_imputed_mar$mech      <- "MAR"
      res_imputed_mar$n         <- n
      res_imputed_mar$pct       <- pct
      res_imputed_mar$rep       <- rep

      all_results <- c(all_results,
                       list(res_complete,
                            res_listwise_mcar, res_imputed_mcar,
                            res_listwise_mar,  res_imputed_mar))
    }
  }
}

# -----------------------------------------------------------------------------
# Save all results
# -----------------------------------------------------------------------------
final <- do.call(rbind, lapply(all_results, as.data.frame))
if (!dir.exists("results")) dir.create("results")
write_xlsx(final, paste0("results/simulation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"))
cat("\nSimulation complete. Results saved in /results\n")

library(readxl)

latest_file <- tail(sort(list.files("results/", pattern = "^simulation_.*\\.xlsx$", full.names = TRUE)), 1)
cat("Reading:", latest_file, "\n")

results <- read_xlsx(latest_file)

# Power — proportion of replications where the GxE effect was significant
power <- results %>%
  group_by(n, pct, mech, condition, PRSxEnvxOut) %>%
  summarise(
    power        = mean(P_value_GxE < 0.05, na.rm = TRUE),
    mean_F       = mean(F_ratio, na.rm = TRUE),
    mean_GxE_est = mean(GxE_estimate, na.rm = TRUE),
    bias_GxE = mean(abs(GxE_estimate), na.rm = TRUE) - 0.5,
    .groups = "drop"
  )

# Model recovery — how often the correct interaction type was identified:
model_recovery <- results %>%
  group_by(n, pct, mech, condition, PRSxEnvxOut) %>%
  summarise(
    pct_correct_model = mean(Type_of_interaction == "diff_suscept_STRONG", na.rm = TRUE),
    .groups = "drop"
  )

# imputed vs listwise vs complete, split by mechanism:
power %>%
  tidyr::pivot_wider(
    names_from  = c(mech, condition), 
    values_from = c(power, mean_F, mean_GxE_est, bias_GxE)
  ) %>%
  print(n = Inf, width = Inf)