library(foreign)
library(dplyr)
library(writexl)
library(parallel)
library(doParallel)
library(foreach)
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

n_values   <- c(1000)
pct_values <- c(0.15, 0.2)
n_reps     <- 50

genes_values  <- c("G1", "G2")
env_values    <- c("E1", "E2")
out_values    <- c("Y1", "Y2")
neg_env       <- c("E2")
analysis_type <- ANALYSIS_TYPE$UNILEVEL


# -----------------------------------------------------------------------------
# Parallelization setup
# -----------------------------------------------------------------------------
n_cores <- max(1, detectCores() - 2)  # leave 2 cores free
cat("Using", n_cores, "cores\n")
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Export everything the workers need
clusterExport(cl, c("genes_values", "env_values", "out_values", "neg_env",
                    "analysis_type", "nimp", "tolerance", "rev_models",
                    "interaction_term", "covariate1", "covariate2", "ANALYSIS_TYPE"))
clusterEvalQ(cl, {
  library(LEGIT)
  library(mice)
  library(dplyr)
  source("config.R")
  source("utils.R")
  source("data_simulation.R")
  source("data_imputation.R")
  source("data_analysis_LEGIT.R")
})

if (!dir.exists("results/checkpoints")) dir.create("results/checkpoints", recursive = TRUE)

# -----------------------------------------------------------------------------
# Simulation loop — checkpointed per n/pct condition
# -----------------------------------------------------------------------------
all_results <- list()

for (n in n_values) {
  for (pct in pct_values) {
    
    checkpoint_file <- paste0("results/checkpoints/checkpoint_n", n, "_pct", round(pct*100), ".xlsx")
    
    # Skip if already done
    if (file.exists(checkpoint_file)) {
      cat("Skipping n =", n, "| pct =", pct, "— checkpoint exists\n")
      all_results <- c(all_results, list(readxl::read_xlsx(checkpoint_file)))
      next
    }
    
    cat("\n=== Running n =", n, "| pct =", pct, "===\n")
    
    condition_results <- foreach(
      rep = 1:n_reps,
      .combine  = "c",
      .packages = c("LEGIT", "mice", "dplyr"),
      .export   = c("genes_values", "env_values", "out_values", "neg_env",
                    "analysis_type", "nimp", "tolerance", "rev_models",
                    "interaction_term", "covariate1", "covariate2", "ANALYSIS_TYPE")
    ) %dopar% {
      
      seed <- rep * 1000 + n + round(pct * 100)
      
      dt           <- simulate_data(n = n, seed = seed)
      dt_miss_mcar <- missings_data_MCAR(dt, percentage = pct, seed = seed)
      dt_miss_mar  <- missings_data_MAR( dt, percentage = pct, seed = seed)
      dt_imp_mcar  <- imputate_data(dt_miss_mcar, seed = seed, verbose = FALSE)
      dt_imp_mar   <- imputate_data(dt_miss_mar,  seed = seed, verbose = FALSE)
      
      res_complete            <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                    dt, analysis_type, imputated = FALSE,
                                                    save_results = FALSE, verbose = FALSE)
      res_complete$condition  <- "complete"
      res_complete$mech       <- "none"
      res_complete$n          <- n
      res_complete$pct        <- pct
      res_complete$rep        <- rep
      
      res_listwise_mcar            <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                         dt_miss_mcar, analysis_type, imputated = FALSE,
                                                         save_results = FALSE, verbose = FALSE)
      res_listwise_mcar$condition  <- "listwise"
      res_listwise_mcar$mech       <- "MCAR"
      res_listwise_mcar$n          <- n
      res_listwise_mcar$pct        <- pct
      res_listwise_mcar$rep        <- rep
      
      res_imputed_mcar            <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                        dt_imp_mcar, analysis_type, imputated = TRUE,
                                                        save_results = FALSE, verbose = FALSE)
      res_imputed_mcar$condition  <- "imputed"
      res_imputed_mcar$mech       <- "MCAR"
      res_imputed_mcar$n          <- n
      res_imputed_mcar$pct        <- pct
      res_imputed_mcar$rep        <- rep
      
      res_listwise_mar            <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                        dt_miss_mar, analysis_type, imputated = FALSE,
                                                        save_results = FALSE, verbose = FALSE)
      res_listwise_mar$condition  <- "listwise"
      res_listwise_mar$mech       <- "MAR"
      res_listwise_mar$n          <- n
      res_listwise_mar$pct        <- pct
      res_listwise_mar$rep        <- rep
      
      res_imputed_mar            <- run_analysis_LEGIT(genes_values, env_values, out_values,
                                                       dt_imp_mar, analysis_type, imputated = TRUE,
                                                       save_results = FALSE, verbose = FALSE)
      res_imputed_mar$condition  <- "imputed"
      res_imputed_mar$mech       <- "MAR"
      res_imputed_mar$n          <- n
      res_imputed_mar$pct        <- pct
      res_imputed_mar$rep        <- rep
      
      list(res_complete, res_listwise_mcar, res_imputed_mcar,
           res_listwise_mar, res_imputed_mar)
    }
    
    # Save checkpoint for this n/pct condition
    condition_df <- do.call(rbind, lapply(condition_results, as.data.frame))
    write_xlsx(condition_df, checkpoint_file)
    cat("Checkpoint saved:", checkpoint_file, "\n")
    
    all_results <- c(all_results, list(condition_df))
  }
}

stopCluster(cl)

# -----------------------------------------------------------------------------
# Combine and save final results
# -----------------------------------------------------------------------------
final <- do.call(rbind, all_results)
if (!dir.exists("results")) dir.create("results")
write_xlsx(final, paste0("results/simulation_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"))
cat("\nSimulation complete. Results saved in /results\n")

