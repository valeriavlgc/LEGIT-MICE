library(foreign)

source("config.R")
source("utils.R")
source("data_simulation.R")
source("data_imputation.R")
source("data_analysis_LEGIT.R")

# -----------------------------------------------------------------------------
# Step 1: Simulate data
#
# simulate_data() returns a data.frame with:
#   - G1, G2            : genetic indicators (LEGIT estimates their weights)
#   - E1, E2            : environmental indicators (E2 is adverse, see neg_env)
#   - C1, C2            : model covariates included in the LEGIT formula
#   - Y1, Y2            : outcomes
#   - Aux1, Aux2, Aux3, Aux4 : auxiliary variables — used in MICE only,
#                         never entered into LEGIT. Always fully observed.
# -----------------------------------------------------------------------------
dt <- simulate_data(n = 1000, seed = 123)

# -----------------------------------------------------------------------------
# Step 2: Introduce MCAR missingness and impute
#
# missings_data_MCAR() applies 15% Missing Completely At Random missingness
# to the analysis variables only (G, E, C, Y). Auxiliary variables remain
# fully observed — their completeness is what makes them useful for MICE.
#
# imputate_data() runs MICE using all available variables as predictors,
# including Y predicting G/E and the auxiliary variables predicting everything.
# -----------------------------------------------------------------------------
dt_miss <- missings_data_MCAR(dt, percentage = 0.15, seed = 123)
dt_imp  <- imputate_data(dt_miss, seed = 123)

# -----------------------------------------------------------------------------
# Step 3: Define analysis variables
# -----------------------------------------------------------------------------
genes_values  <- c("G1", "G2")
env_values    <- c("E1", "E2")
out_values    <- c("Y1", "Y2")
neg_env       <- c("E2")      # E2 was simulated with reversed sign
multilevel_id <- "idtw"       # grouping variable — used only for MULTILEVEL

analysis_type <- ANALYSIS_TYPE$UNILEVEL

# -----------------------------------------------------------------------------
# Step 4: Run the three-way comparison
#
#   1. Complete data   — true GxE signal, upper benchmark
#   2. MCAR + listwise — signal degradation under missingness (no imputation)
#   3. MCAR + MICE     — signal recovery after imputation
#
# Expected pattern: F_complete > F_listwise, F_imputed > F_listwise
# -----------------------------------------------------------------------------

# 1. Complete data — no missingness
cat("\n===== COMPLETE DATA =====\n")
run_analysis_LEGIT(genes_values, env_values, out_values, dt,
                   analysis_type, imputated = FALSE, multilevel_id = multilevel_id)

# 2. MCAR data — listwise deletion (no imputation)
cat("\n===== MCAR — LISTWISE DELETION =====\n")
run_analysis_LEGIT(genes_values, env_values, out_values, dt_miss,
                   analysis_type, imputated = FALSE, multilevel_id = multilevel_id)

# 3. MCAR data — MICE imputed (5 imputations, pooled with Rubin's rules)
cat("\n===== MCAR — MICE IMPUTED =====\n")
run_analysis_LEGIT(genes_values, env_values, out_values, dt_imp,
                   analysis_type, imputated = TRUE, multilevel_id = multilevel_id)

# To run on real data (uncomment and set DATA_PATH in config.R):
# run_analysis_LEGIT(genes_values, env_values, out_values,
#                    read.spss(DATA_PATH, to.data.frame = TRUE),
#                    i, analysis_type, imputated = FALSE, multilevel_id = multilevel_id)






