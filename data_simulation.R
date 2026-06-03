library(data.table)
library(MASS)

# -----------------------------------------------------------------------------
# simulate_data()
#
# Simulates a dataset whose structure matches LEGIT's mathematical model,
# enriched with auxiliary variables that mimic the richness of real datasets
# (e.g. Aux1, Aux2, Aux3, Aux4).
#
# LEGIT assumes:
#   - Observed genetic variables (G1, G2) are indicators of a latent genetic
#     index:  G_latent = w1*G1 + w2*G2
#   - Observed environmental variables (E1, E2) are indicators of a latent
#     environmental index: E_latent = v1*E1 + v2*E2
#   - The outcome is generated as:
#       Y = b0 + bG*G_latent + bE*E_latent + bGE*(G_latent * E_latent) + noise
#
# Auxiliary variables (Aux1, Aux2, Aux3, Aux4) are correlated with G, E,
# and/or Y but are NOT part of the LEGIT model. Their role is purely in the
# imputation step: by including them in the MICE predictor matrix, MICE gains
# indirect information about the missing values of G and E, which improves
# imputation quality and partial recovery of the GxE effect. This mirrors
# what happens with real epidemiological data, where many auxiliary measures
# are available alongside the main analysis variables.
#
# The returned data.frame contains only observable columns — MICE and LEGIT
# never see the latent indices, exactly as with real data.
# -----------------------------------------------------------------------------

simulate_data <- function(n, seed,
                          # True weights defining the latent genetic index (sum to 1)
                          w_genes    = c(0.6, 0.4),
                          # True weights defining the latent environmental index (sum to 1)
                          w_envs     = c(0.7, 0.3),
                          # Means for the genetic indicators (0 = already centred)
                          gene_means = c(0, 0),
                          # Means for the environmental indicators
                          env_means  = c(0, 0)) {
  set.seed(seed)

  # ---------------------------------------------------------------------------
  # 1. Genetic indicators (G1, G2)
  #    Multivariate normal with mild correlation (r = 0.3, mimicking LD).
  #    SD = 1 so all variables are on a common scale from the start.
  # ---------------------------------------------------------------------------
  cor_genes   <- matrix(c(1, 0.3, 0.3, 1), nrow = 2)
  Sigma_genes <- diag(2) %*% cor_genes %*% diag(2)
  genes <- mvrnorm(n, mu = gene_means, Sigma = Sigma_genes)
  colnames(genes) <- c("G1", "G2")

  # ---------------------------------------------------------------------------
  # 2. Environmental indicators (E1, E2)
  #    Moderate correlation between environments (r = 0.4). SD = 1.
  #    E1 is a positive environment; E2 will be treated as adverse (neg_env),
  #    meaning its sign is reversed when generating Y2.
  # ---------------------------------------------------------------------------
  cor_envs   <- matrix(c(1, 0.4, 0.4, 1), nrow = 2)
  Sigma_envs <- diag(2) %*% cor_envs %*% diag(2)
  envs <- mvrnorm(n, mu = env_means, Sigma = Sigma_envs)
  colnames(envs) <- c("E1", "E2")

  # ---------------------------------------------------------------------------
  # 3. Model covariates (C1, C2)
  #    These enter the LEGIT formula directly as nuisance variables.
  #    Independent, SD = 1.
  # ---------------------------------------------------------------------------
  C1 <- rnorm(n, mean = 0, sd = 1)
  C2 <- rnorm(n, mean = 0, sd = 1)

  # ---------------------------------------------------------------------------
  # 4. Auxiliary variables — NOT in the LEGIT model, only in MICE
  #
  #    In real studies you always have more variables than those in the focal
  #    analysis (demographics, other biomarkers, etc.). Including them in the
  #    imputation model allows MICE to better predict missing values of G and E,
  #    because they carry correlated information.
  #
  #    Aux1: continuous, mildly correlated with E (r ≈ 0.2). SD = 10.
  #    Aux2: continuous, slight association with G and with the outcome Y.
  #    Aux3: continuous, correlated with E1 (r ≈ 0.8).
  #    Aux4: continuous, correlated with G1 (r ≈ 0.85).
  #    Aux5: continuous, correlated with G2 (r ≈ 0.8) — directly aids G2 imputation.
  #    Aux6: continuous, correlated with E2 (r ≈ 0.8) — directly aids E2 imputation.
  # ---------------------------------------------------------------------------
  Aux1 <- rnorm(n, mean = 0, sd = 10)

  # Aux2 has a slight association with G and Y (no structural dependency here)
  Aux2 <- rnorm(n, mean = 0, sd = 1)

  # Aux3 is correlated with E1 (r ≈ 0.8) plus independent noise
  Aux3 <- 0.8 * envs[, "E1"] + rnorm(n, mean = 0, sd = sqrt(1 - 0.8^2))

  # Aux4 is correlated with G1 (r ≈ 0.85) plus noise
  Aux4 <- 0.85 * genes[, "G1"] + rnorm(n, mean = 0, sd = sqrt(1 - 0.85^2))

  # Aux5 is correlated with G2 (r ≈ 0.8) plus noise — directly improves G2 imputation
  Aux5 <- 0.8 * genes[, "G2"] + rnorm(n, mean = 0, sd = sqrt(1 - 0.8^2))

  # Aux6 is correlated with E2 (r ≈ 0.8) plus noise — directly improves E2 imputation
  Aux6 <- 0.8 * envs[, "E2"] + rnorm(n, mean = 0, sd = sqrt(1 - 0.8^2))

  # ---------------------------------------------------------------------------
  # 5. Latent indices
  #    Centered raw variables are combined with the true weights.
  #    Centering before the product prevents the interaction term being
  #    confounded by the means of G and E.
  # ---------------------------------------------------------------------------
  G1_c <- scale(genes[, "G1"], center = TRUE, scale = FALSE)
  G2_c <- scale(genes[, "G2"], center = TRUE, scale = FALSE)
  E1_c <- scale(envs[, "E1"],  center = TRUE, scale = FALSE)
  E2_c <- scale(envs[, "E2"],  center = TRUE, scale = FALSE)

  G_latent <- w_genes[1] * G1_c + w_genes[2] * G2_c   # latent genetic index
  E_latent <- w_envs[1]  * E1_c + w_envs[2]  * E2_c   # latent env index

  # ---------------------------------------------------------------------------
  # 6. Outcome regression coefficients
  #    Standardised scale (SD = 1 for all predictors).
  #    bGE = 0.5 gives a medium-large interaction (R² ≈ 15-20% for GxE alone).
  #    Residual SD = 1 keeps SNR high so the interaction is clearly detectable
  #    in complete data but degrades meaningfully under 15% MCAR.
  #
  #    Y1: positive environment (E_latent enters with its natural sign).
  #    Y2: adverse environment (E_latent sign reversed), to test neg_env logic.
  #    Both outcomes are also influenced by Aux1, Aux2, and Aux3 so that the
  #    auxiliary variables carry real predictive information for imputation.
  # ---------------------------------------------------------------------------
  b0       <-  0      # intercept
  bG       <-  0.3    # main effect of latent G
  bE       <-  0.3    # main effect of latent E
  bGE      <-  0.5    # GxE interaction — the signal LEGIT must detect
  bC1      <-  0.1    # nuisance covariate 1
  bC2      <- -0.1    # nuisance covariate 2
  bAux1    <-  0.02   # small Aux1 effect on outcome
  bAux2    <-  0.15   # modest Aux2 effect on outcome
  bAux3    <-  0.10   # modest Aux3 effect on outcome

  lp_Y1 <- b0 + bG * G_latent + bE * E_latent + bGE * (G_latent * E_latent) +
            bC1 * C1 + bC2 * C2 +
            bAux1 * Aux1 + bAux2 * Aux2 + bAux3 * Aux3   # aux vars inform Y1

  # E2 is adverse: negate E_latent so a higher score means worse environment
  lp_Y2 <- b0 + bG * G_latent - bE * E_latent + bGE * (G_latent * (-E_latent)) +
            bC1 * C1 + bC2 * C2 +
            bAux1 * Aux1 + bAux2 * Aux2 + bAux3 * Aux3   # aux vars inform Y2

  # ---------------------------------------------------------------------------
  # 7. Add residual noise and assemble the data.frame
  #    Auxiliary variables are included as columns so MICE can use them.
  #    Latent indices are NOT included — they are unobservable constructs.
  # ---------------------------------------------------------------------------
  Y1 <- lp_Y1 + rnorm(n, mean = 0, sd = 1)
  Y2 <- lp_Y2 + rnorm(n, mean = 0, sd = 1)

  dt.data <- data.frame(
    genes,                             # G1, G2  — LEGIT genetic indicators
    envs,                              # E1, E2  — LEGIT environmental indicators
    C1, C2,                            # covariates in the LEGIT formula
    Y1, Y2,                            # outcomes
    Aux1, Aux2, Aux3, Aux4, Aux5, Aux6 # auxiliary: used in imputation only
  )

  return(dt.data)
}


# -----------------------------------------------------------------------------
# missings_data_MCAR()
#
# Introduces Missing Completely At Random (MCAR) missingness using
# mice::ampute() - the same package used for imputation, keeping the full
# pipeline within a single coherent framework (van Buuren, 2018).
#
# Under MCAR the probability of a value being missing is independent of both
# observed and unobserved data - it is purely random. Listwise deletion is
# also unbiased under MCAR, so the main advantage of MICE here is efficiency
# (retaining the full sample size) rather than bias correction.
#
# Implementation details:
#   - Only the 8 analysis variables (G, E, C, Y) receive missingness.
#     Auxiliary variables (Aux1, Aux2, Aux3, Aux4) are always fully observed.
#   - ampute() is called with mech = "MCAR" so all rows have equal probability
#     of being selected for each missingness pattern.
#   - The default patterns matrix creates one pattern per analysis variable
#     (each pattern makes exactly one variable missing). With bycases = FALSE
#     and prop = percentage, each variable ends up with ~percentage% missing
#     cells, matching the intended 15% per variable.
#   - The incomplete analysis columns are then recombined with the always-
#     observed auxiliary columns to return the full dataset.
# -----------------------------------------------------------------------------
missings_data_MCAR <- function(dt, percentage, seed) {
  set.seed(seed)

  # ---------------------------------------------------------------------------
  # Cell-level MCAR: each analysis variable independently loses `percentage`
  # of its values, selected uniformly at random.
  #
  # We implement this manually rather than using ampute() because ampute works
  # at the ROW level (a row receives a pattern). Row-level patterns cause
  # correlated missingness within groups of variables (e.g. G1 and G2 both
  # missing for the same rows), which removes the natural predictors MICE
  # relies on (G2 predicting G1) and severely degrades imputation quality.
  #
  # Cell-level MCAR is the standard in simulation studies for GxE:
  # missingness in G1 is independent of missingness in G2, E1, Y, etc.
  # For each row, missing values in different variables are uncorrelated,
  # so MICE always has the other variables available to impute each one.
  #
  # Auxiliary variables (Aux1, Aux2, Aux3, Aux4) are excluded — they
  # remain fully observed to serve as predictors in the imputation model.
  # ---------------------------------------------------------------------------
  aux_vars      <- c("Aux1", "Aux2", "Aux3", "Aux4", "Aux5", "Aux6")
  analysis_vars <- setdiff(names(dt), aux_vars)

  n_miss <- round(nrow(dt) * percentage)

  for (var in analysis_vars) {
    # Randomly select n_miss row indices to become missing for this variable.
    # Each variable is sampled independently, so missingness across variables
    # is uncorrelated — the defining property of MCAR.
    miss_idx      <- sample(nrow(dt), n_miss)
    dt[miss_idx, var] <- NA
  }

  return(dt)
}


# Alias for backwards compatibility
missings_data <- missings_data_MCAR


# -----------------------------------------------------------------------------
# missings_data_MAR()
#
# Introduces Missing At Random (MAR) missingness.
#
# Under MAR the probability of a value being missing depends on OTHER observed
# variables but NOT on the missing value itself. Listwise deletion is biased
# under MAR; MICE recovers valid estimates provided the variables that predict
# missingness are included in the imputation model (which they are here, since
# Aux1–Aux4 are always fully observed and entered as predictors).
#
# MAR predictor assignment (all fully observed Aux variables):
#   G1 -> Aux4  (r = 0.85 with G1)
#   G2 -> Aux5  (r = 0.80 with G2)
#   E1 -> Aux3  (r = 0.80 with E1)
#   E2 -> Aux6  (r = 0.80 with E2)
#   C1, C2,
#   Y1, Y2  -> Aux1  (mild correlation with the environment)
#
# Mechanism: for each analysis variable, observations are sampled for
# missingness with probability proportional to the MAR predictor value
# (weighted sampling without replacement). This guarantees exactly n_miss
# missing entries per variable while creating a clear monotone relationship
# between the predictor and missingness probability — the defining feature
# of a MAR mechanism that MICE can exploit.
# -----------------------------------------------------------------------------
missings_data_MAR <- function(dt, percentage, seed) {
  set.seed(seed)

  aux_vars      <- c("Aux1", "Aux2", "Aux3", "Aux4", "Aux5", "Aux6")
  analysis_vars <- setdiff(names(dt), aux_vars)
  n_miss        <- round(nrow(dt) * percentage)

  mar_predictor <- function(var) {
    if      (var == "G1")  dt[["Aux4"]]   # r = 0.85 with G1
    else if (var == "G2")  dt[["Aux5"]]   # r = 0.80 with G2
    else if (var == "E1")  dt[["Aux3"]]   # r = 0.80 with E1
    else if (var == "E2")  dt[["Aux6"]]   # r = 0.80 with E2
    else                   dt[["Aux1"]]   # C1, C2, Y1, Y2
  }

  for (var in analysis_vars) {
    pred     <- mar_predictor(var)
    # Shift to strictly positive so values can serve as sampling weights.
    # Higher Aux value → higher probability of being selected as missing.
    weights  <- pred - min(pred) + 1e-6
    miss_idx <- sample(nrow(dt), n_miss, prob = weights, replace = FALSE)
    dt[miss_idx, var] <- NA
  }

  return(dt)
}
