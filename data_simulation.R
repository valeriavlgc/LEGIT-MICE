library(data.table)
library(missMethods)
library(MASS)

# Explore if this can actually be done or not, I don't think so, if i create missingness in there i cannot get it back inside the 2wayexample, 
# and extracting the values wouldn't work

# simulate_data <- function(n, seed) {
# 
#   Twoway <- example_2way(n,1,logit=FALSE, seed=seed)
#   dt <- Twoway$data
#   dt2 <- Twoway$G
#   dt3 <- Twoway$E
#   dt.combined <- cbind(dt, dt2, dt3)
# 
#   return(dt.combined)
# }


simulate_data <- function(n, seed,
                          gene_means = c(9.25, 7.5),
                          sd_genes = c(5.10, 4.0),
                          env_means = c(2.3, 3.0),
                          sd_envs = c(1.0, 1.2)) {
  set.seed(seed)
  
  # Correlated genetic variables (G1, G2)
  cor_genes <- matrix(c(1, 0.5,
                        0.5, 1), nrow = 2)
  Sigma_genes <- diag(sd_genes) %*% cor_genes %*% diag(sd_genes)
  genes <- mvrnorm(n, mu = gene_means, Sigma = Sigma_genes)
  colnames(genes) <- c("G1", "G2")
  
  # Correlated environmental variables (E1, E2)
  cor_envs <- matrix(c(1, 0.6,
                       0.6, 1), nrow = 2)
  Sigma_envs <- diag(sd_envs) %*% cor_envs %*% diag(sd_envs)
  envs <- mvrnorm(n, mu = env_means, Sigma = Sigma_envs)
  colnames(envs) <- c("E1", "E2")
  
  # Covariates
  covariates <- data.frame(
    C1 = rnorm(n, mean = 0, sd = 0.03),
    C2 = rnorm(n, mean = 0, sd = 0.04)
  )
  
  # Centered variables for interaction
  G1_c <- scale(genes[, "G1"], center = TRUE, scale = FALSE)
  G2_c <- scale(genes[, "G2"], center = TRUE, scale = FALSE)
  E1_c <- scale(envs[, "E1"], center = TRUE, scale = FALSE)
  E2_c <- scale(envs[, "E2"], center = TRUE, scale = FALSE)
  
  interactions <- data.frame(
    G1xE1 = G1_c * E1_c,
    G1xE2 = G1_c * E2_c,
    G2xE1 = G2_c * E1_c,
    G2xE2 = G2_c * E2_c
  )
  
  # Coefficients
  b0 <- 10
  bG1 <- -0.002
  bG2 <- 0.004
  bE1 <- -0.003
  bE2 <- 0.002
  bG1E1 <- 0.6
  bG1E2 <- -0.6
  bG2E1 <- 0.5
  bG2E2 <- -0.3
  bC1 <- 0.02
  bC2 <- -0.01
  
  # Linear predictors for outcomes
  lp_Y1 <- b0 + bG1 * genes[, "G1"] + bE1 * envs[, "E1"] + bG1E1 * interactions$G1xE1 +
    bC1 * covariates$C1 + bC2 * covariates$C2
  
  lp_Y2 <- b0 + bG1 * genes[, "G1"] + bE2 * envs[, "E2"] + bG1E2 * interactions$G1xE2 +
    bC1 * covariates$C1 + bC2 * covariates$C2
  
  lp_Y3 <- b0 + bG2 * genes[, "G2"] + bE1 * envs[, "E1"] + bG2E1 * interactions$G2xE1 +
    bC1 * covariates$C1 + bC2 * covariates$C2
  
  lp_Y4 <- b0 + bG2 * genes[, "G2"] + bE2 * envs[, "E2"] + bG2E2 * interactions$G2xE2 +
    bC1 * covariates$C1 + bC2 * covariates$C2
  
  # Correlated residual errors
  cor_Y <- matrix(c(1, 0.5, 0.4, 0.3,
                    0.5, 1, 0.4, 0.3,
                    0.4, 0.4, 1, 0.6,
                    0.3, 0.3, 0.6, 1), nrow = 4)
  residual_sd <- rep(14, 4)
  Sigma_Y <- diag(residual_sd) %*% cor_Y %*% diag(residual_sd)
  errors <- mvrnorm(n, mu = rep(0, 4), Sigma = Sigma_Y)
  
  # Final outcomes
  outcomes <- data.frame(
    Y1 = lp_Y1 + errors[, 1],
    Y2 = lp_Y2 + errors[, 2],
    Y3 = lp_Y3 + errors[, 3],
    Y4 = lp_Y4 + errors[, 4]
  )
  
  dt.data <- data.frame(
    genes, envs, covariates, outcomes
  )
  
  return(dt.data)
}



missings_data <- function(dt,percentage, seed) {
  set.seed(seed)

  # https://cran.r-project.org/web/packages/missMethods/vignettes/Generating-missing-values.html
  dt_mcar <- delete_MCAR(dt, percentage)

  return(dt_mcar)
}





# missings_data <- function(dt, seed) {
#   set.seed(seed)
# 
#   # 1. First create COMPLETELY OBSERVED anchor variables
#   # These will control MAR in other variables
#   anchor_vars <- c("G_anchor", "E_anchor", "Y_anchor")
# 
#   # Create anchors by sampling existing variables
#   dt$G_anchor <- dt$G1
#   dt$E_anchor <- dt$E1
#   dt$Y_anchor <- dt$Y1
# 
#   # 2. Create missingness in BLOCKS
# 
#   dt <- delete_MAR_rank(
#     ds = dt,
#     p = 0.10,
#     cols_mis = c("C1", "C2", "G1", "G2"),
#     cols_ctrl = rep("G_anchor", 4),
#     n_mis_stochastic = TRUE
#   )
# 
#   # Block 2: Environment variables (E1, E2)
#   # 10% missing, using E_anchor (always observed)
#   dt <- delete_MAR_rank(
#     ds = dt,
#     p = 0.10,
#     cols_mis = c("E1", "E2"),
#     cols_ctrl = rep("E_anchor", 2),
#     n_mis_stochastic = TRUE
#   )
# 
#   # Block 3: Outcome variables (Y1-Y4)
#   # 10% missing, using Y_anchor (always observed)
#   dt <- delete_MAR_rank(
#     ds = dt,
#     p = 0.10,
#     cols_mis = c("Y1", "Y2", "Y3", "Y4"),
#     cols_ctrl = rep("Y_anchor", 4),
#     n_mis_stochastic = TRUE
#   )
# 
#   # Remove temporary anchor variables
#   dt <- dt[, !names(dt) %in% anchor_vars]
# 
#   # Add missingness diagnostics
#   attr(dt, "missing_pattern") <- list(
#     covariates_genes = mean(rowSums(is.na(dt[, c("C1", "C2", "G1", "G2")])) > 0),
#     environment = mean(rowSums(is.na(dt[, c("E1", "E2")])) > 0),
#     outcomes = mean(rowSums(is.na(dt[, paste0("Y", 1:4)])) > 0)
#   )
# 
#   return(dt)
# }





  




# 
# 
# generate_missing <- function(data, 
#                              missing_prop = 0.15,
#                              seed = NULL,
#                              miss_type = "MCAR") {
#   
#   if (!is.null(seed)) set.seed(seed)
#   
#   if (miss_type == "MCAR") {
#     # Missing completely at random
#     miss_data <- delete_MCAR(
#       data, 
#       p = missing_prop,
#       cols_mis = setdiff(colnames(data), c("C1", "C2"))
#     )
#   } else {
#     stop("Currently only MCAR missingness implemented")
#   }
#   
#   # Add missing data diagnostics
#   attr(miss_data, "missing_prop") <- missing_prop
#   attr(miss_data, "missing_type") <- miss_type
#   attr(miss_data, "missing_cols") <- 
#     names(which(colSums(is.na(miss_data)) > 0))
#   
#   return(miss_data)
# }




