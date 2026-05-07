library(LEGIT)
library(foreign)
library(writexl)
library(dplyr)
library(naniar)
library(mice)


data_complete <- complete(dt_imp, action = 1)
params <- list(genes = "G1", env = "E2")
env <- data_complete[, params$env, drop = FALSE]
genes <- data_complete[, params$genes, drop = FALSE]
formula_dynamic <- as.formula(paste("Y1", "~", interaction_term, "+", covariate1, "+", covariate2))
#
fit <- LEGIT(data = data_complete, genes = genes, env = env,
              formula = formula_dynamic)
#
print(fit$fit_main)
#

fit_glm = glm(formula=fit$fit_main$formula, family=fit$fit_main$family, data=fit$fit_main$data)
anova(fit_glm, test = "F")