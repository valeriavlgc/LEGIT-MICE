# Solve the package problem
library(foreign)

source("config.R")
source("utils.R")
source("data_simulation.R")
source("data_imputation.R")
source("data_analysis_LEGIT.R")

## add seed to config ?
# see in which other source files I would need to upload rest of R scripts?
# probaar todo, multilevel y unilevel f ratio both for imputated? is neg env working?

dt <- simulate_data(n=1000, seed=123)
dt_miss <- missings_data(dt = dt, percentage=0.15, seed=123)
dt_imp <- imputate_data(dt= dt_miss, seed=123)


# genes_values <- c("G1")
# env_values <- c("E1", "E2")
# neg_env <- c("E1")
# out_values <- c("Y1", "Y2")

genes_values <- c("G1")
env_values <- c("E2")
out_values <- c("Y2")
neg_env <- c("")


#add amount of data missingnes to fomrula missing.



#check_vars(colnames(complete(dt_imp, 1)), c(genes_values, env_values, out_values, neg_env))
           
analysis_type <- ANALYSIS_TYPE$UNILEVEL

# Original db
run_analysis_LEGIT(genes_values, env_values, out_values, dt, i, analysis_type, imputated = FALSE)
# Db with MCAR data 
run_analysis_LEGIT(genes_values, env_values, out_values, dt_miss, i, analysis_type, imputated = FALSE)
# imputated data
run_analysis_LEGIT(genes_values, env_values, out_values, dt_imp, i, analysis_type, imputated = TRUE)

# To run on real data
#run_analysis_LEGIT(genes_values, env_values, out_values, read.spss(DATA_PATH, to.data.frame=TRUE), i, analysis_type, imputated = FALSE)






