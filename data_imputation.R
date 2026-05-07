library(dplyr)
library(naniar)
library(mice)
library(ggmice)

#max iterations 5?

imputate_data <- function(dt_miss, seed) {
  
  # Number of missings in each variable // so far same amount of missings per var.
  colSums(is.na(dt_miss))
  
  # Pattern of missing data
  print(plot_pattern(dt_miss, vrb = c("Y1", "Y2", "Y3", "Y4"), square = FALSE, rotate = TRUE))

  # Little’s (1988) test statistic: The null hypothesis in this test is that the data is MCAR. p > 0.05 data is likely MCAR
  mcar <- mcar_test(dt_miss)
  
  if (mcar$p.value > 0.05) {
    message("Data is likely MCAR (p = ", round(mcar$p.value, 4), "). Proceeding with imputation...")
  } else {
    message("Data is not MCAR (p = ", round(mcar$p.value, 4), "). Proceeding under assumption that data is MAR...")
  }
  
  
  # Definition of the imputation model
  predictor_matrix <- make.predictorMatrix(dt_miss)
  
  predictor_matrix[,] <- 0
  
  Genes <- c("G1", "G2")
  Covariates <- c("C1", "C2")
  
  Environment <- c("E1", "E2")
  
  Outcomes1 <- c("Y1", "Y2")
  Outcomes2 <- c("Y3", "Y4")
  
  groups <- list(Environment, Outcomes1, Outcomes2, Covariates, Genes)
  
  
  for(group in groups) {
    predictor_matrix[group, group] <- 1
    print(predictor_matrix[group, group])
  }
  
  
  print(plot_pred(predictor_matrix))

  data_imp <- mice(dt_miss, m = nimp, 
                   maxit = 5, method = "pmm", # pmm adequate for continuus data, other methods available polr for ordered categorical data and polyreg unordered. 
                   predictorMatrix = predictor_matrix,
                   seed = seed) 

  summary(data_imp)
  print(xyplot(data_imp, Y2 ~ E2))


  print(stripplot(data_imp))
  
  print(ggmice(data_imp, aes(x = .imp, y = Y1)) +
          geom_jitter(height = 0, width = 0.25) +
          labs(x = "Imputation number"))
  
  # Imputed values vs. observed values. Example on var 1
  print(dt_miss.Y1 <- rbind(data.table(Y1 = unlist(data_imp$imp$Y1), imputed = TRUE),
                            data.table(Y1 = na.omit(dt_miss$Y1), imputed = FALSE)))
  
  # Histogram: Explame on var1
  gg1.Y1 <- ggplot(dt_miss.Y1, aes(Y1, group = imputed, fill = imputed))
  gg1.Y1 <- gg1.Y1 + geom_histogram(aes(y=..count../sum(..count..)),
                                    position = "dodge")
  
  print(gg1.Y1)

  return(data_imp)
}


# imputate_data <- function(dt_miss, seed) {
#   
#   # Number of missings in each variable // so far same amount of missings per var.
#   colSums(is.na(dt_miss))
#   
#   #vis_miss(dt_miss)
#   
#   # Pattern of missing data
#   print(plot_pattern(dt_miss, vrb = c("Y1", "Y2", "Y3", "Y4"), square = FALSE, rotate = TRUE))
#   # npat = 30
#   
#   # print(ggmice(dt_miss, aes(Y2, E2)) +
#   # geom_point())
#   
#   # print(plot_corr(dt_miss))
# 
#   # Little’s (1988) test statistic: The null hypothesis in this test is that the data is MCAR. p > 0.05 data is likely MCAR
#   mcar <- mcar_test(dt_miss)
#   
#   if (mcar$p.value > 0.05) {
#     message("Data is likely MCAR (p = ", round(mcar$p.value, 4), "). Proceeding with imputation...")
#   } else {
#     message("Data is not MCAR (p = ", round(mcar$p.value, 4), "). Proceeding under assumption that data is MAR...")
#   }
#   
# 
#   
#   # Definition of the imputation model
#   
#   predictor_matrix <- make.predictorMatrix(dt_miss)
# 
#   predictor_matrix[,] <- 0
# 
#   Genes <- c("G1", "G2")
#   Covariates <- c("C1", "C2")
#   
#   Environment <- c("E1", "E2")
#   
#   Outcomes1 <- c("Y1", "Y2")
#   Outcomes2 <- c("Y3", "Y4")
# 
#   groups <- list(Environment, Outcomes1, Outcomes2, Covariates, Genes)
# 
#   
#   for(group in groups) {
#     predictor_matrix[group, group] <- 1
#     print(predictor_matrix[group, group])
#   }
# 
# 
#   print(plot_pred(predictor_matrix))
#   #print(predictor_matrix)
#   
#   
#   # Generating imputated datasets
#   
#   # Imputate the data
#   # nimp is the number of imputations in this case set to 5
#   data_imp <- mice(dt_miss, m = nimp, maxit = 5, predictorMatrix = predictor_matrix, seed = seed)
#   # method = methods, "pmm"
#   # polr for ordered categorical data and polyreg unordered. Even with categorical vars: Imputation method polyreg is for categorical data.
#   # print(plot(data_imp))
#   summary(data_imp)
#   print(xyplot(data_imp, Y2 ~ E2))
#   # print(densityplot(data_imp))
#   # colSums(is.na(data_imp))
#   
#   # maybe add if interesting: table(cci(dt_miss.mice)); str(cc(dt_miss.mice)); str(ic(dt_miss.mice)); check first db: str(complete(dt_miss.mice, action = 1))
#   
#   
#   # Plots:
#   #boxplot(data_imp$imp$Y1)
#   
#   print(stripplot(data_imp))
#   
#   print(ggmice(data_imp, aes(x = .imp, y = Y1)) +
#     geom_jitter(height = 0, width = 0.25) +
#     labs(x = "Imputation number"))
#   
#   # Imputed values vs. observed values
#   
#   print(dt_miss.Y1 <- rbind(data.table(Y1 = unlist(data_imp$imp$Y1), imputed = TRUE),
#                       data.table(Y1 = na.omit(dt_miss$Y1), imputed = FALSE)))
#   
#   # Histogram: 
#   
#   gg1.Y1 <- ggplot(dt_miss.Y1, aes(Y1, group = imputed, fill = imputed))
#   gg1.Y1 <- gg1.Y1 + geom_histogram(aes(y=..count../sum(..count..)),
#                                     position = "dodge")
#   
#   print(gg1.Y1)
#   
#   
#   #stripplot(data_imp, Y1~.imp, pch=20, cex=2)
#   
#   return(data_imp)
# }
# 
