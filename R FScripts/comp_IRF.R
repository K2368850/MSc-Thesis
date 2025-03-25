# library(vars)
library(dplyr)
library(zoo)
library(bvartools)
library(factorstochvol)
library(Matrix)

setwd("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code")
# load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code/shrinkTVPVAR-NG-SV.RData")
load("StartingPointNoXForm.RData")
load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code/all_models_new_new.Rdata")
# rm(list=ls())

cc_var_irf_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/Data\ 07-05/Code/IRF/"

actual_vars <- c("FEDFUNDS", "WFRBST01134", "WFRBSB50215", "MCOILWTICO", 
                 "T10Y3MM", "MICH", "T10Y2YM", "M2SL", 
                 "HOUST", "PCECC96", "PSAVERT", "AAAFFM", "AAA10YM", "UMCSENT", "GDPC1",
                 "INDPRO", "PPIACO", "DPCERD3Q086SBEA", 
                 "A822RD3Q086SBEA", "OPHNFB", "NASDAQCOM",
                 "LNS14000001", "LNS14000002","LNS14000003", "LNS14000006", "LNS14000009")
houst_ind <- which(actual_vars=="HOUST")
savings_ind <- which(actual_vars=="PSAVERT")

actual_vars <- actual_vars[-c(houst_ind, savings_ind)]

imp <-c("FEDFUNDS", "GDPC1", "INDPRO", "PCECC96", "PPIACO", "OPHNFB", "WFRBST01134", "WFRBSB50215", "DPCERD3Q086SBEA", "A822RD3Q086SBEA",
        "M2SL", "T10Y2YM", "LNS14000001", "LNS14000002","LNS14000003", "LNS14000006", "LNS14000009")
resp <- actual_vars

split_data <- function(data){
  x.train <- head(data, -20)
  x.test <- tail(data, 20)
  
  return (list(train = x.train, test=x.test))
}

scale_no_demean <- function(data){
  data <- macro %>% select(-all_of("date"))
  
  # Step 1: Compute the mean of each column
  means <- colMeans(data)
  
  # Step 2: Compute the standard deviation of each column
  sds <- apply(data, 2, sd)
  
  # Step 3: Scale each column and add back the mean
  scaled_data <- as.data.frame(
    scale(data, center = FALSE, scale = sds) + means
  )
  return (scaled_data)
}

create_lagged_data <- function(data, max_lag) {
  original_data_trunc <- data[(max_lag+1):nrow(data), ]
  # Generate lagged variables
  lagged_data <- lapply(1:max_lag, function(l) {
    as.data.frame(lapply(data, function(x) lag(x, l)))
  })
  lagged_data <- do.call(cbind, lagged_data)
  colnames(lagged_data) <- paste0(rep(colnames(data), max_lag), "_lag", 
                                  rep(1:max_lag, each=ncol(data)))
  
  # Combine original data with lagged variables, excluding the first max_lag rows
  combined_data <- cbind(data[(max_lag+1):nrow(data), ], 
                         lagged_data[(max_lag+1):nrow(data), ])
  
  # Remove rows with missing values
  combined_data <- na.trim(combined_data)
  
  # Create the formula for the model
  lhs <- paste("cbind(", paste(colnames(data), collapse = ","), ")") # lhs <- paste(colnames(data), collapse = "+")
  rhs <- paste(colnames(lagged_data), collapse = "+")
  formula <- as.formula(paste(lhs, "~", rhs))
  
  # Return the cleaned data and formula
  list(data = combined_data, formula = formula,
       original_var_data=original_data_trunc,
       lag_data = lagged_data[(max_lag+1):nrow(data), ],
       f_rhs = rhs)
}

to_4D <- function(x, remove=TRUE){
  d <- dim(x)
  # if(remove) d[2] <- d[2]-1
  new_mat <- array(0, dim=c(d[1], d))
  for(j in 1:d[2]){
    for(k in 1:d[3]){
      new_mat[ , ,j, k] <- diag(x[, j, k], d[1], d[1])
    }
  }
  return (new_mat)
}

to_3D <- function(x, txpose = TRUE){
  if(txpose) x <- t(x)
  d <- dim(x)
  new_mat <- array(0, dim=c(d[1], d))
  for(j in 1:d[2]){
    new_mat[ , ,j] <- diag(x[, j], d[1], d[1])
  }
  return(new_mat)
}

construct_A0 <- function(mat, dim=24, niter=5000){
  result_array <- array(0, dim = c(dim, dim, niter))
  
  # Define the indices for the lower triangular part of a 24x24 matrix
  lower_tri_indices <- which(lower.tri(matrix(0, dim, dim), diag = FALSE), arr.ind = TRUE)
  
  # Loop over each column (sample) in your original matrix
  for (i in 1:niter) {
    # Create an empty 24x24 matrix for the current sample
    current_matrix <- matrix(0, dim, dim)
    
    # Fill the lower triangular part of the matrix
    current_matrix[lower_tri_indices] <- mat[, i]
    
    diag(current_matrix) <- 1
    
    # Assign the current matrix to the result array
    result_array[, , i] <- current_matrix
  }
  return(result_array)
}

export_irfs <- function(var_obj, impulses, var_type, save_path, n=20){
  for(i in 1:length(impulses)){
    png(filename=paste0(save_path, var_type, "-IRF-", impulses[i], ".png"), width = 1000, height = 2000)
    plot(irf(var_obj, n.ahead=n, ci=0.68, impulse=impulses[i], cumulative=TRUE, boot=TRUE, ortho=FALSE))
    dev.off()
  }
}

export_irfs_1by1 <- function(var_obj, impulses, all_vars, var_type, save_path, n=4){
  for(i in 1:length(impulses)){
    png(filename=paste0(save_path, var_type, "-IRF-", impulses[i], ".png"), width = 1000, height = 2000)
    par(mfrow = c(12, 2))
    for(j in 1:length(all_vars)){
      print(paste(impulses[i], "on", all_vars[j]))
      # png(filename=paste0(save_path, var_type, "-IRF-", impulses[i], ".png"), width = 1000, height = 2000)
      plot(bvartools::irf(var_obj, n.ahead=n, ci=0.68, impulse=impulses[i], response=all_vars[j], cumulative=TRUE, shock=1), #, type="gir"
           main=paste("Responses to", impulses[i]), xlab="Time", ylab=all_vars[j])
      # dev.off()
    }
    dev.off()
  }
}

data <- macro[, -1]

prior_sd <- apply(data, 2, sd)

data <- scale(data, center=FALSE)
p <- 4

data_list <- create_lagged_data(data, p)

lag_data <- split_data(data_list$lag_data)
x.train <- lag_data$train
x.test <- lag_data$test

# data_split <- split_data(data_list$original_var_data) #split_data(scaled_data)
data_split <- split_data(data)
y.train <- data_split$train
y.test <- data_split$test

t <- nrow(y.train) - p
# t <- nrow(y.train)
K <- ncol(y.train)
pK <- p*K
niter <- 5000

# load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code/shrinkTVPVAR-NG-SV.RData")

# stvpvar_bvar <- bvartools::bvar(y=ts(y.train[-c(1:p), ]), #data=y.train,
#                  A0 = list(coeffs=matrix(array(rep(temp_var$pred_objs$final_A, t), dim = c(K, K, t, niter)), nrow=(t*K^2), ncol=niter),
#                            errors = matrix(temp_var$pred_objs$theta_sr_SIGMA, nrow=(K^2), ncol=niter)), # A0 =matrix(temp_var$pred_objs$final_A, nrow=(K^2), ncol=niter),
#                  A = list(coeffs=matrix(temp_var$beta, nrow=(K^2 * p * t), ncol=niter),
#                           errors=matrix(temp_var$theta_sr[, -97, ], nrow=(p*(K^2)), ncol=niter)),
#                  C = list(coeffs=matrix(temp_var$beta_const, nrow=(K*t), ncol=niter),
#                           errors=matrix((temp_var$theta_sr^2)[, pK+1, ], ncol=niter)),
#                  Sigma = list(coeffs = matrix(temp_var$Sigma, nrow=(K^2 * t), ncol=niter),
#                               errors = matrix(to_3D(temp_var$pred_objs$sv_sigma2), nrow=(K^2), ncol=niter))
#                  )

bayesianVARs_names <- names(all_models)[sapply(all_models, function(x) class(x) == "bayesianVARs_bvar")]
bayesianVARs_names <- bayesianVARs_names[-which(bayesianVARs_names == "VAR.SV.SSVS")]
bayesianVARs_names <- "CC.VAR.SV.HS"
# bayesianVARs_names < "VAR.SV.HS"
save_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/FEIR/CC\ VARs/"
setwd(save_path)
p <- 4
gc()
for(i in 1:length(bayesianVARs_names)){
  # load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code/all_models_new_new.Rdata")
  name <- bayesianVARs_names[i]
  temp_var <- all_models[[name]]
  if(!grepl("FSV", name, ignore.case = TRUE)){
  print(name)
  gc()
  # temp_A0 <- construct_A0(temp_var$U, K, niter)
  # temp_A0 <- matrix(rep(temp_A0, t), nrow=(t*K^2), ncol=niter)
  # temp_A <- matrix(rep(coef(temp_var)[-97, ,], t), nrow=(p*t*K*K), ncol=niter)
  # temp_C <- matrix(rep(coef(temp_var)[97, ,], t), nrow=(t*K), ncol=niter)
  # temp_A0 <- 
  temp_bayesianVARs <- bvartools::bvar(y = ts(y.train[-c(1:4),]),
                                       # A0=list(coeffs=matrix(array(rep(temp_var$U, t), dim = c(K, K, t, niter)), nrow=(t*K^2), ncol=niter),
                                       #         errors = matrix(0, nrow=(K*K), ncol=niter)),
                                       A0 = list(coeffs = matrix(rep(construct_A0(temp_var$U, K, niter), t), nrow=(t*K^2), ncol=niter),
                                                 errors = matrix(0, nrow=(K*K), ncol=niter)),
                                       A = list(coeffs = matrix(rep(coef(temp_var)[-97, ,], t), nrow=(p*t*K*K), ncol=niter),
                                                errors=matrix(0, nrow=(K*K), ncol=niter)),
                                       C = list(coeffs=matrix(rep(coef(temp_var)[97, ,], t), nrow=(t*K), ncol=niter),
                                                error=matrix(0, nrow=(K), ncol=niter)),
                                       Sigma = list(coeffs = matrix(vcov(temp_var), nrow=(t*K^2), ncol=niter), 
                                                    errors = matrix(to_3D(temp_var$sv_para[3, ,], txpose=FALSE), nrow=K*K, ncol=niter))
                                       # A0 = matrix(construct_A0(temp_var$U, K, niter), nrow=(K*K), ncol=(niter)),
                                       # A = matrix(coef(temp_var)[-97, ,], nrow=p*K*K, ncol= (niter)),
                                       # A = list(coeffs=matrix(array(rep(coef(temp_var)[-97, ,], t), dim = c(p*K, K, t, niter)), nrow=(p*t*K^2), ncol=niter),
                                       #      errors = matrix(0, nrow=(K*K), ncol=niter)),
                                       # C = list(coeffs = matrix(array(rep(coef(temp_var)[97,,], t), dim=c(K, t, niter)), nrow=K*t, ncol=niter),
                                       #          errors = matrix(0, nrow=(K), ncol=niter)),
                                       # Sigma = list(coeffs = matrix(vcov(temp_var), nrow=(t*K^2), ncol=niter), 
                                       #              errors = matrix(to_3D(temp_var$sv_para[3, ,], txpose=FALSE), nrow=K*K, ncol=niter))
  )
  print("removing")
  rm(temp_A0, temp_A, temp_C)
  gc()
  print("plotting")
  export_irfs_1by1(temp_bayesianVARs, imp, resp, name, save_path, n=20)
  rm(temp_bayesianVARs)
  gc()
  } else {
    print(name)
    gc()
    # temp_A0 <- construct_A0(temp_var$U, K, niter)
    # temp_A0 <- as.matrix(rep(temp_A0, t), nrow=(t*K^2), ncol=niter)
    # temp_A <- matrix(rep(coef(temp_var)[-97, ,], t), nrow=(p*t*K*K), ncol=niter)
    # temp_C <- matrix(rep(coef(temp_var)[97, ,], t), nrow=(t*K), ncol=niter)
    # temp_A0 <- 
    temp_bayesianVARs <- bvartools::bvar(y = ts(y.train[-c(1:4), ]),
                                         # A0 = list(coeffs = temp_A0, errors = matrix(0, nrow=(K*K), ncol=niter)),
                                         A = list(coeffs = matrix(rep(coef(temp_var)[-97, ,], t), nrow=(p*t*K*K), ncol=niter),
                                                  errors=matrix(0, nrow=(K*K), ncol=niter)),
                                         C = list(coeffs=matrix(rep(coef(temp_var)[97, ,], t), nrow=(t*K), ncol=niter),
                                                  error=matrix(0, nrow=(K), ncol=niter)),
                                         Sigma = list(coeffs = matrix(vcov(temp_var), nrow=(t*K^2), ncol=niter), 
                                                      errors = matrix(to_3D(temp_var$sv_para[3, ,], txpose=FALSE), nrow=K*K, ncol=niter))
    )
    gc()
    print("removing")
    rm(temp_A, temp_C)
    gc()
    print("plotting")
    export_irfs_1by1(temp_bayesianVARs, imp, resp, name, save_path, n=20)
    rm(temp_bayesianVARs)
    gc()
  }
}


# temp_bayesianVARs <- bvartools::bvar(y = ts(y.train[-c(1:p), ]),
#                           A0 = matrix(construct_A0(temp_var$U, K, niter), nrow=(K*K), ncol=(niter)),
#                           A = matrix(coef(temp_var)[-97, ,], nrow=p*K*K, ncol= (niter)),
#                           C = matrix(coef(temp_var)[97,,], nrow=1, ncol=niter),
#                           Sigma = list(coeffs = matrix(vcov(temp_var), nrow=t*K^2, ncol=niter), 
#                                        errors = matrix(to_3D(temp_var$sv_para[3, ,], txpose=FALSE), nrow=K*K, ncol=niter))
#                           )
# 
# stvpvar_bvar <- bvartools::bvar(y=ts(y.train[-c(1:p), ]), #data=y.train,
#                                 A0 = list(coeffs=matrix(array(rep(temp_var$pred_objs$final_A, t), dim = c(K, K, t, niter)), nrow=(t*K^2), ncol=niter),
#                                           errors = matrix(temp_var$pred_objs$theta_sr_SIGMA, nrow=(K^2), ncol=niter)), # A0 =matrix(temp_var$pred_objs$final_A, nrow=(K^2), ncol=niter),
#                                 A = list(coeffs=matrix(temp_var$beta, nrow=(K^2 * p * t), ncol=niter),
#                                          errors=matrix(temp_var$theta_sr[, -97, ], nrow=(p*(K^2)), ncol=niter)),
#                                 C = list(coeffs=matrix(temp_var$beta_const, nrow=(K*t), ncol=niter),
#                                          errors=matrix((temp_var$theta_sr^2)[, pK+1, ], ncol=niter)),
#                                 Sigma = list(coeffs = matrix(temp_var$Sigma, nrow=(K^2 * t), ncol=niter),
#                                              errors = matrix(to_3D(temp_var$pred_objs$sv_sigma2), nrow=(K^2), ncol=niter))
# )

construct_lambda1 <- function(hyper, K, n){
  result <- array(0, dim=c(K, K, n))
  
  for (j in 1:n){
    temp <- matrix(0, nrow=24, ncol=24)
    off_diag <- upper.tri(temp) | lower.tri(temp)
    temp[off_diag] <- hyper[, j]
    diag(temp) <- 1
    result[,,j] <- temp
  }
  return(as.matrix(result, nrow=(K*K), ncol=n))
}

construct_lambda2 <- function(hyper){
  i <- nrow(hyper)/2
  # j <- ncol(hyper)
  return(as.matrix(hyper[1:i, ], nrow=(nrow(hyper)/2), ncol=ncol(hyper)))
}


name <- "VAR.SV.SSVS"
temp_var <- all_models[[name]]
print(name)
gc()
temp_A0 <- construct_A0(temp_var$U, K, niter)
temp_A0 <- as.matrix(rep(temp_A0, t), nrow=(t*K^2), ncol=niter)
temp_A0_lambda <- construct_lambda1(temp_var$u_hyperparameter, K, niter)
temp_A <- as.matrix(rep(coef(temp_var)[-97, ,], t), nrow=(p*t*K*K), ncol=niter)
temp_A_lambda <- construct_lambda2(temp_var$phi_hyperparameter)
temp_C <- as.matrix(rep(coef(temp_var)[97, ,], t), nrow=(t*K), ncol=niter)
gc()# temp_A0 <- 
temp_bayesianVARs <- bvartools::bvar(y = ts(y.train),
                                     # A0=list(coeffs=matrix(array(rep(temp_var$U, t), dim = c(K, K, t, niter)), nrow=(t*K^2), ncol=niter),
                                     #         errors = matrix(0, nrow=(K*K), ncol=niter)),
                                     A0 = list(coeffs = temp_A0, lambda =temp_A0_lambda),
                                     A = list(coeffs = temp_A, lambda=temp_A_lambda),#errors=matrix(0, nrow=(K*K), ncol=niter)),
                                     C = list(coeffs=temp_C), #, error=matrix(0, nrow=(K), ncol=niter)),
                                     Sigma = list(coeffs = matrix(vcov(temp_var), nrow=(t*K^2), ncol=niter), 
                                                  errors = matrix(to_3D(temp_var$sv_para[3, ,], txpose=FALSE), nrow=K*K, ncol=niter))
                                     # A0 = matrix(construct_A0(temp_var$U, K, niter), nrow=(K*K), ncol=(niter)),
                                     # A = matrix(coef(temp_var)[-97, ,], nrow=p*K*K, ncol= (niter)),
                                     # A = list(coeffs=matrix(array(rep(coef(temp_var)[-97, ,], t), dim = c(p*K, K, t, niter)), nrow=(p*t*K^2), ncol=niter),
                                     #      errors = matrix(0, nrow=(K*K), ncol=niter)),
                                     # C = list(coeffs = matrix(array(rep(coef(temp_var)[97,,], t), dim=c(K, t, niter)), nrow=K*t, ncol=niter),
                                     #          errors = matrix(0, nrow=(K), ncol=niter)),
                                     # Sigma = list(coeffs = matrix(vcov(temp_var), nrow=(t*K^2), ncol=niter), 
                                     #              errors = matrix(to_3D(temp_var$sv_para[3, ,], txpose=FALSE), nrow=K*K, ncol=niter))
)
print("removing")
rm(temp_A0, temp_A, temp_A0_lambda, temp_A_lambda, temp_C)
gc()
print("plotting")
export_irfs_1by1(temp_bayesianVARs, imp, resp, name, save_path, n=20)
rm(temp_bayesianVARs)
gc()
