# library(vars)
library(dplyr)
library(zoo)
library(bvartools)
library(factorstochvol)
library(Matrix)
library(shrinkTVPVAR)

setwd("~/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code")
load("StartingPointNoXForm.RData")

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


resp_short <- c("WFRBST01134", "WFRBSB50215", "LNS14000001", "LNS14000002","LNS14000003", "LNS14000006", "LNS14000009")

model_files_loc <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/"
setwd(model_files_loc)
files <- grep("-10000", list.files(), value=TRUE, ignore.case=FALSE)
models_file_names <- paste0(model_files_loc, files)

save_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/Data\ 07-05/Code/IRF\ ShrinkTVPVAR/"
save_path2 <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/Data\ 07-05/Code/IRF\ ShrinkTVPVAR\ Short/"
to_3D <- function(x, txpose = TRUE){
  if(txpose) x <- t(x)
  d <- dim(x)
  new_mat <- array(0, dim=c(d[1], d))
  for(j in 1:d[2]){
    new_mat[ , ,j] <- diag(x[, j], d[1], d[1])
  }
  return(new_mat)
}
export_irfs_1by1 <- function(var_obj, impulses, all_vars, var_type, save_path, n=4){
  for(i in 1:length(impulses)){
    png(filename=paste0(save_path, var_type, "-IRF-", impulses[i], ".png"), width = 1000, height = 2000)
    par(mfrow = c(12, 2))
    for(j in 1:length(all_vars)){
      print(paste(impulses[i], "on", all_vars[j]))
      # png(filename=paste0(save_path, var_type, "-IRF-", impulses[i], ".png"), width = 1000, height = 2000)
      plot(bvartools::irf(var_obj, n.ahead=n, ci=0.68, impulse=impulses[i], response=all_vars[j], cumulative=TRUE, shock=1, type="gir"),
           main=paste("Responses to", impulses[i]), xlab="Time", ylab=all_vars[j])
      # dev.off()
    }
    dev.off()
  }
}

export_irfs_short <- function(var_obj, impulses, all_vars, var_type, save_path, n=4){
  for(i in 1:length(impulses)){
    png(filename=paste0(save_path, var_type, "-IRF-", impulses[i], ".png"), width = 1000, height = 2000)
    par(mfrow = c(5, 1))
    for(j in 1:length(all_vars)){
      print(paste(impulses[i], "on", all_vars[j]))
      # png(filename=paste0(save_path, var_type, "-IRF-", impulses[i], ".png"), width = 1000, height = 2000)
      plot(bvartools::irf(var_obj, n.ahead=n, ci=0.68, impulse=impulses[i], response=all_vars[j], cumulative=TRUE, shock=1, type="gir"),
           main=paste("Responses to", impulses[i]), xlab="Time", ylab=all_vars[j])
      # dev.off()
    }
    dev.off()
  }
}



p<-4
t <- nrow(y.train) - p
K <- ncol(y.train)
pK <- p*K
niter <- 5000

model_names <- sapply(files, function(x) substr(x, 1, nchar(x) - 6))
i <- 3
for(i in c(1,2,4)){
  load(models_file_names[i])
  t <- nrow(y.train) - p
  K <- ncol(y.train)
  pK <- p*K
  niter <- 5000
  print(paste("loaded:", model_names[[i]]))
  temp_var <- list(pred_objs=list(final_A = temp_var$pred_objs$final_A,
                                  theta_sr_SIGMA = temp_var$pred_objs$theta_sr_SIGMA,
                                  sv_sigma2 = temp_var$pred_objs$sv_sigma2),
                   beta = temp_var$beta, beta_const = temp_var$beta_const, theta_sr = temp_var$theta_sr, Sigma = temp_var$Sigma)
  gc()
  temp_A0 <- matrix(rep(temp_var$pred_objs$final_A, t), nrow=(t*K^2), ncol=niter)
  temp_A0_err <- matrix(temp_var$pred_objs$theta_sr_SIGMA, nrow=(K^2), ncol=niter)
  print(sum(is.na(temp_A0)))
  print(sum(is.na(temp_A0_err)))
  temp_var$pred_objs$final_A <- NULL
  temp_var$pred_objs$theta_sr_SIGMA <- NULL
  gc()
  print(dim(temp_A0))
  temp_A <- matrix(temp_var$beta, nrow=(K^2 * p * t), ncol=niter)
  temp_A_err <- matrix(temp_var$theta_sr[, -((p*K)+1), ], nrow=(p*(K^2)), ncol=niter)
  print(sum(is.na(temp_A)))
  print(sum(is.na(temp_A_err)))
  temp_var$beta <- NULL
  print(dim(temp_A))
  temp_C <- matrix(temp_var$beta_const, nrow=(K*t), ncol=niter)
  temp_C_err <- matrix((temp_var$theta_sr^2)[, ((p*K) + 1), ], ncol=niter)
  print(sum(is.na(temp_C)))
  print(sum(is.na(temp_C_err)))
  temp_var$beta_const <- NULL
  temp_var$theta_sr <- NULL
  print(dim(temp_C))
  temp_S <- matrix(temp_var$Sigma, nrow=(K^2 * t), ncol=niter)
  temp_S_err <- matrix(to_3D(temp_var$pred_objs$sv_sigma2), nrow=(K^2), ncol=niter)
  print(sum(is.na(temp_C)))
  print(sum(is.na(temp_C_err)))
  # browser
  rm(temp_var)
  stvpvar_bvar <- bvartools::bvar(y=ts(y.train[-c(1:p), ]),
                                  A0 = list(coeffs=temp_A0, errors=temp_A0_err),
                                  # A0 = list(coeffs=matrix(array(rep(temp_var$pred_objs$final_A, t), dim = c(K, K, t, niter)), nrow=(t*K^2), ncol=niter),
                                  #           errors = matrix(temp_var$pred_objs$theta_sr_SIGMA, nrow=(K^2), ncol=niter)), # A0 =matrix(temp_var$pred_objs$final_A, nrow=(K^2), ncol=niter),
                                  A = list(coeffs=temp_A, errors=temp_A_err),
                                  C = list(coeffs=temp_C, errors=temp_C_err),
                                  Sigma = list(coeffs = temp_S,errors = temp_S_err)
  )
  print("removing")
  # rm(temp_var)
  rm(temp_A, temp_A_err, temp_A0, temp_A0_err, temp_C, temp_C_err, temp_S, temp_S_err)
  gc()
  export_irfs_1by1(stvpvar_bvar, imp, actual_vars, model_names[[i]], save_path, n=20)
  export_irfs_short(stvpvar_bvar, imp, resp_short, model_names[[i]], save_path2, n=20) 
  rm(stvpvar_bvar)
  gc()
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

load("/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/shrinkTVPVAR-NG-SV-SDscale.RData")
load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code/all_models_new_new.Rdata")
name <- "VAR.SV.HS"
temp_var <- all_models[[name]]
p <- 4
t <- nrow(y.train) - p
K <- ncol(y.train)
pK <- p*K
niter <- 2500
print(paste("loaded:", model_names[[i]]))

niter <- 5000

A <- construct_A0(temp_var$U)
# A <- temp_var$pred_objs$final_A
A <- apply(A, 3, solve)
A <- array(A, dim=c(24,24,niter))
beta <- array(temp_var$PHI[-97,,], dim = c(24, 4, 24, niter))
beta <- aperm(beta, c(1,3,2, 4))
# beta <- temp_var$beta
# beta_med <- beta[,,,115,]
phi_l <- list()
phi_l[[1]] <- diag(1, 24)
for(i in 1:20){
   # idx <- i %% 4
   # if(idx == 0) idx <- 4
  print(i)
   if(i == 1) {phi_l[[i+1]] <- array(apply(beta[,,1,], 3, function(x) return(x %*% phi_l[[i]])), dim=c(24, 24, niter))}
   else if(i == 2) {
     # print(dim(phi_l[[i]]))
     phi_l[[i+1]] <- array(sapply(1:niter, function(x) beta[,,1, x] %*% phi_l[[i]][,,x]), dim=c(24,24,niter))
     phi_l[[i+1]] <- phi_l[[i+1]] + array(apply(beta[,,2,], 3, function(x) return(x %*% phi_l[[1]])), dim=c(24, 24, niter))
     # phi_l[[i+1]] <- array(sapply(1:niter, function(x) phi_l[[i+1]][,,x] %*% A[,,x]), dim=c(24,24,niter))
     phi_l[[i+1]] <- array(sapply(1:niter, function(x) phi_l[[i+1]][,,x] %*% t(A[,,x])), dim=c(24,24,niter))
   }
   else if(i==3){
     phi_l[[i+1]] <- array(sapply(1:niter, function(x) beta[,,1, x] %*% phi_l[[i]][,,x]), dim=c(24,24,niter))
     phi_l[[i+1]] <- phi_l[[i+1]] + array(sapply(1:niter, function(x) beta[,,2, x] %*% phi_l[[i-1]][,,x]), dim=c(24,24,niter))
     # phi_l[[i+1]] <- phi_l[[i+1]] + array(sapply(1:niter, function(x) beta[,,3, x] %*% phi_l[[i-2]][,,x]), dim=c(24,24,niter))
     phi_l[[i+1]] <- phi_l[[i+1]] + array(apply(beta[,,3,], 3, function(x) return(x %*% phi_l[[1]])), dim=c(24, 24, niter))
     # phi_l[[i+1]] <- array(sapply(1:niter, function(x) phi_l[[i+1]][,,x] %*% A[,,x]), dim=c(24,24,niter))
     phi_l[[i+1]] <- array(sapply(1:niter, function(x) phi_l[[i+1]][,,x] %*% t(A[,,x])), dim=c(24,24,niter))
   }
  else if(i==4){
    phi_l[[i+1]] <- array(sapply(1:niter, function(x) beta[,,1, x] %*% phi_l[[i]][,,x]), dim=c(24,24,niter))
    phi_l[[i+1]] <- phi_l[[i+1]] + array(sapply(1:niter, function(x) beta[,,2, x] %*% phi_l[[i-1]][,,x]), dim=c(24,24,niter))
    phi_l[[i+1]] <- phi_l[[i+1]] + array(sapply(1:niter, function(x) beta[,,3, x] %*% phi_l[[i-2]][,,x]), dim=c(24,24,niter))
    phi_l[[i+1]] <- phi_l[[i+1]] + array(apply(beta[,,4,], 3, function(x) return(x %*% phi_l[[1]])), dim=c(24, 24, niter))
    # phi_l[[i+1]] <- array(sapply(1:niter, function(x) phi_l[[i+1]][,,x] %*% A[,,x]), dim=c(24,24,niter))
    phi_l[[i+1]] <- array(sapply(1:niter, function(x) phi_l[[i+1]][,,x] %*% t(A[,,x])), dim=c(24,24,niter))
    
  } else {
    phi_l[[i+1]] <- array(sapply(1:niter, function(x) beta[,,1, x] %*% phi_l[[i]][,,x]), dim=c(24,24,niter))
    phi_l[[i+1]] <- phi_l[[i+1]] + array(sapply(1:niter, function(x) beta[,,2, x] %*% phi_l[[i-1]][,,x]), dim=c(24,24,niter))
    phi_l[[i+1]] <- phi_l[[i+1]] + array(sapply(1:niter, function(x) beta[,,3, x] %*% phi_l[[i-2]][,,x]), dim=c(24,24,niter))
    phi_l[[i+1]] <- phi_l[[i+1]] + array(sapply(1:niter, function(x) beta[,,4, x] %*% phi_l[[i-3]][,,x]), dim=c(24,24,niter))
    # phi_l[[i+1]] <- array(sapply(1:niter, function(x) phi_l[[i+1]][,,x] %*% A[,,x]), dim=c(24,24,niter))
    phi_l[[i+1]] <- array(sapply(1:niter, function(x) phi_l[[i+1]][,,x] %*% t(A[,,x])), dim=c(24,24,niter))
  }
   
}
phi_l <- lapply(phi_l, abs)
cumulative_phi <- Reduce(function(x, y) x + y, phi_l[-1], accumulate = TRUE)
theta <- lapply(cumulative_phi, function(y) array(sapply(1:niter, function(x) y[,,x] %*% A[,,x]), dim=c(24,24,niter)))

FF <- 1
Unemployments <- c(20,21,22,23,24)
Wealth <- c(2,3)
resps <- c(Unemployments, Wealth)

IRF_list <- list()

retrieve_impulses <- function(obj, imp, resp){
  return(t(sapply(obj, function(x) return(x[imp, resp, ]))))
}
retrieve_impulses_alt <- function(obj, imp, resp){
  return(t(sapply(obj, function(x) return(x[resp, imp, ]))))
}

for(i in 1:length(resps)){
  name <- actual_vars[resps[i]]
  IRF_list[[name]] <- retrieve_impulses(theta, FF, resps[i])
}

IRF_list_alt <- list()
for(i in 1:length(resps)){
  name <- actual_vars[resps[i]]
  IRF_list_alt[[name]] <- retrieve_impulses_alt(theta, FF, resps[i])
}

plot_IRFs <- function(irf_mat, imp_name, resp_name, lower=0.16, upper=0.84){
  medians <- apply(irf_mat, 1, median)
  lower_bound <- apply(irf_mat, 1, function(x) quantile(x, lower))
  upper_bound <- apply(irf_mat, 1, function(x) quantile(x, upper))
  time <- 1:nrow(irf_mat)
  plot(time, medians, type = "l", col = "black", lwd = 2, 
       ylim = range(c(lower_bound, upper_bound)), 
       ylab = resp_name, xlab = "Horizon", main = paste(resp_name, "impulse response to", imp_name, "shock"))
  
  # Add the shaded confidence interval
  polygon(c(time, rev(time)), c(lower_bound, rev(upper_bound)), 
          col = rgb(0, 0, 1, 0.2), border = NA)
  
  # Optionally, add the confidence interval lines
  lines(time, lower_bound, col = "red", lty = 2)
  lines(time, upper_bound, col = "red", lty = 2)
}

resp_names <- actual_vars[resps]


png("NG2500-FEDFUNDS-IRF.png", width=1100, height=1600)
par(mfrow=c(7,1))
for(i in resp_names){
  plot_IRFs(IRF_list[[i]], "FEDFUNDS", i)
}
dev.off()

png("CC-HS-FEDFUNDS-IRF.png", width=1100, height=1600)
par(mfrow=c(7,1))
for(i in resp_names){
  plot_IRFs(IRF_list[[i]], "FEDFUNDS", i)
}
dev.off()

# plot_IRFs(IRF_list[[1]], "FEDFUNDS", names(IRF_list)[1])

load("/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/shrinkTVPVAR-NG-SV-SDscale.RData")


load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code/all_models_new_new.Rdata")
name <- "VAR.SV.HS"
p <- 4
t <- nrow(y.train) - p
K <- ncol(y.train)
pK <- p*K
niter <- 2500
print(paste("loaded:", model_names[[i]]))

# temp_var <- list(pred_objs=list(final_A = temp_var$pred_objs$final_A,
#                                 theta_sr_SIGMA = temp_var$pred_objs$theta_sr_SIGMA,
#                                 sv_sigma2 = temp_var$pred_objs$sv_sigma2),
#                  beta = temp_var$beta, beta_const = temp_var$beta_const, theta_sr = temp_var$theta_sr, Sigma = temp_var$Sigma)
# gc()
# temp_A0 <- matrix(rep(temp_var$pred_objs$final_A, t), nrow=(t*K^2), ncol=niter)
# temp_A0_err <- matrix(temp_var$pred_objs$theta_sr_SIGMA, nrow=(K^2), ncol=niter)
# print(sum(is.na(temp_A0)))
# print(sum(is.na(temp_A0_err)))
# temp_var$pred_objs$final_A <- NULL
# temp_var$pred_objs$theta_sr_SIGMA <- NULL
# gc()
# print(dim(temp_A0))
# temp_A <- matrix(temp_var$beta, nrow=(K^2 * p * t), ncol=niter)
# temp_A_err <- matrix(temp_var$theta_sr[, -((p*K)+1), ], nrow=(p*(K^2)), ncol=niter)
# print(sum(is.na(temp_A)))
# print(sum(is.na(temp_A_err)))
# temp_var$beta <- NULL
# print(dim(temp_A))
# temp_C <- matrix(temp_var$beta_const, nrow=(K*t), ncol=niter)
# temp_C_err <- matrix((temp_var$theta_sr^2)[, ((p*K) + 1), ], ncol=niter)
# print(sum(is.na(temp_C)))
# print(sum(is.na(temp_C_err)))
# temp_var$beta_const <- NULL
# temp_var$theta_sr <- NULL
# print(dim(temp_C))
# temp_S <- matrix(temp_var$Sigma, nrow=(K^2 * t), ncol=niter)
# temp_S_err <- matrix(to_3D(temp_var$pred_objs$sv_sigma2), nrow=(K^2), ncol=niter)
# print(sum(is.na(temp_C)))
# print(sum(is.na(temp_C_err)))
# # browser
# rm(temp_var)
# stvpvar_bvar <- bvartools::bvar(y=ts(y.train[-c(1:p), ]),
#                                 A0 = list(coeffs=temp_A0, errors=temp_A0_err),
#                                 # A0 = list(coeffs=matrix(array(rep(temp_var$pred_objs$final_A, t), dim = c(K, K, t, niter)), nrow=(t*K^2), ncol=niter),
#                                 #           errors = matrix(temp_var$pred_objs$theta_sr_SIGMA, nrow=(K^2), ncol=niter)), # A0 =matrix(temp_var$pred_objs$final_A, nrow=(K^2), ncol=niter),
#                                 A = list(coeffs=temp_A, errors=temp_A_err),
#                                 C = list(coeffs=temp_C, errors=temp_C_err),
#                                 Sigma = list(coeffs = temp_S,errors = temp_S_err)
# )
# print("removing")
# # rm(temp_var)
# rm(temp_A, temp_A_err, temp_A0, temp_A0_err, temp_C, temp_C_err, temp_S, temp_S_err)
# gc()
# export_irfs_1by1(stvpvar_bvar, imp, actual_vars, "NG2500", save_path, n=20)
# export_irfs_short(stvpvar_bvar, imp, resp_short, "NG2500", save_path2, n=20) 
# rm(stvpvar_bvar)
# gc()
# 
# rm(temp_A, temp_A_err, temp_A0, temp_A0_err, temp_C, temp_C_err, temp_S, temp_S_err)
# gc()
# rm(temp_var)
# gc()
# stvpvar_bvar <- bvartools::bvar(y=ts(y.train[-c(1:p), ]),
#                                 A0 = list(coeffs=matrix(array(rep(temp_var$pred_objs$final_A, t), dim = c(K, K, t, niter)), nrow=(t*K^2), ncol=niter),
#                                           errors = matrix(temp_var$pred_objs$theta_sr_SIGMA, nrow=(K^2), ncol=niter)), # A0 =matrix(temp_var$pred_objs$final_A, nrow=(K^2), ncol=niter),
#                                 A = list(coeffs=matrix(temp_var$beta, nrow=(K^2 * p * t), ncol=niter),
#                                          errors=matrix(temp_var$theta_sr[, -97, ], nrow=(p*(K^2)), ncol=niter)),
#                                 C = list(coeffs=matrix(temp_var$beta_const, nrow=(K*t), ncol=niter),
#                                          errors=matrix((temp_var$theta_sr^2)[, pK+1, ], ncol=niter)),
#                                 Sigma = list(coeffs = matrix(temp_var$Sigma, nrow=(K^2 * t), ncol=niter),
#                                              errors = matrix(to_3D(temp_var$pred_objs$sv_sigma2), nrow=(K^2), ncol=niter))
# )
