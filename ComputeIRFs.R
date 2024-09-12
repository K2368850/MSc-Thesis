{# specify_both_irf <- function(coefficients, impulse, response, horizon){
#   # compute_irf <- function(coefficients, impulse, response, horizon){
#   #   n_vars <- dim(coefficients[[1]])[1]
#   #   irf_matrix <- matrix(0, nrow=horizon, ncol=nvars)
#   #   
#   #   shock <- rep(0, n_vars)
#   #   shock[impulse] <- 1
#   #   
#   #   for (h in 1:horizon){
#   #     if (h==1){
#   #       irf_matrix[h, ] <- coefficient[[h]] %*% shock
#   #     } else {
#   #       irf_matrix[h, ] <- coefficients[[h]] %*% irf_matrix[h-1, ]
#   #     }
#   #   }
#   #   return (irf_matrix[, response])
#   # }
#   # 
#   # 
#   # horizon <- 100
#   # impulse <- 29 #index of impulse
#   # response <- 2 #index of response
#   # 
#   # irfs <- lapply(1:length(coefficients), function(t){
#   #   compute_irf(coefficients[[t]], impulse, response, horizon)}
#   #   )
#   # 
#   # # Convert list of IRFs to a matrix
#   # irf_matrix <- do.call(rbind, irfs)
#   # 
#   # # Plotting the IRFs over time
#   # plot_irf <- function(irf_matrix, horizon) {
#   #   matplot(irf_matrix, type = "l", lty = 1, col = 1:ncol(irf_matrix), 
#   #           ylab = "Response", xlab = "Horizon", main = "Impulse Response Functions")
#   #   legend("topright", legend = paste0("Time Period ", 1:horizon), col = 1:horizon, lty = 1)
#   # }
#   # 
#   # plot_irf(irf_matrix, horizon)
# }
# 
# 
# # Function to compute IRFs for all variables given an impulse to one variable
# compute_irf_all <- function(coefficients, impulse, horizon) {
#   n_vars <- dim(coefficients[[1]])[2]  # Adjust based on the structure of coefficients
#   irf_matrix <- matrix(0, nrow = horizon, ncol = n_vars)
#   
#   # Initial shock to one variable
#   shock <- rep(0, n_vars)
#   shock[impulse] <- 1
#   
#   for (h in 1:horizon) {
#     if (h == 1) {
#       irf_matrix[h, ] <- coefficients[[h]] %*% shock
#     } else {
#       irf_matrix[h, ] <- coefficients[[h]] %*% irf_matrix[h-1, ]
#     }
#   }
#   
#   return(irf_matrix)
# }
# 
# # Example usage:
# horizon <- 10  # Specify how many steps ahead you want
# impulse <- 1   # Index of the variable that receives the impulse
# 
# # Initialize a list to store IRFs for all time periods
# irf_results <- lapply(1:length(coefficients), function(t) {
#   compute_irf_all(coefficients[[t]], impulse, horizon)
# })
# 
# # Convert list of IRFs to a 3D array for easier handling
# irf_array <- array(unlist(irf_results), dim = c(horizon, n_vars, length(coefficients)))
# 
# # Plot the IRFs for all variables in response to the impulse
# plot_irf_all <- function(irf_array, impulse_var_name) {
#   par(mfrow = c(n_vars, 1))  # Create a grid of plots
#   for (var_idx in 1:n_vars) {
#     matplot(1:horizon, irf_array[, var_idx, ], type = "l", lty = 1, col = 1:length(coefficients),
#             ylab = paste0("Response of Variable ", var_idx), xlab = "Horizon",
#             main = paste0("Impulse Response to a Shock in ", impulse_var_name))
#     legend("topright", legend = paste0("Time Period ", 1:length(coefficients)), col = 1:length(coefficients), lty = 1)
#   }
# }
# 
# # Plot the results (assuming variable names are available)
# plot_irf_all(irf_array, impulse_var_name = "Variable 1")
}

### Better IRF

coefficients_matrix <- as.matrix(big_boy_model_sv$beta_mean)

### Function to compute IRFs for all variables given an impulse to one variable
compute_irf_all <- function(coefficients, impulse, horizon) {
  n_vars <- dim(coefficients[[1]])[1]
  irf_matrix <- matrix(0, nrow = horizon, ncol = n_vars)
  
  # Initial shock to one variable
  shock <- rep(0, n_vars)
  shock[impulse] <- 1
  
  for (h in 1:horizon) {
    if (h == 1) {
      irf_matrix[h, ] <- coefficients[[h]] %*% shock
    } else {
      irf_matrix[h, ] <- coefficients[[h]] %*% irf_matrix[h-1, ]
    }
  }
  
  return(irf_matrix)
}

### Example usage:
### Assuming `coefficients` is the list of time-varying coefficients
horizon <- 10  # Specify how many steps ahead you want
impulse <- 29   # Index of the variable that receives the impulse

### Initialize a list to store IRFs for all time periods
irf_results <- lapply(1:length(coefficients_matrix), function(t) {
  compute_irf_all(coefficients_matrix[[t]], impulse, horizon)
})

### Convert list of IRFs to a 3D array for easier handling
irf_array <- array(unlist(irf_results), dim = c(horizon, n_vars, length(coefficients)))

### Plot the IRFs for all variables in response to the impulse
plot_irf_all <- function(irf_array, impulse_var_name) {
  par(mfrow = c(n_vars, 1))  # Create a grid of plots
  for (var_idx in 1:n_vars) {
    matplot(1:horizon, irf_array[, var_idx, ], type = "l", lty = 1, col = 1:length(coefficients),
            ylab = paste0("Response of Variable ", var_idx), xlab = "Horizon",
            main = paste0("Impulse Response to a Shock in ", impulse_var_name))
    legend("topright", legend = paste0("Time Period ", 1:length(coefficients)), col = 1:length(coefficients), lty = 1)
  }
}

### Plot the results (assuming variable names are available)
plot_irf_all(irf_array, impulse_var_name = "Variable 1")

# model <- shrinkTVP(PCEPI ~ ., data=data_macro)
# model_sv <- shrinkTVP(PCEPI ~ ., data=data_macro, sv=TRUE)
# model_ngg_sv <- shrinkTVP(PCEPI + UNRATE~. )
# 
# rm(model, model_sv, big_model, big_model_sv)


