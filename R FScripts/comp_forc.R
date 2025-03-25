# library(vars)
library(forecast)
# library(bvartools)
library(Metrics)
library(bayesianVARs)
library(shrinkTVPVAR)
library(dplyr)

setwd("/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/Forecasts/")
load("/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/macro.Rdata")
load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code/all_models_new_new.Rdata")
# load("macro.Rdata")

multiply_by_sd <- function(values, sd){
  return(sweep(values, 2, sd, "*"))
}
add_by_sd <- function(values, sd){
  return(sweep(values, 2, sd, "+"))
}
split_data <- function(data){
  x.train <- head(data, -20)
  x.test <- tail(data, 20)
  
  return (list(train = x.train, test=x.test))
}

data <- macro[, -1]

prior_sd <- apply(data, 2, sd)

# scaled_uncentered <- scale_w_center(data)

data_split <- split_data(scale(data, center=FALSE))
y.train <- data_split$train
y.test <- data_split$test

results <- list()
horizons <- c(2, 3, 4, 8, 12, 16, 20)

cc <- all_models[["Constant"]]
cc_fit <- fitted(cc)
cc_rmse <- accuracy(as.numeric(multiply_by_sd(cc_fit, prior_sd)), multiply_by_sd(y.train[-c(1:4), ], prior_sd))
for(h in horizons){
  cc_pred <- predict(cc, n.ahead=h, ci=0.68)
  cc_pred <- sapply(names(cc_pred$fcst), function(x) return(cc_pred$fcst[[x]][, 1]))
  cc_rmfse <- forecast::accuracy(as.numeric(multiply_by_sd(cc_pred, prior_sd)), multiply_by_sd(y.test[1:h, ], prior_sd))
  print(cc_rmfse)
  results[[paste0("CC_jrmfse", h)]] <- cc_rmfse
}
for(h in c(1, horizons)){
  cc_pred <- predict(cc, n.ahead=h, ci=0.68)
  cc_pred <- sapply(names(cc_pred$fcst), function(x) return(cc_pred$fcst[[x]][h, 1]))
  # cc_rmfse <- 0
  # if (h==1){
  cc_rmfse <- rmse(cc_pred*prior_sd, y.test[h, ]*prior_sd)
  print(cc_rmfse)
  # } else {
    # cc_rmfse <- rmse(cc)
  # }
  results[[paste0("CC_rmsfe", h)]] <- cc_rmfse
}
# cc_pred <- predict(cc, n.ahead=20, ci=0.68)
# cc_pred <- sapply(names(cc_pred$fcst), function(x) return(cc_pred$fcst[[x]][, 1]))

# cc_rmfse <- accuracy(as.numeric(multiply_by_sd(cc_pred, prior_sd)), multiply_by_sd(y.test, prior_sd))

results[["CC_rmse"]] <- cc_rmse
results[["CC_rmfse"]] <- cc_rmfse


all_models[["VAR.FSV.N"]] <- NULL
all_models[["VAR.SV.N"]] <- NULL
bayesianVARs_names <- names(all_models)[sapply(all_models, function(x) class(x) == "bayesianVARs_bvar")]
bayesianVARs_names <- bayesianVARs_names[-c(1, 5)]
for(i in 1:length(bayesianVARs_names)){
  name <- bayesianVARs_names[i]
  print(name)
  name_rmse <- paste0(name, "_rmse")
  name_rmsfe <- paste0(name, "_rmsfe")
  name_jrmsfe <- paste0(name, "_jrmsfe")
  name_LPL <- paste0(name, "_LPL")
  m <- all_models[[name]]
  m_fit <- fitted(m, error_term=FALSE)
  m_point_fit <- apply(m_fit$fitted, c(1, 2), mean)
  m_rmse <- accuracy(as.numeric(multiply_by_sd(m_point_fit, prior_sd)), multiply_by_sd(y.train[-c(1:4), ], prior_sd))
  
  tryCatch({
    m_LPL <- rowMeans(predict(m, ahead=1:20, Y_obs=y.test, LPL=TRUE, stable=FALSE)$LPL_draws)
    results[[name_LPL]] <- m_LPL}, error=function(e){print("chol decomp fail")}, finally={print("ok")})
  
  m_forc <- predict(m, ahead=1:20, Y_obs=y.test, stable=FALSE)
  m_forc <- apply(m_forc$predictions, c(1,2), mean)
  
  for(h in horizons){
    results[[paste0(name_jrmsfe, h)]] <- forecast::accuracy(as.numeric(multiply_by_sd(m_forc[1:h, ], prior_sd)), multiply_by_sd(y.test[1:h, ], prior_sd))
  }
  for(h in c(1, horizons)){
    results[[paste0(name_rmsfe, h)]] <- rmse(prior_sd*m_forc[h, ], prior_sd*y.test[h, ])
  }
  
  # m_rmsfe <- accuracy(as.numeric(m_forc), y.test)
  
  results[[name_rmse]] <- m_rmse
  # results[[name_rmsfe]] <- m_rmsfe
  
}#, finally = {print("oops?")})

for(i in 1:length(bayesianVARs_names)){
  name <- bayesianVARs_names[i]
  name_rmse <- paste0(name, "_rmse")
  m <- all_models[[name]]
  m_fit <- fitted(m, error_term=FALSE)$fitted
  m_fit <- sapply(1:(dim(m_fit))[3], function(x) rowMeans(multiply_by_sd(m_fit[,,x]-y.train[-c(1:4), ], prior_sd)^2))
  m_rmse <- mean(m_fit)
  results[[name_rmse]] <- m_rmse
}

# for(i in 1:length(bayesianVARs_names)){
#   name <- bayesianVARs_names[i]
name_rmse <- paste0(name, "_rmse")
m <- all_models[[name]]
m_fit <- fitted(m, error_term=FALSE)$fitted
m_fit <- sapply(1:(dim(m_fit))[3], function(x) rowMeans(multiply_by_sd(m_fit[,,x]-y.train[-c(1:4), ], prior_sd)^2))
m_rmse <- mean(m_fit)
results[[name_rmse]] <- m_rmse
# }

setwd("/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/Forecasts/")
# load()
save.image("forecast_res_w_scale.Rdata")

save(results, file="/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/results_list.Rdata")
save.image("halfway_done_forc.Rdata")
# load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/R Results/Forecasts/halfway_done_forc.Rdata")

load("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/Data 07-05/Code/all_models_new_new.Rdata")


name <- "VAR.SV.SSVS"
name_rmse <- paste0(name, "_rmse")
name_rmsfe <- paste0(name, "_rmsfe")
name_jrmsfe <- paste0(name, "_jrmsfe")
name_LPL <- paste0(name, "_LPL")
m <- all_models[[name]]
m_fit <- fitted(m, error_term=TRUE)
m_point_fit <- apply(m_fit$fitted, c(1, 2), mean)
m_rmse <- accuracy(as.numeric(m_point_fit), y.train[-c(1:4), ])

m_LPL <- predict(m, ahead=1:20, Y_obs=y.test, LPL=TRUE, stable=FALSE)$LPL_draws
m_LPL <- rowMeans(m_LPL)

m_forc <- predict(m, ahead=1:20, Y_obs=y.test, stable=FALSE)
m_forc <- apply(m_forc$predictions, c(1,2), mean)

for(h in horizons){
  results[[paste0(name_jrmsfe, h)]] <- forecast::accuracy(as.numeric(multiply_by_sd(m_forc[1:h, ], prior_sd)), multiply_by_sd(y.test[1:h, ], prior_sd))
}
for(h in c(1, horizons)){
  results[[paste0(name_rmsfe, h)]] <- rmse(prior_sd*m_forc[h, ], prior_sd*y.test[h, ])
}

# m_rmsfe <- accuracy(as.numeric(m_forc), y.test)

results[[name_rmse]] <- m_rmse
# results[[name_rmsfe]] <- m_rmsfe
results[[name_LPL]] <- m_LPL


name <- "VAR.FSV.N"
name_rmse <- paste0(name, "_rmse")
name_rmsfe <- paste0(name, "_rmsfe")
name_jrmsfe <- paste0(name, "_jrmsfe")
name_LPL <- paste0(name, "_LPL")
m <- all_models[[name]]
m_fit <- fitted(m, error_term=FALSE)
m_point_fit <- apply(m_fit$fitted, c(1, 2), mean)
m_rmse <- accuracy(as.numeric(m_point_fit), y.train[-c(1:4), ])

m_LPL <- predict(m, ahead=1:20, Y_obs=y.test, LPL=TRUE, stable=FALSE)$LPL_draws
m_LPL <- rowMeans(m_LPL)

m_forc <- predict(m, ahead=1:20, Y_obs=y.test, stable=FALSE)
m_forc <- apply(m_forc$predictions, c(1,2), mean)

for(h in horizons){
  results[[paste0(name_jrmsfe, h)]] <- forecast::accuracy(as.numeric(multiply_by_sd(m_forc[1:h, ], prior_sd)), multiply_by_sd(y.test[1:h, ], prior_sd))
}
for(h in c(1, horizons)){
  results[[paste0(name_rmsfe, h)]] <- rmse(prior_sd*m_forc[h, ], prior_sd*y.test[h, ])
}

# m_rmsfe <- accuracy(as.numeric(m_forc), y.test)

results[[name_rmse]] <- m_rmse
# results[[name_rmsfe]] <- m_rmsfe
results[[name_LPL]] <- m_LPL

save(results, file="forecast_res_only_sddata_newest.Rdata")
save.image("forecast_environment_newest.Rdata")

model_files_loc <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/"
# setwd(model_files_loc)
setwd("/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/R Results/Forecast/")
files <- grep("10000", list.files(), value=TRUE, ignore.case=FALSE)
models_file_names <- paste0(model_files_loc, files)

save_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King'sCollegeLondon/Documents/KCL Stats KCL/MSc Thesis-Dissertation/R Results/Forecast/"
density_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/Images/Density/"
state_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/Images/State/"
# heatmap_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/Images/Heatmap/"
heatmap_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/Images\ 10000/Heatmap/"
state_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/Images\ 10000/State/"
density_path <- "/Users/geraldpress/Library/CloudStorage/OneDrive-King\'sCollegeLondon/Documents/KCL\ Stats\ KCL/MSc\ Thesis-Dissertation/R\ Results/Images\ 10000/Density/"
compute_LPDS_from_forecast <- function(actual, all_forecasts, old_sd, rescale=FALSE){
  prior_means <-sapply(names(all_forecasts), function(x) return(colMeans(all_forecasts[[x]]$y_pred)))
  if(rescale & !is.null(old_sd)) prior_means <- multiply_by_sd(prior_means, old_sd)
  prior_variances <- sapply(names(all_forecasts), function(x) return(mat2var(all_forecasts[[x]]$y_pred)))
  if(rescale & !is.null(old_sd)) prior_variances <- add_by_sd(prior_variances, old_sd^2)#prior_variances <- multiply_by_sd(prior_variances, old_sd^2)
  return (approx_func(prior_means, prior_variances, actual))
}

approx_func <- function(mean, var, actual){
  result <- log(2*pi) + log(var) + (((actual-mean)^2)/(var))
  return (-0.5*result)
}
mat2var <- function(mat){
  return(apply(mat, 2, sd)^2)
}

model_names <- sapply(files, function(x) substr(x, 1, nchar(x) - 6))

to_plot <- c("a_tau", "a_xi", "beta", "beta_const", "beta_mean", "c_tau", "c_xi",
             "consts", "kappa2", "kappa2_B", "lambda2", "lambda2_B", "tau2", "theta_sr", "xi2")
t1 <- Sys.time()
for(i in 1:length(models_file_names)){
  name <- model_names[[i]]
  name_rmse <- paste0(name, "_rmse")
  name_rmsfe <- paste0(name, "_rmsfe")
  name_jrmsfe <- paste0(name, "_jrmsfe")
  load(models_file_names[i])
  for(p in to_plot){
    png(filename=paste0(density_path, name, "density_plot_", p, ".png"), width = 4000, height = 3000)
    density_plotter(temp_var, to_plot=p)
    dev.off()
  }
  png(filename=paste0(state_path, name, "state_plot", ".png"), width=4000, height=3000)
  state_plotter(temp_var, Sigma=FALSE)
  dev.off()
  png(filename=paste0(state_path, name, "state_plot_Sigma", ".png"), width=4000, height=3000)
  state_plotter(temp_var, Sigma=TRUE)
  dev.off()
  for(p in to_plot){
    png(filename=paste0(heatmap_path, name, "heatmap_plot_", p, ".png"), width = 4000, height = 3000)
    TV_heatmap(temp_var, to_plot=p)
    dev.off()
  }
  
  fcasts <- forecast_shrinkTVPVAR(temp_var, n.ahead=20)
  temp_fit <- fitted(temp_var)
  results[[paste0(name, "_LPL")]] <- compute_LPDS_from_forecast(y.test, fcasts, prior_sd, rescale=FALSE)
  rm(temp_var)
  gc()
  fcasts <- sapply(names(fcasts), function(x) return(colMeans(fcasts[[x]]$y_pred)))
  for(h in horizons){
    results[[paste0(name_jrmsfe, h)]] <- forecast::accuracy(as.numeric(multiply_by_sd(fcasts[1:h, ], prior_sd)), multiply_by_sd(y.test[1:h, ], prior_sd))
  }
  for(h in c(1, horizons)){
    results[[paste0(name_rmsfe, h)]] <- rmse(fcasts[h, ]*prior_sd, y.test[h, ]*prior_sd)
  }
  results[[paste0(name_rmse, h)]] <- forecast::accuracy(as.numeric(multiply_by_sd(temp_fit, prior_sd)), multiply_by_sd(y.train, prior_sd))
  
}
t2 <- Sys.time()
save(results, file="forecast_res_all_except_svss.Rdata")


rmse_names <- grep("rmse", names(results), value=TRUE)
jrmsfe_names <- grep("jrmsfe", names(results), value=TRUE)

jrmsfe_vals <- sapply(jrmsfe_names, function(x) return(results[[x]][2]))
# load("shrinkTVPVAR-Hier-Lasso-SDscale.Rdata")
# setwd(model_files_loc)

library(corrplot)
fsv_names <- grep(".FSV.", names(all_models), value=TRUE)[-4]
real_names <- c("CC-FSV-VAR-HS", "CC-FSV-VAR", "CC-FSV-VAR-NG")
fsv12_list <- list()
fsv4_list <- list()
for(i in 1:length(fsv_names)){
  m <- all_models[[fsv_names[i]]]
  # temp_fit_noerr <- fitted(m, error_term=FALSE)
  # temp_fit_err <- fitted(m, error_term=TRUE)
  # temp_fit <- temp_fit_err$fitted - temp_fit_noerr$fitted
  # avg_fit <- apply(temp_fit, c(1,2), mean)
  # fsv12_list[[fsv_names[i]]] <- fsvsample(avg_fit, factors=12)
  # fsv4_list[[fsv_names[i]]] <- fsvsample(avg_fit, factors=4)
  median_fac <- apply(m$facload, c(1,2), median)
  rownames(median_fac) <- actual_vars
  standard_median_fac <- scale(median_fac)
  corrplot(standard_median_fac, method="circle", is.corr=FALSE) #title=paste0("Standardized Posterior Median of Factor Loadings for ", real_names[i])
}

for(i in 1:length(fsv_names)){
  m <- all_models[[fsv_names[i]]]
  print(fsv_names[i])
  # temp_fit_noerr <- fitted(m, error_term=FALSE)
  # temp_fit_err <- fitted(m, error_term=TRUE)
  # temp_fit <- temp_fit_err$fitted - temp_fit_noerr$fitted
  # avg_fit <- apply(temp_fit, c(1,2), mean)
  # fsv12_list[[fsv_names[i]]] <- fsvsample(avg_fit, factors=12)
  # fsv4_list[[fsv_names[i]]] <- fsvsample(avg_fit, factors=4)
  median_fac <- apply(m$facload, c(1,2), median)
  rownames(median_fac) <- actual_vars
  # standard_median_fac <- scale(median_fac)
  standard_median_fac <- apply(median_fac, 2, function(x) x / sum(x))
  corrplot(standard_median_fac, method="circle", is.corr=FALSE) #title=paste0("Standardized Posterior Median of Factor Loadings for ", real_names[i])
}
i <- 1

for(i in 1:length(fsv_names)){
  # browser()
  m <- all_models[[fsv_names[i]]]
  class(m) <- "fsvdraws"
  m[["y"]] <- m[["Y"]]
  factorstochvol::plotalot(m)
}

l <- fsv12_list
for(i in 1:3){
  temp <- l[[i]]
  facload <- temp$facload
  median_fac <- apply(facload, c(1,2), median)
  rownames(median_fac) <- actual_vars
  std_fac <- scale(median_fac)
  corrplot(std_fac, method="circle", is.corr=FALSE, title=paste0("Posterior Median of Standardized 12-Factor SV for", real_names[i]))
}

l <- fsv4_list
for(i in 1:3){
  temp <- l[[i]]
  facload <- temp$facload
  median_fac <- apply(facload, c(1,2), median)
  rownames(median_fac) <- actual_vars
  std_fac <- scale(median_fac)
  corrplot(std_fac, method="circle", is.corr=FALSE, title=paste0("Posterior Median of 4-Factor SV for", real_names[i]))
}

plotalot(fsv12_list[[1]])
plotalot(fsv4_list[[1]])
plotalot(fsv12_list[[2]])
plotalot(fsv4_list[[2]])
plotalot(fsv12_list[[3]])
plotalot(fsv4_list[[3]])


cc <- all_models[["Constant"]]

cc_resid <- residuals(cc)

cc_fsv12 <- fsvsample(cc_resid, factors=4)


LPL_names <- grep("LPL", names(results), value=TRUE)
LPL_tvpvar_names <- grep("SDscale", LPL_names, value=TRUE)
LPL_cc_names <- setdiff(LPL_names, LPL_tvpvar_names)

LPL_avg <- list()

for(i in LPL_tvpvar_names){
  LPL_avg[[i]] <- rowMeans(results[[i]])
}

for(i in LPL_cc_names){
  LPL_avg[[i]] <- results[[i]]
}
