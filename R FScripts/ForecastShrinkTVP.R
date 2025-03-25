forecast_shrinkTVP <- function(mod, newdata, n.ahead){

  # Check if mod is of class shrinkTVP
  if (!inherits(mod, "shrinkTVP")) {
    stop("mod has to be an object of class shrinkTVP")
  }

  if (!missing(n.ahead)) {
    if (int_input_bad(n.ahead) || n.ahead == 0) {
      stop("n.ahead has to be a single, positive integer")
    }
  }

  # Check x is data frame and that either newdata or n.ahead is supplied
  if (!missing(newdata)){
    if (is.data.frame(newdata) == FALSE){
      stop("newdata has to be a data frame")
    }

    if (missing(n.ahead)){
      n.ahead <- nrow(newdata)
    } else {
      if (n.ahead > nrow(newdata)) {
        stop("n.ahead can not be larger than the number of rows in newdata")
      }
    }
  } else {
    # Create a dummy newdata data frame, so the formula interface knows the dimensions
    newdata <- data.frame(dummy = rep(NA, n.ahead))
  }

  terms <- delete.response(mod$model$terms)
  m <- model.frame(terms, data = newdata, xlev = mod$model$xlevels)
  x_test <- model.matrix(terms, m)

  # Number of saved posterior draws, number of covariates, number of observations, number of ar coefficients
  nsamp <- nrow(mod$beta_mean)
  N <- length(mod$model$y)
  d <- ncol(mod$beta_mean)
  p <- attr(mod, "p")


  # Get final beta
  if(is.list(mod$beta)){
    beta_prev <- sapply(mod$beta, function(x){
      return(x[, ncol(x)])
    })
  } else {
    beta_prev <- as.matrix(mod$beta[, ncol(mod$beta)])
  }

  # If dynamic (and not iid), get final lambda and other hyperparameters
  if ("shrinkTVP_dyn" %in% class(mod)){

    a_psi <- attr(mod, "a_psi")
    c_psi <- attr(mod, "c_psi")

    if (!attr(mod, "iid")) {
      if(is.list(mod$lambda_p)){
        lambda_prev <- sapply(mod$lambda_p, function(x){
          return(x[, ncol(x)])
        })
      } else {
        lambda_prev <- as.matrix(mod$lambda_p[, ncol(mod$lambda_p)])
      }
    }
  }

  # Storage mods
  beta_pred <- matrix(NA, nrow = nsamp, ncol = d)
  y_pred <- as.data.frame(matrix(0, ncol = n.ahead, nrow = nsamp))

  if (attr(mod, "sv") == TRUE){
    ht_prev <- log(mod$sigma2[, ncol(mod$sigma2)])
    ht_store <- data.frame(matrix(NA, ncol = n.ahead, nrow = nsamp))
  } else {
    sigma2 <- mod$sigma2
  }

  for (ti in 1:n.ahead){

    for (j in 1:d){

      if ("shrinkTVP_dyn" %in% class(mod)){

        # Forecast psi into the future (except for inter_column if it is not shrunken)
        if (attr(mod, "shrink_inter") == FALSE & attr(mod, "inter_column") == j) {
          psi_pred <- rep(1, nsamp)
        } else {
          if (attr(mod, "iid")) {
            psi_pred <- rf(nsamp, 2*a_psi[j], 2*c_psi[j])
          } else {
            kappa_pred <- rpois(nsamp, lambda_prev[,j] * a_psi[j]/c_psi[j] * mod$rho_p[,j]/(1 - mod$rho_p[,j]))
            lambda_pred <- rgamma(nsamp, a_psi[j] + kappa_pred, a_psi[j]/c_psi[j] * 1/(1 - mod$rho_p[,j]))
            psi_pred <- 1/rgamma(nsamp, c_psi[j], lambda_pred)
            lambda_prev <- lambda_pred
          }
        }
        beta_pred[, j] <- beta_prev[, j] + rnorm(nsamp, 0, abs(mod$theta_sr[, j]) * sqrt(psi_pred))

      } else {
        beta_pred[, j] <- beta_prev[, j] + rnorm(nsamp, 0, abs(mod$theta_sr[, j]))
      }
    }


    if (attr(mod, "sv") == TRUE){
      mean_ht <- mod$sv_mu + mod$sv_phi * (ht_prev - mod$sv_mu)
      ht_pred <- rnorm(nsamp, mean_ht, sqrt(mod$sv_sigma2))
      sigma2 <- exp(ht_pred)

      ht_store[, ti] <- ht_pred
      ht_prev <- ht_pred
    }

    # Add regular regressors (excluding AR components)
    if (! 0 %in% dim(x_test)) {
      y_pred[, ti] <- t(x_test[ti, ] %*% t(beta_pred[, 1:(d - p)]))
    }

    # Add AR components
    if (p != 0) {
      for (curr_p in 1:p){
        # Differentiate between case where the previous y is unknown and where it is known
        curr_prev_t <- ti - curr_p

        if (curr_prev_t <= 0) {
          curr_prev_y <- mod$model$y[length(mod$model$y) + curr_prev_t]
        } else {
          curr_prev_y <- y_pred[, ti - curr_p]
        }

        y_pred[, ti] <- y_pred[, ti] + curr_prev_y * beta_pred[, (d - p) + curr_p]
      }
    }

    # Add error term
    y_pred[, ti] <- y_pred[, ti] + rnorm(nsamp, 0, sqrt(sigma2))

    beta_prev <- beta_pred
  }

  res <- list(y_pred = y_pred, y_orig = mod$model$y)

  if (attr(mod, "sv") == TRUE){
    res$n.ahead <- ht_store
  }

  class(res) <- c("shrinkTVP_forc")
  attr(res, "index") <- attr(mod, "index")

  return(res)
}