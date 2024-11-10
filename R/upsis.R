#' Update a model using PSIS
#' @param model a model object (currently only brmsfit supported)
#' @param data_add additional data collected after the model was fit
#' @param data_remove data included in model fitting whose influence to remove
#' @export
upsis <- function (model, data_add = NULL, data_remove = NULL) {
  UseMethod("upsis")
}

#' Update a brms model using PSIS
#' @param model a brmsfit object
#' @param data_add additional data collected after the model was fit
#' @param data_remove data included in model fitting whose influence to remove
upsis.brmsfit <- function (model, data_add = NULL, data_remove = NULL) {
  if(is.null(data_add) & is.null(data_remove)){
    warning("no data added or removed; returning original model object")
    return(model)
  }
  ll_add <- ll_remove <- rep(0, brms::ndraws(model))
  if(!is.null(data_add)) {
    ll_add <- brms::log_lik(model, newdata = data_add) |>
      rowSums()
  }
  if(!is.null(data_remove)) {
    ll_remove <- brms::log_lik(model, newdata = data_remove) |>
      rowSums()
  }
  ll_change <- ll_add - ll_remove
  
  message("performing PSIS")
  psis_out <- loo::psis(ll_change)
  pareto_k <- psis_out$diagnostics$pareto_k
  print(paste0("Pareto k is ", pareto_k))
  if(pareto_k > .7) {
    warning(
      paste0(
        "Pareto k is greater than 0.7, indicating that the PSIS",
        "approximation to the updated posterior may be unreliable"
        )
    )
  }
  draw_weights <- loo::weights.importance_sampling(psis_out, log = FALSE)
  
  updated_model <- model
  updated_model$fit <- reweight_stanfit_draws(model$fit, draw_weights)
  
  out <- list(updated_model = updated_model, psis = psis_out)
  out
}

#' Update post-warmup draws of a stanfit based on a vector of weights
#' @param stanfit a stanfit
#' @param weights a vector of weights
reweight_stanfit_draws <- function(stanfit, draw_weights) {
  stanargs <- stanfit@stan_args
  n_chains <- length(stanargs)
  draws_per_chain <- length(draw_weights) / n_chains
  assertthat::assert_that(draws_per_chain == round(draws_per_chain))
  
  draw_weight_mat <- matrix(draw_weights, nrow = draws_per_chain)
  
  iter <- lapply(stanargs, function(x){x$iter}) |>
    unlist() |>
    unique()
  assertthat::assert_that(length(iter) == 1)
  
  warmup <- lapply(stanargs, function(x){x$warmup}) |>
    unlist() |>
    unique()
  assertthat::assert_that(length(warmup) == 1)
  
  thin <- lapply(stanargs, function(x){x$thin}) |>
    unlist() |>
    unique()  
  assertthat::assert_that(length(thin) == 1)
  
  draws_per_chain_2 <- (iter - warmup) / thin
  assertthat::assert_that(draws_per_chain == draws_per_chain_2)
  
  samples <- stanfit@sim$samples
  assertthat::assert_that(length(samples) == n_chains)
  
  for (i in 1:n_chains) {
    si <- samples[[i]]
    wi <- draw_weight_mat[,i]
    ids <- sample(
      1:draws_per_chain, 
      draws_per_chain, 
      replace = TRUE,
      prob = wi
      )
    for(j in 1:length(si)) {
      si_length <- length(si[[j]])
      si[[j]][(1 + si_length - draws_per_chain) : si_length] <- 
        si[[j]][(1 + si_length - draws_per_chain) : si_length][ids]
      }
    stanfit@sim$samples[[i]] <- si
  }
  
  stanfit
}
