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
  draws_per_chain <- brms::ndraws(model) / brms::nchains(model)
  draw_weights_list <- matrix(draw_weights, nrow = draws_per_chain) |>
    as.data.frame() |>
    as.list()
  
  updated_model <- model
  draws_list <- posterior::as_draws_list(model$fit@sim$samples)
  n_warmup <- (posterior::ndraws(draws_list) - brms::ndraws(model)) / brms::nchains(model)

  updated_model$fit@sim$samples <-
    resample_by_chain(
      draws_list, n_warmup, draw_weights_list
    )
  
  for(i in seq_along(updated_model$fit@sim$samples)) {
    attributes(updated_model$fit@sim$samples[[i]]) <-
      attributes(model$fit@sim$samples[[i]])
  }
  
  out <- list(updated_model = updated_model, psis = psis_out)
  out
}

resample_by_chain <- function(draws_list, n_warmup, weights) {
  if(length(n_warmup) == 1) {
    n_warmup <- rep(n_warmup, posterior::nchains(draws_list))
  }
  assertthat::assert_that("draws_list" %in% class(draws_list))
  mapply(resample_one_chain, draws_list, n_warmup, weights)
}

resample_one_chain <- function(draws_list, n_warmup_draws, weights) {
  draws_list <- posterior::as_draws_list(draws_list)
  assertthat::assert_that(posterior::nchains(draws_list) == 1)
  assertthat::assert_that(length(n_warmup_draws) == 1)
  n_draws <- posterior::ndraws(draws_list)
  assertthat::assert_that(n_draws > n_warmup_draws)
  assertthat::assert_that((n_draws - n_warmup_draws) == length(weights))
  draws_df <- posterior::as_draws_df(draws_list)
  post_warmup_draws <- draws_df[(n_warmup_draws + 1) : n_draws, ]
  resampled_draws <- posterior::resample_draws(
    post_warmup_draws, 
    weights = weights, 
    method = "stratified"
    )
  if(n_warmup_draws > 0) {
    warmup_draws <- draws_df[1 : n_warmup_draws, ]
    resampled_draws <- rbind(warmup_draws, resampled_draws)
  }
  posterior::as_draws_list(resampled_draws)
  
}






