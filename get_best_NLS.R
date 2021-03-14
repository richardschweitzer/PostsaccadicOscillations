# run an optimizer through a grid of different starting parameters and return the best model
# by Richard Schweitzer
get_best_NLS <- function(df, model_form, model_control, 
                         start_params_low, start_params_high, 
                         start_params_n, 
                         absolute_lowest, 
                         use_robust=FALSE, debug_mode=FALSE) {
  require(assertthat)
  require(data.table)
  # check whether params_low and _high are correctly specified
  assert_that(length(start_params_high)==length(start_params_low))
  assert_that(all(names(start_params_high)==names(start_params_low)))
  assert_that(length(start_params_high)==length(absolute_lowest))
  assert_that(all(names(start_params_high)==names(absolute_lowest)))
  # what are the different levels of values for each parameter?
  start_param_list = vector(mode = 'list', length = length(start_params_low))
  for (param_i in 1:length(start_params_high)) {
    start_param_list[[param_i]] <- seq(start_params_low[param_i], start_params_high[param_i], 
                                       length.out = start_params_n)
  }
  # get all combinations for a grid-search approach
  start_params <- expand.grid(start_param_list)
  colnames(start_params) <- names(start_params_high)
  setDT(start_params)
  if (debug_mode) {
    print("Starting parameters:")
    print(start_params)
  }
  # here we'll save the results
  result_params <- copy(start_params) 
  result_params[ , BIC := NaN]
  result_params[ , AIC := NaN]
  n_col <- ncol(result_params)
  # here we'll save all the fitted models
  model_list <- vector(mode = "list", length = nrow(result_params))
  # get the model formula
  model_formula <- as.formula(model_form)
  if (debug_mode) {
    print("Model:")
    print(model_formula)
  }
  # run through
  for (r in as.integer(1:nrow(start_params))) {
    if (debug_mode) {
      print(paste0("r=", r, " , starting params: ", paste(round(start_params[r, ],2), collapse = " ")))
    }
    this_row <- as.list(rep(NaN, n_col))
    this_fit <- NULL
    # go fit
    tryCatch({
      if (!use_robust) {
        require(minpack.lm)
        suppressWarnings(this_fit <- nlsLM(data = df, 
                                           formula = model_formula, 
                                           start = start_params[r, ], 
                                           control = model_control, 
                                           lower = absolute_lowest))
      } else {
        require(robustbase)
        suppressWarnings(this_fit <- nlrob(data = df, 
                                           formula = model_formula, 
                                           start = start_params[r, ], 
                                           lower = absolute_lowest, 
                                           method = "MM"))
      }
      this_row <- as.list(c(as.numeric(coef(this_fit)), BIC(this_fit), AIC(this_fit)))
      if (debug_mode) {
        print(c(as.numeric(coef(this_fit)), BIC(this_fit), AIC(this_fit)))
      }
      model_list[[r]] <- this_fit
    }, error = function(ex) {
      this_row <- as.list(rep(NaN, n_col))
      if (debug_mode) {
        print(ex)
      }
    })
    # put into data.table
    set(result_params, r, names(result_params), this_row)
  }
  # get best BIC
  best_index <- which(result_params$BIC==min(result_params$BIC, na.rm = TRUE))[1]
  # return values
  if (length(best_index)==1) {
    return(list(best_fit = model_list[[best_index]], 
                best_params = result_params[best_index], 
                start_params = start_params,
                result_params = result_params))
  } else {
    return(list(best_fit = NULL, 
                best_params = rep(NaN, length(start_params_low)+2), 
                start_params = start_params,
                result_params = result_params))
  }
  
}