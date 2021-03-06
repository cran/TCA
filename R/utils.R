#' @importFrom matrixStats rowVars
#' @importFrom futile.logger flog.appender
#' @importFrom futile.logger appender.tee
#' @importFrom futile.logger appender.console
#' @importFrom futile.logger flog.threshold
#' @importFrom futile.logger flog.debug
#' @importFrom futile.logger flog.info
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom data.table fread
#' @importFrom data.table fwrite
#' @importFrom parallel stopCluster
#' @importFrom parallel detectCores
#' @importFrom stats lm
#' @importFrom stats logLik
#' @importFrom stats p.adjust
#' @importFrom stats pchisq
#' @importFrom stats residuals
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats sd
start_logger <- function(log_file, debug, verbose = TRUE){
  config_level <- if (debug) "debug" else "default"
  Sys.setenv(R_CONFIG_ACTIVE = config_level)
  if (is.null(log_file)){
    invisible(flog.appender(appender.console()))
  }else{
    invisible(flog.appender(appender.tee(log_file)))
  }
  invisible(flog.threshold(if(debug) "DEBUG" else "INFO"))
  if (!verbose) (flog.threshold("ERROR"))
}


init_cluster <- function(num_cores = NULL){
  flog.debug("Initiate cluster...")
  cl <- makeCluster(get_num_cores(num_cores))
  flog.debug("Parallel is on with %s nodes.",get_num_cores(num_cores))
  invisible(clusterEvalQ(cl, c(library("pracma"), library("matrixcalc"), library("nloptr") )))
  flog.debug("Packages were loaded into the cluster nodes.")
  return(cl)
}


stop_cluster <- function(cl){
  flog.debug("Stop cluster")
  stopCluster(cl)
}


get_num_cores <- function(num_cores){
  if (is.null(num_cores)){
    return (max(1,detectCores() - 1))
  }else{
    return (num_cores)
  }
}


assert <- function (expr, error) {
  if (! expr) stop(error, call. = FALSE)
}

tensor.validate_input <- function(X, scale, parallel, num_cores, log_file, debug){

  flog.debug("Validating input types...")

  assert(is.matrix(X), "X must be of class 'matrix'")
  assert(!is.null(rownames(X)) & !is.null(colnames(X)), "X must have row names and column names")
  assert(is.logical(scale), "scale must be of class 'logical'")
  assert(is.logical(debug), "debug must be of class 'logical'")
  assert(is.logical(parallel), "parallel must be of class 'logical'")
  assert(is.null(num_cores) | is.numeric(num_cores), "argument num_cores must take a numric value or NULL")
  assert(is.character(log_file) | is.null(log_file), "log_file must be of class 'character' or NULL")

}

tca.validate_input <- function(X, W, C1, C1.map, C2, refit_W, refit_W.features, refit_W.sparsity, refit_W.sd_threshold, tau, constrain_mu, parallel, num_cores, max_iters, log_file, debug){

  flog.debug("Validating input types...")

  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)

  assert(is.matrix(X), "X must be of class 'matrix'")
  assert(is.matrix(W), "W must be of class 'matrix'")
  assert(is.null(C1) | is.matrix(C1), "C1 must be of class 'matrix' or NULL")
  assert(is.null(C1.map) | is.matrix(C1.map), "C1.map must be of class 'matrix' or NULL")
  assert(is.null(C2) | is.matrix(C2), "C2 must be of class 'matrix' or NULL")
  assert(is.null(refit_W.features) | is.character(refit_W.features), "refit_W.features must be of class 'character' or NULL")
  assert(is.null(tau) | is.numeric(tau), "tau must be of class 'numeric' or NULL")

  assert(is.numeric(refit_W.sparsity), "refit_W.sparsity must be of class 'numeric'")
  assert(is.numeric(refit_W.sd_threshold), "refit_W.sparsity must be of class 'numeric'")
  assert(is.numeric(max_iters), "max_iters must be of class 'numeric'")
  assert(is.logical(refit_W), "refit_W must be of class 'logical'")
  assert(is.logical(parallel), "parallel must be of class 'logical'")
  assert(is.logical(debug), "debug must be of class 'logical'")
  assert(is.character(log_file) | is.null(log_file), "log_file must be of class 'character' or NULL")

  flog.debug("Validating input stucture and values...")
  assert(!is.null(rownames(X)) & !is.null(colnames(X)), "X must have row names and column names")
  assert(!is.null(rownames(W)) & !is.null(colnames(W)), "W must have row names and column names")
  if (!is.null(C1)) assert(!is.null(rownames(C1)) & !is.null(colnames(C1)), "C1 must have row names and column names")
  if (!is.null(C2)) assert(!is.null(rownames(C2)) & !is.null(colnames(C2)), "C2 must have row names and column names")

  flog.debug("Validating input conditions...")
  if (refit_W) assert(refit_W.sparsity <= nrow(X) , "argument refit_W.sparsity must satisfy refit_W.sparsity <= nrow(X)")
  if (refit_W) assert(refit_W.sd_threshold >= 0 , "argument refit_W.sd_threshold must satisfy refit_W.sd_threshold >= 0")
  if (!(is.null(C1.map))) assert(constrain_mu , "argument C1.map cannot be useed with constrain_mu set to FALSE")

  flog.debug("Validating matrix dimensions...")
  assert(dim(X)[2] == dim(W)[1] , "The number of columns in X is inconsistent with the number of rows in W")
  if (!is.null(C1)) assert(dim(X)[2] == dim(C1)[1] , "The number of columns in X is inconsistent with the number of rows in C1")
  if (!is.null(C2)) assert(dim(X)[2] == dim(C2)[1] , "The number of columns in X is inconsistent with the number of rows in C2")

  flog.debug("Validating the order of observations across matrices...")
  assert(all(colnames(X) == rownames(W)), "The order of observations in W (in the rows) must match the order of the observations in X (in the columns).")
  if (!is.null(C1)) assert(all(colnames(X) == rownames(C1)), "The order of observations in C1 (in the rows) must match the order of the observations in X (in the columns).")
  if (!is.null(C2)) assert(all(colnames(X) == rownames(C2)), "The order of observations in C2 (in the rows) must match the order of the observations in X (in the columns).")

  flog.debug("Validating that W is non-negative and each row sums up to 1...")
  assert(all(W >= 0), "The entries of W must be non-negative.")
  assert(all(abs(rowSums(W) - 1) < 0.0001), "Each row in W must sum up to 1.")

  if (!is.null(C1.map)) assert(sum(C1.map == 1 | C1.map == 0) == ncol(C1.map)*nrow(C1.map), "The entries of C1.map must all be 0 or 1")

  if (!is.null(tau)) assert(tau >= 0, "tau must be non-negative")

  th <- config[["min_sd"]]**2
  assert(sum(rowVars(X) < th) == 0, paste("X must not include features with variance less than ",as.character(th),sep=""))

}


tcasub.validate_input <- function(tca.mdl, features, log_file, debug){
  flog.debug("Validating input...")
  assert(is.list(tca.mdl) & is.matrix(tca.mdl[["W"]]) &
    is.matrix(tca.mdl[["mus_hat"]]) &
    is.matrix(tca.mdl[["sigmas_hat"]]) &
    is.numeric(tca.mdl[["tau_hat"]]) &
    is.matrix(tca.mdl[["deltas_hat"]]) &
    is.matrix(tca.mdl[["gammas_hat"]]) &
    is.matrix(tca.mdl[["C1"]]) &
    is.matrix(tca.mdl[["C2"]]), "'tca.mdl' must be a returned value of the function 'tca'")
  assert((!is.null(features)) & length(intersect(rownames(tca.mdl$mus), features)) == length(features), "'features' must be a subset of the features in 'tca.mdl'")
  assert(is.logical(debug), "debug must be of class 'logical'")
  assert(is.character(log_file) | is.null(log_file), "log_file must be of class 'character' or NULL")
}


tcareg.validate_input <- function(X, W, y, C3, test, null_model, alternative_model, save_results, fast_mode, output, sort_results, parallel, num_cores, log_file, features_metadata ,debug){

  flog.debug("Validating input type for tcareg...")

  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)

  options <- c("marginal","marginal_conditional","joint","single_effect","custom")
  flog.debug("Validating input type...")
  assert(is.matrix(X), "argument X must be of class 'matrix'")
  assert(is.matrix(y), "argument y must be of class 'matrix'")
  assert(is.null(C3) | is.matrix(C3), "argument C3 must be of class 'matrix' or NULL")

  assert(is.element(test,options), paste("argument test must be one of the following options: ", paste(options, collapse=', ' ),sep = ""))

  assert(is.null(num_cores) | is.numeric(num_cores), "argument num_cores must take a numric value or NULL")
  assert(is.null(null_model) | is.character(null_model), "argument null_model must take a character value or NULL")
  assert(is.null(alternative_model) | is.character(alternative_model), "argument alternative_model must take a character value or NULL")
  assert(is.null(features_metadata) | is.character(features_metadata), "argument features_metadata must take a character value or NULL")
  assert(is.character(output), "argument output must take a character value")

  assert(is.logical(save_results), "argument save_results must take a logical value")
  assert(is.logical(fast_mode), "argument fast_mode must take a logical value")
  assert(is.logical(sort_results), "argument sort_results must take a logical value")
  assert(is.logical(parallel), "argument parallel must take a logical value")
  assert(is.logical(debug), "argument debug must take a logical value")
  assert(is.character(log_file) | is.null(log_file), "argument log_file must take a character value or NULL")

  flog.debug("Validating input stucture and values...")
  assert(!is.null(rownames(X)) & !is.null(colnames(X)), "X must have row names and column names")
  assert(!is.null(rownames(W)) & !is.null(colnames(W)), "W must have row names and column names")
  assert(!is.null(rownames(y)), "y must have row names")
  if (!is.null(C3)) assert(!is.null(rownames(C3)) & !is.null(colnames(C3)), "C3 must have row names and column names")
  if (test == "custom") assert( (!is.null(alternative_model)), "argument alternative_model cannot be NULL when test=custom")
  if ( (!is.null(alternative_model))){
    assert(test == "custom", "argument test must be set to 'custom' if argument alternative_model is not NULL")
    # make sure that null_model and alternative_model include labels that exist in the column names of W
    assert(all(is.element(null_model, colnames(W))) & all(is.element(alternative_model, colnames(W))), "null_model and alternative_model must include values that exist in the column names of W")
    if ( (!is.null(null_model))){
      # make sure the two models are nested
      assert(length(setdiff(alternative_model, null_model)) > 0 & length(setdiff(null_model, alternative_model)) == 0, "null_model must be nested within alternative_model")
    }
  }
  assert(!( (sort_results & test == "marginal_conditional") & fast_mode), "sort_results cannot be set to true under fast_mode == TRUE and test = 'marginal_conditional'")

  flog.debug("Validating matrix dimensions...")
  if (!is.null(C3)) assert(dim(X)[2] == dim(C3)[1] , "the number of columns in X is inconsistent with the number of rows in C3")
  assert(dim(X)[2] == length(y) , "the number of columns in X is inconsistent with the number of rows in y")

  flog.debug("Validating the order of observations across matrices...")
  assert(all(colnames(X) == rownames(W)), "The order of observations in W (in the rows) must match the order of the observations in X (in the columns).")
  if (!is.null(C3)) assert(all(colnames(X) == rownames(C3)), "The order of observations in C3 (in the rows) must match the order of the observations in X (in the columns).")
  assert(all(colnames(X) == rownames(y)), "The order of observations in y (in the rows) must match the order of the observations in X (in the columns).")

  th <- config[["min_sd"]]**2
  assert(sum(rowVars(X) < th) == 0, paste("X must not include features with variance less than ",as.character(th),sep=""))

}


# log likelihood ratio test
# ll0 - vector with the log likelihood values for m tests under the null model
# ll1 - vector with the log likelihood values for m tests under the alternative model
# df - vector with the degrees of freedom for each of the m tests
lrt.test <- function(ll0, ll1, df){
  flog.debug("Running lrt.test...")
  m <- length(ll1)
  stats = numeric(m)
  pvals = numeric(m)
  for (j in 1:m){
    stats[j] <- -2*(ll0[j]-ll1[j])
    pvals[j] <- pchisq(stats[j], df = df, lower.tail=FALSE)
  }
  return(list("stats" = stats, "pvals" = pvals))
}


save_association_results <- function(res, output, test, fast_mode, alternative_model, feature_ids, W_names, C3_names, sort_results, features_metadata){
  flog.debug("Running save_association_results...")
  m <- length(feature_ids)
  metadata <- if (is.null(features_metadata)) matrix(0,m,0) else parse_features_metadata(feature_ids, features_metadata)
  if (test == "marginal" | (!fast_mode & test == "marginal_conditional") ){
    for (i in 1:length(res)){
      filename <- paste(output, ".", test, ".", W_names[i], ".txt", sep ="")
      wnames <- if (test == "marginal") W_names[i] else W_names
      save_association_results.save(res[[i]], test, fast_mode, filename, feature_ids, sort_results, wnames, C3_names, metadata)
    }
  }else{
    filename <- paste(output, ".", test, ".txt", sep ="")
    keep <- 1:length(W_names)
    if (test == "custom") keep <- alternative_model
    if (test == "single_effect") keep <- 1
    save_association_results.save(res, test, fast_mode, filename, feature_ids, sort_results, W_names[keep], C3_names, metadata)
  }
}


save_association_results.save <- function(res, test, fast_mode, filename, feature_ids, sort_results, W_names, C3_names, metadata){
  flog.debug("save_association_results.save...")
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  m <- length(feature_ids)
  betas <- if(length(W_names) == 1) "beta" else unlist(lapply(1:length(W_names), function(i) paste("beta.",W_names[i],sep="")))
  pvals <- if(test == "marginal_conditional" & fast_mode) unlist(lapply(1:length(W_names), function(i) paste("pval.",W_names[i],sep=""))) else "pval"
  qvals <- if(test == "marginal_conditional" & fast_mode) unlist(lapply(1:length(W_names), function(i) paste("qval.",W_names[i],sep=""))) else "qval"
  alpha <- if (length(C3_names)) unlist(lapply(1:length(C3_names), function(i) paste("alpha.",C3_names[i],sep=""))) else c()
  df <- if(is.null(res$df)) c() else "df"
  df.values <- if(is.null(res$df)) matrix(0,m,0) else matrix(0,m,1)+res$df

  if (fast_mode){
    data <- data.frame(feature_ids, metadata, res$pval, res$qval, res$beta, res$alternative_ll, res$stat, df.values, res$intercept, res$alpha, res$phi)
    colnames(data) <- c("ID", colnames(metadata), pvals, qvals, betas, "alternative_ll", "stat", df, "intercept", alpha, "phi")
  }else{
    data <- data.frame(feature_ids, metadata, res$pval, res$qval, res$beta, res$null_ll, res$alternative_ll, res$stat, res$df, res$intercept, res$alpha, res$phi)
    colnames(data) <- c("ID", colnames(metadata), "pval", "qval", betas, "null_ll", "alternative_ll", "chi_squared", "df", "intercept", alpha, "phi")
  }
  if (sort_results){
    # sort by p-value
    data <- data[order(data$pval),]
  }
  fwrite(data, file = filename, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = config$sep, showProgress = config$show_progress, verbose = config$verbose)
}


parse_features_metadata <- function(feature_ids, features_metadata_file){
  flog.debug("parse_features_metadata...")
  metadata <- data.frame(fread(file = features_metadata_file, header=TRUE, showProgress=FALSE))
  rownames(metadata) <- metadata[,1]
  return(metadata[feature_ids,2:ncol(metadata)])
}
