#General Functions for Quantile Regression
library(quantreg)
library(foreach)
library(doParallel)
library(data.table)


rq.fit.sfn_start_val <- function(a,y,tau=.5,
                                 rhs = (1-tau)*c(t(a) %*% rep(1,length(y))), control, sv, weight_vec = NULL)
{
  y <- -y
  n <- length(y)
  m <- a@dimension[2]
  if(n != a@dimension[1])
    stop("Dimensions of design matrix and the response vector not compatible")

  # additional syntax to incorporate weights is included here
  if (!is.null(weight_vec)){
    if(n != dim(as.matrix(weight_vec))[1]){
      stop("Dimensions of design matrix and the weight vector not compatible")
    }
    # multiplying y by the weights
    y <- y * weight_vec

    # pre-multiplying the a matrix by a diagonal matrix of weights
    #a <- sweep(a,MARGIN=1,weight_vec,`*`)
    a <- as(as.vector(weight_vec), "matrix.diag.csr") %*% a
  }

  u <- rep(1,length=n)
  x <- rep((1-tau),length=n)
  nnzdmax <- nnza <- a@ia[n+1]-1
  iwmax <- 7*m+3
  ao <- t(a)
  e <- ao %*% a
  nnzemax <- e@ia[m+1]-1
  ctrl <- sfn.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  nsubmax <- ctrl$nsubmax
  tmpmax <- ctrl$tmpmax
  nnzlmax <- ctrl$nnzlmax
  if (is.null(ctrl$nsubmax)) nsubmax <- nnzemax
  if (is.null(ctrl$tmpmax)) tmpmax <- 6 * m
  if (is.null(ctrl$nnzlmax)) nnzlmax <- 4 * nnzdmax
  wwm <- vector("numeric",3*m)
  s <- u - x
  if(missing(sv)){
    b1 <- solve(e, ao %*% y, tmpmax=tmpmax,nnzlmax=nnzlmax,nsubmax=nsubmax)
  }
  else {
    # note: LDWS flipped the sign here, since formula above yields OLS coeff * -1
    b1 = -sv
  }

  r <- y - a %*% b1
  z <- ifelse(abs(r)<ctrl$small,(r*(r>0)+ctrl$small),r*(r>0))
  w <- z - r
  wwn <- matrix(0,n,14)
  wwn[,1] <- r
  wwn[,2] <- z
  wwn[,3] <- w
  fit <- .Fortran("srqfn",
                  n = as.integer(n),
                  m = as.integer(m),
                  nnza = as.integer(nnza),
                  a = as.double(a@ra),
                  ja = as.integer(a@ja),
                  ia = as.integer(a@ia),
                  ao = as.double(ao@ra),
                  jao = as.integer(ao@ja),
                  iao = as.integer(ao@ia),
                  nnzdmax = as.integer(nnzdmax),
                  d = double(nnzdmax),
                  jd = integer(nnzdmax),
                  id = integer(m+1),
                  dsub = double(nnzemax+1),
                  jdsub = integer(nnzemax+1),
                  nnzemax = as.integer(nnzemax),
                  e = as.double(e@ra),
                  je = as.integer(e@ja),
                  ie = as.integer(e@ia),
                  nsubmax = as.integer(nsubmax),
                  lindx = integer(nsubmax),
                  xlindx = integer(m+1),
                  nnzlmax = as.integer(nnzlmax),
                  lnz = double(nnzlmax),
                  xlnz = integer(m+1),
                  iw = integer(m*5),
                  iwmax = as.integer(iwmax),
                  iwork = integer(iwmax),
                  xsuper = integer(m+1),
                  tmpmax = as.integer(tmpmax),
                  tmpvec = double(tmpmax),
                  wwm = as.double(wwm),
                  wwn = as.double(wwn),
                  cachsz = as.integer(ctrl$cachsz),
                  level = as.integer( 8 ),
                  x = as.double(x),
                  s = as.double(s),
                  u = as.double(u),
                  c = as.double(y),
                  sol = as.double(b1),
                  rhs = as.double(rhs),
                  small = as.double(ctrl$small),
                  ierr = integer(1),
                  maxiter = as.integer(ctrl$maxiter),
                  time = double(7),
                  PACKAGE = "quantreg")[c("sol","ierr","maxiter","time")]
  ierr <- fit$ierr
  if(!(ierr==0) && ctrl$warn.mesg)
    warning(sfnMessage(ierr))
  coefficients <- -fit$sol

  residuals <- -y - a %*% coefficients
  if (!is.null(weight_vec)){
    residuals <- residuals / weight_vec
  }

  list(coefficients = coefficients,
       residuals = residuals,
       control = ctrl,
       ierr = ierr,
       it = fit$maxiter,
       weight_vec = weight_vec)
}

getRank = function(m) {
  # Computes the rank of a matrix. Outputs consistent results across
  # various matrix types.
  #
  # Args:
  #   A matrix m.
  #
  # Returns:
  #   The rank of m.
  #
  transp_prod <- as.matrix(t(m) %*% (m))
  return(sum(abs(diag(qr.R(qr(transp_prod)))) > TOL))
}

findRedundantCols = function(m) {
  # Given a rank-deficient M by N sparse matrix, where M >= N, in dgCMatrix format,
  # returns the indices of columns to remove from the original matrix so that the
  # resulting matrix is full rank.
  #
  # Args:
  #   m: Rank-deficient M by N sparse matrix, where M >= N, in dgCMatrix format.
  #
  # Returns:
  #   Indices of columns to remove from m so that the remaining matrix is full rank.
  #
  decomp <- qr(m)
  orig_col_names <- colnames(m)
  R <- qr.R(decomp)
  R_col_names <- colnames(R)
  return(setdiff(
    1:length(orig_col_names) ,
    match(
      R_col_names[which(abs(diag(R)) > TOL)],
      orig_col_names)))
}

ensureSpecFullRank = function(spec_mat, col_names) {
  # Verifies if a regression specification is full-rank. If the input
  # is rank-deficient, identifies and drops columns so that the remaining
  # matrix is full-rank.
  #
  # Args:
  #   spec_mat: An M by N regression specification matrix, where M > N.
  #   col_names: Column names of the regression specification.
  #
  # Returns:
  #   A full-rank specification matrix. If the input is full-rank, returns the
  #   input unmodified. Otherwise, returns a matrix with a subset of the columns
  #   from the input.
  #

  # Check if input is already matrix full rank


  if (getRank(spec_mat) == ncol(spec_mat)) {
    return(list(
      "spec_matrix" = spec_mat,
      "var_names" = col_names))
  }

  transp_prod <- t(spec_mat) %*% spec_mat
  # First check if any dummy columns are all zero
  zero_cols <- which(diag(transp_prod) == 0)
  if (length(zero_cols) > 0) {
    spec_mat <- spec_mat[,-zero_cols]
    col_names <- col_names[-zero_cols]
  }

  # Check if updated matrix is full rank
  if (getRank(spec_mat) == ncol(spec_mat)) {
    return(list(
      "spec_matrix" = spec_mat,
      "var_names" = col_names))
  }

  # Use a generic routine to identify columns to drop
  dgc_spec_mat <- as(spec_mat, "dgCMatrix")
  colnames(dgc_spec_mat) <- col_names
  cols_to_drop = findRedundantCols(dgc_spec_mat)
  return(list(
    "spec_matrix" = spec_mat[,-cols_to_drop],
    "var_names" = col_names[-cols_to_drop]))
}

denseMatrixToSparse = function(m) {
  # Utility to convert a matrix to a SparseM csr matrix.
  #
  # Args:
  #   m: A matrix in a standard, dense format.
  #
  # Returns:
  #   A matrix in SparseM csr format.
  #

  if (length(m) <= .Machine$integer.max) {
    return(as.matrix.csr(m))
  }

  # as.matrix.csr cannot coerce a long vector to csr,
  # so break the input into maximally-sized chunks with
  # are coercible and then bind them into one sparse matrix

  chunk_size <- floor(.Machine$integer.max / ncol(m)) * ncol(m)
  num_chunks <- floor(length(m) / chunk_size)

  items_to_cut <- chunk_size * num_chunks
  sparse_chunks <- lapply(
    split(t(m)[1:items_to_cut], sort(rep(1:num_chunks, chunk_size))),
    FUN = function(x) {
      return(as.matrix.csr(matrix(x, ncol=ncol(m), byrow=TRUE)))
    })

  combined <- do.call(rbind, sparse_chunks)

  if (length(m) > items_to_cut) {
    tail <- as.matrix.csr(
      matrix(
        t(m)[(items_to_cut + 1):length(m)],
        ncol=ncol(m),
        byrow=TRUE))
    combined <- rbind(combined, tail)
  }

  return(combined)
}

addMissingSpecColumns = function(df, names) {
  missing_cols <- setdiff(names, colnames(df))
  df[missing_cols] <- 0
  return(df[names])
}

#For Rho function for calculating Pseudo-R^2 values
rho <- function(u,tau=.5,weight_vec = NULL){
  # Evaluates the check objective function (possibly weighted) for QREG
  #
  #Args:
  #   u: vector of residuals
  #   tau: probability index parameter (quantile of interest)
  #   weight_vec: optional vector of weights
  #
  #Returns:
  #  sum of absolute residuals
  #

  if (is.null(weight_vec)){
    ( u*(tau - (u < 0)) )
  } else {
    ( weight_vec*u*(tau - (u < 0)) )
  }

}

getColNums = function(start_list,
                      reg_spec_data,
                      alpha,
                      j){
  #
  #Args:
  #   start_list: starting values (can be NA's) to be fd into sfn_start_val function
  #   reg_spec_data: result of ensureSpecRank function; regression matrix with full rank
  #   alpha: column vector of quantiles to be estimated
  #   j: index of quantile currently being calculated
  #
  #Returns:
  #  col_nums: if start_list is supplied, then returns the correct column numbers to be used in regression;
  #            else, returns NULL
  #
  if(!any(is.na(start_list))){
    #get the columns appropriate for starting values
    cols = reg_spec_data$var_names
    cols = paste(alpha[j], cols, sep = '_')
    col_nums = colnames(start_list) %in% cols
    return(col_nums)
  }
  return(NULL)
}

quantRegress = function(reg_spec_data,
                        ehat,
                        sv,
                        ind_hat,
                        tau,
                        trunc,
                        small,
                        control,
                        weight_vec = NULL) {
  #Runs quantile regression on residuals of the model (calculates spaces around jstar quantile)
  #
  #Args:
  #   reg_spec_data: result of ensureSpecRank function; regression matrix with full rank
  #   ehat: current residuals; subset of which to be used as dependent column
  #   sv: starting values (can be NA's) to be fd into sfn_start_val function
  #   ind_hat: column vector indicating which rows to be used in quantile regression
  #   tau: estimated quantile
  #   trunc: Boolean value; if true, replace those dependent values less than small with small itself;
  #          else, only use rows with residuals greater than small
  #   small: Value used with trunc; values less than small 'blow up' too greatly when logged
  #   control: control list to be fed to sfn_start_val function
  #   weight_vec: vector of optional weights
  #
  #Returns:
  #   j_model: list of estimated coefficients, warnings, iterations, and controls as in
  #            standard quantile regression function
  #
  if (trunc) { # then replace values using pmax
    #     if(missing(sv)){
    #       j_model <- rq.fit.sfn_start_val(
    #         a = reg_spec_data$spec_matrix,
    #         y = log(pmax(ehat[ind_hat],small)),
    #         tau = tau,
    #         control = control) # Model the quantile
    #     }
    if (!is.null(weight_vec)){
      weight_vec = as.matrix(weight_vec[ind_hat])
    }

    j_model <- rq.fit.sfn_start_val(
      a = reg_spec_data$spec_matrix,
      y = log(pmax(ehat[ind_hat],small)),
      tau = tau,
      sv = sv,
      control = control,
      weight_vec = weight_vec) # Model the quantile
  }
  else {
    resids = ehat[ind_hat]

    if (!is.null(weight_vec)){
      weight_vec = as.matrix(weight_vec[ind_hat])
      weight_vec = weight_vec[resids > small]
    }
    if(dim(reg_spec_data$spec_matrix)[1] == sum(resids > small)) {
      # if the dimensions already match, then we assume that the matrix
      # was already spliced for the correct rows outside of the function
      j_model <- rq.fit.sfn_start_val(
        a = reg_spec_data$spec_matrix,
        y = log(resids[resids > small]),
        tau = tau,
        sv = sv,
        control = control,
        weight_vec = weight_vec) # Model the quantile
    }
    else{
      # else, we splice the correct rows here
      j_model <- rq.fit.sfn_start_val(
        a = reg_spec_data$spec_matrix[resids > small,],
        y = log(resids[resids > small]),
        tau = tau,
        sv = sv,
        control = control,
        weight_vec = weight_vec) # Model the quantile
    }
  }
  return(j_model)
}

printWarnings = function(model){
  #handles warning messages; print output for better logging
  #
  #Args:
  #   model: result of quantRegress call
  #
  #Returns:
  #   only prints statements
  #

  if(model$ierr != 0){
    print(paste('Quantile Regresion ran into warning', model$ierr, ' and had ',model$it,' iterations.'))
    if(model$it < 12){
      print('Iterations failed to go past threshold')
    }
  }
}

quantRegSpacing = function(
  dep_col,
  data,
  var_names,
  alpha,
  jstar,
  small = 1e-3,
  trunc = FALSE,
  start_list,
  weight_vec = NULL) {
  # Computes coefficients for the quantile regression spacing method.
  #
  # Args:
  #   dep_col: Column of response variable.
  #   data: Regression specification matrix.
  #   var_names: RHS regression variable names.
  #   alpha: Quantiles to be estimated.
  #   jstar: First quantile to be estimated (usually the center one)
  #   p: Length of alpha.
  #   small: Minimum size of residuals for computational accuracy.
  #   trunc: Boolean value; if true, replace those dependent values less than small with small itself;
  #          else, only use rows with residuals greater than small
  #   start_list: Starting values for regression optimization.
  #   weight_vec: vector of optional weights
  #
  # Returns:
  #   list of $coef: num_betas x p matrix of estimated parameters for each supplied quantiles,
  #           $pseudo_r: 1 x p matrix of psuedo R^2 values for each quantile estimate,
  #           $warnings: 1 x p matrix of warnings produced by each quantile regression call,
  #           $iter: 1 x p matrix of iterations ran by each quantile regression call
  #

  width = dim(data)[2]
  tau = alpha[jstar]
  p = length(alpha)

  #create logs for output
  count_log = list()
  length(count_log) = p
  warnings_log = list()
  length(warnings_log) = p
  iter_log = list()
  length(iter_log) = p
  pseudo_r = list()
  length(pseudo_r) = p
  model = list() # Collect the quantile regressions into a list
  length(model) = p

  # check to see if regression matrix is sparse. If not, then turn into CSR matrix
  if(!is(data, 'matrix.csr')) data = denseMatrixToSparse(data)

  tmpmax <- floor(1e5 + exp(-12.1)*(data@ia[width+1]-1)^2.35)

  # Ensure matrix is not rank deficient
  reg_spec_starting_data <- ensureSpecFullRank(spec_mat = data, col_names = var_names)

  # Calculate initial fit
  ##print(jstar)
  ptm <- proc.time()
  if(!missing(start_list)){
    col_nums = getColNums(start_list, reg_spec_starting_data, alpha, jstar)
    sv = as.numeric(start_list[col_nums])
    star_model = rq.fit.sfn_start_val(
      a = reg_spec_starting_data$spec_matrix,
      y = dep_col,
      tau = tau,
      control = list(tmpmax= tmpmax),
      sv = sv,
      weight_vec = weight_vec)
  }
  else{
    star_model = rq.fit.sfn_start_val(
      a = reg_spec_starting_data$spec_matrix,
      y = dep_col,
      tau = tau,
      control = list(tmpmax= tmpmax),
      weight_vec = weight_vec)
  }
  ##printWarnings(star_model)

  ehat0 = star_model$residuals

  #Calculate R^2
  V <- sum(rho(u = ehat0, tau = tau, weight_vec = weight_vec))
  V0 <- rq.fit.sfn_start_val(a = as.matrix.csr(rep(1, length(dep_col))),
                             y = dep_col,
                             tau = tau,
                             weight_vec = weight_vec)$residuals
  V0 <- sum(rho(u = V0, tau = tau,weight_vec = weight_vec))

  ##print(paste('Pseudo-R^2 Value of', (1 - V/V0), 'for', alpha[jstar], 'quantile'))
  ##print(proc.time() - ptm)

  #set column names
  coef_df <- as.data.frame(t(star_model$coefficients))
  colnames(coef_df) <- reg_spec_starting_data$var_names
  coef_df <- addMissingSpecColumns(
    coef_df,
    var_names)
  colnames(coef_df) <- paste(alpha[jstar], colnames(coef_df), sep="_")

  #log output for return
  pseudo_r[[jstar]] = (1 - V/V0)
  model[[jstar]] = coef_df
  warnings_log[[jstar]] = star_model$ierr
  iter_log[[jstar]] = star_model$it
  count_log[[jstar]] = dim(reg_spec_starting_data$spec_matrix)[1]

  rm(star_model)

  # Estimate upper quantiles sequentially
  ehat = ehat0
  for (j in (jstar+1):p) {
    ind_hat = which(ehat > 0)

    # Determine quantile to estimate
    tau.t = (alpha[j] - alpha[j-1])/(1 - alpha[j-1])

    # Ensure the cut of the starting data that we take for
    # current spacing is not rank-deficient
    if(!trunc) reg_spec_data <- ensureSpecFullRank(spec_mat = reg_spec_starting_data$spec_matrix[which(ehat > small),],
                                                  col_names = reg_spec_starting_data$var_names)
    # else, handle rank specification for typically-sized matrix
    else reg_spec_data <- ensureSpecFullRank(spec_mat = reg_spec_starting_data$spec_matrix[ind_hat,],
                                             col_names = reg_spec_starting_data$var_names)

    #run quantile regression
    coef <- NULL
    ##print(j)
    ptm <- proc.time()
    if(!missing(start_list)){
      col_nums = getColNums(start_list, reg_spec_data, alpha, j)
      sv = as.numeric(start_list[col_nums])
      j_model <- quantRegress(reg_spec_data = reg_spec_data, ehat = ehat,
                              sv = sv, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                              small = small, control = list(tmpmax = tmpmax), weight_vec = weight_vec)
    }
    else{
      j_model <- quantRegress(reg_spec_data = reg_spec_data, ehat = ehat,
                              ind_hat = ind_hat, tau = tau.t, trunc =  trunc,
                              small = small, control = list(tmpmax = tmpmax), weight_vec = weight_vec)
    }
    #printWarnings(j_model)

    #Calculate R^2
    V <- sum(rho( u = j_model$residuals, tau = tau.t, weight_vec = j_model$weight_vec))
    V0 <- quantRegress(reg_spec_data = list('spec_matrix' = as.matrix.csr(rep(1, length(ind_hat)))),
                       ehat = ehat, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                       small = small, weight_vec = weight_vec)
    V0 <- sum(rho(u = V0$residuals, tau = tau.t, weight_vec = V0$weight_vec))
    #print(paste('Pseudo-R^2 Value of', (1 - V/V0), 'for', alpha[j], 'quantile'))
    #print(proc.time() - ptm)

    # Update residuals
    coef = j_model$coefficients
    coef_df <- as.data.frame(t(coef))


    #get column names
    colnames(coef_df) <- reg_spec_data$var_names
    coef_df <- addMissingSpecColumns(
      coef_df,
      var_names)
    colnames(coef_df) <- paste(alpha[j], colnames(coef_df), sep="_")

    #log results
    model[[j]] = coef_df
    pseudo_r[[j]] = (1 - V/V0)
    warnings_log[[j]] = j_model$ierr
    iter_log[[j]] = j_model$it
    count_log[[j]] = dim(reg_spec_data$spec_matrix)[1]

    # Update residuals
    ehat = ehat - exp(
      as.matrix(
        data %*%
          unname(t(as.matrix(model[[j]])))))

  }

  # Estimate lower quantiles sequentially
  ehat = ehat0
  for (j in (jstar-1):1) {
    ind_hat = which(ehat < 0)

    # Determine quantile to estimate
    tau.t = (alpha[j + 1] - alpha[j])/(alpha[j + 1])

    # Ensure the cut of the starting data that we take for
    # current spacing is not rank-deficient
    # if truncating, then ensure exact rows of the regression matrix is included
    if(!trunc) reg_spec_data <- ensureSpecFullRank(spec_mat = reg_spec_starting_data$spec_matrix[which(-ehat > small),],
                                                  col_names = reg_spec_starting_data$var_names)
    # else, handle rank specification for typically-sized matrix
    else reg_spec_data <- ensureSpecFullRank(spec_mat = reg_spec_starting_data$spec_matrix[ind_hat,],
                                             col_names = reg_spec_starting_data$var_names)
    #run quantile regression

    coef <- NULL
    #print(j)
    ptm <- proc.time()
    if(!missing(start_list)){
      col_nums = getColNums(start_list, reg_spec_data, alpha, j)
      sv = as.numeric(start_list[col_nums])
      j_model <- quantRegress(reg_spec_data = reg_spec_data,
                              ehat = -ehat,
                              sv = sv,
                              ind_hat = ind_hat,
                              tau = tau.t, trunc = trunc, small = small,
                              control = list(tmpmax = tmpmax),
                              weight_vec = weight_vec)
    }
    else{
      j_model <- quantRegress(reg_spec_data = reg_spec_data,
                              ehat = -ehat,
                              ind_hat = ind_hat,
                              tau = tau.t, trunc = trunc, small = small,
                              control = list(tmpmax = tmpmax),
                              weight_vec = weight_vec)
    }

    #printWarnings(j_model)

    #Calculate pseudo-R^2
    V <- sum(rho(u = j_model$residuals, tau = tau.t, weight_vec = j_model$weight_vec))
    V0 <- quantRegress(reg_spec_data = list('spec_matrix' = as.matrix.csr(rep(1, length(ind_hat)))),
                       ehat = -ehat, ind_hat = ind_hat, tau = tau.t, trunc = trunc,
                       small = small, weight_vec = weight_vec)
    V0 <- sum(rho(u = V0$residuals, tau = tau.t, weight_vec = V0$weight_vec))
    ##print(paste('Pseudo-R^2 Value of', (1 - V/V0), 'for', alpha[j], 'quantile'))
    ##print(proc.time() - ptm)

    # Update residuals
    coef = j_model$coefficients
    coef_df <- as.data.frame(t(coef))

    #get column names
    colnames(coef_df) <- reg_spec_data$var_names
    coef_df <- addMissingSpecColumns(
      coef_df,
      var_names)
    colnames(coef_df) <- paste(alpha[j], colnames(coef_df), sep="_")

    #log results
    model[[j]] = coef_df
    pseudo_r[[j]] = (1 - V/V0)
    warnings_log[[j]] = j_model$ierr
    iter_log[[j]] = j_model$it
    count_log[[j]] = dim(reg_spec_data$spec_matrix)[1]

    # Update residuals
    ehat = ehat + exp(as.matrix(
      data %*%
        unname(t(as.matrix(model[[j]])))))
  }

  return(list('coef' = do.call(cbind, model),
              'pseudo_r' = do.call(cbind, pseudo_r),
              'warnings' = do.call(cbind, warnings_log),
              'iter' = do.call(cbind, iter_log),
              'counts' = do.call(cbind, count_log)))
}






clusterSample = function(cluster_indices,
                         stratum_indices = NULL,
                         M = 1,
                         draw_weights = FALSE) {
  #Draws a subsample of clusters and (optionally) assigns
  #random, exponential weights to each cluster
  #
  #Args:
  #   cluster_indices: N_cluster x 2 vector with starting
  #                    and ending positions of each cluster.
  #                    NOTE that this requires that data are
  #                    sorted in this dimension
  #   stratum_indices: N_cluster x 1 vector with stratum indices
  #                    (factor variables). This enables the
  #                    researcher to stratify (e.g., on cluster size)
  #                    Set to NULL to take a random subsample
  #   M: greater than 0 but less equal to 1, percentage of
  #      clusters to sample.
  #   draw_weights : Boolean value; if true, draw a vector of exponential
  #                  weights to use in subsample
  #Returns:
  #   subsample_rows: Boolean vector of selected rows of sampled clusters
  #   subsample_weights: either NULL, or vector of cluster weights

  # TODO: probably should add some logic to validate inputs here...

  # checking conditions on M
  if(M > 1 | 0 >= M){
    stop("Subsampling percentage must be between 0 and 1")
  }

  # draw an M% stratified sample of indices
  if (!is.null(stratum_indices)){
    sampled_indices <- do.call(rbind,
                               lapply(split(as.data.frame(cluster_indices), stratum_indices),
                                      function(x) x[sample(nrow(x), floor(nrow(x)*M)),]))
  } else {
    sampled_indices <- as.data.frame(cluster_indices[
      sample(nrow(cluster_indices),
             floor(nrow(cluster_indices)*M) ),
      ])
  }

  # now, re-sort the vector
  sampled_indices <- sampled_indices[order(sampled_indices[,1]),]
  n_sampled_rows <- sum((sampled_indices$V2-sampled_indices[,1]+1))
  # types seem to get messed up easily, so I'll just initialize with the sample function
  subsample_rows <- matrix(FALSE,max(cluster_indices[,2]),1)

  # drawing weights (if applicable)
  if (!draw_weights){
    subsample_weights <- NULL
  } else {
    cluster_weights <- pmax(rexp(dim(sampled_indices)[1]),5e-3)
    subsample_weights <- matrix(0,n_sampled_rows,1)
  }

  # looping over sampled clusters to populate output vectors
  row_s_ct <- 0
  for (j in 1:dim(sampled_indices)[1] ){
    sidx <- sampled_indices[j,1]
    eidx <- sampled_indices[j,2]

    subsample_rows[sidx:eidx] <- TRUE

    if (draw_weights){
      subsample_weights[(row_s_ct+1):(row_s_ct+(eidx-sidx+1)),] <- cluster_weights[j]
    }
    row_s_ct <- row_s_ct + (eidx-sidx+1)
  }
  subsample_rows <- which(subsample_rows)

  return(list('subsample_rows' = subsample_rows,
              'subsample_weights' = subsample_weights))

}













subsampleStandardErrors = function(
  dep_col,
  data,
  var_names,
  alpha,
  jstar,
  cluster_indices = NULL,
  stratum_indices = NULL,
  M = 0.2,
  draw_weights = FALSE,
  num_bs = 100, parallelize = F, num_cores = 1, small = 1e-6,
  trunc = FALSE, save_rows = FALSE, start_model, weight_vec = NULL,
  square_ols_weights = FALSE) {
  # Computes standard errors for the quantile regression spacing method using
  # subsampling.
  #
  # Inputs:
  #   data: Regression specification matrix.
  #   dep_col: Column of response variable.
  #   var_names: RHS regression variable names.
  #   alpha: Quantiles to be estimated.
  #   jstar: First quantile to be estimated (usually the center one)
  #   p: Length of alpha.
  #   cluster_indices: N_cluster x 2 vector with starting
  #                    and ending positions of each cluster.
  #                    NOTE that this requires that data are
  #                    sorted in this dimension
  #   stratum_indices: N_cluster x 1 vector with stratum indices
  #                    (factor variables). This enables the
  #                    researcher to stratify (e.g., on cluster size)
  #                    Set to NULL to take a random subsample
  #   M: greater than 0 but less equal to 1, percentage of
  #      clusters to sample.
  #   draw_weights : Boolean value; if true, draw a vector of exponential
  #                  weights to use in subsample
  #   num_bs: Number of subsample draws (must be greater than 1).
  #   small: Minimum size of residuals for computational accuracy.
  #   trunc: Boolean value; if true, replace those dependent values less than small with small itself;
  #          else, only use rows with residuals greater than small
  #   start_model: Starting values for regression's optimization.
  #   save_rows: If set to TRUE, function will save vector of row numbers
  #              that triggered warnings by quantile regression
  #   weight_vec: vector of same length and order as dependent column, to be used as weights for estimation
  #               (note, if draw weights is set to TRUE, this variable will be the element-wise product
  #                of itself and a random vector of weights)
  #   square_ols_weights : boolean value; if true, square user-supplied weight_vec for OLS
  #
  # Returns:
  #   list of $cov: num_betas x num_betas covariance matrix using
  #                 bootstrapped subsampling covariances
  #           $warnings: num_bs x p matrix of warning messages
  #                      produced in each bootstrap's quantile regression call
  #           $iter: num_bs x p matrix of iterations ran by each
  #                  bootstrapped quantile regression call

  num_betas = dim(data)[2]

  # checking conditions on M
  if(M > 1 | 0 >= M){
    stop("Subsampling percentage must be between 0 and 1")
  }

  if(parallelize){

    cl = makeCluster(num_cores, outfile = '')

    #parallelize and create status bar
    #pb <- tkProgressBar(max=num_bs)
    #progress <- function(n) setTkProgressBar(pb, n)
    #opts <- list(progress=progress)
    registerDoParallel(cl)
    #registerDoSNOW(cl)


    print(paste0(
      "Number of workers for subsampling is: ",
      getDoParWorkers()))
    print(paste0(
      "Registered parallel backend is: ",
      getDoParName()))
    print(paste0(
      "Do parallel version is: ",
      getDoParVersion()))

    fit = foreach(
      bs=1:num_bs,
      .packages = c('quantreg'),
      .export = c('quantRegSpacing', 'rq.fit.sfn_start_val','findRedundantCols',
                  'ensureSpecFullRank', 'getRank', 'addMissingSpecColumns', 'clusterSample',
                  'rho', 'getColNums', 'quantRegress', 'printWarnings', 'ols_sparse_fit'),
      .combine = function(...) mapply('rbind',...),
      .errorhandling = 'remove' #, .options.snow = opts
    ) %dopar% {
      TOL <- 0.000000001 #make sure this matches earlier global parameter
      SMALL <- 1e-3

      # if no clustering is specified, draw sample here.
      # Need to provide a cluster_id variable to stratify
      if (is.null(cluster_indices)) {
        rows = sample(dim(data)[1],floor(M*dim(data)[1]))
        if (draw_weights){
          if(any(is.null(weight_vec))) rand_weight_vec <- pmax(rexp(length(rows)),5e-3)
          else {
            rand_weight_vec <- pmax(rexp(length(rows)),5e-3) * weight_vec[rows]
            if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
            else ols_rand_weight_vec <- rand_weight_vec
          }
        }
        else {
          rand_weight_vec <- weight_vec[rows]
          if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
          else ols_rand_weight_vec <- rand_weight_vec
        }
      } else{
        # call clusterSample here to generate clustered standard errors
        subsample_outputs <- clusterSample(cluster_indices = cluster_indices,
                                           stratum_indices = stratum_indices,
                                           M = M,
                                           draw_weights = draw_weights)
        rows <- subsample_outputs$subsample_rows
        if(any(is.null(weight_vec))) rand_weight_vec <- subsample_outputs$subsample_weights
        else{
          if(draw_weights) {
            rand_weight_vec <- subsample_outputs$subsample_weights * weight_vec[rows]
            if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
            else ols_rand_weight_vec <- rand_weight_vec
          }
          else {
            rand_weight_vec <- weight_vec[rows]
            if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
            else ols_rand_weight_vec <- rand_weight_vec
          }
        }
      }
      if (!exists("ols_rand_weight_vec")) ols_rand_weight_vec <- NULL
      if (!exists("rand_weight_vec"))     rand_weight_vec <- NULL


      cur_fit <- quantRegSpacing(
        data = data[rows,],
        dep_col = dep_col[rows],
        var_names = var_names,
        alpha = alpha,
        jstar = jstar,
        small = small,
        trunc = trunc,
        start_list = start_model,
        weight_vec = rand_weight_vec)
      cur_fit$OLS <- as.data.frame(t(ols_sparse_fit(a = data[rows,],
                                    y = dep_col[rows],
                                    weight_vec = ols_rand_weight_vec)))

      return(cur_fit)
    }
    stopCluster(cl)
  }
  else{
    # if no clustering is specified, draw sample here.
    # Need to provide a cluster_id variable to stratify
    if (is.null(cluster_indices)) {
      rows = sample(dim(data)[1],floor(M*dim(data)[1]))
      if (draw_weights){
        if(any(is.null(weight_vec))) rand_weight_vec <- pmax(rexp(length(rows)),5e-3)
        else {
          rand_weight_vec <- pmax(rexp(length(rows)),5e-3) * weight_vec[rows]
          if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
          else ols_rand_weight_vec <- rand_weight_vec
        }
      }
      else {
        rand_weight_vec <- weight_vec[rows]
        if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
        else ols_rand_weight_vec <- rand_weight_vec
      }
    } else{
      # call clusterSample here to generate clustered standard errors
      subsample_outputs <- clusterSample(cluster_indices = cluster_indices,
                                         stratum_indices = stratum_indices,
                                         M = M,
                                         draw_weights = draw_weights)
      rows <- subsample_outputs$subsample_rows
      if(any(is.null(weight_vec))) rand_weight_vec <- subsample_outputs$subsample_weights
      else{
        if(draw_weights) {
          rand_weight_vec <- subsample_outputs$subsample_weights * weight_vec[rows]
          if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
          else ols_rand_weight_vec <- rand_weight_vec
        }
        else {
          rand_weight_vec <- weight_vec[rows]
          if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
          else ols_rand_weight_vec <- rand_weight_vec
        }
      }
    }
    if (!exists("ols_rand_weight_vec")) ols_rand_weight_vec <- NULL
    if (!exists("rand_weight_vec"))     rand_weight_vec <- NULL
    fit <- quantRegSpacing(
      data = data[rows,],
      dep_col = dep_col[rows],
      var_names = var_names,
      alpha = alpha,
      jstar = jstar,
      small = small,
      trunc = trunc,
      start_list = start_model,
      weight_vec = rand_weight_vec)
    fit$OLS <- as.data.frame(t(ols_sparse_fit(a = data[rows,],
                                              y = dep_col[rows],
                                              weight_vec = ols_rand_weight_vec)))

    for(bs in 2:num_bs){
      # if no clustering is specified, draw sample here.
      # Need to provide a cluster_id variable to stratify
      if (is.null(cluster_indices)) {
        rows = sample(dim(data)[1],floor(M*dim(data)[1]))
        if (draw_weights){
          if(any(is.null(weight_vec))) rand_weight_vec <- pmax(rexp(length(rows)),5e-3)
          else {
            rand_weight_vec <- pmax(rexp(length(rows)),5e-3) * weight_vec[rows]
            if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
            else ols_rand_weight_vec <- rand_weight_vec
          }
        }
        else {
          rand_weight_vec <- weight_vec[rows]
          if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
          else ols_rand_weight_vec <- rand_weight_vec
        }
      } else{
        # call clusterSample here to generate clustered standard errors
        subsample_outputs <- clusterSample(cluster_indices = cluster_indices,
                                           stratum_indices = stratum_indices,
                                           M = M,
                                           draw_weights = draw_weights)
        rows <- subsample_outputs$subsample_rows
        if(any(is.null(weight_vec))) rand_weight_vec <- subsample_outputs$subsample_weights
        else{
          if(draw_weights) {
            rand_weight_vec <- subsample_outputs$subsample_weights * weight_vec[rows]
            if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
            else ols_rand_weight_vec <- rand_weight_vec
          }
          else {
            rand_weight_vec <- weight_vec[rows]
            if (square_ols_weights) ols_rand_weight_vec <- rand_weight_vec * weight_vec[rows]
            else ols_rand_weight_vec <- rand_weight_vec
          }
        }
      }
      if (!exists("ols_rand_weight_vec")) ols_rand_weight_vec <- NULL
      if (!exists("rand_weight_vec"))     rand_weight_vec <- NULL
      cur_fit <- quantRegSpacing(
        data = data[rows,],
        dep_col = dep_col[rows],
        var_names = var_names,
        alpha = alpha,
        jstar = jstar,
        small = small,
        trunc = trunc,
        start_list = start_model,
        weight_vec = rand_weight_vec)
      cur_fit$OLS <- as.data.frame(t(ols_sparse_fit(a = data[rows,],
                                                    y = dep_col[rows],
                                                    weight_vec = ols_rand_weight_vec)))
      fit = mapply('rbind', fit, cur_fit)
    }
  }

  quant_cov_mat <- cov(fit$coef)
  quant_cov_mat <- quant_cov_mat * (M)

  colnames(fit$OLS) <- var_names
  ols_cov_mat <- cov(fit$OLS)
  ols_cov_mat <- ols_cov_mat * (M)

  return(list('quant_cov' = quant_cov_mat,
              'ols_cov' = ols_cov_mat,
              'warnings' = fit$warnings,
              'iter' = fit$iter,
              'OLS' = fit$OLS,
              'counts' = fit$counts,
              'coef_boot' = fit$coef))
}




summary_output = function(
  quant_fit,
  ols_fit,
  se,
  alpha,
  varnames,
  ols_r_squared) {
  # Writes quantile and OLS output to a file.
  #
  # Args:
  #  quant_fit: Named rows of quantile beta estimates.
  #  ols_fit: Named rows of OLS beta estimates.
  #  se: Vector containing covariance matrices for both
  #      quantile and OLS fits.
  #  alpha: Column vector representing quantiles being estimated.
  #  varnames: vector of variable names for output
  #  ols_r_squared: R-squared for the OLS fit.
  #
  # Returns: Summary output matrix

  quant_betas <- quant_fit$coef
  pseudo_r <- quant_fit$pseudo_r

  quant_se <- diag(se$quant_cov)^(0.5)
  quant_out_vec <- c()
  quant_out_vec[seq(1,length(quant_betas)*2,2)] <- unlist(quant_betas)
  quant_out_vec[seq(2,length(quant_betas)*2,2)] <- unlist(quant_se)

  num_betas <- length(ols_fit)
  quant_out <- matrix(quant_out_vec, (num_betas*2), length(alpha))
  quant_out <- rbind(quant_out, pseudo_r)

  ols_se <- diag(se$ols_cov)^(0.5)
  ols_out <- c()
  ols_out[seq(1,length(ols_fit)*2,2)] <- unlist(ols_fit)
  ols_out[seq(2,length(ols_fit)*2,2)] <- unlist(ols_se)
  ols_out <- append(ols_out, ols_r_squared)

  final_output <- cbind(quant_out, ols_out)
  colnames(final_output) <- c(alpha, "OLS")

  r_names <- c()
  r_names[seq(1,num_betas*2,2)] <- varnames
  r_names[seq(2,num_betas*2,2)] <- paste(varnames, 'SE', sep = '_')
  r_names <- append(r_names, 'Pseudo-R^2 (Quantile); R^2 (OLS)')
  rownames(final_output) = r_names

  return(final_output)
}
