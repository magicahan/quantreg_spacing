# Set working directory to project root
this.dir <- paste(dirname(parent.frame(2)$ofile), "/../../..", sep="")
setwd(this.dir)

# Load configurable settings
source("code/R/utils/setup.R")

###################### Global configurations ######################

TEMP_DIR <-"temp"
LOG_DIR = "logs"
OUT_DIR = "output"
N_CLUSTER = 10
TOL <- 0.000000001
SMALL <- 1e-3
# Quantiles of interest
if (!exists("ALPHA")){
  ALPHA <- c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
}
if (!exists("JSTAR")){
  JSTAR <- 4
}
if (!exists("ALT_MIN_IND_DUMMY_SUM")){
  ALT_MIN_IND_DUMMY_SUM <- NULL
}
OLS_OUT <- TRUE

# HORIZON <- 3
# SPECIFICATION <- "ind_x_wage"

LOG_FILENAME = paste(
  LOG_DIR,
  "/",
  file = format(Sys.time(), "%m_%d_%y_"),
  "_Rquantlog_",SPECIFICATION,"_LOG",HORIZON,".txt",
  sep = '')

log_file <- file(LOG_FILENAME)
sink(file = log_file, append = TRUE, type = "output", split = TRUE)
sink(file = log_file, append = TRUE, type = "message")

###################### Libraries ######################

library(quantreg)
library(foreach)
library(doParallel)
library(data.table)
library(Matrix)
source("/Users/luxihan/Downloads/quantile_regression_utils.R")
source("/Users/luxihan/Downloads/ols_regression_utils.R")

source("code/R/utils/get_relevant_vcv_rows.R")
source("code/R/utils/quantile_regression_utils.R")
source("code/R/utils/ols_regression_utils.R")
source("/Users/luxihan/Downloads/quantile_regression_utils.R")
source("/Users/luxihan/Downloads/ols_regression_utils.R")

###################### Utilities ######################

getBetaMeans = function(x, betas, alpha, group_id_col=NULL, TRIM=.01) {
  # Computes means for various slices of a regression_spec by betas product.
  #
  # Args:
  #  x: Regression specification.
  #  betas: Fitted quantile coefficients for specification.
  #  alpha: Quantiles targeted.
  #  group_id_col: Optional vector to compute means for specific groups in
  #   addition to the whole specification.
  #  TRIM: Fraction (0 to 0.5) of observations to be trimmed from each end of
  #   a given column in the regression_spec by betas product before the mean is computed.
  #
  # Returns: Matrix with rows containing means for each quantile for various slices
  #  of the specification.
  #

  # First reshape the betas into wide format
  num_quantiles <- length(alpha)
  stopifnot(num_quantiles %% 2 == 1)

  num_betas <- length(betas) / num_quantiles
  betas_wide <- matrix(betas, num_betas, num_quantiles)

  x_beta_prod <- as.matrix(x %*% betas_wide)

  colnames(x_beta_prod) <- as.character(alpha)
  median_col_idx <- ((num_quantiles - 1) / 2) + 1
  x_beta_prod[,-median_col_idx] <- exp(x_beta_prod[,-median_col_idx])

  beta_means <- apply(x_beta_prod, 2, mean)
  beta_trimmed_means <- apply(x_beta_prod, 2, FUN = function(k) { mean(k, trim = TRIM) })

  means_out <- rbind(beta_means, beta_trimmed_means)
  rownames(means_out) <- c("FitQ_or_S_mean_full_spec", "FitQ_or_S_mean_full_spec_trimmed")

  if (!is.null(group_id_col)) {
    x_beta_prod <- as.data.table(x_beta_prod)
    cols_to_mean <- colnames(x_beta_prod)
    x_beta_prod <- cbind(x_beta_prod, group_id_col)
    colnames(x_beta_prod) <- c(cols_to_mean, "group_id")

    group_means <- x_beta_prod[,lapply(.SD, mean, na.rm=TRUE), by="group_id", .SDcols=cols_to_mean]
    group_trimmed_means <-
      x_beta_prod[,lapply(.SD, function(k) { mean(k, trim = TRIM, na.rm=TRUE) }), by="group_id", .SDcols=cols_to_mean]

    group_means_no_group_col <- as.matrix(group_means[,-1])
    group_trimmed_means_no_group_col <- as.matrix(group_trimmed_means[,-1])

    rownames(group_means_no_group_col) <- paste(
      "FitQ_or_S_mean_group_id_",
      unlist(group_means[,1]), sep="")
    rownames(group_trimmed_means_no_group_col) <- paste(
      "FitQ_or_S_mean_trimmed_group_id_",
      unlist(group_trimmed_means[,1]), sep="")

    means_out <- rbind(
      means_out,
      group_means_no_group_col,
      group_trimmed_means_no_group_col)
  }

  return(means_out)
}




getIndMeans = function(x, betas, varnames, alpha, group_id_col=NULL, TRIM=.01) {
  # Computes means for portion of linear index coming from industry dummies, overall and for subgroups
  #
  # Args:
  #  x: Regression specification.
  #  betas: Fitted quantile coefficients for specification.
  #  varnames: Vector of variable names
  #  alpha: Quantiles targeted.
  #  group_id_col: Optional vector to compute means for specific groups in
  #   addition to the whole specification.
  #  TRIM: Fraction (0 to 0.5) of observations to be trimmed from each end of
  #   a given column in the regression_spec by betas product before the mean is computed.
  #
  # Returns: Matrix with rows containing means for each quantile for various slices
  #  of the specification.
  #

  # First reshape the betas into wide format
  num_quantiles <- length(alpha)
  stopifnot(num_quantiles %% 2 == 1)

  num_betas <- length(betas) / num_quantiles
  betas_wide <- matrix(betas, num_betas, num_quantiles)

  # get industry columns & calculate industry fitted values
  ind_cols <- non_ind_dum_rows <- substr(varnames,1,12) == "factor(sicc)"
  ind_linidx_matrix <- as.matrix(x[,ind_cols] %*% betas_wide[ind_cols,])

  # add column names to ind_linidx_matrix
  colnames(ind_linidx_matrix) <- as.character(alpha)

  ind_means <- apply(ind_linidx_matrix, 2, mean)
  ind_trimmed_means <- apply(ind_linidx_matrix, 2, FUN = function(k) { mean(k, trim = TRIM) })

  means_out <- rbind(ind_means,ind_trimmed_means)
  rownames(means_out) <- c("IndFE_mean_full_spec", "IndFE_mean_full_spec_trimmed")

  if (!is.null(group_id_col)) {
    ind_linidx_matrix <- as.data.table(ind_linidx_matrix)
    cols_to_mean <- colnames(ind_linidx_matrix)
    ind_linidx_matrix <- cbind(ind_linidx_matrix, group_id_col)
    colnames(ind_linidx_matrix) <- c(cols_to_mean, "group_id")

    group_means <- ind_linidx_matrix[,lapply(.SD, mean, na.rm=TRUE), by="group_id", .SDcols=cols_to_mean]
    group_trimmed_means <-
      ind_linidx_matrix[,lapply(.SD, function(k) { mean(k, trim = TRIM, na.rm=TRUE) }), by="group_id", .SDcols=cols_to_mean]

    group_means_no_group_col <- as.matrix(group_means[,-1])
    group_trimmed_means_no_group_col <- as.matrix(group_trimmed_means[,-1])

    rownames(group_means_no_group_col) <- paste(
      "IndFE_mean_group_id_",
      unlist(group_means[,1]), sep="")
    rownames(group_trimmed_means_no_group_col) <- paste(
      "IndFE_mean_trimmed_group_id_",
      unlist(group_trimmed_means[,1]), sep="")

    means_out <- rbind(
      means_out,
      group_means_no_group_col,
      group_trimmed_means_no_group_col)
  }

  return(means_out)
}


out_to_csv = function(
  quant_fit,
  ols_fit,
  se,
  alpha,
  varnames,
  ols_r_squared,
  quantile_means,
  indfe_means,
  csv_filename) {
  # Writes quantile and OLS output to a file.
  #
  # Args:
  #  quant_fit: Named rows of quantile beta estimates.
  #  ols_fit: Named rows of OLS beta estimates.
  #  se: Vector containing covariance matrices for both
  #      quantile and OLS fits.
  #  alpha: Column vector representing quantiles being estimated.
  #  ols_r_squared: R-squared for the OLS fit.
  #  varnames: vector of variable names for output
  #  quantile_means: Means for various slices of regression_spec by betas product
  #  indfe_means: Means of industry FE by betas product
  #  csv_filename: Filename to use for output.
  #
  # Returns: Nothing. Writes OLS and quantile regression output to
  #          specified CSV file.
  #
  write.csv(se$warnings, file = paste(
    csv_filename,
    "_Warnings.csv",
    sep = ''))

  write.csv(se$iter, file = paste(
    csv_filename,
    "_Iters.csv",
    sep = ''))



  quant_betas <- quant_fit$coef
  pseudo_r <- quant_fit$pseudo_r

  quant_se <- diag(se$quant_cov)^(0.5)
  quant_out_vec <- c()
  quant_out_vec[seq(1,length(quant_betas)*2,2)] <- unlist(quant_betas)
  quant_out_vec[seq(2,length(quant_betas)*2,2)] <- unlist(quant_se)

  num_betas <- length(ols_fit)
  quant_out <- matrix(quant_out_vec, (num_betas*2), length(alpha))
  quant_out <- rbind(quant_out, pseudo_r)
  quant_out <- rbind(quant_out, quant_fit$counts)

  ols_se <- diag(se$ols_cov)^(0.5)
  ols_out <- c()
  ols_out[seq(1,length(ols_fit)*2,2)] <- unlist(ols_fit)
  ols_out[seq(2,length(ols_fit)*2,2)] <- unlist(ols_se)
  ols_out <- append(ols_out, ols_r_squared)
  ols_out <- append(ols_out, max(quant_fit$counts)) # assume that the jstar-th number of observations is the same as
                                                    # the OLS number of observations

  final_output <- cbind(quant_out, ols_out)
  colnames(final_output) <- c(alpha, "OLS")

  r_names <- c()
  r_names[seq(1,num_betas*2,2)] <- varnames
  r_names[seq(2,num_betas*2,2)] <- paste(varnames, 'SE', sep = '_')
  r_names <- append(r_names, c('Pseudo-R^2 (Quantile); R^2 (OLS)', 'Observation Count'))
  rownames(final_output) = r_names

  # for disclosure purposes, suppressing industry dummies for now
  non_ind_dum_rows <- (substr(r_names,1,12) != "factor(sicc)") & (substr(r_names,1,12) != "factor(year)")
  final_output <- final_output[non_ind_dum_rows,]

  #final_output <- rbind(final_output, cbind(quantile_means, 0), cbind(indfe_means, 0))
  final_output <- rbind(final_output, cbind(quantile_means, 0))

  write.csv(final_output, paste(csv_filename, "_Summary_Stats.csv", sep = ''),
            row.names = TRUE)
}

###################### Main script to run quantile regression ######################

# read in data
input_filename <- paste(TEMP_DIR, file="/reg_specV3.Rda", sep="")
reg_spec_df <- readRDS(input_filename)

orig_num_obs <- dim(reg_spec_df)[1]
print(paste0("Original obs count after loading reg_spec.Rda ", orig_num_obs))

# sorting by clusterEIN variable here
reg_spec_df = reg_spec_df[order(reg_spec_df$clusterEIN),]

if (!is.null(FILTER_COL)) {
  print(paste0(
    "Keeping all observations where ",
    FILTER_COL,
    " is ",
    paste(FILTER_VALUES, collapse=" ")))
  reg_spec_df <- reg_spec_df[reg_spec_df[,FILTER_COL] %in% FILTER_VALUES,]
  cur_num_obs <- dim(reg_spec_df)[1]
  print(paste0("Obs count after filtering ", cur_num_obs))
}

# selecting dependent variable here
if (HORIZON == 3) {
  reg_spec_df$reswagech <- reg_spec_df$reswagech_log3
}
if (HORIZON == 5) {
  reg_spec_df$reswagech <- reg_spec_df$reswagech_log5
}
if (HORIZON == 10) {
  reg_spec_df$reswagech <- reg_spec_df$reswagech_log10
}

print(paste("Initial obs count: sorting variables only = ",dim(reg_spec_df)[1]))

reg_spec_df <- reg_spec_df[complete.cases(reg_spec_df$reswagech),] #make sure all dependent variables have values
print(paste("Obs count: sorting + dependent variables = ",dim(reg_spec_df)[1]))

# handling weights here (to prevent very large clusters from dominating results)
if(WEIGHT_BY_CLUSTER_SIZE){
  # change to data table first
  reg_spec_df = as.data.table(reg_spec_df)
  # get average firm size by year
  weights = reg_spec_df[, pmax(.N,50),by = c('clusterEIN','year')]
  avgFirmSize = weights[,mean(V1), by = 'year']
  colnames(avgFirmSize)[2] = 'avgFirmSize'

  weights = merge(weights, avgFirmSize, by = 'year')
  weights[,'weight'] = sqrt(pmax(weights[,'avgFirmSize'] / weights[,'V1']*10,0.01))
  # square root is handled within the OLS function
  weights[,'weight_ols'] = weights[,'weight']^2

  reg_spec_df = merge(reg_spec_df, weights[,c('year', 'clusterEIN', 'weight','weight_ols')], by = c('year', 'clusterEIN'), all.x = TRUE)
  # turn back into a dataframe
  reg_spec_df = as.data.frame(reg_spec_df)
}



# Constrain SIC codes to 3 digits
reg_spec_df$sicc <- trunc(reg_spec_df$sicc / 10)

# Drop observations with NA values for age / lagged wages
reg_spec_df <- reg_spec_df[complete.cases(reg_spec_df$age),]
print(paste("Obs count after dropping missing age measures = ",dim(reg_spec_df)[1]))
reg_spec_df <- reg_spec_df[complete.cases(reg_spec_df$avgreswage),]
print(paste("Obs count after dropping missing lagged inc measures = ",dim(reg_spec_df)[1]))

# Adding an additional check for the income growth rate control, when it is present
if (any(grepl('incgth', FIRM_CONTROLS))) {
  reg_spec_df <- reg_spec_df[complete.cases(reg_spec_df$reswagech_log3_lagged_3_years),]
  print(paste("Obs count after dropping missing lagged inc growth measures = ",dim(reg_spec_df)[1]))
}



# building regression specification formula
# Age X industry wage rank controls here
age_wage_ctrl_list <- paste("V",1:((HERMITE_POL_DEGREE+1)*length(AGE_GRID)),sep="")
age_wage_ctrl_list <- paste(age_wage_ctrl_list, collapse = " + ")

# firm controls next
firm_ctrl_list <- paste(FIRM_CONTROLS, collapse = " + ")

# main coefficients of interest next
if (GROUP_ID == "") {
  main_coeff_list <- MAIN_COEFFS
  group_id_col <- NULL
} else {
  for (j in 1:length(MAIN_COEFFS))  {
    if (j == 1) {
      main_coeff_list <- paste("factor(", GROUP_ID, "):", MAIN_COEFFS[j], sep="")
    } else {
      main_coeff_list <- c(main_coeff_list,
                           paste("factor(", GROUP_ID, "):", MAIN_COEFFS[j], sep=""))
    }
  }
}
main_coeff_list <- paste(main_coeff_list, collapse = " + ")


# this is a counter variable that will be used to identify the complete rows of the data
reg_spec_df$nmiss_ind <- seq(1,dim(reg_spec_df)[1],1)


reg_formula <-as.formula(paste("reswagech ~ nmiss_ind + factor(sicc) + factor(year) + ",
                               main_coeff_list," + ",
                               firm_ctrl_list, " + ", age_wage_ctrl_list))

reg_formula_noind <-as.formula(paste("reswagech ~ nmiss_ind + factor(year) + ",
                                     main_coeff_list," + ",
                                     firm_ctrl_list, " + ", age_wage_ctrl_list))


# recoding SIC3 code with largest number of clusters as -1
# this will make it the base category in the regressions

# getting list of clusters here
ind_cluster_list <- as.data.frame(table(reg_spec_df[, c("sicc","clusterEIN")]))
ind_cluster_list <- ind_cluster_list[ind_cluster_list$Freq > 0,]

sic3_cts <- as.data.frame(table(ind_cluster_list$sicc))
max_sic <- min(as.numeric(as.character(sic3_cts[
  sic3_cts[,2]==max(sic3_cts[,2]),1])))
reg_spec_df[reg_spec_df$sicc == max_sic,'sicc'] <- -1

print('Regression formula is (minus nmiss_ind, which is a helper variable):')
print(reg_formula)

reg_spec_with_dummies <- model.matrix(reg_formula_noind, reg_spec_df[!is.na(reg_spec_df$sicc),])

# triggering garbage collection here
for (i in 1:10) {gc()}
# print("current memory usage")
# mem_used()


# extracting list of non-missing rows here
final_nmiss_ind <- reg_spec_with_dummies[,2]

indFormula = as.formula("reswagech ~  factor(sicc)")
regspecIndustry = sparse.model.matrix(object = indFormula, data = reg_spec_df[final_nmiss_ind,])
# drop constant term from industry dummy matrix
regspecIndustry <- regspecIndustry[,-1]

if (GROUP_ID != "") {
  # setting up group_id_col, which is needed for computing avg marginal effects
  group_id_col <- reg_spec_df[final_nmiss_ind, GROUP_ID]
  # sanity check - make sure that number of obs line up
  stopifnot(dim(reg_spec_with_dummies)[1] == length(group_id_col))
}
reg_spec_df <- reg_spec_df[final_nmiss_ind,]

# Drop intercept to avoid collinearity with age lagged income variables
# also, dropping final_nmiss_ind column from regression matrix
reg_spec_with_dummies <- reg_spec_with_dummies[,-c(1,2)]


# triggering garbage collection here
for (i in 1:10) {gc()}


print(paste("Obs count after dropping other missing control variables = ",dim(reg_spec_with_dummies)[1]))

num_ind_dummies <- length(unique(reg_spec_df$sicc)) - 1

# Drop dummies with column sums lower than MIN_IND_DUMMY_SUM
# (to avoid near-singularity issues in quantreg's sfn implementation)
ind_col_sums <- colSums(regspecIndustry[, 1:num_ind_dummies])
print("Industry dummy column summary stats:")
print(summary(ind_col_sums))
print(paste("SD: ", sd(ind_col_sums), sep=""))

# Combine and convert to sparse matrix
reg_spec_sparse <- cbind(as(regspecIndustry,"matrix.csr"),denseMatrixToSparse(reg_spec_with_dummies))
reg_spec_var_names <- c(colnames(regspecIndustry),colnames(reg_spec_with_dummies))
rm(reg_spec_with_dummies)
rm(regspecIndustry)




# print("current memory usage - after dropping dense matrix")
# mem_used()


if (is.null(ALT_MIN_IND_DUMMY_SUM)){
  dummies_to_drop <- which(ind_col_sums < MIN_IND_DUMMY_SUM)
} else {
  dummies_to_drop <- which(ind_col_sums < ALT_MIN_IND_DUMMY_SUM)
}

#in case no dummies are to be dropped
if(!is.null(dummies_to_drop)){
  #Here, we drop ROWS of industries that are also dropped.
  rows_to_drop <- as.matrix(reg_spec_sparse[,dummies_to_drop] %*%
                              matrix(1,length(dummies_to_drop),1)) == 1
  # Note that we could probably use an 'Other' catergory instead.

  # Here, we drop the columns
  reg_spec_var_names = reg_spec_var_names[-dummies_to_drop]
  reg_spec_sparse <- reg_spec_sparse[,-dummies_to_drop]


  print(paste0("Note: ",sum(rows_to_drop)," rows will be dropped from small industries."))

  reg_spec_sparse <- reg_spec_sparse[which(!rows_to_drop),]
  reg_spec_df <- reg_spec_df[which(!rows_to_drop),]
  if (GROUP_ID != "") {
    group_id_col <- group_id_col[which(!rows_to_drop)]
    # sanity check - make sure that number of obs line up
    stopifnot(dim(reg_spec_sparse)[1] == length(group_id_col))
  }
  print(paste("Obs count after dropping small industries = ",dim(reg_spec_sparse)[1]))

}

# dropping lots of extraneous columns from reg spec df
if (is.null(reg_spec_df$weight)) {
  reg_spec_df <- reg_spec_df[,c("reswagech","clusterEIN","cluster_size_bin")]
} else {
  reg_spec_df <- reg_spec_df[,c("reswagech","weight","weight_ols","clusterEIN","cluster_size_bin")]
}

# triggering garbage collection here
for (i in 1:10) {gc()}
# print("current memory usage - after dropping dense matrix")
# mem_used()


#### create cluster_indices and stratum indices matrices

uniqClusters <- table(reg_spec_df[, c("clusterEIN","cluster_size_bin")])
uniqClusters <- as.data.frame(uniqClusters)
# drop zeros and NAs
uniqClusters <- uniqClusters[uniqClusters$Freq > 0 & !is.na(uniqClusters$Freq) ,]

# re-ordering by clusterEIN to match order of data
uniqClusters = uniqClusters[order(uniqClusters$clusterEIN),]

# now, populating list of start and end indices
cluster_indices <- matrix(0,dim(uniqClusters)[1],2)
for (jj in 1:dim(uniqClusters)[1]){
  if (jj==1){
    cluster_indices[jj,1] <- 1
    cluster_indices[jj,2] <- uniqClusters[jj,"Freq"]
  } else {
    cluster_indices[jj,1]<- 1 + sum(uniqClusters[1:(jj-1),"Freq"])
    cluster_indices[jj,2]<- sum(uniqClusters[1:(jj),"Freq"])
  }
}
# peeling off stratum_indices vector
stratum_indices <- as.vector(uniqClusters$cluster_size_bin)

# test: call clusterSample function here...
# clusterSample(cluster_indices = cluster_indices,
#               stratum_indices = stratum_indices,
#               M = SUBSAMPLE_PCT*0.01,
#               draw_weights = FALSE)

# Winsorize the data before gathering point estimates and standard errors

percent = 0.001
lowerBound = quantile(reg_spec_df$reswagech, percent)
upperBound = quantile(reg_spec_df$reswagech, 1 - percent)

# replace values less than lower 0.1 % of the data with the lower bound
reg_spec_df$reswagech[which(reg_spec_df$reswagech < lowerBound)] = lowerBound
# replace values greater than upper 99.9 % of the data with the upper bound
reg_spec_df$reswagech[which(reg_spec_df$reswagech > upperBound)] = upperBound

# Prints warnings as they occur
options(warn=1)
print("Calculating initial quantile fit")
ptm <- proc.time()
quantreg_fit <- quantRegSpacing(
  dep_col = reg_spec_df$reswagech,
  data = reg_spec_sparse,
  var_names = reg_spec_var_names,
  alpha = ALPHA,
  jstar = JSTAR,
  trunc = T, small = 1e-3, weight_vec = reg_spec_df$weight)
print(proc.time() - ptm)
print("Initial quantile fit complete")

print("Calculating initial OLS fit")
ptm <- proc.time()
ols_fit <- ols_sparse_fit(
  a = reg_spec_sparse,
  y = reg_spec_df$reswagech,
  weight_vec = reg_spec_df$weight)
print(proc.time() - ptm)
print("Initial OLS fit complete")

ols_r_squared <- get_ols_r_squared(
  a = reg_spec_sparse,
  y = reg_spec_df$reswagech,
  betas = ols_fit,
  weight_vec = reg_spec_df$weight_ols)

print("Starting subsampling")
ptm <- proc.time()
se = subsampleStandardErrors(
  dep_col = reg_spec_df$reswagech,
  data = reg_spec_sparse,
  var_names = reg_spec_var_names,
  alpha = ALPHA,
  jstar = JSTAR,
  M = SUBSAMPLE_PCT * 0.01,
  cluster_indices = cluster_indices,
  stratum_indices = factor(stratum_indices),
  draw_weights = FALSE,
  num_bs = NUM_BOOTSTRAP,
  parallelize = PARALLELIZE_BOOTSTRAP,
  num_cores = NUM_CORES_BOOTSTRAP,
  trunc = TRUE,
  start_model = quantreg_fit$coef,
  small = 1e-3,
  weight_vec = reg_spec_df$weight,
  square_ols_weights = TRUE)
print(proc.time() - ptm)
print("Subsampling complete")

print("Calculating average fitted quantiles for average marginal effects")
quantile_means <- getBetaMeans(
  reg_spec_sparse,
  t(quantreg_fit$coef),
  ALPHA,
  group_id_col=group_id_col)

print("Calculating average value of industry dummies")
indfe_means <- getIndMeans(
  reg_spec_sparse,
  t(quantreg_fit$coef),
  varnames = reg_spec_var_names,
  ALPHA,
  group_id_col=group_id_col)

vcv_relevant_subset <- getRelevantVCVRows(as.data.frame(se$quant_cov))

# Write output
csv_filename = paste(
  OUT_DIR,
  "/",
  file = format(Sys.time(), "%m_%d_%y_"),SPECIFICATION,"_LOG",HORIZON, sep = '')

write.csv(
  vcv_relevant_subset,
  file=paste0(csv_filename, "_vcv_relevant_rows_only.csv"))

out_to_csv(
  quant_fit = quantreg_fit,
  ols_fit = ols_fit,
  se = se,
  alpha = ALPHA,
  varnames = reg_spec_var_names,
  ols_r_squared = ols_r_squared,
  quantile_means = quantile_means,
  indfe_means = indfe_means,
  csv_filename = csv_filename)

# writing coefficients and VCV matrix to a temp directory

csv_backup_filename = paste(
  OUT_DIR,
  "/",
  file = format(Sys.time(),"DO_NOT_EMAIL_", "%m_%d_%y_"),SPECIFICATION,"_LOG",HORIZON, sep = '')

write.csv(quantreg_fit$coef, file = paste(
  csv_backup_filename,
  "_AllQREGCoeffs.csv",
  sep = ''))

write.csv(se$quant_cov, file = paste(
  csv_backup_filename,
  "_QREGVCVMatrix.csv",
  sep = ''))

print(paste('Finished Analysis on', dim(reg_spec_sparse)[1], 'Observations'))

# return output back to the console
sink()
sink(type = "message")
