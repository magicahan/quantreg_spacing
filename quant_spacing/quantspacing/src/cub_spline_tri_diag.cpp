#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_qrng.h>
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <string>
#include <Rmath.h>
#include <math.h>
#include<cmath>
#include <stdlib.h>
#include <ctime>
#include <fstream>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppGSL)]]
#define MIN(x,y) x < y ? x : y

void quad_form(double a, double b, double discriminant, double * root1, double * root2) {
  // Calculate the roots of a 2nd degree polynomial using the quadratic formula and the already calculated discriminant.
  *root1 = (-b + sqrt(discriminant)) / (2*a);
  *root2 = (-b - sqrt(discriminant)) / (2*a);
}


arma::mat init_matrix(double num, int m, int n){
  arma::mat result(m ,n);
  result.zeros();
  for(int i = 0; i < m; i++ ){
    for(int j = 0; j < n; j++){
      result(i, j) = num;
    }
  }
  return result;
}

int find(double data, Rcpp::NumericVector array, int max) {
  // Currently an O(N) search algorithm, but since N will usually be less than 10, it might be more trouble than it's worth to implement something asymptotically faster.
  if (data <= array[0]) {
    return -1;
  }
  if (array[max-1] <= data) {
    return max-1;
  }
  for (size_t i = 0; i < max-1; i++) {
    if (array[i] <= data && data <= array[i+1]) {
      return i;
    }
  }
  return -2;
}


int find_arma(double data, arma::vec array, int max) {
  // Currently an O(N) search algorithm, but since N will usually be less than 10, it might be more trouble than it's worth to implement something asymptotically faster.
  if (data <= array(0)) {
    return -1;
  }
  if (array(max-1) <= data) {
    return max-1;
  }
  for (size_t i = 0; i < max-1; i++) {
    if (array(i) <= data && data <= array(i+1)) {
      return i;
    }
  }
  return -2;
}


// [[Rcpp::export]]
Rcpp::NumericMatrix c_spline(Rcpp::NumericMatrix quantiles, Rcpp::NumericVector alphas, Rcpp::NumericVector yp1, Rcpp::NumericVector ypp) {
  // Perform the tri-diagonalization two-time pass through to compute the second derivatives of the cubic spline

  int N = quantiles.nrow();
  int p = quantiles.ncol();
  Rcpp::NumericMatrix y2(N, p);

  // Loop through the rows (no need to use vectors in C!), row designated by i. 0-indexed!
  for (size_t i = 0; i < N; i++) {
    int j;
    double d, qp, sig, up;
    double * u = (double *)calloc(sizeof(double), p-1);
    Rcpp::NumericMatrix::Row quantiles_i = quantiles(i, Rcpp::_);
    Rcpp::NumericMatrix::Row y2_i = y2(i, Rcpp::_);

    // Use a specific first derivative at the first quantile
    y2_i[0] = -0.5;
    u[0] = (3.0 / (quantiles_i[1] - quantiles_i[0])) * ((alphas[1]-alphas[0]) / (quantiles_i[1]-quantiles_i[0]) - yp1[i]);

    for (j = 1; j <= p-2; j++) {
      sig = (quantiles_i[j] - quantiles_i[j-1]) / (quantiles_i[j+1] - quantiles_i[j-1]);
      d = sig * y2_i[j-1] + 2.0;
      y2_i[j] = (sig - 1.0) / d;
      u[j] = (alphas[j+1]-alphas[j]) / (quantiles_i[j+1] - quantiles_i[j]) - (alphas[j]-alphas[j-1]) / (quantiles_i[j] - quantiles_i[j-1]);
      u[j] = (6.0 * u[j] / (quantiles_i[j+1] - quantiles_i[j-1]) - sig * u[j-1]) / d;
      // printf("%f ", u[i]);
    }
    // printf("\n");

    // Use a specific first derivative at the last quantile
    qp = 0.5;
    up = (3.0 / (quantiles_i[p-1] - quantiles_i[p-2])) * (ypp[i] - ((alphas[p-1]-alphas[p-2]) / (quantiles_i[p-1]-quantiles_i[p-2])));

    y2_i[p-1] = (up - (qp * u[p-2])) / ((qp * y2_i[p-2]) + 1.0);

    for (j = p-2; j >= 0; j--) {
      y2_i[j] = (y2_i[j] * y2_i[j+1]) + u[j];
      // printf("%f ", y2[i]);
    }
    // printf("\n");

    free(u);
    // memcpy(y2[j], y2_j, sizeof(double)*p);
  }
  return y2;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix matrix_mult(Rcpp::NumericMatrix A, Rcpp::NumericMatrix B){
  int m = A.nrow();
  int n = A.ncol();
  int p = B.ncol();
  double temp;
  Rcpp::NumericMatrix C(m, p);
  for (int i = 0; i < m; i++){
    for(int j = 0; j < p; j++){
      temp = 0;
      for(int k = 0; k < n; k ++){
        temp = temp + A(i, k) * B(k, j);
      }
      C(i, j) = temp;
    }
  }
  return C;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix transpose(Rcpp::NumericMatrix A){
  int m = A.nrow();
  int n = A.ncol();
  Rcpp::NumericMatrix B(n, m);
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      B(j, i) = A(i, j);
    }
  }
  return B;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix invert_two(Rcpp::NumericMatrix A){
  double det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
  Rcpp::NumericMatrix B(2, 2);
  B(0, 0) = A(1, 1);
  B(0, 1) = -A(0, 1);
  B(1, 0) = -A(1, 0);
  B(1, 1) = A(0, 0);
  B = 1/det * B;
  return B;
}

int non_neg(double a, double b, double c, double q_lo, double q_hi) {
  double discriminant = 4*b*b - 12 * a * c;
  // printf("disc: %f\n", discriminant);
  if (discriminant >= 0) {
    // The slope was negative at some point, now we check if it is negative inside the interval.

    // Find the roots of the cubic function stored in the roots array.
    double * roots = (double *)malloc(2 * sizeof(double));
    gsl_poly_solve_quadratic(3*a, 2*b, c, &roots[0], &roots[1]);
    // printf("%f %f\n", roots[0], roots[1]);

    // int leading = (a > 0) ? 1 : 0;
    // printf("%d\n", leading);

    // If neither of the roots lie in the interval, the cubic must be positive on the interval because the interval is outside the zeros of the quadratic and so the cubic is increasing on the interval because the leading coefficient is always positive.
    if ((roots[0] >= q_lo && roots[0] <= q_hi) || (roots[1] >= q_lo && roots[1] <= q_hi)) {
      // The function was decreasing on the interval, so linearly-interpolate between the quantiles.
      // printf("in int\n");
      return 1;
    }
  }
  return 0;
}

// [[Rcpp::export]]
Rcpp::NumericVector c_splint(Rcpp::NumericVector y, Rcpp::NumericMatrix quantiles, Rcpp::NumericVector alphas, Rcpp::NumericMatrix y2, Rcpp::NumericMatrix tail_param_u, Rcpp::NumericMatrix tail_param_l, std::string tails, Rcpp::IntegerVector distn) {
  int N = quantiles.nrow();
  int p = quantiles.ncol();
  Rcpp::NumericVector y_hat(N);

  // Loop through the rows (no need to use vectors in C!), row designated by i. 0-indexed!
  for (size_t i = 0; i < N; i++) {
    Rcpp::NumericMatrix::Row quantiles_i = quantiles(i, Rcpp::_);
    Rcpp::NumericMatrix::Row y2_i = y2(i, Rcpp::_);
    int j;

    // Quantile function:
    if (distn[i] == 3) {
      // Search for the correct interval (based on alphas, not quantiles, since we are passed a probability)
      j = find(y[i], alphas, p);
      if (j == -1) {
        // Probability was lower than all the quantiles, so use the tail distribution
        if (tails == "gaussian") {
          y_hat[i] = gsl_cdf_gaussian_Pinv(y[i], tail_param_l(i, 1)) + tail_param_l(i, 0);
          // printf("low: %f %f\n", y[i], y_hat[i]);
        } else if (tails == "exponential") {
          y_hat[i] = tail_param_l(i, 0) + tail_param_l(i, 1) * log(y[i]);
        }
      } else if (j == (p-1)) {
        // Probability was higher than all the quantiles, so use the tail distribution
        if (tails == "gaussian") {
          y_hat[i] = gsl_cdf_gaussian_Pinv(y[i], tail_param_u(i, 1)) + tail_param_u(i, 0);
          // printf("hi: %f %f\n", y[i], y_hat[i]);
        } else if (tails == "exponential") {
          y_hat[i] = tail_param_u(i, 0) - tail_param_u(i, 1) * log(1 - y[i]);
        }
      } else {
        // Probability was in an interval, designated by j. 0-indexed!
        // printf("mid\n");
        // h is a helper quantity that is used frequently.
        // a, b, c, and d are the coefficients of the cubic function on the interval in question, with a the 3rd degree and d the 0th.
        double h = quantiles_i[j+1] - quantiles_i[j];
        double a = (y2_i[j+1] - y2_i[j]) / (6*h);
        double b = (quantiles_i[j+1]*y2_i[j] - quantiles_i[j]*y2_i[j+1]) / (2*h);
        double c = (alphas[j+1] - alphas[j]) / h + h * (y2_i[j] - y2_i[j+1]) / 6 + (quantiles_i[j] * quantiles_i[j] * y2_i[j+1] - quantiles_i[j+1] * quantiles_i[j+1] * y2_i[j]) / (2*h);
        double d = -((alphas[j+1] * quantiles_i[j])/h) + (alphas[j] * quantiles_i[j+1])/h - (h*quantiles_i[j+1] * y2_i[j])/6 + (quantiles_i[j+1]*quantiles_i[j+1]*quantiles_i[j+1] * y2_i[j])/(6*h) + (h*quantiles_i[j]*y2_i[j+1])/6 - (quantiles_i[j]*quantiles_i[j]*quantiles_i[j] * y2_i[j+1])/(6*h);
        // Divide the coefficients through by a to be able to pass to the solver. Subtract d by the given y to make the roots of the function the solutions.
        double b2 = b / a;
        double c2 = c / a;
        double d2 = (d - y[i]) / a;
        // printf("b2: %f, c2: %f, d2: %f\n", b2, c2, d2);
        double * roots = (double *)malloc(3 * sizeof(double));
        // Find the roots of the cubic function stored in the roots array.
        int num_roots = gsl_poly_solve_cubic(b2, c2, d2, &roots[0], &roots[1], &roots[2]);
        // printf("num_roots: %d\n", num_roots);
        // Cubic functions either have 1 root or 3 roots, test for both cases.
        if (num_roots == 1) {
          // If there is 1 root, it must lie in the interval, so return this root.
          y_hat[i] = roots[0];
          free(roots);
        } else {
          // If there are 3 roots, check which one lies in the interval.
          int boundaries = 1;
          for (size_t k = 0; k < 3; k++) {
            // printf("root%lu: %f\n", k, roots[k]);
            if (quantiles_i[j] <= roots[k] && roots[k] <= quantiles_i[j+1]) {
              y_hat[i] = roots[k];
              boundaries = 0;
              free(roots);
              break;
            }
          }
          // Makes sure edge cases near quantiles don't fall outside the interval due to slight approximations (just take the closest root to the interval).
          if (boundaries) {
            size_t min_index;
            double min_dist = DBL_MAX;
            double dist;
            for (size_t k = 1; k < 3; k++) {
              dist = MIN(fabs(roots[k] - quantiles_i[j]), fabs(roots[k] - quantiles_i[j+1]));
              if (dist < min_dist) {
                min_dist = dist;
                min_index = k;
              }
            }
            y_hat[i] = roots[min_index];
            free(roots);
          }
        }
      }
      continue;
    }
    // Search for the correct interval on which to approximate the distribution.
    j = find(y[i], quantiles_i, p);
    // printf("%d %f\n", j, y[i]);
    if (j == -1) {
      // Value was lower than all the quantiles, so use the tails distribution to approximate.
      if (tails == "gaussian") {
        if (distn[i] == 1) { // CDF
          y_hat[i] = gsl_cdf_gaussian_P(y[i] - tail_param_l(i, 0), tail_param_l(i, 1));
        } else { // PDF
          y_hat[i] = gsl_sf_erf_Z((y[i] - tail_param_l(i, 0)) / tail_param_l(i, 1)) / tail_param_l(i, 1);
        }
      } else if (tails == "exponential") {
        if (distn[i] == 1) { // CDF
          y_hat[i] = pow(M_E, ((y[i] - tail_param_l(i, 0)) / tail_param_l(i, 1)));
        } else { // PDF
          y_hat[i] = pow(M_E, ((y[i] - tail_param_l(i, 0)) / tail_param_l(i, 1))) / tail_param_l(i, 1);
        }
      }
    } else if (j == (p-1)) {
      // Value was higher than all the quantiles, so use the tails distribution to approximate.
      if (tails == "gaussian") {
        if (distn[i] == 1) {
          y_hat[i] = gsl_cdf_gaussian_P(y[i] - tail_param_u(i, 0), tail_param_u(i, 1));
        } else {
          y_hat[i] = gsl_sf_erf_Z((y[i] - tail_param_u(i, 0)) / tail_param_u(i, 1)) / tail_param_u(i, 1);
        }
      } else if (tails == "exponential") {
        if (distn[i] == 1) {
          y_hat[i] = 1 - pow(M_E, (-(y[i] - tail_param_u(i, 0)) / tail_param_u(i, 1)));
        } else {
          y_hat[i] = pow(M_E, (-(y[i] - tail_param_u(i, 0)) / tail_param_u(i, 1))) / tail_param_u(i, 1);
        }
      }
    } else {
      // Value was in an interval, designated by j. 0-indexed!
      double h = quantiles_i[j+1] - quantiles_i[j];
      double a = (quantiles_i[j+1] - y[i]) / h;
      double b = (y[i] - quantiles_i[j]) / h;
      // printf("h: %f %f %f\n", h, quantiles_i[j+1], quantiles_i[j]);
      // printf("a: %f\n", a);
      // printf("b: %f\n", b);
      // Calculate the discriminant to make sure the cubic function isn't decreasing.
      double discriminant = (1/(3*h*h))*(-6*alphas[j]*y2_i[j] + 6*alphas[j+1]*y2_i[j] + h*h*y2_i[j]*y2_i[j] + 6*alphas[j]*y2_i[j+1] - 6*alphas[j+1]*y2_i[j+1] - 2*h*h*y2_i[j]*y2_i[j+1] + 3*quantiles_i[j]*quantiles_i[j]*y2_i[j]*y2_i[j+1] - 6*quantiles_i[j]*quantiles_i[j+1]*y2_i[j]*y2_i[j+1] + 3*quantiles_i[j+1]*quantiles_i[j+1]*y2_i[j]*y2_i[j+1] + h*h*y2_i[j+1]*y2_i[j+1]);
      // printf("disc: %f\n", discriminant);
      if (discriminant >= 0) {
        // The slope was negative at some point, now we check if it is negative inside the interval.
        double qa = y2_i[j+1]/(2*h) - (y2_i[j]/(2*h));
        double qb = (quantiles_i[j+1]*y2_i[j])/h - (quantiles_i[j]*y2_i[j+1])/h;
        double root1 = (-qb + sqrt(discriminant)) / (2*qa);
        double root2 = (-qb - sqrt(discriminant)) / (2*qa);
        // printf("roots: %f %f\n", *root1, *root2);
        if ((root1 >= quantiles_i[j] && root1 <= quantiles_i[j+1]) || (root2 >= quantiles_i[j] && root2 <= quantiles_i[j+1])) {
          // The function was decreasing on the interval, so linearly-interpolate between the quantiles.
          // printf("in int\n");
          if (distn[i]) { // CDF
            y_hat[i] = (y[i] - quantiles_i[j]) * (alphas[j+1] - alphas[j]) / h + alphas[j];
          } else { // PDF
            y_hat[i] = (alphas[j+1] - alphas[j]) / h;
          }
          continue;
        }
      }
      // printf("pos slope\n");
      // Everything checked out, so splint together the cubic function using the second derivatives and the location of y inside the interval.
      if (distn[i] == 1) { // CDF
        // printf("distn\n");
        y_hat[i] = a*alphas[j] + b*alphas[j+1] + ((a*a*a-a)*y2_i[j] + (b*b*b-b)*y2_i[j+1]) * (h*h) / 6;
      } else { // PDF (just the derivative of the CDF)
        // printf("pdf\n");
        double b_prime = 1 / h;
        double a_prime = -b_prime;
        y_hat[i] = a_prime*alphas[j] + b_prime*alphas[j+1] + ((3*(a*a)*a_prime-a_prime)*y2_i[j] + (3*(b*b)*b_prime-b_prime)*y2_i[j+1]) * (h*h) / 6;
      }
    }
  }
  return y_hat;
}


// [[Rcpp::export]]
arma::mat spacingsToQuantiles_C(arma::mat spacingCoef, arma::mat data, int jstar){
	int p = spacingCoef.n_cols;
	int N = data.n_rows;
	arma::mat quantiles(N, p);
	arma::mat starResids = data*(arma::mat)spacingCoef(arma::span::all, arma::span(jstar-1,jstar-1));
	quantiles(arma::span::all, jstar-1) = starResids(arma::span::all, 0);

	arma::mat resids(N, 1);

	resids(arma::span::all, arma::span(0, 0)) = starResids(arma::span::all, arma::span(0, 0));

	for (int j = jstar; j < p; j++) {
	arma::mat spacing = data*spacingCoef(arma::span::all, arma::span(j,j));
	spacing(arma::span::all, arma::span(0,0)) = exp(spacing(arma::span::all, arma::span(0,0)));
	resids(arma::span::all, arma::span(0, 0)) = resids(arma::span::all, arma::span(0,0)) + spacing(arma::span::all, arma::span(0,0));
	quantiles(arma::span::all, arma::span(j, j)) = resids(arma::span::all, arma::span(0,0));
	}


	resids(arma::span::all, arma::span(0,0)) = starResids(arma::span::all, arma::span(0,0));

	for (int j = jstar - 2; j >= 0; j--) {
	  arma::mat spacing = data*spacingCoef(arma::span::all, arma::span(j,j));
	  spacing(arma::span::all, arma::span(0,0)) = exp(spacing(arma::span::all, arma::span(0,0)));
	  resids(arma::span::all, arma::span(0, 0)) = resids(arma::span::all, arma::span(0,0)) - spacing(arma::span::all, arma::span(0,0));
	  quantiles(arma::span::all, arma::span(j, j)) = resids(arma::span::all, arma::span(0,0));
	}

	return quantiles;
}


// [[Rcpp::export]]
arma::mat c_spline_quant_comb(arma::vec y, arma::mat quantiles, arma::vec alphas,  std::string tails, arma::vec distn) {
  // Perform the tri-diagonalization two-time pass through to compute the second derivatives of the cubic spline

  int N = quantiles.n_rows;
  int p = quantiles.n_cols;

  arma::vec q_shift = quantiles.col(0);
  arma::vec q_stretch = quantiles.col(p-1) - quantiles.col(0);
  // q_stretch.fill(1.0);
  // for(int i = 0; i < N; i++){
  //   q_stretch(i) = 1;
  // }
  for (size_t j = 0; j < p; j++) {
    quantiles.col(j) = (quantiles.col(j) - q_shift) / q_stretch;
  }

  arma::mat temp_l(2, 2, arma::fill::ones);
  if (tails == "gaussian") {
    temp_l(0, 1) = gsl_cdf_gaussian_Pinv(alphas(0), 1);
    temp_l(1, 1) = gsl_cdf_gaussian_Pinv(alphas(1), 1);
  } else if (tails == "exponential") {
    temp_l(0, 1) = gsl_cdf_exponential_Pinv(alphas(0), 1);
    temp_l(1, 1) = gsl_cdf_exponential_Pinv(alphas(1), 1);
  }
  arma::mat tail_param_l = quantiles.cols(0, 1) * inv(temp_l).t();
  arma::mat temp_u(2, 2, arma::fill::ones);
  if (tails == "gaussian") {
    temp_u(0, 1) = gsl_cdf_gaussian_Pinv(alphas(p - 2), 1);
    temp_u(1, 1) = gsl_cdf_gaussian_Pinv(alphas(p - 1), 1);
  } else if (tails == "exponential") {
    temp_u(0, 1) = gsl_cdf_exponential_Pinv(alphas(p - 2), 1);
    temp_u(1, 1) = gsl_cdf_exponential_Pinv(alphas(p - 1), 1);
  }
  arma::mat tail_param_u = quantiles.cols(p - 2, p - 1) * inv(temp_u).t();
  arma::vec yp1(N);
  arma::vec ypp(N);
  if (tails == "gaussian") {
    for (size_t i = 0; i < N; i++) {
      yp1(i) = gsl_ran_gaussian_pdf(quantiles(i, 1) - tail_param_l(i, 0), tail_param_l(i, 1));
      ypp(i) = gsl_ran_gaussian_pdf(quantiles(i, p - 2) - tail_param_u(i, 0), tail_param_u(i, 1));
    }
  } else if (tails == "exponential") {
    for (size_t i = 0; i < N; i++) {
      yp1(i) = gsl_ran_exponential_pdf(quantiles(i, 1) - tail_param_l(i, 0), tail_param_l(i, 1));
      ypp(i) = gsl_ran_exponential_pdf(quantiles(i, p - 2) - tail_param_u(i, 0), tail_param_u(i, 1));
    }
  }

  arma::mat A(4*(p-3), 4*(p-3), arma::fill::zeros);

  // interpolation restrictions
  int i;
  for (i = 1; i <= p-3; i++) {
    A(i-1, (4*(i-1)))   = 1;
    A(i-1, (4*(i-1)+1)) = alphas[i];
    A(i-1, (4*(i-1)+2)) = alphas[i] * alphas[i];
    A(i-1, (4*(i-1)+3)) = alphas[i] * alphas[i] * alphas[i];
  }

  // add last one onto the end
  A(p-3, 4*(p-3)-4) = 1;
  A(p-3, 4*(p-3)-3) = alphas[p-2];
  A(p-3, 4*(p-3)-2) = alphas[p-2] * alphas[p-2];
  A(p-3, 4*(p-3)-1) = alphas[p-2] * alphas[p-2] * alphas[p-2];

  // known derivatives
  A(p-2, 0) = 0;
  A(p-2, 1) = 1;
  A(p-2, 2) = 2 * alphas[1];
  A(p-2, 3) = 3 * alphas[1] * alphas[1];
  A(p-1, 4*(p-3)-4) = 0;
  A(p-1, 4*(p-3)-3) = 1;
  A(p-1, 4*(p-3)-2) = 2 * alphas[p-2];
  A(p-1, 4*(p-3)-1) = 3 * alphas[p-2] * alphas[p-2];

  // continuity, first, and second derivative restrictions
  for (i = 2; i <= p-3; i++) {
    A(p-2+i,    4*(i-2))   =  1;
    A(p-2+i,    4*(i-2)+1) =  alphas[i];
    A(p-2+i,    4*(i-2)+2) =  alphas[i] * alphas[i];
    A(p-2+i,    4*(i-2)+3) =  alphas[i] * alphas[i] * alphas[i];
    A(p-2+i,    4*(i-2)+4) = -1;
    A(p-2+i,    4*(i-2)+5) = -alphas[i];
    A(p-2+i,    4*(i-2)+6) = -alphas[i] * alphas[i];
    A(p-2+i,    4*(i-2)+7) = -alphas[i] * alphas[i] * alphas[i];
    A(2*p-6+i,  4*(i-2))   =  0;
    A(2*p-6+i,  4*(i-2)+1) =  1;
    A(2*p-6+i,  4*(i-2)+2) =  2 * alphas[i];
    A(2*p-6+i,  4*(i-2)+3) =  3 * alphas[i] * alphas[i];
    A(2*p-6+i,  4*(i-2)+4) =  0;
    A(2*p-6+i,  4*(i-2)+5) = -1;
    A(2*p-6+i,  4*(i-2)+6) = -2 * alphas[i];
    A(2*p-6+i,  4*(i-2)+7) = -3 * alphas[i] * alphas[i];
    A(3*p-10+i, 4*(i-2))   =  0;
    A(3*p-10+i, 4*(i-2)+1) =  0;
    A(3*p-10+i, 4*(i-2)+2) =  2;
    A(3*p-10+i, 4*(i-2)+3) =  6 * alphas[i];
    A(3*p-10+i, 4*(i-2)+4) =  0;
    A(3*p-10+i, 4*(i-2)+5) =  0;
    A(3*p-10+i, 4*(i-2)+6) = -2;
    A(3*p-10+i, 4*(i-2)+7) = -6 * alphas[i];
  }
  // A.print();


  // quantile matrix (all 0 other than these)
  arma::mat Q(4*(p-3), N, arma::fill::zeros);
  for (i = 0; i <= p-3; i++) {
    Q(i, arma::span::all) = quantiles(arma::span::all, i+1).t();
  }
  Q(p-2, arma::span::all) = 1 / yp1.t();
  Q(p-1, arma::span::all) = 1 / ypp.t();
  // Q.print();
  // inv(A).print();

  arma::mat theta = inv(A) * Q;

  // theta.print();
  // arma::cout << arma::size(alphas) << arma::endl;
  // arma::cout << arma::size(quantiles) << arma::endl;
  alphas.shed_row(0);
  alphas.shed_row(alphas.n_rows-1);
  quantiles.shed_col(0);
  quantiles.shed_col(quantiles.n_cols-1);
  // arma::cout << arma::size(alphas) << arma::endl;
  // arma::cout << arma::size(quantiles) << arma::endl;
  p -= 2;

  arma::vec y_hat(N);

  // Loop through the rows (no need to use vectors in C!), row designated by i. 0-indexed!
  for (size_t i = 0; i < N; i++) {
    int j;
    if (distn(i) == 3) {
      j = find_arma(y(i), alphas, p);
      // printf("%f %d\n", y(i), j);
      if (j == -1) {
        // Probability was lower than all the alphas, so use the tail distribution
        if (tails == "gaussian") {
          y_hat(i) = gsl_cdf_gaussian_Pinv(y(i), tail_param_l(i, 1)) + tail_param_l(i, 0);
          // printf("low: %f %f\n", y(i), y_hat(i));
        } else if (tails == "exponential") {
          y_hat(i) = gsl_cdf_exponential_Pinv(y(i), tail_param_l(i, 1)) + tail_param_l(i, 0);
        }
      } else if (j == (p-1)) {
        // Probability was higher than all the alphas, so use the tail distribution
        if (tails == "gaussian") {
          y_hat(i) = gsl_cdf_gaussian_Pinv(y(i), tail_param_u(i, 1)) + tail_param_u(i, 0);
          // printf("hi: %f %f\n", y(i), y_hat(i));
        } else if (tails == "exponential") {
          y_hat(i) = gsl_cdf_exponential_Pinv(y(i), tail_param_u(i, 1)) + tail_param_u(i, 0);
        }
      } else {
        if (non_neg(theta(4*j+3, i), theta(4*j+2, i), theta(4*j+1, i), alphas(j), alphas(j+1))) {
          y_hat(i) = (y(i) - alphas(j)) * (quantiles(i, j+1) - quantiles(i, j)) / (alphas(j+1) - alphas(j)) + quantiles(i, j);
        } else {
          y_hat(i) = theta(4*j, i) + y(i) * theta(4*j+1, i) + y(i) * y(i) * theta(4*j+2, i) + y(i) * y(i) * y(i) * theta(4*j+3, i);
        }
      }
      y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
    } else {
      y(i) = (y(i) - q_shift(i)) / q_stretch(i);
      j = find_arma(y(i), quantiles.row(i).t(), p);
      if (j == -1) {
        // Value was lower than all the quantiles, so use the tails distribution to approximate.
        if (tails == "gaussian") {
          double x = (y(i) - tail_param_l(i, 0))/tail_param_l(i, 1);
          if (distn(i) == 1) { // CDF
            if(fabs(x) <= 35.0){
              y_hat(i) = gsl_cdf_gaussian_P(x*tail_param_l(i, 1), tail_param_l(i, 1));
            }
            else{
              y_hat(i) = R::pnorm(x*tail_param_l(i, 1), 0.0, tail_param_l(i,1), 1, 0);
            }
          } else { // PDF
            if(fabs(x) <= 35.0){
              y_hat(i) = gsl_sf_erf_Z(x)/tail_param_l(i, 1);
            }
            else{
              y_hat(i) = R::dnorm(x, 0.0, 1.0, 0)/tail_param_l(i, 1);
            }
            y_hat(i) /= q_stretch(i);
          }
        } else if (tails == "exponential") {
          if (distn(i) == 1) { // CDF
            y_hat(i) = gsl_cdf_exponential_P(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
          } else { // PDF
            y_hat(i) = gsl_ran_exponential_pdf(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
            y_hat(i) /= q_stretch(i);
          }
        }
      } else if (j == (p-1)) {
        // Value was higher than all the quantiles, so use the tails distribution to approximate.
        if (tails == "gaussian") {
          double x = (y(i) - tail_param_u(i, 0))/tail_param_u(i, 1);
          if (distn(i) == 1) {
            if(fabs(x) <= 35.0){
              y_hat(i) = gsl_cdf_gaussian_P(x * tail_param_u(i, 1), tail_param_u(i, 1));
            } else{
              y_hat(i) = R::pnorm(x*tail_param_u(i, 1), 0.0, tail_param_u(i, 1), 1, 0);
            }
          } else {
            //y_hat(i) = gsl_sf_erf_Z((y(i) - tail_param_u(i, 0)) / tail_param_u(i, 1)) / tail_param_u(i, 1);
            if(fabs(x) <= 35.0){
              //y_hat[i] = gsl_sf_erf_Z(x) / tail_param_u(i, 1);
              y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_u(i, 1) ;
            } else{
              y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_u(i, 1);
            }
            y_hat(i) /= q_stretch(i);
          }
        } else if (tails == "exponential") {
          if (distn(i) == 1) {
            y_hat(i) = gsl_cdf_exponential_P(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
          } else {
            y_hat(i) = gsl_ran_exponential_pdf(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
            y_hat(i) /= q_stretch(i);
          }
        }
      } else {
        // Find the roots of the cubic function stored in the roots array.
        double * roots = (double *)malloc(3 * sizeof(double));
        int num_roots = gsl_poly_solve_cubic(theta(4*j+2, i)/theta(4*j+3, i), theta(4*j+1, i)/theta(4*j+3, i), (theta(4*j, i)-y(i))/theta(4*j+3, i), &roots[0], &roots[1], &roots[2]);
        double y_hat_temp;

        // Cubic functions either have 1 root or 3 roots, test for both cases.
        if (num_roots == 1) {
          // If there is 1 root, it must lie in the interval, so return this root.
          y_hat_temp = roots[0];
          free(roots);
        } else {
          // If there are 3 roots, check which one lies in the interval.
          int boundaries = 1;
          for (size_t k = 0; k < 3; k++) {
            // printf("root%lu: %f\n", k, roots[k]);
            if (quantiles(i, j) <= roots[k] && roots[k] <= quantiles(i, j+1)) {
              y_hat_temp = roots[k];
              boundaries = 0;
              free(roots);
              break;
            }
          }
          // Makes sure edge cases near quantiles don't fall outside the interval due to slight approximations (just take the closest root to the interval).
          if (boundaries) {
            size_t min_index;
            double min_dist = DBL_MAX;
            double dist;
            for (size_t k = 1; k < 3; k++) {
              dist = MIN(fabs(roots[k] - quantiles(i, j)), fabs(roots[k] - quantiles(i, j+1)));
              if (dist < min_dist) {
                min_dist = dist;
                min_index = k;
              }
            }
            y_hat_temp = roots[min_index];
            free(roots);
          }
        }
        if (distn(i) == 1) {
          y_hat(i) = y_hat_temp;
        } else {
          y_hat(i) = 1 / (theta(4*j+1, i) + 2 * y_hat_temp * theta(4*j+2, i) + 3 * y_hat_temp * y_hat_temp * theta(4*j+3, i));
          y_hat(i) /= q_stretch(i);
        }
      }
    }
  }
  return y_hat;
}


// [[Rcpp::export]]
arma::mat c_spline_comb(arma::vec y, arma::mat quantiles, arma::vec alphas, std::string tails, arma::vec distn) {
  // Perform the tri-diagonalization two-time pass through to compute the second derivatives of the cubic spline
  int N = quantiles.n_rows;
  int p = quantiles.n_cols;

  arma::vec q_shift = quantiles.col(0);
  arma::vec q_stretch = quantiles.col(p-1) - quantiles.col(0);
  for (size_t j = 0; j < p; j++) {
    quantiles.col(j) = (quantiles.col(j) - q_shift) / q_stretch;
  }

  arma::mat temp_l(2, 2, arma::fill::ones);
  if (tails == "gaussian") {
    temp_l(0, 1) = gsl_cdf_gaussian_Pinv(alphas(0), 1);
    temp_l(1, 1) = gsl_cdf_gaussian_Pinv(alphas(1), 1);
  } else if (tails == "exponential") {
    temp_l(0, 1) = gsl_cdf_exponential_Pinv(alphas(0), 1);
    temp_l(1, 1) = gsl_cdf_exponential_Pinv(alphas(1), 1);
  }
  arma::mat tail_param_l = quantiles.cols(0, 1) * inv(temp_l).t();
  arma::mat temp_u(2, 2, arma::fill::ones);
  if (tails == "gaussian") {
    temp_u(0, 1) = gsl_cdf_gaussian_Pinv(alphas(p - 2), 1);
    temp_u(1, 1) = gsl_cdf_gaussian_Pinv(alphas(p - 1), 1);
  } else if (tails == "exponential") {
    temp_u(0, 1) = gsl_cdf_exponential_Pinv(alphas(p - 2), 1);
    temp_u(1, 1) = gsl_cdf_exponential_Pinv(alphas(p - 1), 1);
  }
  arma::mat tail_param_u = quantiles.cols(p - 2, p - 1) * inv(temp_u).t();
  arma::vec yp1(N);
  arma::vec ypp(N);
  if (tails == "gaussian") {
    for (size_t i = 0; i < N; i++) {
      yp1(i) = gsl_ran_gaussian_pdf(quantiles(i, 1) - tail_param_l(i, 0), tail_param_l(i, 1));
      ypp(i) = gsl_ran_gaussian_pdf(quantiles(i, p - 2) - tail_param_u(i, 0), tail_param_u(i, 1));
    }
  } else if (tails == "exponential") {
    for (size_t i = 0; i < N; i++) {
      yp1(i) = gsl_ran_exponential_pdf(quantiles(i, 1) - tail_param_l(i, 0), tail_param_l(i, 1));
      ypp(i) = gsl_ran_exponential_pdf(quantiles(i, p - 2) - tail_param_u(i, 0), tail_param_u(i, 1));
    }
  }

  alphas.shed_row(0);
  alphas.shed_row(alphas.n_rows-1);
  quantiles.shed_col(0);
  quantiles.shed_col(quantiles.n_cols-1);
  p -= 2;

  arma::mat y2(N, p);

  // Loop through the rows (no need to use vectors in C!), row designated by i. 0-indexed!
  for (size_t i = 0; i < N; i++) {
    int j;
    double d, qp, sig, up;
    // double * u = (double *)calloc(sizeof(double), p-1);
    arma::vec u = arma::vec(p-1);
    arma::rowvec quantiles_i = quantiles.row(i);
    arma::rowvec y2_i = y2.row(i);

    // Use a specific first derivative at the first quantile
    y2(i, 0) = -0.5;
    u(0) = (3.0 / (quantiles_i(1) - quantiles_i(0))) * ((alphas(1)-alphas(0)) / (quantiles_i(1)-quantiles_i(0)) - yp1(i));

    for (j = 1; j <= p-2; j++) {
      sig = (quantiles_i(j) - quantiles_i(j-1)) / (quantiles_i(j+1) - quantiles_i(j-1));
      d = sig * y2(i, j-1) + 2.0;
      y2(i, j) = (sig - 1.0) / d;
      u(j) = (alphas(j+1)-alphas(j)) / (quantiles_i(j+1) - quantiles_i(j)) - (alphas(j)-alphas(j-1)) / (quantiles_i(j) - quantiles_i(j-1));
      u(j) = (6.0 * u(j) / (quantiles_i(j+1) - quantiles_i(j-1)) - sig * u(j-1)) / d;
      // printf("%f ", u(j));
    }
    // printf("\n");

    // Use a specific first derivative at the last quantile
    qp = 0.5;
    up = (3.0 / (quantiles_i(p-1) - quantiles_i(p-2))) * (ypp(i) - ((alphas(p-1)-alphas(p-2)) / (quantiles_i(p-1)-quantiles_i(p-2))));
    // printf("%f %f\n", qp, up);

    y2(i, p-1) = (up - (qp * u(p-2))) / ((qp * y2(i, p-2)) + 1.0);
    // printf("%f\n", y2(i, p-1));

    for (j = p-2; j >= 0; j--) {
      y2(i, j) = (y2(i, j) * y2(i, j+1)) + u(j);
      // printf("%f ", y2(i));
    }
    // printf("\n");

    // free(u);
    // memcpy(y2[j], y2_j, sizeof(double)*p);
  }

  // alphas.shed_row(0);
  // alphas.shed_row(alphas.n_rows-1);
  // quantiles.shed_col(0);
  // quantiles.shed_col(quantiles.n_cols-1);
  // arma::cout << arma::size(alphas) << arma::endl;
  // arma::cout << arma::size(quantiles) << arma::endl;
  // p -= 2;

  arma::vec y_hat(N);

  // Loop through the rows (no need to use vectors in C!), row designated by i. 0-indexed!
  for (size_t i = 0; i < N; i++) {
    // printf("i: %lu y: %f\n", i, y(i));
    // Search for the correct interval on which to approximate the distribution.
    int j;
    if (distn(i) == 3) {
      // Based on alphas, not quantiles, since we are passed a probability.
      j = find_arma(y(i), alphas, p);
    } else {
      y(i) = (y(i) - q_shift(i)) / q_stretch(i);
      j = find_arma(y(i), quantiles.row(i).t(), p);
    }
    // printf("j: %d\n", j);

    if (j == -1) {
      // Value was lower than all the quantiles, so use the tails distribution to approximate.
      if (tails == "gaussian") {
        double x = (y(i) - tail_param_l(i, 0))/tail_param_l(i, 1);
        if (distn(i) == 1) {
          //y_hat(i) = gsl_cdf_gaussian_P(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
          if(fabs(x) <= 35.0){
            y_hat(i) = gsl_cdf_gaussian_P(x * tail_param_l(i, 1), tail_param_l(i, 1));
          } else{
            y_hat(i) = R::pnorm(x*tail_param_l(i, 1), 0.0, tail_param_l(i, 1), 1, 0);
          }
        } else if (distn(i) == 2) {
          //y_hat(i) = gsl_ran_gaussian_pdf(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
          if(fabs(x) <= 25.0){
            y_hat(i) = gsl_sf_erf_Z(x)/tail_param_l(i, 1);
          }
          else{
            y_hat(i) = R::dnorm(x, 0.0, 1.0, 0)/tail_param_l(i, 1);
          }
        } else if (distn(i) == 3) {
          y_hat(i) = gsl_cdf_gaussian_Pinv(y(i), tail_param_l(i, 1)) + tail_param_l(i, 0);
        }
      } else if (tails == "exponential") {
        if (distn(i) == 1) {
          y_hat(i) = gsl_cdf_exponential_P(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
        } else if (distn(i) == 2) {
          y_hat(i) = gsl_ran_exponential_pdf(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
        } else if (distn(i) == 3) {
          y_hat(i) = gsl_cdf_exponential_Pinv(y(i), tail_param_l(i, 1)) + tail_param_l(i, 0);
        }
      }
      if (distn(i) == 2) {
        y_hat(i) = y_hat(i) / q_stretch(i);
      }
      if (distn(i) == 3) {
        y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
      }
      continue;
    } else if (j == (p-1)) {
      double x = (y(i) - tail_param_u(i, 0))/tail_param_u(i, 1);
      // Value was higher than all the quantiles, so use the tails distribution to approximate.
      if (tails == "gaussian") {
        if (distn(i) == 1) {
          //y_hat(i) = gsl_cdf_gaussian_P(y(i) - tail_param_u(i, 0), tail_param_u(i, 1));
          if(fabs(x) <= 35.0){
            y_hat(i) = gsl_cdf_gaussian_P(x * tail_param_u(i, 1), tail_param_u(i, 1));
          } else{
            y_hat(i) = R::pnorm(x*tail_param_u(i, 1), 0.0, tail_param_u(i, 1), 1, 0);
          }
        } else if (distn(i) == 2) {
          if(fabs(x) <= 35.0){
            //y_hat[i] = gsl_sf_erf_Z(x) / tail_param_u(i, 1);
            y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_u(i, 1) ;
          } else{
            y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_u(i, 1);
          }
        } else if (distn(i) == 3) {
          y_hat(i) = gsl_cdf_gaussian_Pinv(y(i), tail_param_u(i, 1)) + tail_param_u(i, 0);
        }
      } else if (tails == "exponential") {
        if (distn(i) == 1) {
          y_hat(i) = gsl_cdf_exponential_P(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
        } else if (distn(i) == 2) {
          y_hat(i) = gsl_ran_exponential_pdf(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
        } else if (distn(i) == 3) {
          y_hat(i) = gsl_cdf_exponential_Pinv(y(i), tail_param_u(i, 1)) + tail_param_u(i, 0);
        }
      }
      if (distn(i) == 2) {
        y_hat(i) = y_hat(i) / q_stretch(i);
      }
      if (distn(i) == 3) {
        y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
      }
      continue;
    }

    double h = quantiles(i, j+1) - quantiles(i, j);
    double a = (quantiles(i, j+1) - y(i)) / h;
    double b = (y(i) - quantiles(i, j)) / h;
    // printf("h: %f %f %f\n", h, quantiles(i, j+1), quantiles(i, j));
    // printf("a: %f\n", a);
    // printf("b: %f\n", b);

    double a1 = (y2(i, j+1) - y2(i, j)) / (6 * h);
    double b1 = (quantiles(i, j+1) * y2(i, j) - quantiles(i, j) * y2(i, j+1)) / (2 * h);
    double c1 = (alphas(j+1) - alphas(j)) / h + (h * (y2(i, j) - y2(i, j+1))) / 6 + (quantiles(i, j) * quantiles(i, j) * y2(i, j+1) - quantiles(i, j+1) * quantiles(i, j+1) * y2(i, j)) / (2 * h);
    double d1 = (alphas(j) * quantiles(i, j+1) - alphas(j+1) * quantiles(i, j)) / h + (quantiles(i, j+1) * quantiles(i, j+1) * quantiles(i, j+1) * y2(i, j) - quantiles(i, j) * quantiles(i, j) * quantiles(i, j) * y2(i, j+1)) / (6 * h) + (h * (quantiles(i, j) * y2(i, j+1) - quantiles(i, j+1) * y2(i, j))) / 6;
    // printf("a1: %f, b1: %f, c1: %f, d1: %f\n", a1, b1, c1, d1);

    // Calculate the discriminant to make sure the cubic function isn't decreasing.
    // printf("dist: %d\n", distn(i));
    if (non_neg(a1, b1, c1, quantiles(i, j), quantiles(i, j+1))) {
      if (distn(i) == 1) { // CDF
        // printf("cdf non_neg\n");
        y_hat(i) = (y(i) - quantiles(i, j)) * (alphas(j+1) - alphas(j)) / h + alphas(j);
      } else if (distn(i) == 2) { // PDF
        // printf("pdf non_neg\n");
        y_hat(i) = (alphas(j+1) - alphas(j)) / h;
        y_hat(i) /= q_stretch(i);
      } else if (distn(i) == 3) {
        y_hat(i) = (y(i) - alphas(j)) * h / (alphas(j+1) - alphas(j)) + quantiles(i, j);
        y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
      }
      continue;
    }

    // Quantile function:
    if (distn(i) == 1) { // CDF
      // printf("distn\n");
      y_hat(i) = a*alphas(j) + b*alphas(j+1) + ((a*a*a-a)*y2(i, j) + (b*b*b-b)*y2(i, j+1)) * (h*h) / 6;
    } else if (distn(i) == 2) { // PDF (just the derivative of the CDF)
      // printf("pdf\n");
      double b_prime = 1 / h;
      double a_prime = -b_prime;
      // printf("%f %f %f %f %f %f %f %f\n", y2(i, j), a_prime, b_prime, a, b, h, alphas(j), alphas(j + 1));
      // printf("%f %f %f %f\n", a*a, 3*(a*a)*a_prime, 3*(a*a)*a_prime - a_prime, (3*(a*a)*a_prime - a_prime) * y2(i, j));
      double temp = a_prime*alphas(j) + b_prime*alphas(j+1) + ((3*(a*a)*a_prime - a_prime) * y2(i, j) + (3*(b*b)*b_prime - b_prime) * y2(i, j+1)) * (h*h) / 6.0;
      // printf("temp: %f\n", temp);
      y_hat(i) = temp / q_stretch(i);
    } else if (distn(i) == 3) {
      // Probability was in an interval, designated by j. 0-indexed!
      // printf("mid\n");

      // Divide the coefficients through by a to be able to pass to the solver. Subtract d by the given y to make the roots of the function the solutions.
      double b2 = b1 / a1;
      double c2 = c1 / a1;
      double d2 = (d1 - y(i)) / a1;
      // printf("b2: %f, c2: %f, d2: %f\n", b2, c2, d2);

      // Find the roots of the cubic function stored in the roots array.
      double * roots = (double *)malloc(3 * sizeof(double));
      int num_roots = gsl_poly_solve_cubic(b2, c2, d2, &roots[0], &roots[1], &roots[2]);
      // printf("num_roots: %d\n", num_roots);

      // Cubic functions either have 1 root or 3 roots, test for both cases.
      if (num_roots == 1) {
        // If there is 1 root, it must lie in the interval, so return this root.
        y_hat(i) = roots[0];
        free(roots);
      } else {
        // If there are 3 roots, check which one lies in the interval.
        int boundaries = 1;
        for (size_t k = 0; k < 3; k++) {
          // printf("root%lu: %f\n", k, roots[k]);
          if (quantiles(i, j) <= roots[k] && roots[k] <= quantiles(i, j+1)) {
            y_hat(i) = roots[k];
            boundaries = 0;
            free(roots);
            break;
          }
        }
        // Makes sure edge cases near quantiles don't fall outside the interval due to slight approximations (just take the closest root to the interval).
        if (boundaries) {
          size_t min_index;
          double min_dist = DBL_MAX;
          double dist;
          for (size_t k = 1; k < 3; k++) {
            dist = MIN(fabs(roots[k] - quantiles(i, j)), fabs(roots[k] - quantiles(i, j+1)));
            if (dist < min_dist) {
              min_dist = dist;
              min_index = k;
            }
          }
          y_hat(i) = roots[min_index];
          free(roots);
        }
      }
      // printf("y_hat: %f\n", y_hat(i));
      y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
    }
  }
  return y_hat;
}

// [[Rcpp::export]]
arma::mat c_spline_comb_track(arma::vec y, arma::mat quantiles, arma::vec alphas, std::string tails, arma::vec distn) {
  // Perform the tri-diagonalization two-time pass through to compute the second derivatives of the cubic spline

  int N = quantiles.n_rows;
  int p = quantiles.n_cols;

  arma::vec q_shift = quantiles.col(0);
  arma::vec q_stretch = quantiles.col(p-1) - quantiles.col(0);
  for (size_t j = 0; j < p; j++) {
    quantiles.col(j) = (quantiles.col(j) - q_shift) / q_stretch;
  }

  arma::mat temp_l(2, 2, arma::fill::ones);
  if (tails == "gaussian") {
    temp_l(0, 1) = gsl_cdf_gaussian_Pinv(alphas(0), 1);
    temp_l(1, 1) = gsl_cdf_gaussian_Pinv(alphas(1), 1);
  } else if (tails == "exponential") {
    temp_l(0, 1) = gsl_cdf_exponential_Pinv(alphas(0), 1);
    temp_l(1, 1) = gsl_cdf_exponential_Pinv(alphas(1), 1);
  }
  arma::mat tail_param_l = quantiles.cols(0, 1) * inv(temp_l).t();
  arma::mat temp_u(2, 2, arma::fill::ones);
  if (tails == "gaussian") {
    temp_u(0, 1) = gsl_cdf_gaussian_Pinv(alphas(p - 2), 1);
    temp_u(1, 1) = gsl_cdf_gaussian_Pinv(alphas(p - 1), 1);
  } else if (tails == "exponential") {
    temp_u(0, 1) = gsl_cdf_exponential_Pinv(alphas(p - 2), 1);
    temp_u(1, 1) = gsl_cdf_exponential_Pinv(alphas(p - 1), 1);
  }
  arma::mat tail_param_u = quantiles.cols(p - 2, p - 1) * inv(temp_u).t();
  arma::vec yp1(N);
  arma::vec ypp(N);
  if (tails == "gaussian") {
    for (size_t i = 0; i < N; i++) {
      yp1(i) = gsl_ran_gaussian_pdf(quantiles(i, 1) - tail_param_l(i, 0), tail_param_l(i, 1));
      ypp(i) = gsl_ran_gaussian_pdf(quantiles(i, p - 2) - tail_param_u(i, 0), tail_param_u(i, 1));
    }
  } else if (tails == "exponential") {
    for (size_t i = 0; i < N; i++) {
      yp1(i) = gsl_ran_exponential_pdf(quantiles(i, 1) - tail_param_l(i, 0), tail_param_l(i, 1));
      ypp(i) = gsl_ran_exponential_pdf(quantiles(i, p - 2) - tail_param_u(i, 0), tail_param_u(i, 1));
    }
  }
  alphas.shed_row(0);
  alphas.shed_row(alphas.n_rows-1);
  quantiles.shed_col(0);
  quantiles.shed_col(quantiles.n_cols-1);
  p -= 2;

  arma::mat y2(N, p);

  // Loop through the rows (no need to use vectors in C!), row designated by i. 0-indexed!
  for (size_t i = 0; i < N; i++) {
    int j;
    double d, qp, sig, up;
    // double * u = (double *)calloc(sizeof(double), p-1);
    arma::vec u = arma::vec(p-1);
    arma::rowvec quantiles_i = quantiles.row(i);
    arma::rowvec y2_i = y2.row(i);

    // Use a specific first derivative at the first quantile
    y2(i, 0) = -0.5;
    u(0) = (3.0 / (quantiles_i(1) - quantiles_i(0))) * ((alphas(1)-alphas(0)) / (quantiles_i(1)-quantiles_i(0)) - yp1(i));

    for (j = 1; j <= p-2; j++) {
      sig = (quantiles_i(j) - quantiles_i(j-1)) / (quantiles_i(j+1) - quantiles_i(j-1));
      d = sig * y2(i, j-1) + 2.0;
      y2(i, j) = (sig - 1.0) / d;
      u(j) = (alphas(j+1)-alphas(j)) / (quantiles_i(j+1) - quantiles_i(j)) - (alphas(j)-alphas(j-1)) / (quantiles_i(j) - quantiles_i(j-1));
      u(j) = (6.0 * u(j) / (quantiles_i(j+1) - quantiles_i(j-1)) - sig * u(j-1)) / d;
      // printf("%f ", u(j));
    }
    // printf("\n");

    // Use a specific first derivative at the last quantile
    qp = 0.5;
    up = (3.0 / (quantiles_i(p-1) - quantiles_i(p-2))) * (ypp(i) - ((alphas(p-1)-alphas(p-2)) / (quantiles_i(p-1)-quantiles_i(p-2))));
    // printf("%f %f\n", qp, up);

    y2(i, p-1) = (up - (qp * u(p-2))) / ((qp * y2(i, p-2)) + 1.0);
    // printf("%f\n", y2(i, p-1));

    for (j = p-2; j >= 0; j--) {
      y2(i, j) = (y2(i, j) * y2(i, j+1)) + u(j);
      // printf("%f ", y2(i));
    }
    // printf("\n");

    // free(u);
    // memcpy(y2[j], y2_j, sizeof(double)*p);
  }

  // alphas.shed_row(0);
  // alphas.shed_row(alphas.n_rows-1);
  // quantiles.shed_col(0);
  // quantiles.shed_col(quantiles.n_cols-1);
  // arma::cout << arma::size(alphas) << arma::endl;
  // arma::cout << arma::size(quantiles) << arma::endl;
  // p -= 2;

  arma::vec y_hat(N);

  // Loop through the rows (no need to use vectors in C!), row designated by i. 0-indexed!
  for (size_t i = 0; i < N; i++) {
    // printf("i: %lu y: %f\n", i, y(i));
    // Search for the correct interval on which to approximate the distribution.
    int j;
    if (distn(i) == 3) {
      // Based on alphas, not quantiles, since we are passed a probability.
      j = find_arma(y(i), alphas, p);
    } else {
      y(i) = (y(i) - q_shift(i)) / q_stretch(i);
      j = find_arma(y(i), quantiles.row(i).t(), p);
    }
    // printf("j: %d\n", j);

    if (j == -1) {
      // Value was lower than all the quantiles, so use the tails distribution to approximate.
      if (tails == "gaussian") {
        double x = (y(i) - tail_param_l(i, 0))/tail_param_l(i, 1);
        if (distn(i) == 1) {
          //y_hat(i) = gsl_cdf_gaussian_P(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
          if(fabs(x) <= 35.0){
            y_hat(i) = gsl_cdf_gaussian_P(x * tail_param_l(i, 1), tail_param_l(i, 1));
          } else{
            y_hat(i) = R::pnorm(x*tail_param_l(i, 1), 0.0, tail_param_l(i, 1), 1, 0);
          }
        } else if (distn(i) == 2) {
          //y_hat(i) = gsl_ran_gaussian_pdf(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
          if(fabs(x) <= 25.0){
            y_hat(i) = gsl_sf_erf_Z(x)/tail_param_l(i, 1);
          }
          else{
            y_hat(i) = R::dnorm(x, 0.0, 1.0, 0)/tail_param_l(i, 1);
          }
        } else if (distn(i) == 3) {
          y_hat(i) = gsl_cdf_gaussian_Pinv(y(i), tail_param_l(i, 1)) + tail_param_l(i, 0);
        }
      } else if (tails == "exponential") {
        if (distn(i) == 1) {
          y_hat(i) = gsl_cdf_exponential_P(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
        } else if (distn(i) == 2) {
          y_hat(i) = gsl_ran_exponential_pdf(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
        } else if (distn(i) == 3) {
          y_hat(i) = gsl_cdf_exponential_Pinv(y(i), tail_param_l(i, 1)) + tail_param_l(i, 0);
        }
      }
      if (distn(i) == 2) {
        y_hat(i) = y_hat(i) / q_stretch(i);
      }
      if (distn(i) == 3) {
        y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
      }
      continue;
    } else if (j == (p-1)) {
      double x = (y(i) - tail_param_u(i, 0))/tail_param_u(i, 1);
      // Value was higher than all the quantiles, so use the tails distribution to approximate.
      if (tails == "gaussian") {
        if (distn(i) == 1) {
          //y_hat(i) = gsl_cdf_gaussian_P(y(i) - tail_param_u(i, 0), tail_param_u(i, 1));
          if(fabs(x) <= 35.0){
            y_hat(i) = gsl_cdf_gaussian_P(x * tail_param_u(i, 1), tail_param_u(i, 1));
          } else{
            y_hat(i) = R::pnorm(x*tail_param_u(i, 1), 0.0, tail_param_u(i, 1), 1, 0);
          }
        } else if (distn(i) == 2) {
          if(fabs(x) <= 35.0){
            //y_hat[i] = gsl_sf_erf_Z(x) / tail_param_u(i, 1);
            y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_u(i, 1) ;
          } else{
            y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_u(i, 1);
          }
        } else if (distn(i) == 3) {
          y_hat(i) = gsl_cdf_gaussian_Pinv(y(i), tail_param_u(i, 1)) + tail_param_u(i, 0);
        }
      } else if (tails == "exponential") {
        if (distn(i) == 1) {
          y_hat(i) = gsl_cdf_exponential_P(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
        } else if (distn(i) == 2) {
          y_hat(i) = gsl_ran_exponential_pdf(y(i) - tail_param_l(i, 0), tail_param_l(i, 1));
        } else if (distn(i) == 3) {
          y_hat(i) = gsl_cdf_exponential_Pinv(y(i), tail_param_u(i, 1)) + tail_param_u(i, 0);
        }
      }
      if (distn(i) == 2) {
        y_hat(i) = y_hat(i) / q_stretch(i);
      }
      if (distn(i) == 3) {
        y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
      }
      continue;
    }


    double h = quantiles(i, j+1) - quantiles(i, j);
    double a = (quantiles(i, j+1) - y(i)) / h;
    double b = (y(i) - quantiles(i, j)) / h;
    // printf("h: %f %f %f\n", h, quantiles(i, j+1), quantiles(i, j));
    // printf("a: %f\n", a);
    // printf("b: %f\n", b);

    double a1 = (y2(i, j+1) - y2(i, j)) / (6 * h);
    double b1 = (quantiles(i, j+1) * y2(i, j) - quantiles(i, j) * y2(i, j+1)) / (2 * h);
    double c1 = (alphas(j+1) - alphas(j)) / h + (h * (y2(i, j) - y2(i, j+1))) / 6 + (quantiles(i, j) * quantiles(i, j) * y2(i, j+1) - quantiles(i, j+1) * quantiles(i, j+1) * y2(i, j)) / (2 * h);
    double d1 = (alphas(j) * quantiles(i, j+1) - alphas(j+1) * quantiles(i, j)) / h + (quantiles(i, j+1) * quantiles(i, j+1) * quantiles(i, j+1) * y2(i, j) - quantiles(i, j) * quantiles(i, j) * quantiles(i, j) * y2(i, j+1)) / (6 * h) + (h * (quantiles(i, j) * y2(i, j+1) - quantiles(i, j+1) * y2(i, j))) / 6;
    // printf("a1: %f, b1: %f, c1: %f, d1: %f\n", a1, b1, c1, d1);

    // Calculate the discriminant to make sure the cubic function isn't decreasing.
    // printf("dist: %d\n", distn(i));
    if (non_neg(a1, b1, c1, quantiles(i, j), quantiles(i, j+1))) {
      if (distn(i) == 1) { // CDF
        // printf("cdf non_neg\n");
        y_hat(i) = (y(i) - quantiles(i, j)) * (alphas(j+1) - alphas(j)) / h + alphas(j);
      } else if (distn(i) == 2) { // PDF
        // printf("pdf non_neg\n");
        y_hat(i) = (alphas(j+1) - alphas(j)) / h;
        y_hat(i) /= q_stretch(i);
      } else if (distn(i) == 3) {
        y_hat(i) = (y(i) - alphas(j)) * h / (alphas(j+1) - alphas(j)) + quantiles(i, j);
        y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
      }
      continue;
    }

    // Quantile function:
    if (distn(i) == 1) { // CDF
      // printf("distn\n");
      y_hat(i) = a*alphas(j) + b*alphas(j+1) + ((a*a*a-a)*y2(i, j) + (b*b*b-b)*y2(i, j+1)) * (h*h) / 6;
    } else if (distn(i) == 2) { // PDF (just the derivative of the CDF)
      // printf("pdf\n");
      double b_prime = 1 / h;
      double a_prime = -b_prime;
      // printf("%f %f %f %f %f %f %f %f\n", y2(i, j), a_prime, b_prime, a, b, h, alphas(j), alphas(j + 1));
      // printf("%f %f %f %f\n", a*a, 3*(a*a)*a_prime, 3*(a*a)*a_prime - a_prime, (3*(a*a)*a_prime - a_prime) * y2(i, j));
      double temp = a_prime*alphas(j) + b_prime*alphas(j+1) + ((3*(a*a)*a_prime - a_prime) * y2(i, j) + (3*(b*b)*b_prime - b_prime) * y2(i, j+1)) * (h*h) / 6.0;
      // printf("temp: %f\n", temp);
      y_hat(i) = temp / q_stretch(i);
    } else if (distn(i) == 3) {
      // Probability was in an interval, designated by j. 0-indexed!
      // printf("mid\n");

      // Divide the coefficients through by a to be able to pass to the solver. Subtract d by the given y to make the roots of the function the solutions.
      double b2 = b1 / a1;
      double c2 = c1 / a1;
      double d2 = (d1 - y(i)) / a1;
      // printf("b2: %f, c2: %f, d2: %f\n", b2, c2, d2);

      // Find the roots of the cubic function stored in the roots array.
      double * roots = (double *)malloc(3 * sizeof(double));
      int num_roots = gsl_poly_solve_cubic(b2, c2, d2, &roots[0], &roots[1], &roots[2]);
      // printf("num_roots: %d\n", num_roots);

      // Cubic functions either have 1 root or 3 roots, test for both cases.
      if (num_roots == 1) {
        // If there is 1 root, it must lie in the interval, so return this root.
        y_hat(i) = roots[0];
        free(roots);
      } else {
        // If there are 3 roots, check which one lies in the interval.
        int boundaries = 1;
        for (size_t k = 0; k < 3; k++) {
          // printf("root%lu: %f\n", k, roots[k]);
          if (quantiles(i, j) <= roots[k] && roots[k] <= quantiles(i, j+1)) {
            y_hat(i) = roots[k];
            boundaries = 0;
            free(roots);
            break;
          }
        }
        // Makes sure edge cases near quantiles don't fall outside the interval due to slight approximations (just take the closest root to the interval).
        if (boundaries) {
          size_t min_index;
          double min_dist = DBL_MAX;
          double dist;
          for (size_t k = 1; k < 3; k++) {
            dist = MIN(fabs(roots[k] - quantiles(i, j)), fabs(roots[k] - quantiles(i, j+1)));
            if (dist < min_dist) {
              min_dist = dist;
              min_index = k;
            }
          }
          y_hat(i) = roots[min_index];
          free(roots);
        }
      }
      // printf("y_hat: %f\n", y_hat(i));
      y_hat(i) = y_hat(i) * q_stretch(i) + q_shift(i);
    }
  }
  return y_hat;
}



int find2(double data, Rcpp::NumericVector mat, int max) {
  // Currently an O(N) search algorithm, but since N will usually be less than 10, it might be more trouble than it's worth to implement something asymptotically faster.
  if (data <= mat[1]) {
    return -1;
  }
  if (mat[max] <= data) {
    return max-1;
  }
  for (size_t i = 0; i < max-1; i++) {
    if (mat[i + 1] <= data && data <= mat[i + 2]) {
      return i;
    }
  }
  return -2;
}


// [[Rcpp::export]]
Rcpp::NumericVector c_splincomb(Rcpp::NumericVector y, Rcpp::NumericMatrix quantiles_full, Rcpp::NumericVector alphas_full, std::string tails, Rcpp::IntegerVector distn) {
  int N = quantiles_full.nrow();
  int p = quantiles_full.ncol() - 2;
  //int p_full = quantiles_full.ncol();
  Rcpp::NumericVector alphas(p);
  alphas = alphas_full[Rcpp::Range(1, p)];
  Rcpp::NumericVector y_hat(N);
  Rcpp::NumericMatrix y2_i(1, p);
  Rcpp::NumericMatrix vcov_l(2, 2);
  Rcpp::NumericMatrix vcov_u(2,2);
  vcov_l(0, 0) = 1.0;
  vcov_l(1, 0) = 1.0;
  vcov_u(0, 0) = 1.0;
  vcov_u(1, 0) = 1.0;
  if (tails == "gaussian"){
    vcov_l(0, 1) = gsl_cdf_gaussian_Pinv(alphas_full[0], 1.0);
    vcov_l(1, 1) = gsl_cdf_gaussian_Pinv(alphas_full[1], 1.0);
    vcov_u(0, 1) = gsl_cdf_gaussian_Pinv(alphas_full[p], 1.0);
    vcov_u(1, 1) = gsl_cdf_gaussian_Pinv(alphas_full[p + 1], 1.0);
    vcov_l = transpose(invert_two(vcov_l));
    vcov_u = transpose(invert_two(vcov_u));

  }
  else if(tails == "exponential"){
    vcov_l(0, 1) = log(alphas_full[0]);
    vcov_l(1, 1) = log(alphas_full[1]);
    vcov_u(0, 1) = -log(alphas_full[p]);
    vcov_u(1, 1) = -log(alphas_full[p + 1]);
    vcov_l = transpose(invert_two(vcov_l));
    vcov_u = transpose(invert_two(vcov_u));
  }
  double tail_param_l[2];
  double tail_param_u[2];
  double yp1;
  double ypp;
  // Rcpp::NumericMatrix quantiles_sub = quantiles_full(Rcpp::_, Rcpp::Range(1, p));
  //Rcpp::NumericMatrix quantiles_i(1, p);
  //Rcpp::NumericMatrix quantiles_full_i(1, p + 2);

  // Loop through the rows (no need to use vectors in C!), row designated by i. 0-indexed!
  for (size_t i = 0; i < N; i++) {
    int j;
    double d, qp, sig, up;
    double * u = (double *)calloc(sizeof(double), p-1);
    //quantiles_i(0, Rcpp::_) = quantiles_full(Rcpp::Range(i, i), Rcpp::Range(1, p));
    //quantiles_full_i(0, Rcpp::_) = quantiles_full(i, Rcpp::_);
    //quantiles_i = quantiles_full_i(Rcpp::_, Rcpp::Range(1, p));
    Rcpp::NumericMatrix::Row quantiles_full_i = quantiles_full(i, Rcpp::_);
    Rcpp::NumericMatrix quantiles_i = quantiles_full(Rcpp::Range(i, i), Rcpp::Range(1, p));
    // Rcpp::NumericMatrix::Row quantiles_i = quantiles_sub(i, Rcpp::_);
    //Rcpp::SubMatrix<REALSXP> quantiles_i = quantiles_full(Rcpp::Range(i, i), Rcpp::Range(1, p));

    if(tails == "gaussian"){
      //compute tail_param_l
      tail_param_l[0] = quantiles_full_i[0] * vcov_l(0, 0) + quantiles_full_i[1] * vcov_l(1, 0);
      tail_param_l[1] = quantiles_full_i[0] * vcov_l(0, 1) + quantiles_full_i[1] * vcov_l(1, 1);
      // tail_param_l = matrix_mult(quantiles_full_i(Rcpp::_, Rcpp::Range(0, 1)), vcov_l);

      //compute tail_param_u
      tail_param_u[0] = quantiles_full_i[p] * vcov_u(0, 0) + quantiles_full_i[p+1] * vcov_u(1, 0);
      tail_param_u[1] = quantiles_full_i[p] * vcov_u(0, 1) + quantiles_full_i[p+1] * vcov_u(1, 1);
      // tail_param_u = matrix_mult(quantiles_full_i(Rcpp::_, Rcpp::Range(p , p + 1)), vcov_u);
      yp1 = gsl_ran_gaussian_pdf(quantiles_full_i[1]  - tail_param_l[0], tail_param_l[1]);
      //yp1 = R::dnorm(quantiles_full_i[1]  - tail_param_l[0], 0, tail_param_l[1],  0);
      ypp = gsl_ran_gaussian_pdf(quantiles_full_i[p] - tail_param_u[0], tail_param_u[1]);
      //ypp = R::dnorm(quantiles_full_i[p] - tail_param_u[0],0, tail_param_u[1] , 0);
    }
    else if(tails == "exponential"){
      // tail_param_l = matrix_mult(quantiles_full_i(Rcpp::_, Rcpp::Range(0, 1)), vcov_l);
      tail_param_l[0] = quantiles_full_i[0] * vcov_l(0, 0) + quantiles_full_i[1] * vcov_l(1, 0);
      tail_param_l[1] = quantiles_full_i[0] * vcov_l(0, 1) + quantiles_full_i[1] * vcov_l(1, 1);

      //compute tail_param_u
      // tail_param_u = matrix_mult(quantiles_full_i(Rcpp::_, Rcpp::Range(p , p + 1)), vcov_u);
      tail_param_u[0] = quantiles_full_i[p] * vcov_u(0, 0) + quantiles_full_i[p+1] * vcov_u(1, 0);
      tail_param_u[1] = quantiles_full_i[p] * vcov_u(0, 1) + quantiles_full_i[p+1] * vcov_u(1, 1);
      yp1 = pow(exp(1),  (quantiles_full_i[1]   - tail_param_l[0])/tail_param_l[1]) / tail_param_l[1];
      ypp = pow(exp(1), -(quantiles_full_i[p-2] - tail_param_u[0])/tail_param_u[1]) / tail_param_u[1];
    }
    // Use a specific first derivative at the first quantile
    int index;
    index = find(y[i], alphas, p);
    if (index != -1 & index != p - 1){
      y2_i[0] = -0.5;
      u[0] = (3.0 / (quantiles_full(i, 2) - quantiles_full(i, 1))) * ((alphas[1]-alphas[0]) / (quantiles_full(i, 2)-quantiles_full(i, 1)) - yp1);

      for (int j = 1; j <= p-2; j++) {
        sig = (quantiles_full(i, j + 1) - quantiles_full(i, j)) / (quantiles_full(i, j + 2) - quantiles_full(i, j));
        d = sig * y2_i[j-1] + 2.0;
        y2_i[j] = (sig - 1.0) / d;
        u[j] = (alphas[j+1]-alphas[j]) / (quantiles_full(i, j + 2) - quantiles_full(i, j + 1)) - (alphas[j]-alphas[j-1]) / (quantiles_full(i, j + 1) - quantiles_full(i, j));
        u[j] = (6.0 * u[j] / (quantiles_full(i, j+2) - quantiles_full(i, j)) - sig * u[j-1]) / d;
        //printf("%f ", u[i]);
      }
      // printf("\n");

      // Use a specific first derivative at the last quantile
      qp = 0.5;
      up = (3.0 / (quantiles_full(i, p) - quantiles_full(i, p - 1))) * (ypp - ((alphas[p-1]-alphas[p-2]) / (quantiles_full(i, p)-quantiles_full(i, p - 1))));

      y2_i[p-1] = (up - (qp * u[p-2])) / ((qp * y2_i[p-2]) + 1.0);

      for (int j = p-2; j >= 0; j--) {
        y2_i[j] = (y2_i[j] * y2_i[j+1]) + u[j];
        // printf("%f ", y2[i]);
      }}
    // printf("\n");

    free(u);

    // Quantile function:
    if (distn[i] == 3) {
      // Search for the correct interval (based on alphas, not quantiles, since we are passed a probability)
      j = index;
      if (j == -1) {
        // Probability was lower than all the quantiles, so use the tail distribution
        if (tails == "gaussian") {
          y_hat[i] = gsl_cdf_gaussian_Pinv(y[i], tail_param_l[1]) + tail_param_l[0];
          //y_hat[i] = R::qnorm(y[i], 0, tail_param_l[1], 1, 0) + tail_param_l[0];
          // printf("low: %f %f\n", y[i], y_hat[i]);
        } else if (tails == "exponential") {
          y_hat[i] = tail_param_l[0] + tail_param_l[1] * log(y[i]);
        }
      } else if (j == (p-1)) {
        // Probability was higher than all the quantiles, so use the tail distribution
        if (tails == "gaussian") {
          y_hat[i] = gsl_cdf_gaussian_Pinv(y[i], tail_param_u[1]) + tail_param_u[0];
          //y_hat[i] = R::qnorm(y[i], 0, tail_param_u[1], 1, 0) + tail_param_u[0];
          // printf("hi: %f %f\n", y[i], y_hat[i]);
        } else if (tails == "exponential") {
          y_hat[i] = tail_param_u[0] - tail_param_u[1] * log(1 - y[i]);
        }
      } else {
        // Probability was in an interval, designated by j. 0-indexed!
        // printf("mid\n");
        // h is a helper quantity that is used frequently.
        // a, b, c, and d are the coefficients of the cubic function on the interval in question, with a the 3rd degree and d the 0th.
        double h = quantiles_full(i, j + 2) - quantiles_full(i, j + 1);
        double a = (y2_i[j+1] - y2_i[j]) / (6*h);
        double b = (quantiles_full(i, j + 2)*y2_i[j] - quantiles_full(i, j + 1)*y2_i[j+1]) / (2*h);
        double c = (alphas[j+1] - alphas[j]) / h + h * (y2_i[j] - y2_i[j+1]) / 6 + (quantiles_full(i, j + 1) * quantiles_full(i, j + 1) * y2_i[j+1] - quantiles_full(i, j + 2) * quantiles_full(i, j + 2) * y2_i[j]) / (2*h);
        double d = -((alphas[j+1] * quantiles_full(i, j + 1))/h) + (alphas[j] * quantiles_full(i, j + 2))/h - (h*quantiles_full(i, j + 2) * y2_i[j])/6 + (quantiles_full(i, j + 2)*quantiles_full(i, j + 2)*quantiles_full(i, j + 2) * y2_i[j])/(6*h) + (h*quantiles_full(i, j + 1)*y2_i[j+1])/6 - (quantiles_full(i, j + 1)*quantiles_full(i, j + 1)*quantiles_full(i, j + 1) * y2_i[j+1])/(6*h);
        // Divide the coefficients through by a to be able to pass to the solver. Subtract d by the given y to make the roots of the function the solutions.
        double b2 = b / a;
        double c2 = c / a;
        double d2 = (d - y[i]) / a;
        // printf("b2: %f, c2: %f, d2: %f\n", b2, c2, d2);
        double * roots = (double *)malloc(3 * sizeof(double));
        // Find the roots of the cubic function stored in the roots array.
        int num_roots = gsl_poly_solve_cubic(b2, c2, d2, &roots[0], &roots[1], &roots[2]);
        // printf("num_roots: %d\n", num_roots);
        // Cubic functions either have 1 root or 3 roots, test for both cases.
        if (num_roots == 1) {
          // If there is 1 root, it must lie in the interval, so return this root.
          y_hat[i] = roots[0];
          free(roots);
        } else {
          // If there are 3 roots, check which one lies in the interval.
          int boundaries = 1;
          for (size_t k = 0; k < 3; k++) {
            // printf("root%lu: %f\n", k, roots[k]);
            if (quantiles_full(i, j + 1) <= roots[k] && roots[k] <= quantiles_full(i, j + 2)) {
              y_hat[i] = roots[k];
              boundaries = 0;
              free(roots);
              break;
            }
          }
          // Makes sure edge cases near quantiles don't fall outside the interval due to slight approximations (just take the closest root to the interval).
          if (boundaries) {
            size_t min_index;
            double min_dist = DBL_MAX;
            double dist;
            for (size_t k = 1; k < 3; k++) {
              dist = MIN(fabs(roots[k] - quantiles_full(i, j + 1)), fabs(roots[k] - quantiles_full(i, j + 2)));
              if (dist < min_dist) {
                min_dist = dist;
                min_index = k;
              }
            }
            y_hat[i] = roots[min_index];
            free(roots);
          }
        }
      }
      continue;
    }
    // Search for the correct interval on which to approximate the distribution.
    j = find(y[i], quantiles_i, p);
    //printf("d is %d \n", j);
    // // printf("%d %f\n", j, y[i]);
    if (j == -1) {
      double x = (y[i] - tail_param_l[0])/tail_param_l[1];
      // Value was lower than all the quantiles, so use the tails distribution to approximate.
      if (tails == "gaussian") {
        if (distn[i] == 1) { // CDF
          if(fabs(x) <= 35.0){
            y_hat[i] = gsl_cdf_gaussian_P(x * tail_param_l[1], tail_param_l[1]);
          } else{
            y_hat[i] = R::pnorm(x*tail_param_l[1], 0.0, tail_param_l[1], 1, 0);
          }
        } else { // PDF
          // y_hat[i] = gsl_sf_erf_Z((y[i] - tail_param_l[0]) / tail_param_l[1]) / tail_param_l[1];
          // printf("%f", gsl_sf_erf_Z(-32.0) / tail_param_l[1]);
          if(fabs(x) <= 35.0){
            y_hat[i] = gsl_sf_erf_Z(x) / tail_param_l[1];
            y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_l[1];
          } else{
            y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_l[1];
          }
        }
      } else if (tails == "exponential") {
        if (distn[i] == 1) { // CDF
          y_hat[i] = pow(M_E, ((y[i] - tail_param_l[0]) / tail_param_l[1]));
        } else { // PDF
          y_hat[i] = pow(M_E, ((y[i] - tail_param_l[0]) / tail_param_l[1])) / tail_param_l[1];
        }
      }
    } else if (j == (p-1)) {
      double x = (y[i] - tail_param_u[0])/tail_param_u[1];
      // Value was higher than all the quantiles, so use the tails distribution to approximate.
      if (tails == "gaussian") {
        if (distn[i] == 1) {
          // y_hat[i] = gsl_cdf_gaussian_P(y[i] - tail_param_u[0], tail_param_u[1]);
          if(fabs(x) <= 35.0){
            y_hat[i] = gsl_cdf_gaussian_P(x * tail_param_u[1], tail_param_u[1]);
          } else{
            y_hat[i] = R::pnorm(x*tail_param_u[1], 0.0, tail_param_u[1], 1, 0);

          }
        } else {
          // y_hat[i] = gsl_sf_erf_Z((y[i] - tail_param_u[0]) / tail_param_u[1]) / tail_param_u[1];
          if(fabs(x) <= 35.0){
            y_hat[i] = gsl_sf_erf_Z(x) / tail_param_u[1];
            y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_u[1] ;
          } else{
            y_hat[i] = R::dnorm(x, 0.0, 1.0, 0)/tail_param_u[1] ;
          }

        }
      } else if (tails == "exponential") {
        if (distn[i] == 1) {
          y_hat[i] = 1 - pow(M_E, (-(y[i] - tail_param_u[0]) / tail_param_u[1]));
        } else {
          y_hat[i] = pow(M_E, (-(y[i] - tail_param_u[0]) / tail_param_u[1])) / tail_param_u[1];
        }
      }
    }
    else {
      // Value was in an interval, designated by j. 0-indexed!
      double h = quantiles_full(i, j+2) - quantiles_full(i, j + 1);
      double a = (quantiles_full(i, j + 2) - y[i]) / h;
      double b = (y[i] - quantiles_full(i, j + 1)) / h;
      // printf("h: %f %f %f\n", h, quantiles_i[j+1], quantiles_i[j]);
      // printf("a: %f\n", a);
      // printf("b: %f\n", b);
      // Calculate the discriminant to make sure the cubic function isn't decreasing.
      double discriminant = (1/(3*h*h))*(-6*alphas[j]*y2_i[j] + 6*alphas[j+1]*y2_i[j] + h*h*y2_i[j]*y2_i[j] + 6*alphas[j]*y2_i[j+1] - 6*alphas[j+1]*y2_i[j+1] - 2*h*h*y2_i[j]*y2_i[j+1] + 3*quantiles_full(i, j + 1)*quantiles_full(i, j + 1)*y2_i[j]*y2_i[j+1] - 6*quantiles_full(i, j + 1)*quantiles_full(i, j + 2)*y2_i[j]*y2_i[j+1] + 3*quantiles_full(i, j + 2)*quantiles_full(i, j + 2)*y2_i[j]*y2_i[j+1] + h*h*y2_i[j+1]*y2_i[j+1]);
      // printf("disc: %f\n", discriminant);
      if (discriminant >= 0) {
        // The slope was negative at some point, now we check if it is negative inside the interval.
        double qa = y2_i[j+1]/(2*h) - (y2_i[j]/(2*h));
        double qb = (quantiles_full(i, j + 2)*y2_i[j])/h - (quantiles_full(i, j + 1)*y2_i[j+1])/h;
        double root1 = (-qb + sqrt(discriminant)) / (2*qa);
        double root2 = (-qb - sqrt(discriminant)) / (2*qa);
        // printf("roots: %f %f\n", *root1, *root2);
        if ((root1 >= quantiles_full(i, j + 1) && root1 <= quantiles_full(i, j + 2)) || (root2 >= quantiles_full(i, j + 1) && root2 <= quantiles_full(i, j + 2))) {
          // The function was decreasing on the interval, so linearly-interpolate between the quantiles.
          // printf("in int\n");
          if (distn[i]) { // CDF
            y_hat[i] = (y[i] - quantiles_full(i, j + 1)) * (alphas[j+1] - alphas[j]) / h + alphas[j];
          } else { // PDF
            y_hat[i] = (alphas[j+1] - alphas[j]) / h;
          }
          continue;
        }
      }
      // printf("pos slope\n");
      // Everything checked out, so splint together the cubic function using the second derivatives and the location of y inside the interval.
      if (distn[i] == 1) { // CDF
        // printf("distn\n");
        y_hat[i] = a*alphas[j] + b*alphas[j+1] + ((a*a*a-a)*y2_i[j] + (b*b*b-b)*y2_i[j+1]) * (h*h) / 6;
      } else { // PDF (just the derivative of the CDF)
        // printf("pdf\n");
        double b_prime = 1 / h;
        double a_prime = -b_prime;
        y_hat[i] = a_prime*alphas[j] + b_prime*alphas[j+1] + ((3*(a*a)*a_prime-a_prime)*y2_i[j] + (3*(b*b)*b_prime-b_prime)*y2_i[j+1]) * (h*h) / 6;
      }
    }
  }
  return y_hat;
}


