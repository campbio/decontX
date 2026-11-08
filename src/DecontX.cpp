#include <RcppEigen.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List decontXEM(const Eigen::MappedSparseMatrix<double> &counts,
                     const NumericVector &counts_colsums,
                     const NumericVector &theta,
                     const bool &estimate_eta,
                     const NumericMatrix &eta,
                     const NumericMatrix &phi,
                     const IntegerVector &z,
                     const bool &estimate_delta,
                     const NumericVector &delta,
                     const double &pseudocount)
{

  // Perform error checking. These checks also keep the unchecked indexing
  // in the loop below in bounds, so keep them complete.
  if (counts.cols() != theta.size())
  {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if (counts.cols() != z.size())
  {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }
  if (counts.cols() != counts_colsums.size())
  {
    stop("Length of 'counts_colsums' must be equal to the number of columns in 'counts'.");
  }
  if (counts.rows() != phi.nrow())
  {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if (counts.rows() != eta.nrow())
  {
    stop("The number of rows in 'eta' must be equal to the number of rows in 'counts'.");
  }
  if (phi.ncol() != eta.ncol())
  {
    stop("The number of columns in 'eta' must be equal to the number of columns in 'phi'.");
  }
  if (min(z) < 1 || max(z) > eta.ncol())
  {
    stop("The entries in 'z' need to be between 1 and the number of columns in eta and phi.");
  }
  if (delta.size() != 2 || sum(delta < 0) > 0)
  {
    stop("'delta' must be a numeric vector of length 2 with positive integers.");
  }

  // Declare variables and functions
  NumericVector new_theta(theta.size());
  NumericVector native_total(theta.size());
  NumericMatrix new_phi(phi.nrow(), phi.ncol());
  // Instantiate new_eta but only update it when estimate_eta is TRUE
  NumericMatrix new_eta(eta.nrow(), eta.ncol());


  // Obtaining 'fit_dirichlet' function from MCMCprecision package
  Environment pkg = Environment::namespace_env("MCMCprecision");
  Function f = pkg["fit_dirichlet"];

  int i;
  int j;
  int nr = phi.nrow();
  double x;
  double pcontamin;
  double pnative;
  double normp;
  double px;

  // The loop below runs for every nonzero count on every EM iteration, so
  // it reads and writes through raw pointers. Rcpp's checked accessors
  // (operator() looks up the dim attribute on every call) more than doubled
  // its time. The guards above already check every size and the range of z,
  // and the sparse iterator of a valid dgCMatrix only yields rows below
  // counts.rows().
  const double *phi_p = phi.begin();
  const double *eta_p = eta.begin();
  const double *theta_p = theta.begin();
  const int *z_p = z.begin();
  double *new_phi_p = new_phi.begin();
  double *native_total_p = native_total.begin();
  for (int j = 0; j < counts.cols(); ++j)
  {
    const R_xlen_t col = static_cast<R_xlen_t>(nr) * (z_p[j] - 1);
    const double *phi_k = phi_p + col;
    const double *eta_k = eta_p + col;
    double *new_phi_k = new_phi_p + col;
    const double theta_j = theta_p[j];
    // Same terms in the same order as adding to native_total[j] directly
    double native_j = 0;
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_)
    {
      i = i_.index();
      x = i_.value();

      // Calculate variational probabilities
      // Removing the log/exp speeds it up and produces the same result since
      // there are only 2 probabilities being multiplied

      //pnative = log(phi(i,k) + pseudocount) + log(theta(j) + pseudocount);
      //pcontamin = log(eta(i,k) + pseudocount) + log(1 - theta(j) + pseudocount);
      pnative = (phi_k[i] + pseudocount) * (theta_j + pseudocount);
      pcontamin = (eta_k[i] + pseudocount) * (1 - theta_j + pseudocount);

      // Normalize probabilities and add to proper components
      //normp = exp(pnative) / (exp(pcontamin) + exp(pnative));
      normp = pnative / (pcontamin + pnative);
      px = normp * x;
      new_phi_k[i] += px;
      native_j += px;
    }
    native_total_p[j] = native_j;
  }

  // Calculate Eta using Weights from Phi
  if (estimate_eta == TRUE)
  {
    NumericVector phi_rowsum = rowSums(new_phi);
    for (i = 0; i < new_eta.ncol(); i++)
    {
      for (j = 0; j < new_eta.nrow(); j++)
      {
        new_eta(j, i) = phi_rowsum[j] - new_phi(j, i);
      }
    }
  }

  // Normalize Phi
  NumericVector phi_colsum = colSums(new_phi);
  for (i = 0; i < new_phi.ncol(); i++)
  {
    new_phi(_, i) = new_phi(_, i) / phi_colsum[i];
  }

  // Normalize Eta
  if (estimate_eta == TRUE)
  {
    NumericVector eta_colsum = colSums(new_eta);
    for (i = 0; i < new_phi.ncol(); i++)
    {
      new_eta(_, i) = new_eta(_, i) / eta_colsum[i];
    }
  }

  // Update Theta
  NumericVector contamination_prop = (counts_colsums - native_total) / counts_colsums;
  NumericVector native_prop = 1 - contamination_prop;
  NumericMatrix theta_raw = cbind(native_prop, contamination_prop);

  NumericVector new_delta = delta;
  if (estimate_delta == TRUE)
  {
    Rcpp::List result = f(Named("x", theta_raw));
    new_delta = result["alpha"];
  }

  // Estimate new theta
  new_theta = (native_total + new_delta[0]) / (counts_colsums + sum(new_delta));

  // If not to estimate new_eta, just return input eta
  if (estimate_eta == FALSE)
  {
    new_eta = eta;
  }

  return Rcpp::List::create(Rcpp::Named("phi") = new_phi,
                            Rcpp::Named("eta") = new_eta,
                            Rcpp::Named("theta") = new_theta,
                            Rcpp::Named("delta") = new_delta,
                            Rcpp::Named("contamination") = contamination_prop);
}

// [[Rcpp::export]]
double decontXLogLik(const Eigen::MappedSparseMatrix<double> &counts,
                     const NumericVector &theta,
                     const NumericMatrix &eta,
                     const NumericMatrix &phi,
                     const IntegerVector &z,
                     const double &pseudocount)
{

  // Perform error checking. These checks also keep the unchecked indexing
  // in the loop below in bounds, so keep them complete.
  if (counts.cols() != theta.size())
  {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if (counts.cols() != z.size())
  {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }
  if (counts.rows() != phi.nrow())
  {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if (counts.rows() != eta.nrow())
  {
    stop("The number of rows in 'eta' must be equal to the number of rows in 'counts'.");
  }
  if (phi.ncol() != eta.ncol())
  {
    stop("The number of columns in 'eta' must be equal to the number of columns in 'phi'.");
  }
  if (min(z) < 1 || max(z) > eta.ncol())
  {
    stop("The entries in 'z' need to be between 1 and the number of columns in eta and phi.");
  }

  // Declare variables and functions
  double loglik = 0;

  int i;
  double x;
  int nr = phi.nrow();

  // Original R code:
  // ll <- sum(Matrix::t(counts) * log(theta * t(phi)[z, ] +
  //       (1 - theta) * t(eta)[z, ] + 1e-20))

  // Raw pointers as in decontXEM(); the guards above check every size and
  // the range of z
  const double *phi_p = phi.begin();
  const double *eta_p = eta.begin();
  const double *theta_p = theta.begin();
  const int *z_p = z.begin();
  for (int j = 0; j < counts.cols(); ++j)
  {
    const R_xlen_t col = static_cast<R_xlen_t>(nr) * (z_p[j] - 1);
    const double *phi_k = phi_p + col;
    const double *eta_k = eta_p + col;
    const double theta_j = theta_p[j];
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_)
    {
      i = i_.index();
      x = i_.value();

      loglik += x * log((phi_k[i] * theta_j) + (eta_k[i] * (1 - theta_j)) + pseudocount);
    }
  }

  return loglik;
}

// [[Rcpp::export]]
Rcpp::List decontXInitialize(const Eigen::MappedSparseMatrix<double> &counts,
                             const NumericVector &theta,
                             const IntegerVector &z,
                             const double &pseudocount)
{

  // Perform error checking
  if (counts.cols() != theta.size())
  {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if (counts.cols() != z.size())
  {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }
  // NA is INT_MIN, so this check catches it too
  if (min(z) < 1)
  {
    stop("The entries in 'z' need to be at least 1.");
  }

  // Declare variables and functions
  NumericMatrix new_phi(counts.rows(), max(z));
  NumericMatrix new_eta(counts.rows(), max(z));
  std::fill(new_phi.begin(), new_phi.end(), pseudocount);
  std::fill(new_eta.begin(), new_eta.end(), pseudocount);

  int k;
  int i;
  double x;
  for (int j = 0; j < counts.cols(); ++j)
  {
    for (Eigen::MappedSparseMatrix<double>::InnerIterator i_(counts, j); i_; ++i_)
    {
      i = i_.index();
      x = i_.value();
      k = z[j] - 1;

      new_phi(i, k) += x * theta(j);
    }
  }

  // Calculate Eta using Weights from Phi
  NumericVector phi_rowsum = rowSums(new_phi);
  int j;
  for (i = 0; i < new_eta.ncol(); i++)
  {
    for (j = 0; j < new_eta.nrow(); j++)
    {
      new_eta(j, i) = phi_rowsum[j] - new_phi(j, i);
    }
  }

  // Normalize Phi and Eta
  NumericVector phi_colsum = colSums(new_phi);
  NumericVector eta_colsum = colSums(new_eta);
  for (i = 0; i < new_phi.ncol(); i++)
  {
    new_phi(_, i) = new_phi(_, i) / phi_colsum[i];
    new_eta(_, i) = new_eta(_, i) / eta_colsum[i];
  }

  return Rcpp::List::create(Rcpp::Named("phi") = new_phi,
                            Rcpp::Named("eta") = new_eta);
}

// The native counts: each nonzero count times its native probability.
// Returns only the new values, in the order of the counts' x slot; the
// caller puts them in a copy of the counts, which shares the i and p slots.
// [[Rcpp::export]]
NumericVector calculateNativeCounts(const Eigen::MappedSparseMatrix<double> &counts,
                                    const NumericVector &theta,
                                    const NumericMatrix &eta,
                                    const NumericMatrix &phi,
                                    const IntegerVector &z,
                                    const double &pseudocount)
{

  // Perform error checking. These checks also keep the unchecked indexing
  // in the loop below in bounds, so keep them complete.
  if (counts.cols() != theta.size())
  {
    stop("Length of 'theta' must be equal to the number of columns in 'counts'.");
  }
  if (counts.cols() != z.size())
  {
    stop("Length of 'z' must be equal to the number of columns in 'counts'.");
  }
  if (counts.rows() != phi.nrow())
  {
    stop("The number of rows in 'phi' must be equal to the number of rows in 'counts'.");
  }
  if (counts.rows() != eta.nrow())
  {
    stop("The number of rows in 'eta' must be equal to the number of rows in 'counts'.");
  }
  if (phi.ncol() != eta.ncol())
  {
    stop("The number of columns in 'eta' must be equal to the number of columns in 'phi'.");
  }
  if (min(z) < 1 || max(z) > eta.ncol())
  {
    stop("The entries in 'z' need to be between 1 and the number of columns in eta and phi.");
  }

  // A dgCMatrix is always compressed, so the column starts (p), row indices
  // (i) and values (x) are R's own vectors, read here without copying
  const int *col_start = counts.outerIndexPtr();
  const int *row_index = counts.innerIndexPtr();
  const double *count_value = counts.valuePtr();
  NumericVector native_value(col_start[counts.cols()]);
  double *native_p = native_value.begin();

  // Raw pointers as in decontXEM()
  const double *phi_p = phi.begin();
  const double *eta_p = eta.begin();
  const double *theta_p = theta.begin();
  const int *z_p = z.begin();
  int nr = phi.nrow();
  int i;
  double pcontamin;
  double pnative;
  double normp;
  for (int j = 0; j < counts.cols(); ++j)
  {
    const R_xlen_t col = static_cast<R_xlen_t>(nr) * (z_p[j] - 1);
    const double *phi_k = phi_p + col;
    const double *eta_k = eta_p + col;
    const double theta_j = theta_p[j];
    for (int idx = col_start[j]; idx < col_start[j + 1]; ++idx)
    {
      i = row_index[idx];

      // Calculate variational probabilities, in the same direct form as
      // decontXEM(). This used to go through log() and exp(), which agrees
      // to within about 1e-14 relative and is much slower.
      pnative = (phi_k[i] + pseudocount) * (theta_j + pseudocount);
      pcontamin = (eta_k[i] + pseudocount) * (1 - theta_j + pseudocount);

      // Normalize probabilities and scale the count
      normp = pnative / (pcontamin + pnative);
      native_p[idx] = count_value[idx] * normp;
    }
  }

  return native_value;
}
