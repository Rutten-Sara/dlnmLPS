// Helpfunctions model comparison criteria


#include <RcppEigen.h>
#include <cmath>

// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
double Hutchinson_pD(
    const Eigen::SparseMatrix<double> &Ipost,
    const Eigen::SparseMatrix<double> &Ilike,
    int R = 5000)  // number of random probes
{
  const int p = Ipost.rows();
  double trace_est = 0.0;

  // Sparse Cholesky factorization of Ipost
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Ipost);
  if (chol.info() != Eigen::Success) Rcpp::stop("Cholesky decomposition failed");

  Rcpp::RNGScope scope;
  for (int r = 0; r < R; r++) {
    // ±1 Rademacher random vector
    Eigen::VectorXd z(p);
    for (int i=0; i<p; i++) z[i] = (R::runif(0,1) < 0.5 ? -1.0 : 1.0);

    // b = I_like * z
    Eigen::VectorXd b = Ilike * z;

    // Solve I_post * x = b
    Eigen::VectorXd x = chol.solve(b);

    trace_est += z.dot(x);
  }

  trace_est /= R;
  return trace_est;
}




// [[Rcpp::export]]
double Hutchinson_pD_con(
    const Eigen::SparseMatrix<double> &Ipost,
    const Eigen::SparseMatrix<double> &Ilike,
    const Eigen::SparseMatrix<double> &A,
    int R = 5000)  // number of random probes
{
  const int p = Ipost.rows();
  double trace_est = 0.0;
  const int m = A.rows();

  // Sparse Cholesky factorization of Ipost
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Ipost);
  if (chol.info() != Eigen::Success) Rcpp::stop("Cholesky decomposition failed");


  // Compute A H^{-1} A'
  Eigen::MatrixXd At_dense = A.transpose();
  Eigen::MatrixXd HinvAt = chol.solve(At_dense);
  Eigen::MatrixXd S = A * HinvAt;

  Eigen::LDLT<Eigen::MatrixXd> S_solver(S);

  if (S_solver.info() != Eigen::Success)
    Rcpp::stop("Constraint system failed");


  Rcpp::RNGScope scope;
  for (int r = 0; r < R; r++) {
    // ±1 Rademacher random vector
    Eigen::VectorXd z(p);
    for (int i=0; i<p; i++) z[i] = (R::runif(0,1) < 0.5 ? -1.0 : 1.0);

    // b = I_like * z
    Eigen::VectorXd b = Ilike * z;

    // Solve I_post * x = b
    Eigen::VectorXd Hinv_b = chol.solve(b);

    // Constraint correction:
    // H^-1 A' (A H^-1 A')^-1 A H^-1 b
    Eigen::VectorXd AHinv_b = A * Hinv_b;

    Eigen::VectorXd correction =
      HinvAt * S_solver.solve(AHinv_b);

    Eigen::VectorXd Vc_b =
      Hinv_b - correction;


    trace_est += z.dot(Vc_b);
  }

  trace_est /= R;
  return trace_est;
}










// [[Rcpp::export]]
Eigen::MatrixXd simulate_from_precision(const Eigen::SparseMatrix<double> &Prec,
                                                const int nsim = 500) {
  const int n = Prec.rows();
  Eigen::MatrixXd draws(n, nsim);

  // 1) Sparse Cholesky decomposition: Prec = P^T L L^T P
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(Prec);
  if (solver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed.");

  // 2) Generate standard normals
  Eigen::MatrixXd Z = Eigen::MatrixXd::NullaryExpr(n, nsim, [](){ return R::rnorm(0.0, 1.0); });

  // 3) Solve L^T Y = Z (upper-triangular solve)
  Eigen::MatrixXd Y = solver.matrixU().triangularView<Eigen::Upper>().solve(Z);

  // 4) Apply permutation P^T to get correct order
  draws = solver.permutationP().transpose() * Y;

  return draws;
}





// [[Rcpp::export]]
Eigen::MatrixXd simulate_from_precision_con(const Eigen::SparseMatrix<double> &Prec,
                                                    const Eigen::SparseMatrix<double> &A,
                                                    const int nsim = 500) {
  const int n = Prec.rows();
  const int m = A.rows();

  Eigen::MatrixXd draws(n, nsim);

  // 1) Sparse Cholesky decomposition: Prec = P^T L L^T P
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(Prec);
  if (solver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed.");

  // 2) Generate standard normals
  Eigen::MatrixXd Z = Eigen::MatrixXd::NullaryExpr(n, nsim, [](){ return R::rnorm(0.0, 1.0); });

  // 3) Solve L^T Y = Z (upper-triangular solve)
  Eigen::MatrixXd Y = solver.matrixU().triangularView<Eigen::Upper>().solve(Z);

  // 4) Apply permutation P^T to get correct order
  draws = solver.permutationP().transpose() * Y;


  // ----- Build Schur complement -----

  Eigen::MatrixXd At_dense = A.transpose();

  Eigen::MatrixXd HinvAt = solver.solve(At_dense);

  Eigen::MatrixXd S = A * HinvAt;

  double eps = 1e-5;
  S += eps * Eigen::MatrixXd::Identity(m, m);

  Eigen::LDLT<Eigen::MatrixXd> S_solver(S);

  if (S_solver.info() != Eigen::Success)
    Rcpp::stop("Schur complement decomposition failed.");

  // ----- Project every draw -----

  for (int s = 0; s < nsim; ++s) {

    Eigen::VectorXd Ax = A * draws.col(s);

    Eigen::VectorXd correction =
      HinvAt * S_solver.solve(Ax);

    draws.col(s) -= correction;
  }

  return draws;
}





// [[Rcpp::export]]
Rcpp::NumericMatrix poisson_loglik_samples_cpp(
    Rcpp::NumericVector xi_mode,
    Rcpp::NumericMatrix betas,
    Rcpp::S4 Xv,
    Rcpp::NumericVector y,
    Rcpp::NumericVector offset,
    SEXP v_mode,
    SEXP Ntrials
) {
  // Extract sparse matrix slots
  Rcpp::IntegerVector i = Xv.slot("i");   // row indices
  Rcpp::IntegerVector p = Xv.slot("p");   // column pointers
  Rcpp::NumericVector x = Xv.slot("x");   // values
  Rcpp::IntegerVector dims = Xv.slot("Dim");

  int n = dims[0];
  int p_dim = dims[1];
  int S = betas.ncol();

  Rcpp::NumericMatrix result(n, S);

  // Precompute eta_base
  Rcpp::NumericVector eta_base = Rcpp::clone(offset);

  for (int col = 0; col < p_dim; col++) {
    for (int idx = p[col]; idx < p[col + 1]; idx++) {
      int row = i[idx];
      eta_base[row] += x[idx] * xi_mode[col];
    }
  }

  // Precompute log-factorial
  Rcpp::NumericVector lfy(n);
  for (int r = 0; r < n; r++) {
    lfy[r] = std::lgamma(y[r] + 1.0);
  }

  // Main computation
  for (int s = 0; s < S; s++) {

    // Start from eta_base
    Rcpp::NumericVector eta = Rcpp::clone(eta_base);

    // Add Xv %*% betas[, s] using sparse structure
    for (int col = 0; col < p_dim; col++) {
      double beta_val = betas(col, s);
      if (beta_val == 0.0) continue;

      for (int idx = p[col]; idx < p[col + 1]; idx++) {
        int row = i[idx];
        eta[row] += x[idx] * beta_val;
      }
    }

    // Final likelihood
    for (int r = 0; r < n; r++) {
      double mu = std::exp(eta[r]);
      result(r, s) = y[r] * eta[r] - mu - lfy[r];
    }
  }

  return result;
}







// [[Rcpp::export]]
Rcpp::NumericMatrix gaussian_loglik_samples_cpp(
    Rcpp::NumericVector xi_mode,
    Rcpp::NumericMatrix betas,
    Rcpp::S4 Xv,
    Rcpp::NumericVector y,
    Rcpp::NumericVector offset,
    SEXP v_mode,
    SEXP Ntrials
) {
  // Extract sparse matrix slots
  Rcpp::IntegerVector i = Xv.slot("i");   // row indices
  Rcpp::IntegerVector p = Xv.slot("p");   // column pointers
  Rcpp::NumericVector x = Xv.slot("x");   // values
  Rcpp::IntegerVector dims = Xv.slot("Dim");

  int n = dims[0];
  int p_dim = dims[1];
  int S = betas.ncol();
  const double pi = std::acos(-1.0);

  Rcpp::NumericVector v(v_mode);
  double tau = std::exp(v[0]);

  Rcpp::NumericMatrix result(n, S);

  // Precompute eta_base
  Rcpp::NumericVector eta_base = Rcpp::clone(offset);

  for (int col = 0; col < p_dim; col++) {
    for (int idx = p[col]; idx < p[col + 1]; idx++) {
      int row = i[idx];
      eta_base[row] += x[idx] * xi_mode[col];
    }
  }


  // Main computation
  for (int s = 0; s < S; s++) {

    // Start from eta_base
    Rcpp::NumericVector eta = Rcpp::clone(eta_base);

    // Add Xv %*% betas[, s] using sparse structure
    for (int col = 0; col < p_dim; col++) {
      double beta_val = betas(col, s);
      if (beta_val == 0.0) continue;

      for (int idx = p[col]; idx < p[col + 1]; idx++) {
        int row = i[idx];
        eta[row] += x[idx] * beta_val;
      }
    }

    // Final likelihood
    for (int r = 0; r < n; r++) {
      double resid = y[r] - eta[r];

      result(r, s) =
        -0.5 * log(2*pi)
        +0.5 * log(tau)
        -0.5 * resid * resid *tau;

    }
  }

  return result;
}




// [[Rcpp::export]]
Rcpp::NumericMatrix nb_loglik_samples_cpp(
    Rcpp::NumericVector xi_mode,
    Rcpp::NumericMatrix betas,
    Rcpp::S4 Xv,
    Rcpp::NumericVector y,
    Rcpp::NumericVector offset,
    SEXP v_mode,
    SEXP Ntrials
) {
  // Extract sparse matrix slots
  Rcpp::IntegerVector i = Xv.slot("i");   // row indices
  Rcpp::IntegerVector p = Xv.slot("p");   // column pointers
  Rcpp::NumericVector x = Xv.slot("x");   // values
  Rcpp::IntegerVector dims = Xv.slot("Dim");

  int n = dims[0];
  int p_dim = dims[1];
  int S = betas.ncol();

  Rcpp::NumericMatrix result(n, S);

  // Extract phi
  Rcpp::NumericVector v(v_mode);
  double phi = std::exp(v[0]);
  double inv_phi = 1.0 / phi;


  // Precompute eta_base
  Rcpp::NumericVector eta_base = Rcpp::clone(offset);

  for (int col = 0; col < p_dim; col++) {
    for (int idx = p[col]; idx < p[col + 1]; idx++) {
      int row = i[idx];
      eta_base[row] += x[idx] * xi_mode[col];
    }
  }

  // Precompute log-factorial
  Rcpp::NumericVector lfy(n);
  for (int r = 0; r < n; r++) {
    lfy[r] = std::lgamma(y[r] + 1.0);
  }

  // Main computation
  for (int s = 0; s < S; s++) {

    // Start from eta_base
    Rcpp::NumericVector eta = Rcpp::clone(eta_base);

    // Add Xv %*% betas[, s] using sparse structure
    for (int col = 0; col < p_dim; col++) {
      double beta_val = betas(col, s);
      if (beta_val == 0.0) continue;

      for (int idx = p[col]; idx < p[col + 1]; idx++) {
        int row = i[idx];
        eta[row] += x[idx] * beta_val;
      }
    }

    // Final likelihood
    for (int r = 0; r < n; r++) {
      double mu = std::exp(eta[r]);
      double gamma_nb =
        phi * std::log(mu / (mu + phi));

      double beta_gamma_nb =
        -phi * phi * std::log(phi / (phi + mu));

        result(r, s) =
        inv_phi * (y[r] * gamma_nb - beta_gamma_nb)
          + std::lgamma(y[r] + phi)
          - std::lgamma(phi)
          - lfy[r];
    }
  }

  return result;
}



// [[Rcpp::export]]
Rcpp::NumericMatrix logistic_loglik_samples_cpp(
    Rcpp::NumericVector xi_mode,
    Rcpp::NumericMatrix betas,
    Rcpp::S4 Xv,
    Rcpp::NumericVector y,
    Rcpp::NumericVector offset,
    SEXP v_mode,
    SEXP Ntrials
) {
  // Extract sparse matrix slots
  Rcpp::IntegerVector i = Xv.slot("i");   // row indices
  Rcpp::IntegerVector p = Xv.slot("p");   // column pointers
  Rcpp::NumericVector x = Xv.slot("x");   // values
  Rcpp::IntegerVector dims = Xv.slot("Dim");

  // Ntrials
  Rcpp::NumericVector ntrials(Ntrials);

  int n = dims[0];
  int p_dim = dims[1];
  int S = betas.ncol();

  Rcpp::NumericMatrix result(n, S);

  // Precompute eta_base
  Rcpp::NumericVector eta_base = Rcpp::clone(offset);

  for (int col = 0; col < p_dim; col++) {
    for (int idx = p[col]; idx < p[col + 1]; idx++) {
      int row = i[idx];
      eta_base[row] += x[idx] * xi_mode[col];
    }
  }

  // Precompute log-factorial
  Rcpp::NumericVector lfy(n);
  for (int r = 0; r < n; r++) {
    lfy[r] = std::lgamma(y[r] + 1.0);
  }

  // Main computation
  for (int s = 0; s < S; s++) {

    // Start from eta_base
    Rcpp::NumericVector eta = Rcpp::clone(eta_base);

    // Add Xv %*% betas[, s] using sparse structure
    for (int col = 0; col < p_dim; col++) {
      double beta_val = betas(col, s);
      if (beta_val == 0.0) continue;

      for (int idx = p[col]; idx < p[col + 1]; idx++) {
        int row = i[idx];
        eta[row] += x[idx] * beta_val;
      }
    }

    // Final likelihood
    for (int r = 0; r < n; r++) {
      double mu = std::exp(eta[r]);
      double p_mu = mu/(1.0+mu);
      result(r, s) = std::pow(p_mu, y[r]) *
          std::pow(1.0 - p_mu, ntrials[r] - y[r]);
    }
  }

  return result;
}









