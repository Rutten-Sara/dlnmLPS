// Helpfunctions uncertainty calculation


#include <RcppEigen.h>
#include <cmath>
#include <algorithm>
#include <vector>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::VectorXd Predict_se_con(const Eigen::SparseMatrix<double> &Prec,
                               const Eigen::MatrixXd &Xpred,
                               const Rcpp::IntegerVector &ind,
                               const Eigen::SparseMatrix<double> &A) {

  const int n = Xpred.rows();
  const int full_dim = Prec.rows();
  const int m = A.rows();

  Eigen::VectorXd se(n);

  // H = -Prec (Make sure Prec is passed matching your Newton sign convention)
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Prec);
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition of H failed");
  }

  // -------------------------------------------------------------
  // 1. Build S = A * H^{-1} * A^T efficiently in a single batch
  // -------------------------------------------------------------
  Eigen::MatrixXd At_dense = A.transpose();
  Eigen::MatrixXd H_inv_At = chol.solve(At_dense); // Multi-solve (n x m) done instantly
  Eigen::MatrixXd S = A * H_inv_At;                // (m x m) dense matrix multiplication

  // Set up a dense solver for S instead of explicitly inverting it
  Eigen::LDLT<Eigen::MatrixXd> s_solver(S);
  if (s_solver.info() != Eigen::Success) {
    Rcpp::stop("Schur complement decomposition failed");
  }

  // -------------------------------------------------------------
  // 2. Compute Predictions Standard Errors
  // -------------------------------------------------------------
  for (int i = 0; i < n; ++i) {

    Eigen::VectorXd x = Eigen::VectorXd::Zero(full_dim);
    for (int j = 0; j < ind.size(); ++j) {
      x[ind[j] - 1] = Xpred(i, j);
    }

    // t = H^{-1} x
    Eigen::VectorXd Hx = chol.solve(x);

    // Unconstrained variance: x^T * H^{-1} * x
    double v1 = x.dot(Hx);

    // Constraint term correction: (A * H^{-1} * x)^T * S^{-1} * (A * H^{-1} * x)
    Eigen::VectorXd AHx = A * Hx;
    double v2 = AHx.dot(s_solver.solve(AHx)); // Uses LDLT solve instead of explicit S_inv multiplication

    double var = v1 - v2;

    // Numerical guard for values close to zero
    se[i] = std::sqrt(std::max(var, 0.0));
  }

  return se;
}



// [[Rcpp::export]]
Eigen::VectorXd Predict_se(const Eigen::SparseMatrix<double> &Prec,
                           const Eigen::MatrixXd &Xpred,
                           const Rcpp::IntegerVector &ind) {

  const int n = Xpred.rows();
  const int full_dim = Prec.rows();
  Eigen::VectorXd se(n);


  // Cholesky factorization
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Prec);
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed");
  }


  for (int i = 0; i < n; ++i) {
    // Zero-padded vector of full_dim
    Eigen::VectorXd xfull = Eigen::VectorXd::Zero(full_dim);
    for (int j = 0; j < ind.size(); ++j) {
      xfull[ind[j]-1] = Xpred(i, j);  // put row i of Xpred into full-size vector
    }

    Eigen::VectorXd zi = chol.solve(xfull);
    se[i] = std::sqrt(xfull.dot(zi));
  }

  return se;
}



// [[Rcpp::export]]
Eigen::MatrixXd Predict_se_all_con(const Eigen::SparseMatrix<double> &Prec,
                                   const Eigen::MatrixXd &Xpred,
                                   const Rcpp::IntegerVector &ind_cb,
                                   const Rcpp::List &ind_cb_spat_list,
                                   const Eigen::SparseMatrix<double> &A) {

  const int n_areas = ind_cb_spat_list.size();
  const int n_preds = Xpred.rows();
  const int full_dim = Prec.rows();
  const int m = A.rows();

  Eigen::MatrixXd se_matrix(n_areas, n_preds);

  // H = -Prec (Make sure Prec is passed matching your Newton sign convention)
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Prec);
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition of H failed");
  }

  // -------------------------------------------------------------
  // 1. Build S = A * H^{-1} * A^T efficiently in a single batch
  // -------------------------------------------------------------
  Eigen::MatrixXd At_dense = A.transpose();
  Eigen::MatrixXd H_inv_At = chol.solve(At_dense); // Multi-solve (n x m) done instantly
  Eigen::MatrixXd S = A * H_inv_At;                // (m x m) dense matrix multiplication

  // Set up a dense solver for S instead of explicitly inverting it
  Eigen::LDLT<Eigen::MatrixXd> s_solver(S);
  if (s_solver.info() != Eigen::Success) {
    Rcpp::stop("Schur complement decomposition failed");
  }

  // -------------------------------------------------------------
  // 2. Compute Predictions Standard Errors
  // -------------------------------------------------------------
  for (int a = 0; a < n_areas; ++a) {
    Rcpp::IntegerVector ind_cb_spat = ind_cb_spat_list[a];

    // Combine indices
    std::vector<int> ind_combined(ind_cb.size() + ind_cb_spat.size());
    std::copy(ind_cb.begin(), ind_cb.end(), ind_combined.begin());
    std::copy(ind_cb_spat.begin(), ind_cb_spat.end(), ind_combined.begin() + ind_cb.size());
    for (int i = 0; i < n_preds; ++i) {
      Eigen::VectorXd xfull = Eigen::VectorXd::Zero(full_dim);
      for (int j = 0; j < ind_combined.size(); ++j) {
        xfull[ind_combined[j]-1] = Xpred(i, j);
      }

      Eigen::VectorXd Hx = chol.solve(xfull);


      // Unconstrained variance: x^T * H^{-1} * x
      double v1 = xfull.dot(Hx);

      // Constraint term correction: (A * H^{-1} * x)^T * S^{-1} * (A * H^{-1} * x)
      Eigen::VectorXd AHx = A * Hx;
      double v2 = AHx.dot(s_solver.solve(AHx)); // Uses LDLT solve instead of explicit S_inv multiplication

      double var = v1 - v2;

      // Numerical guard for values close to zero
      se_matrix(a, i) = std::sqrt(std::max(var, 0.0));
    }
  }
  return se_matrix;
}



// [[Rcpp::export]]
Eigen::MatrixXd Predict_se_all(
    const Eigen::SparseMatrix<double> &Prec,
    const Eigen::MatrixXd &Xpred,
    const Rcpp::IntegerVector &ind_cb,
    const Rcpp::List &ind_cb_spat_list) {

  const int n_areas = ind_cb_spat_list.size();
  const int n_preds = Xpred.rows();
  const int full_dim = Prec.rows();

  Eigen::MatrixXd se_matrix(n_areas, n_preds);

  // Cholesky decomposition (once!)
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> chol(Prec);
  if (chol.info() != Eigen::Success) {
    Rcpp::stop("Cholesky decomposition failed");
  }

  for (int a = 0; a < n_areas; ++a) {
    Rcpp::IntegerVector ind_cb_spat = ind_cb_spat_list[a];

    // Combine indices
    std::vector<int> ind_combined(ind_cb.size() + ind_cb_spat.size());
    std::copy(ind_cb.begin(), ind_cb.end(), ind_combined.begin());
    std::copy(ind_cb_spat.begin(), ind_cb_spat.end(), ind_combined.begin() + ind_cb.size());

    for (int i = 0; i < n_preds; ++i) {
      Eigen::VectorXd xfull = Eigen::VectorXd::Zero(full_dim);
      for (int j = 0; j < ind_combined.size(); ++j) {
        xfull[ind_combined[j]-1] = Xpred(i, j);
      }

      Eigen::VectorXd zi = chol.solve(xfull);
      se_matrix(a, i) = std::sqrt(xfull.dot(zi));
    }
  }

  return se_matrix;
}

