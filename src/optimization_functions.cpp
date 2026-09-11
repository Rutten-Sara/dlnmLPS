// Helpfunctions parameter optimization



#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::SparseMatrix<double> Hess_logpxi(const Eigen::SparseMatrix<double> &Qv,
                                        const Eigen::VectorXd &Cvxi,
                                        const Eigen::SparseMatrix<double> &Xv) {
  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;

  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }

  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;

  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_W_Xv;
  result += Qv;
  result *= -1;

  return result;
}



// [[Rcpp::export]]


Eigen::SparseMatrix<double> Hess_logpxi_NB(const Eigen::SparseMatrix<double> &Qv,
                                           const Eigen::SparseMatrix<double> & Wnb,
                                           const Eigen::SparseMatrix<double> &MVnb,
                                           const Eigen::SparseMatrix<double> &Xv) {

  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * Wnb*Xv;
  Eigen::SparseMatrix<double> Xt_M_Xv = Xv.transpose() * MVnb*Xv;

  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_M_Xv;
  result -= Xt_W_Xv;
  result -= Qv;

  return result;

}




// [[Rcpp::export]]
Eigen::MatrixXd solve_sparse_cholesky(const Eigen::SparseMatrix<double> &Qv,
                                      const Eigen::VectorXd &Cvxi,
                                      const Eigen::SparseMatrix<double> &Xv,
                                      const Eigen::VectorXd &b) {

  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;

  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }

  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;

  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_W_Xv;
  result += Qv;
  result *= -1;

  result += 1e-5 * Eigen::MatrixXd::Identity(result.rows(), result.cols()).sparseView();


  // Step 4: Solve result*x = b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(result);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed");

  return solver.solve(b);

}




// [[Rcpp::export]]

Eigen::MatrixXd solve_sparse_cholesky_NB(const Eigen::SparseMatrix<double> &Qv,
                                         const Eigen::SparseMatrix<double> & Wnb,
                                         const Eigen::SparseMatrix<double> &MVnb,
                                         const Eigen::SparseMatrix<double> &Xv,
                                         const Eigen::VectorXd &b) {

  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * Wnb*Xv;
  Eigen::SparseMatrix<double> Xt_M_Xv = Xv.transpose() * MVnb*Xv;

  // Result = -Xt_W_Xv - Qv
  Eigen::SparseMatrix<double> result = Xt_M_Xv;
  result -= Xt_W_Xv;
  result -= Qv;

  result += 1e-5 * Eigen::MatrixXd::Identity(result.rows(), result.cols()).sparseView();


  // Step 4: Solve result*x = b
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(result);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Cholesky decomposition failed");

  return solver.solve(b);

}



// [[Rcpp::export]]
Eigen::VectorXd solve_gaussian_mode(
    const Eigen::SparseMatrix<double> &Qv,
    const Eigen::SparseMatrix<double> &XtX,
    const Eigen::VectorXd &Xty,
    const double v1
) {

  Eigen::SparseMatrix<double> A = v1 * XtX;
  A += Qv;

  A += 1e-5 *
    Eigen::MatrixXd::Identity(A.rows(), A.cols()).sparseView();

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Sparse LDLT decomposition failed");

  Eigen::VectorXd xi = solver.solve(Xty);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Sparse LDLT solve failed");

  return v1 * xi;
}




// [[Rcpp::export]]
Eigen::MatrixXd solve_sparse_cholesky_constr(const Eigen::SparseMatrix<double> &Qv,
                                             const Eigen::VectorXd &Cvxi,
                                             const Eigen::SparseMatrix<double> &Xv,
                                             const Eigen::SparseMatrix<double>& A,
                                             const Eigen::VectorXd &b) {

  const int n = Qv.rows();
  const int m = A.rows();

  // -------------------------
  // 1. Build H
  // -------------------------
  Eigen::SparseMatrix<double> WX = Xv;

  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }

  Eigen::SparseMatrix<double> H = Qv + Xv.transpose() * WX;
  H *= -1.0;

  H += 1e-5 *
    Eigen::MatrixXd::Identity(n, n).sparseView();

  if (H.rows() != H.cols())
    Rcpp::stop("H is not square");

  if (A.cols() != H.rows())
    Rcpp::stop(
      "Dimension mismatch: A has %d columns but H has %d rows",
      A.cols(), H.rows()
    );

  if (b.size() != H.rows())
    Rcpp::stop(
      "Dimension mismatch: b has length %d but H has %d rows",
      b.size(), H.rows()
    );

  // -------------------------
  // 2. Factorize H once
  // -------------------------
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(H);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("H factorization failed");

  // -------------------------
  // 3. Solve block systems:
  //    Y = H^{-1} A^T   (n x m)
  // -------------------------
  Eigen::MatrixXd At_dense = A.transpose();   // sparse → dense (only once)
  Eigen::MatrixXd Y = solver.solve(At_dense);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("H solve for A^T failed");

  // -------------------------
  // 4. Schur matrix: S = A * Y
  // -------------------------
  Eigen::MatrixXd S = A * Y;   // m x m

  S += 1e-5 *
    Eigen::MatrixXd::Identity(m, m);

  // -------------------------
  // 5. RHS: rhs = A * H^{-1} b
  // -------------------------
  Eigen::VectorXd Hb = solver.solve(b);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("H solve for b failed");

  Eigen::VectorXd rhs = A * Hb;

  // -------------------------
  // 6. Solve small system
  // -------------------------
  Eigen::LDLT<Eigen::MatrixXd> schur_solver;
  schur_solver.compute(S);

  if (schur_solver.info() != Eigen::Success)
    Rcpp::stop("Schur factorization failed");

  Eigen::VectorXd lambda = schur_solver.solve(rhs);

  if (schur_solver.info() != Eigen::Success)
    Rcpp::stop("Schur solve failed");

  // -------------------------
  // 7. Recover primal step
  // -------------------------
  Eigen::VectorXd dxi = Hb - solver.solve(A.transpose() * lambda);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Final H solve failed");

  return dxi;
}



// [[Rcpp::export]]


Eigen::MatrixXd solve_sparse_cholesky_NB_constr(const Eigen::SparseMatrix<double> &Qv,
                                                const Eigen::SparseMatrix<double> & Wnb,
                                                const Eigen::SparseMatrix<double> &MVnb,
                                                const Eigen::SparseMatrix<double> &Xv,
                                                const Eigen::SparseMatrix<double>& A,
                                                const Eigen::VectorXd &b) {

  const int n = Qv.rows();
  const int m = A.rows();

  // -------------------------
  // 1. Build H
  // -------------------------
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * Wnb*Xv;
  Eigen::SparseMatrix<double> Xt_M_Xv = Xv.transpose() * MVnb*Xv;

  Eigen::SparseMatrix<double> H =  Xt_M_Xv - Xt_W_Xv - Qv;

  H += 1e-5 *
    Eigen::MatrixXd::Identity(n, n).sparseView();


  if (H.rows() != H.cols())
    Rcpp::stop("H is not square");

  if (A.cols() != H.rows())
    Rcpp::stop(
      "Dimension mismatch: A has %d columns but H has %d rows",
      A.cols(), H.rows()
    );

  if (b.size() != H.rows())
    Rcpp::stop(
      "Dimension mismatch: b has length %d but H has %d rows",
      b.size(), H.rows()
    );

  // -------------------------
  // 2. Factorize H once
  // -------------------------
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(H);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("H factorization failed");


  // -------------------------
  // 3. Solve block systems:
  //    Y = H^{-1} A^T   (n x m)
  // -------------------------
  Eigen::MatrixXd At_dense = A.transpose();
  Eigen::MatrixXd Y = solver.solve(At_dense);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("H solve for A^T failed");

  // -------------------------
  // 4. Schur matrix: S = A * Y
  // -------------------------
  Eigen::MatrixXd S = A * Y;   // m x m

  S += 1e-5 *
    Eigen::MatrixXd::Identity(m, m);

  // -------------------------
  // 5. RHS: rhs = A * H^{-1} b
  // -------------------------
  Eigen::VectorXd Hb = solver.solve(b);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("H solve for b failed");

  Eigen::VectorXd rhs = A * Hb;

  // -------------------------
  // 6. Solve small system
  // -------------------------
  Eigen::LDLT<Eigen::MatrixXd> schur_solver;
  schur_solver.compute(S);

  if (schur_solver.info() != Eigen::Success)
    Rcpp::stop("Schur factorization failed");

  Eigen::VectorXd lambda = schur_solver.solve(rhs);

  if (schur_solver.info() != Eigen::Success)
    Rcpp::stop("Schur solve failed");

  // -------------------------
  // 7. Recover primal step
  // -------------------------
  Eigen::VectorXd dxi = Hb - solver.solve(A.transpose() * lambda);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Final H solve failed");


  return dxi;
}




// [[Rcpp::export]]
Eigen::VectorXd solve_gaussian_mode_constr(
    const Eigen::SparseMatrix<double> &Qv,
    const Eigen::SparseMatrix<double> &XtX,
    const Eigen::VectorXd &Xty,
    const double v1,
    const Eigen::SparseMatrix<double> &C
) {
  const int n = Qv.rows();
  const int m = C.rows();

  if (Qv.cols() != n)
    Rcpp::stop("Qv must be square");

  if (XtX.rows() != n || XtX.cols() != n)
    Rcpp::stop("XtX dimensions do not match Qv");

  if (Xty.size() != n)
    Rcpp::stop("Xty length does not match Qv");

  if (C.cols() != n)
    Rcpp::stop(
      "C has %d columns but the parameter dimension is %d",
      C.cols(), n
    );

  // ------------------------------------------------------------
  // 1. Gaussian Hessian
  // ------------------------------------------------------------
  Eigen::SparseMatrix<double> G = v1 * XtX;
  G += Qv;

  G += 1e-5 *
    Eigen::MatrixXd::Identity(n, n).sparseView();

  // ------------------------------------------------------------
  // 2. Factorize G ONCE
  // ------------------------------------------------------------
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
  solver.compute(G);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Sparse LDLT decomposition failed");

  // ------------------------------------------------------------
  // 3. Unconstrained solution
  // ------------------------------------------------------------
  Eigen::VectorXd x0 = v1*solver.solve(Xty);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Sparse LDLT solve failed");

  // ------------------------------------------------------------
  // 4. Solve
  //
  //     Y = G^{-1} C^T
  // ------------------------------------------------------------
  Eigen::MatrixXd Ct_dense = Eigen::MatrixXd(C.transpose());

  Eigen::MatrixXd Y = solver.solve(Ct_dense);

  if (solver.info() != Eigen::Success)
    Rcpp::stop("Sparse LDLT solve for C^T failed");

  // ------------------------------------------------------------
  // 5. Schur complement
  //
  //     S = C G^{-1} C^T = C Y
  // ------------------------------------------------------------
  Eigen::MatrixXd S = C * Y;

  // ------------------------------------------------------------
  // 6. Constraint RHS
  //
  //     r = C x0
  // ------------------------------------------------------------
  Eigen::VectorXd r = C * x0;

  // ------------------------------------------------------------
  // 7. Solve
  //
  //     S lambda = r
  // ------------------------------------------------------------
  Eigen::LDLT<Eigen::MatrixXd> schur_solver;
  schur_solver.compute(S);

  if (schur_solver.info() != Eigen::Success)
    Rcpp::stop("Schur factorization failed");

  Eigen::VectorXd lambda = schur_solver.solve(r);

  if (schur_solver.info() != Eigen::Success)
    Rcpp::stop("Schur solve failed");

  // ------------------------------------------------------------
  // 8. Correct the unconstrained solution
  //
  //     xi = x0 - G^{-1} C^T lambda
  // ------------------------------------------------------------
  Eigen::VectorXd xi = x0 - Y * lambda;

  // ------------------------------------------------------------
  // 9. Preserve the convention of the unconstrained function
  // ------------------------------------------------------------
  return xi;
}




// [[Rcpp::export]]
double determinant_Pv(const Eigen::SparseMatrix<double> &Pv) {
  // Perform LDLT decomposition
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(Pv);

  if (ldlt.info() != Eigen::Success) {
    Rcpp::stop("LDLT decomposition failed");
  }

  // log-determinant = sum(log(diagonal of D))
  Eigen::VectorXd D_diag = ldlt.vectorD();  // for SimplicialLDLT, this is the diagonal of D
  if ((D_diag.array() <= 0).any()) {
    Rcpp::stop("Non-positive values in D; log-determinant is not defined");
  }

  double logdet = D_diag.array().log().sum();
  return logdet;

}





// [[Rcpp::export]]
Eigen::SparseMatrix<double> XWX_func (const Eigen::VectorXd &Cvxi,
                                      const Eigen::SparseMatrix<double> &Xv,
                                      const Eigen::SparseMatrix<double> &tXv) {

  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;

  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }

  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = tXv * WX;

  // Result = Xt_W_Xv
  Eigen::SparseMatrix<double> result = Xt_W_Xv;

  return result;
}




// [[Rcpp::export]]
Eigen::SparseMatrix<double> XWX_func_approx (const Eigen::VectorXd &Cvxi,
                                       const Eigen::SparseMatrix<double> &Xv) {

  // Scale rows of Xv by Cvxi without creating DiagonalMatrix
  Eigen::SparseMatrix<double> WX = Xv;

  for (int k = 0; k < WX.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(WX, k); it; ++it) {
      it.valueRef() *= Cvxi[it.row()];
    }
  }

  // Compute Xt * (W * X)
  Eigen::SparseMatrix<double> Xt_W_Xv = Xv.transpose() * WX;

  // Result = Xt_W_Xv
  Eigen::SparseMatrix<double> result = Xt_W_Xv;

  return result;
}




// [[Rcpp::export]]
Eigen::SparseMatrix<double> XX_func(
    const Eigen::SparseMatrix<double> &Xv
) {
  Eigen::SparseMatrix<double> XtX =
    Xv.transpose() * Xv;

  return XtX;
}
