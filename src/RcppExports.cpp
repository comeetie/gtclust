// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ldirpost
double ldirpost(const NumericVector& beta, int n, const NumericVector& lambda, const NumericVector& lpi);
RcppExport SEXP _gtclust_ldirpost(SEXP betaSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    rcpp_result_gen = Rcpp::wrap(ldirpost(beta, n, lambda, lpi));
    return rcpp_result_gen;
END_RCPP
}
// grad_ldirpost
NumericVector grad_ldirpost(const NumericVector& beta, int n, const NumericVector& lambda, const NumericVector& lpi);
RcppExport SEXP _gtclust_grad_ldirpost(SEXP betaSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_ldirpost(beta, n, lambda, lpi));
    return rcpp_result_gen;
END_RCPP
}
// dH_ldirpost
NumericVector dH_ldirpost(const NumericVector& beta, int n, const NumericVector& lambda, const NumericVector& lpi);
RcppExport SEXP _gtclust_dH_ldirpost(SEXP betaSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    rcpp_result_gen = Rcpp::wrap(dH_ldirpost(beta, n, lambda, lpi));
    return rcpp_result_gen;
END_RCPP
}
// z_ldirpost
double z_ldirpost(const NumericVector& beta, int n, const NumericVector& lambda, const NumericVector& lpi);
RcppExport SEXP _gtclust_z_ldirpost(SEXP betaSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    rcpp_result_gen = Rcpp::wrap(z_ldirpost(beta, n, lambda, lpi));
    return rcpp_result_gen;
END_RCPP
}
// logdetH_ldirpost
double logdetH_ldirpost(const NumericVector& beta, int n, const NumericVector& lambda, const NumericVector& lpi);
RcppExport SEXP _gtclust_logdetH_ldirpost(SEXP betaSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    rcpp_result_gen = Rcpp::wrap(logdetH_ldirpost(beta, n, lambda, lpi));
    return rcpp_result_gen;
END_RCPP
}
// iH_ldirpost
NumericMatrix iH_ldirpost(const NumericVector& beta, int n, const NumericVector& lambda, const NumericVector& lpi);
RcppExport SEXP _gtclust_iH_ldirpost(SEXP betaSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    rcpp_result_gen = Rcpp::wrap(iH_ldirpost(beta, n, lambda, lpi));
    return rcpp_result_gen;
END_RCPP
}
// H_ldirpost
NumericMatrix H_ldirpost(const NumericVector& beta, int n, const NumericVector& lambda, const NumericVector& lpi);
RcppExport SEXP _gtclust_H_ldirpost(SEXP betaSEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    rcpp_result_gen = Rcpp::wrap(H_ldirpost(beta, n, lambda, lpi));
    return rcpp_result_gen;
END_RCPP
}
// ldirpost_norm
List ldirpost_norm(const NumericVector beta0, int n, const NumericVector& lambda, const NumericVector& lpi);
RcppExport SEXP _gtclust_ldirpost_norm(SEXP beta0SEXP, SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    rcpp_result_gen = Rcpp::wrap(ldirpost_norm(beta0, n, lambda, lpi));
    return rcpp_result_gen;
END_RCPP
}
// dirichlet_evidence
double dirichlet_evidence(int n, const NumericVector& lambda, const NumericVector& lpi, const NumericVector& pi);
RcppExport SEXP _gtclust_dirichlet_evidence(SEXP nSEXP, SEXP lambdaSEXP, SEXP lpiSEXP, SEXP piSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lpi(lpiSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type pi(piSEXP);
    rcpp_result_gen = Rcpp::wrap(dirichlet_evidence(n, lambda, lpi, pi));
    return rcpp_result_gen;
END_RCPP
}
// sparseBdiag
Eigen::SparseMatrix<double> sparseBdiag(Rcpp::List B_list);
RcppExport SEXP _gtclust_sparseBdiag(SEXP B_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type B_list(B_listSEXP);
    rcpp_result_gen = Rcpp::wrap(sparseBdiag(B_list));
    return rcpp_result_gen;
END_RCPP
}
// remove1row1col
Eigen::SparseMatrix<double> remove1row1col(Eigen::SparseMatrix<double> L);
RcppExport SEXP _gtclust_remove1row1col(SEXP LSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type L(LSEXP);
    rcpp_result_gen = Rcpp::wrap(remove1row1col(L));
    return rcpp_result_gen;
END_RCPP
}
// buildLaplacianR
Eigen::SparseMatrix<double> buildLaplacianR(Eigen::SparseMatrix<double> Lg, Eigen::SparseMatrix<double> Lh, NumericMatrix cutset, NumericVector permutation);
RcppExport SEXP _gtclust_buildLaplacianR(SEXP LgSEXP, SEXP LhSEXP, SEXP cutsetSEXP, SEXP permutationSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type Lg(LgSEXP);
    Rcpp::traits::input_parameter< Eigen::SparseMatrix<double> >::type Lh(LhSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type cutset(cutsetSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type permutation(permutationSEXP);
    rcpp_result_gen = Rcpp::wrap(buildLaplacianR(Lg, Lh, cutset, permutation));
    return rcpp_result_gen;
END_RCPP
}
// hclustcc_cpp
List hclustcc_cpp(const List nb, const NumericMatrix& X, List method_obj);
RcppExport SEXP _gtclust_hclustcc_cpp(SEXP nbSEXP, SEXP XSEXP, SEXP method_objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type method_obj(method_objSEXP);
    rcpp_result_gen = Rcpp::wrap(hclustcc_cpp(nb, X, method_obj));
    return rcpp_result_gen;
END_RCPP
}
// bayesian_hclustcc_cpp
List bayesian_hclustcc_cpp(const List nb, const NumericMatrix& X, List method_obj);
RcppExport SEXP _gtclust_bayesian_hclustcc_cpp(SEXP nbSEXP, SEXP XSEXP, SEXP method_objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type nb(nbSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type method_obj(method_objSEXP);
    rcpp_result_gen = Rcpp::wrap(bayesian_hclustcc_cpp(nb, X, method_obj));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gtclust_ldirpost", (DL_FUNC) &_gtclust_ldirpost, 4},
    {"_gtclust_grad_ldirpost", (DL_FUNC) &_gtclust_grad_ldirpost, 4},
    {"_gtclust_dH_ldirpost", (DL_FUNC) &_gtclust_dH_ldirpost, 4},
    {"_gtclust_z_ldirpost", (DL_FUNC) &_gtclust_z_ldirpost, 4},
    {"_gtclust_logdetH_ldirpost", (DL_FUNC) &_gtclust_logdetH_ldirpost, 4},
    {"_gtclust_iH_ldirpost", (DL_FUNC) &_gtclust_iH_ldirpost, 4},
    {"_gtclust_H_ldirpost", (DL_FUNC) &_gtclust_H_ldirpost, 4},
    {"_gtclust_ldirpost_norm", (DL_FUNC) &_gtclust_ldirpost_norm, 4},
    {"_gtclust_dirichlet_evidence", (DL_FUNC) &_gtclust_dirichlet_evidence, 4},
    {"_gtclust_sparseBdiag", (DL_FUNC) &_gtclust_sparseBdiag, 1},
    {"_gtclust_remove1row1col", (DL_FUNC) &_gtclust_remove1row1col, 1},
    {"_gtclust_buildLaplacianR", (DL_FUNC) &_gtclust_buildLaplacianR, 4},
    {"_gtclust_hclustcc_cpp", (DL_FUNC) &_gtclust_hclustcc_cpp, 3},
    {"_gtclust_bayesian_hclustcc_cpp", (DL_FUNC) &_gtclust_bayesian_hclustcc_cpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_gtclust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
