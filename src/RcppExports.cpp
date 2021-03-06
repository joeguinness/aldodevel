// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// getLinvEntries
NumericMatrix getLinvEntries(NumericVector covparms, StringVector covfun_name, NumericMatrix locs, IntegerMatrix NNarray);
RcppExport SEXP _aldodevel_getLinvEntries(SEXP covparmsSEXP, SEXP covfun_nameSEXP, SEXP locsSEXP, SEXP NNarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type covfun_name(covfun_nameSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(getLinvEntries(covparms, covfun_name, locs, NNarray));
    return rcpp_result_gen;
END_RCPP
}
// MaternFun
Rcpp::NumericMatrix MaternFun(Rcpp::NumericMatrix& distmat, Rcpp::NumericVector covparms);
RcppExport SEXP _aldodevel_MaternFun(SEXP distmatSEXP, SEXP covparmsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type distmat(distmatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type covparms(covparmsSEXP);
    rcpp_result_gen = Rcpp::wrap(MaternFun(distmat, covparms));
    return rcpp_result_gen;
END_RCPP
}
// OrderedCompLik
NumericVector OrderedCompLik(NumericVector covparms, StringVector covfun_name, NumericVector y, NumericMatrix locs, IntegerMatrix NNarray);
RcppExport SEXP _aldodevel_OrderedCompLik(SEXP covparmsSEXP, SEXP covfun_nameSEXP, SEXP ySEXP, SEXP locsSEXP, SEXP NNarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type covfun_name(covfun_nameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(OrderedCompLik(covparms, covfun_name, y, locs, NNarray));
    return rcpp_result_gen;
END_RCPP
}
// OrderedGroupCompLik
NumericVector OrderedGroupCompLik(NumericVector covparms, NumericVector y, NumericMatrix locs, List NNlist);
RcppExport SEXP _aldodevel_OrderedGroupCompLik(SEXP covparmsSEXP, SEXP ySEXP, SEXP locsSEXP, SEXP NNlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< List >::type NNlist(NNlistSEXP);
    rcpp_result_gen = Rcpp::wrap(OrderedGroupCompLik(covparms, y, locs, NNlist));
    return rcpp_result_gen;
END_RCPP
}
// vecchiaLik_function
NumericVector vecchiaLik_function(NumericVector covparms, StringVector covfun_name, NumericVector y, NumericMatrix locs, IntegerMatrix NNarray);
RcppExport SEXP _aldodevel_vecchiaLik_function(SEXP covparmsSEXP, SEXP covfun_nameSEXP, SEXP ySEXP, SEXP locsSEXP, SEXP NNarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type covfun_name(covfun_nameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(vecchiaLik_function(covparms, covfun_name, y, locs, NNarray));
    return rcpp_result_gen;
END_RCPP
}
// vecchiaLoglik
NumericVector vecchiaLoglik(NumericVector covparms, StringVector covfun_name, NumericVector y, NumericMatrix locs, IntegerMatrix NNarray);
RcppExport SEXP _aldodevel_vecchiaLoglik(SEXP covparmsSEXP, SEXP covfun_nameSEXP, SEXP ySEXP, SEXP locsSEXP, SEXP NNarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type covfun_name(covfun_nameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(vecchiaLoglik(covparms, covfun_name, y, locs, NNarray));
    return rcpp_result_gen;
END_RCPP
}
// vecchiaLinverse
NumericVector vecchiaLinverse(NumericVector covparms, StringVector covfun_name, NumericMatrix locs, IntegerMatrix NNarray);
RcppExport SEXP _aldodevel_vecchiaLinverse(SEXP covparmsSEXP, SEXP covfun_nameSEXP, SEXP locsSEXP, SEXP NNarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type covfun_name(covfun_nameSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(vecchiaLinverse(covparms, covfun_name, locs, NNarray));
    return rcpp_result_gen;
END_RCPP
}
// LinvMultFromEntries
NumericVector LinvMultFromEntries(NumericMatrix LinvEntries, NumericVector z, IntegerMatrix NNarray);
RcppExport SEXP _aldodevel_LinvMultFromEntries(SEXP LinvEntriesSEXP, SEXP zSEXP, SEXP NNarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type LinvEntries(LinvEntriesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(LinvMultFromEntries(LinvEntries, z, NNarray));
    return rcpp_result_gen;
END_RCPP
}
// LMultFromEntries
NumericVector LMultFromEntries(NumericMatrix LinvEntries, NumericVector z, IntegerMatrix NNarray);
RcppExport SEXP _aldodevel_LMultFromEntries(SEXP LinvEntriesSEXP, SEXP zSEXP, SEXP NNarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type LinvEntries(LinvEntriesSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type z(zSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(LMultFromEntries(LinvEntries, z, NNarray));
    return rcpp_result_gen;
END_RCPP
}
// vecchiaLik_grouped_function
NumericVector vecchiaLik_grouped_function(NumericVector covparms, StringVector covfun_name, NumericVector y, NumericMatrix locs, List NNlist);
RcppExport SEXP _aldodevel_vecchiaLik_grouped_function(SEXP covparmsSEXP, SEXP covfun_nameSEXP, SEXP ySEXP, SEXP locsSEXP, SEXP NNlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type covfun_name(covfun_nameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< List >::type NNlist(NNlistSEXP);
    rcpp_result_gen = Rcpp::wrap(vecchiaLik_grouped_function(covparms, covfun_name, y, locs, NNlist));
    return rcpp_result_gen;
END_RCPP
}
// vecchiaLik_grouped
NumericVector vecchiaLik_grouped(NumericVector covparms, NumericVector y, NumericMatrix locs, List NNlist);
RcppExport SEXP _aldodevel_vecchiaLik_grouped(SEXP covparmsSEXP, SEXP ySEXP, SEXP locsSEXP, SEXP NNlistSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< List >::type NNlist(NNlistSEXP);
    rcpp_result_gen = Rcpp::wrap(vecchiaLik_grouped(covparms, y, locs, NNlist));
    return rcpp_result_gen;
END_RCPP
}
// vecchiaLik
NumericVector vecchiaLik(NumericVector covparms, StringVector covfun_name, NumericVector y, NumericMatrix locs, IntegerMatrix NNarray);
RcppExport SEXP _aldodevel_vecchiaLik(SEXP covparmsSEXP, SEXP covfun_nameSEXP, SEXP ySEXP, SEXP locsSEXP, SEXP NNarraySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type covparms(covparmsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type covfun_name(covfun_nameSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type NNarray(NNarraySEXP);
    rcpp_result_gen = Rcpp::wrap(vecchiaLik(covparms, covfun_name, y, locs, NNarray));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_aldodevel_getLinvEntries", (DL_FUNC) &_aldodevel_getLinvEntries, 4},
    {"_aldodevel_MaternFun", (DL_FUNC) &_aldodevel_MaternFun, 2},
    {"_aldodevel_OrderedCompLik", (DL_FUNC) &_aldodevel_OrderedCompLik, 5},
    {"_aldodevel_OrderedGroupCompLik", (DL_FUNC) &_aldodevel_OrderedGroupCompLik, 4},
    {"_aldodevel_vecchiaLik_function", (DL_FUNC) &_aldodevel_vecchiaLik_function, 5},
    {"_aldodevel_vecchiaLoglik", (DL_FUNC) &_aldodevel_vecchiaLoglik, 5},
    {"_aldodevel_vecchiaLinverse", (DL_FUNC) &_aldodevel_vecchiaLinverse, 4},
    {"_aldodevel_LinvMultFromEntries", (DL_FUNC) &_aldodevel_LinvMultFromEntries, 3},
    {"_aldodevel_LMultFromEntries", (DL_FUNC) &_aldodevel_LMultFromEntries, 3},
    {"_aldodevel_vecchiaLik_grouped_function", (DL_FUNC) &_aldodevel_vecchiaLik_grouped_function, 5},
    {"_aldodevel_vecchiaLik_grouped", (DL_FUNC) &_aldodevel_vecchiaLik_grouped, 4},
    {"_aldodevel_vecchiaLik", (DL_FUNC) &_aldodevel_vecchiaLik, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_aldodevel(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
