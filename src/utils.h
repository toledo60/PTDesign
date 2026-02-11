#pragma once

#include <Rcpp.h>
#include <vector>

Rcpp::NumericMatrix initialize_numMat(const std::vector<double>& sched,
                                      int n, int k, int index);

Rcpp::NumericMatrix initialize_numMat_from_vec(const std::vector<double>& block, int n, int k);