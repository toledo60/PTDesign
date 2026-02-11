#pragma once

#include <Rcpp.h>
#include <vector>

Rcpp::NumericMatrix initialize_numMat(const std::vector<double>& sched,
                                      int n, int k, int index);