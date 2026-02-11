#include "utils.h"

using std::vector;
using Rcpp::NumericMatrix;

NumericMatrix initialize_numMat(const vector<double>& sched, int n, int k, int index) {
    NumericMatrix result(n, k);
    const double* source_block = sched.data() + std::size_t(index) * n * k;
    std::memcpy( REAL(result),source_block, sizeof(double) * n * k );
    return result;
}