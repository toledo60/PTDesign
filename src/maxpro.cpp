#include <Rcpp.h>
#include <random>
#include <limits>
#include <cstring>
#include <numeric>
#include "utils.h"

using namespace std;
using namespace Rcpp;


double psi(const double* D, int n, int k){
  double summation = 0.0;
  for(int i = 0 ; i < n - 1 ; i++ ){
    for(int j = i + 1 ; j < n ; j++){
      double prod = 1.0;
      for(int l = 0; l < k; l++){
        double diff = (D[l * n + i] - D[l * n + j]);
        prod *= diff * diff;
      }
      // Optional safety: if prod hits 0, psi is infinite
      if (prod <= 0.0) return std::numeric_limits<double>::infinity();
      summation += 1.0 / prod;
    }
  }
  return std::pow(2.0/(n * (n - 1)) * summation, 1.0 / (double)k);
}

// [[Rcpp::export]]
List MaxProLHD(int n, int k, int M, int Nmax, int Nswap, int tolerance, NumericVector temp_sched){

  mt19937_64 RNG_engine{random_device{}()};
  uniform_int_distribution<int> distSwap(0, M-2);
  uniform_int_distribution<int> distCol(0, k-1);
  uniform_int_distribution<int> distRow(0, n-1);
  uniform_real_distribution<> dist(0.0, 1.0);

  vector<int> klist(n);
  iota(klist.begin(), klist.end(), 1);

  vector<double> D0(size_t(n) * k);

  for(int j = 0; j < k; j++){
    for(int i = n - 1; i >= 0; --i){
      uniform_int_distribution<int> shuffler(0, i);
      swap(klist[i], klist[shuffler(RNG_engine)]);
      D0[size_t(j) * n + i] = (klist[i] - 0.5) / n;
    }
  }

  double psi0 = psi(D0.data(), n, k);

  vector<double> design_sched;
  design_sched.resize(size_t(M) * n * k);

  vector<double> psi_sched;
  psi_sched.resize(M);

  for(int ii = 0; ii < M; ii++){
    double* predestination = design_sched.data() + size_t(ii) * n * k;
    std::memcpy(predestination, D0.data(), sizeof(double) * n * k);
    psi_sched[ii] = psi0;
  }

  double best_psi = psi0;
  vector<double> best_design(size_t(n) * k);
  std::memcpy(best_design.data(), D0.data(), sizeof(double) * n * k);

  int since_best = 0;
  int iters_done = Nmax;

  for(int iter = 0; iter < Nmax; iter++){

    bool improved_best_this_iter = false;

    ///////////// Metropolis within each replica ///////////////////////
    for(int ii = 0; ii < M; ii++){

      int h = distCol(RNG_engine);
      int w = distRow(RNG_engine);
      int v;
      do { v = distRow(RNG_engine); } while (v == w);

      double psi_original = psi_sched[ii];
      const int block = (ii * n * k);

      //////////////////// Fast psi update ///////////////////////////////
      double summation = 0.0;

      for(int i = 0; i < n; i++){
        if(!(i == w || i == v)){
          double prod1 = 1.0;
          double prod2 = 1.0;

          for(int l = 0; l < k ; l++){
            double variable = design_sched[block + l * n + i];
            double diff1 = (design_sched[block + l * n + w] - variable);
            double diff2 = (design_sched[block + l * n + v] - variable);

            prod1 *= diff1 * diff1;
            prod2 *= diff2 * diff2;
          }

          if (prod1 <= 0.0 || prod2 <= 0.0) {
            summation = std::numeric_limits<double>::infinity();
            break;
          }

          double affected_col = design_sched[block + h * n + i];
          double delta_v = design_sched[block + h * n + v] - affected_col;
          double delta_w = design_sched[block + h * n + w] - affected_col;

          // Avoid division by zero if delta_w == 0
          if (delta_w == 0.0) {
            summation = std::numeric_limits<double>::infinity();
            break;
          }

          double replacement = (delta_v * delta_v) / (delta_w * delta_w);

          summation += ( 1.0/(prod1 * replacement) - 1.0/prod1
                       + 1.0/(prod2 * (1.0/replacement)) - 1.0/prod2 );
        }
      }

      double psi_try;
      if (!std::isfinite(summation)) {
        psi_try = std::numeric_limits<double>::infinity();
      } else {
        double base = (n * (n - 1) * 0.5) * std::pow(psi_original, (double)k) + summation;
        double scaled = 2.0 / (n * (n - 1)) * base;
        psi_try = (scaled > 0.0) ? std::pow(scaled, 1.0 / (double)k)
                                : std::numeric_limits<double>::infinity();
      }
      /////////////////// End fast psi update ////////////////////////////////

      bool accept = false;
      if (psi_try < psi_original) {
        accept = true;
      } else {
        double a = std::exp(-(psi_try - psi_original) / temp_sched[ii]);
        if (dist(RNG_engine) <= a) accept = true;
      }

      if (accept) {
        // apply swap in the chosen column
        std::swap(design_sched[block + h * n + w], design_sched[block + h * n + v]);
        psi_sched[ii] = psi_try;

        // BEST-EVER check (across all replicas)
        if (psi_try < best_psi) {
          best_psi = psi_try;
          improved_best_this_iter = true;

          const double* src = design_sched.data() + size_t(ii) * n * k;
          std::memcpy(best_design.data(), src, sizeof(double) * n * k);
        }
      }
    }
    ////////////// Swap temperatures //////////////////
    if((iter + 1) % Nswap == 0){
      int q = distSwap(RNG_engine);
      double psi1 = psi_sched[q];
      double psi2 = psi_sched[q+1];
      double t1 = temp_sched[q];
      double t2 = temp_sched[q+1];

      if(dist(RNG_engine) <= std::exp((psi2 - psi1) * (1.0/t2 - 1.0/t1))){
        std::swap(temp_sched[q], temp_sched[q+1]);
      }
    }
    //////////////////////////////////////////////////////////////////////////////

    if (improved_best_this_iter) since_best = 0;
    else ++since_best;

    if (iter > 5000 && since_best > tolerance) {
      iters_done = iter + 1;
      break;
    }
  }

  return List::create(
    Named("design")  = initialize_numMat_from_vec(best_design, n, k),
    Named("measure") = best_psi,
    Named("t0")      = temp_sched,
    Named("ntotal")  = iters_done
  );
}