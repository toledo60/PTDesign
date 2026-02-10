#include <Rcpp.h>
#include <random>
#include <cmath>
#include <limits>
#include <vector>

using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix random_oofa(int n, int k) {
    mt19937_64 RNG_engine{random_device{}()};
    vector<int> klist(k);

    iota(klist.begin(), klist.end(), 1);     // fill 1:k

    // Column major matrix: element i,j (row/column) are accessed as D0[j * n + i]
    double D0[n * k];

    for (int i = 0; i < n; i++) {
        for (int j = k - 1; j >= 0; --j) {
            // shuffle the integers
            uniform_int_distribution<int> shuffler(0, j);
            swap(klist[j], klist[shuffler(RNG_engine)]);
            D0[j * n + i] = klist[j];
        }
    }
    NumericMatrix result(n, k);
    memcpy(REAL(result), D0, sizeof(double) * n * k);
    return result;
}

// [[Rcpp::export]]
double phi_p(NumericMatrix x, int n, int k, double p) {	// phip with p=15
	double sum_phip = 0.0;
    double minus_p_halves = -p *0.5;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double dist = 0.0;

            for (int l = 0; l < k; ++l) {
                double diff = x(i, l) - x(j, l);
                dist += diff * diff;
            }
			sum_phip += pow(dist, minus_p_halves);
        }
    }
    return pow((2.0 / (n * (n-1))) * sum_phip, 1.0/p);
}

struct PhiPCache {
    int n, k;
    double p = 15.0;
    double minus_p_halves = -0.5 * 15.0;

    NumericMatrix X;          
    std::vector<double> D2;    // packed upper triangle of pairwise squared distances
    double sum_phip;           // sum_{i<j} D2[i,j]^(-p/2)

    // index into packed storage (i<j)
    inline int idx(int i, int j) const {
        // number of pairs before row i is i*n - i*(i+1)/2
        return i * n - i * (i + 1) / 2 + (j - i - 1);
    }

    PhiPCache(const NumericMatrix& X0) : n(X0.nrow()), k(X0.ncol()), X(clone(X0)) {
        const int m = n * (n - 1) / 2;
        D2.resize(m);
        sum_phip = 0.0;

        for (int i = 0; i < n - 1; ++i) {
            for (int j = i + 1; j < n; ++j) {
                double d2 = 0.0;
                for (int l = 0; l < k; ++l) {
                    double diff = X(i, l) - X(j, l);
                    d2 += diff * diff;
                }
                D2[idx(i, j)] = d2;
                sum_phip += std::pow(d2, minus_p_halves);
            }
        }
    }
    // Current phi_p value
    inline double phi() const {
        double avg = (2.0 / (n * (n - 1))) * sum_phip;
        return std::pow(avg, 1.0 / p);
    }

    // Preview the effect of swapping columns (a,b) in row r (no commit)
    double preview_swap(int r, int a, int b) const {
        double delta_sum = 0.0;
        double xra = X(r, a), xrb = X(r, b);

        for (int j = 0; j < n; ++j) {
            if (j == r) continue;
            int i = std::min(r, j), h = std::max(r, j);
            int id = idx(i, h);

            double xja = X(j, a), xjb = X(j, b);

            // old contribution in the two touched columns
            double cold = (xra - xja) * (xra - xja) + (xrb - xjb) * (xrb - xjb);
            // new after swapping row r's entries in a,b
            double cnew = (xrb - xja) * (xrb - xja) + (xra - xjb) * (xra - xjb);

            double oldD2 = D2[id];
            double newD2 = oldD2 + (cnew - cold);

            delta_sum += std::pow(newD2, minus_p_halves) - std::pow(oldD2, minus_p_halves);
        }

        double avg = (2.0 / (n * (n - 1))) * (sum_phip + delta_sum);
        return std::pow(avg, 1.0 / p);
    }
    // Commit the swap (apply changes to D2, sum_phip, and X)
    void commit_swap(int r, int a, int b) {
        double xra = X(r, a), xrb = X(r, b);
        double delta_sum = 0.0;

        for (int j = 0; j < n; ++j) {
            if (j == r) continue;
            int i = std::min(r, j), h = std::max(r, j);
            int id = idx(i, h);

            double xja = X(j, a), xjb = X(j, b);

            double cold = (xra - xja) * (xra - xja) + (xrb - xjb) * (xrb - xjb);
            double cnew = (xrb - xja) * (xrb - xja) + (xra - xjb) * (xra - xjb);

            double oldD2 = D2[id];
            double newD2 = oldD2 + (cnew - cold);

            D2[id] = newD2;
            delta_sum += std::pow(newD2, minus_p_halves) - std::pow(oldD2, minus_p_halves);
        }
        sum_phip += delta_sum;
        std::swap(X(r, a), X(r, b));
    }
};

// Helper function to swap two elements in a matrix row
void swap_elements(NumericMatrix& matrix, int row, int col1, int col2) {
    std::swap(matrix(row, col1), matrix(row, col2));
}

NumericMatrix permute_levels(const NumericMatrix &X, int k, std::mt19937_64 &rng) {
  std::vector<int> perm(k);
  for (int i = 0; i < k; i++) perm[i] = i + 1;
  std::shuffle(perm.begin(), perm.end(), rng);

  NumericMatrix Xnew(X.nrow(), X.ncol());
  for (int i = 0; i < X.nrow(); i++) {
    for (int j = 0; j < X.ncol(); j++) {
      int val = X(i, j);
      Xnew(i, j) = perm[val - 1];
    }
  }
  return Xnew;
}

// [[Rcpp::export]]
List pt_optim(NumericMatrix design, int M, int Nmax, 
    int Nswap, double alpha, double T0, int tolerance) {
    int n = design.nrow();
    int k = design.ncol();

    std::mt19937_64 RNG_engine{std::random_device{}()};
    std::uniform_int_distribution<int> rand_row(0, n - 1);
    std::uniform_int_distribution<int> rand_col(0, k - 1);
    std::uniform_int_distribution<int> rand_replica(0, M - 2);
    std::uniform_real_distribution<> unif(0.0, 1.0);

    // Initialize replicas (states) via caches
    std::vector<PhiPCache> caches;
    caches.reserve(M);
    for (int i = 0; i < M; ++i) caches.emplace_back(design);

    std::vector<double> temperatures(M);
    std::vector<double> criteria(M);

    if (T0 < 0) {
        double d2 = k/6.0 * n * (n+1);
        double d0 = 1.0/std::sqrt(d2 - k) - 1.0/std::sqrt(d2);
        T0 = -d0 * (1.0/std::log(0.99));
        Rcout << "T0=" << T0 << "\n";
    }

    for (int i = 0; i < M; ++i) {
        temperatures[i] = T0 * std::pow(alpha, i);
        criteria[i] = caches[i].phi(); 
    }

    NumericMatrix design_optim = clone(caches[0].X);
    double global_min_crit = *std::min_element(criteria.begin(), criteria.end());
    int global_min_index = 0;
    for (int i = 0; i < M; ++i) {
        if (criteria[i] == global_min_crit) {
            design_optim = clone(caches[i].X);
            global_min_index = i;
            break;
        }
    }

    int unsuccessful_global_min = 0;

    for (int iter = 0; iter < Nmax; ++iter) {
        bool improvement = false;

        for (int i = 0; i < M; ++i) {
            int row = rand_row(RNG_engine);
            int col1 = rand_col(RNG_engine);
            int col2;
            do { col2 = rand_col(RNG_engine); } while (col1 == col2);

            double new_crit = caches[i].preview_swap(row, col1, col2);

            double dE = new_crit - criteria[i];
            if (dE <= 0.0 || std::log(unif(RNG_engine)) < -dE / temperatures[i]) {
                // Accept: commit the swap and update criterion
                caches[i].commit_swap(row, col1, col2);
                criteria[i] = new_crit;

                if (new_crit < global_min_crit) {
                    global_min_crit = new_crit;
                    global_min_index = i;
                    improvement = true;
                    design_optim = clone(caches[i].X);
                }
            }
            // else: reject (do nothing; no swap to revert because we didn't commit)
        }

        if ((iter + 1) % Nswap == 0) {
            int i = rand_replica(RNG_engine) % (M - 1);
            double delta_E = (criteria[i + 1] - criteria[i]);
            double delta_B = (1.0 / temperatures[i + 1] - 1.0 / temperatures[i]);
            if (unif(RNG_engine) < std::exp(delta_E * delta_B)) {
                std::swap(temperatures[i], temperatures[i + 1]);
            }
        }

        if (improvement) unsuccessful_global_min = 0;
        else ++unsuccessful_global_min;

        if (unsuccessful_global_min > tolerance && iter > 2500) {
            Nmax = iter + 1;
            break;
        }
    }
    return List::create(
        Named("design") = design_optim,
        Named("measure") = global_min_crit,
        Named("t0") = wrap(temperatures),
        Named("ntotal") = Nmax
    );
}