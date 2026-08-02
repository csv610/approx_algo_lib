// =====================================================================
// Chapter 27: Shortest Vector in a Lattice
// Vazirani, "Approximation Algorithms" (Springer 2001)
//
// Implements:
//   - Euclid's algorithm for GCD (1D shortest vector)
//   - Gauss' algorithm for 2D lattice shortest vector
//   - Gram-Schmidt orthogonalization (nD)
//   - LLL basis reduction
//
// Usage: echo '{"basis":[[...]],"dimension":...}' | ./shortest_vector
// =====================================================================

#include "json.hpp"
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <string>
#include <sstream>
#include <iomanip>

using json = nlohmann::json;
using namespace std;

namespace aal {

long long gcd(long long a, long long b) {
    a = llabs(a);
    b = llabs(b);
    while (b) {
        a %= b;
        swap(a, b);
    }
    return a;
}

long long shortest_vector_1d(const vector<long long>& basis) {
    long long g = 0;
    for (long long x : basis) g = gcd(g, x);
    return g;
}

double dot(const vector<double>& a, const vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

double norm_sq(const vector<double>& v) {
    return dot(v, v);
}

vector<double> sub(const vector<double>& a, const vector<double>& b) {
    vector<double> r(a.size());
    for (size_t i = 0; i < a.size(); ++i) r[i] = a[i] - b[i];
    return r;
}

vector<double> scale(double c, const vector<double>& v) {
    vector<double> r(v.size());
    for (size_t i = 0; i < v.size(); ++i) r[i] = c * v[i];
    return r;
}

pair<vector<vector<double>>, vector<vector<double>>> gram_schmidt(
    const vector<vector<double>>& basis
) {
    int n = static_cast<int>(basis.size());
    int d = static_cast<int>(basis[0].size());

    vector<vector<double>> star(n, vector<double>(d, 0.0));
    vector<vector<double>> mu(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        star[i] = basis[i];
        for (int j = 0; j < i; ++j) {
            double sj_sq = norm_sq(star[j]);
            if (sj_sq < 1e-15) continue;
            mu[i][j] = dot(basis[i], star[j]) / sj_sq;
            star[i] = sub(star[i], scale(mu[i][j], star[j]));
        }
    }

    return {star, mu};
}

pair<vector<long long>, double> gauss_2d(vector<vector<long long>> basis) {
    auto& b1 = basis[0];
    auto& b2 = basis[1];

    auto v1 = vector<double>(b1.begin(), b1.end());
    auto v2 = vector<double>(b2.begin(), b2.end());
    if (norm_sq(v1) > norm_sq(v2)) swap(b1, b2);

    while (true) {
        v1 = vector<double>(b1.begin(), b1.end());
        v2 = vector<double>(b2.begin(), b2.end());
        if (norm_sq(v2) >= norm_sq(v1)) break;

        swap(b1, b2);
        v1 = vector<double>(b1.begin(), b1.end());
        v2 = vector<double>(b2.begin(), b2.end());

        long long r = static_cast<long long>(round(dot(v1, v2) / norm_sq(v2)));
        for (size_t k = 0; k < b1.size(); ++k)
            b1[k] -= r * b2[k];
    }

    v1 = vector<double>(b1.begin(), b1.end());
    return {b1, sqrt(norm_sq(v1))};
}

vector<vector<double>> lll_reduce(vector<vector<double>> basis, double delta) {
    int n = static_cast<int>(basis.size());
    int d = static_cast<int>(basis[0].size());

    auto [star, mu] = gram_schmidt(basis);
    vector<double> b_star_sq(n);
    for (int i = 0; i < n; ++i) b_star_sq[i] = norm_sq(star[i]);

    int k = 1;
    while (k < n) {
        for (int j = k - 1; j >= 0; --j) {
            if (fabs(mu[k][j]) > 0.5) {
                long long r = static_cast<long long>(round(mu[k][j]));
                for (int t = 0; t < d; ++t)
                    basis[k][t] -= r * basis[j][t];
                auto [new_star, new_mu] = gram_schmidt(basis);
                star = new_star;
                mu = new_mu;
                for (int i = 0; i < n; ++i) b_star_sq[i] = norm_sq(star[i]);
            }
        }

        if (b_star_sq[k] >= (delta - mu[k][k - 1] * mu[k][k - 1]) * b_star_sq[k - 1]) {
            ++k;
        } else {
            swap(basis[k], basis[k - 1]);
            auto [new_star, new_mu] = gram_schmidt(basis);
            star = new_star;
            mu = new_mu;
            for (int i = 0; i < n; ++i) b_star_sq[i] = norm_sq(star[i]);
            k = max(k - 1, 1);
        }
    }

    return basis;
}

double orthogonality_defect(
    const vector<vector<double>>& basis,
    const vector<vector<double>>& gs_star
) {
    double basis_prod = 1.0;
    for (auto& b : basis) basis_prod *= sqrt(norm_sq(b));

    double gs_prod = 1.0;
    for (auto& s : gs_star) {
        double ns = sqrt(norm_sq(s));
        if (ns > 1e-15) gs_prod *= ns;
    }

    return (gs_prod < 1e-15) ? 0.0 : basis_prod / gs_prod;
}

json solve_shortest_vector(const json& input) {
    int dim = input["dimension"];
    auto basis_raw = input["basis"];

    json result;
    result["algorithm"] = "shortest_vector";
    result["dimension"] = dim;

    if (dim == 1) {
        vector<long long> b1d;
        for (auto& v : basis_raw) b1d.push_back(v[0]);
        result["gcd_1d"] = shortest_vector_1d(b1d);
    }

    if (dim == 2) {
        vector<vector<long long>> basis(2, vector<long long>(2));
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                basis[i][j] = basis_raw[i][j];

        auto [sv, norm] = gauss_2d(basis);
        json sv_arr = json::array();
        for (auto x : sv) sv_arr.push_back(static_cast<int>(x));
        result["shortest_vector_2d"] = sv_arr;
        result["shortest_norm_2d"] = norm;
    }

    vector<vector<double>> basis_d(dim, vector<double>(dim));
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            basis_d[i][j] = static_cast<double>(basis_raw[i][j]);

    auto [gs_star, mu] = gram_schmidt(basis_d);
    result["gram_schmidt"] = gs_star;
    result["orthogonality_defect"] = orthogonality_defect(basis_d, gs_star);

    auto reduced = lll_reduce(basis_d);
    result["reduced_basis"] = reduced;

    return result;
}

} // namespace aal

void solve(istream& in, ostream& out) {
    json input = json::parse(in);
    json result = aal::solve_shortest_vector(input);
    out << setw(2) << result << endl;
}

int main() {
    solve(cin, cout);
    return 0;
}
