#include "CG_solver.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace {

std::vector<double> multiply(const SparseMatrix& A, const std::vector<double>& x) {
    std::vector<double> result(A.size(), 0.0);
    for (int i = 0; i < static_cast<int>(A.size()); ++i) {
        double sum = 0.0;
        for (int j = 0; j < static_cast<int>(A[i].indexes.size()); ++j)
            sum += A[i].values[j] * x[A[i].indexes[j]];
        result[i] = sum;
    }
    return result;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double sum = 0.0;
    for (int i = 0; i < static_cast<int>(a.size()); ++i)
        sum += a[i] * b[i];
    return sum;
}


struct IC0 {
    std::vector<int> rp, col;
    std::vector<double> val;

    void build(const SparseMatrix& A) {
        const int n = static_cast<int>(A.size());
        rp.assign(n + 1, 0); col.clear(); val.clear();
        std::vector<std::pair<int,double>> row;
        for (int i = 0; i < n; ++i) {
            row.clear();
            for (std::size_t k = 0; k < A[i].indexes.size(); ++k)
                if (A[i].indexes[k] <= i)
                row.emplace_back(A[i].indexes[k], A[i].values[k]);
            std::sort(row.begin(), row.end());
            for (auto& e : row)
            {
            	col.push_back(e.first);
            	val.push_back(e.second);
            }
            rp[i + 1] = static_cast<int>(col.size());
        }

        for (int i = 0; i < n; ++i) {
            const int bi = rp[i], ei = rp[i + 1];
            for (int p = bi; p < ei; ++p) {
                const int j = col[p];

                double s = val[p];
                int pi = bi, pj = rp[j];
                const int ej = rp[j + 1] - 1;
                while (pi < p && pj < ej) {
                    if (col[pi] == col[pj]) { s -= val[pi] * val[pj]; ++pi; ++pj; }
                    else if (col[pi] < col[pj]) ++pi;
                    else ++pj;
                }
                if (j < i) {
                    val[p] = s / val[rp[j + 1] - 1];
                } else {
                    if (s <= 1e-14 * std::abs(val[p]) || s <= 0.0)
                        s = std::abs(val[p]);
                    val[p] = std::sqrt(s);
                }
            }
        }
    }

    // z = (L L^T)^{-1} r
    void apply(const std::vector<double>& r, std::vector<double>& z) const {
        const int n = static_cast<int>(r.size());
        z = r;
        for (int i = 0; i < n; ++i) {                   // L y = r
            double s = z[i];
            const int d = rp[i + 1] - 1;
            for (int p = rp[i]; p < d; ++p) s -= val[p] * z[col[p]];
            z[i] = s / val[d];
        }
        for (int i = n - 1; i >= 0; --i) {              // L^T z = y
            const int d = rp[i + 1] - 1;
            z[i] /= val[d];
            for (int p = rp[i]; p < d; ++p) z[col[p]] -= val[p] * z[i];
        }
    }
};

}

std::vector<double> CgSolver::solve(const SparseMatrix& A,
                                     const std::vector<double>& b) {
    return solve(A, b, std::vector<double>(A.size(), 0.0));
}

std::vector<double> CgSolver::solve(const SparseMatrix& A,
                                     const std::vector<double>& b,
                                     std::vector<double> x) {
    if (x.size() != A.size())
        x.assign(A.size(), 0.0);

    errTol = constants::calcErrorTolerance(A.size());

	solveLS(A, b, x);

    return x;
}

bool CgSolver::solveLS(const SparseMatrix& A,
					   const std::vector<double>& b,
					   std::vector<double>& x) {
	const int n = static_cast<int>(A.size());
	if (x.size() != A.size()) x.assign(n, 0.0);

	IC0 M;
	M.build(A);

	std::vector<double> r = b;
	{
		std::vector<double> Ax = multiply(A, x);
		for (int i = 0; i < n; ++i) r[i] -= Ax[i];
	}
	const double r0 = std::sqrt(dot(r, r));
	if (r0 == 0.0) return true;

	std::vector<double> z(n), p(n), Ap(n);
	M.apply(r, z);
	p = z;
	double rz = dot(r, z);

	for (int k = 1; k <= maxIterCount; ++k) {
		Ap = multiply(A, p);
		const double pAp = dot(p, Ap);

		const double alpha = rz / pAp;
		for (int i = 0; i < n; ++i) { x[i] += alpha * p[i]; r[i] -= alpha * Ap[i]; }

		if (std::sqrt(dot(r, r)) < relResTol * r0) {
			return true;
		}
		M.apply(r, z);
		const double rz_new = dot(r, z);
		const double beta = rz_new / rz;
		for (int i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
		rz = rz_new;
	}
	return false;
}

