#include "iter_ls_solver.hpp"


double IterLsSolver::calcBackwardError(const SparseMatrix& A, const std::vector<double>& b, 
    const std::vector<double>& x
) {
    std::vector<double> residual = calcResidual(A, b, x);
    double A_norm = 0;
    for (const auto& row : A) {
        double s = 0.0;
        for (double v : row.values) s += std::abs(v);
        A_norm = std::max(A_norm, s);
    }
    double b_norm = 0, x_norm = 0, r_norm = 0;
    for (double v : b) b_norm = std::max(b_norm, std::abs(v));
    for (double v : x) x_norm = std::max(x_norm, std::abs(v));
    for (double v : residual) r_norm = std::max(r_norm, std::abs(v));
    double scale = A_norm * x_norm + b_norm;

    //только в случае если A = 0 и b = 0 (вроде как не возможно)
    if (scale == 0.0) return 0.0;
    return r_norm / scale;
}


std::vector<double> IterLsSolver::calcResidual(const SparseMatrix& A, const std::vector<double>& b, 
    const std::vector<double>& x
) {
    std::vector<double> r;
    r.reserve(A.size());
    long double Ax_i;
    for(int i = 0; i < A.size(); i++) {
        Ax_i = 0;
        for(int j = 0; j < A[i].indexes.size(); j++)
            Ax_i += (long double)A[i].values[j] * (long double)x[A[i].indexes[j]];
        r.push_back((long double) b[i] - Ax_i);
    }
    return r;
}