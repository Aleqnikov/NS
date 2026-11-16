#pragma once
#include "iter_ls_solver.hpp"
#include "util/solver_constants.hpp"
#include <vector>

class CgSolver : public IterLsSolver {
public:
	CgSolver(
		double absDxTol = constants::ABSOLUTE_DX_TOLERANCE,
		double relDxTol = constants::RELATIVE_DX_TOLERANCE,
		int maxIterCount = constants::MAX_ITERATION_COUNT,
		int errCheckInter = constants::ERROR_CHECK_ITER_INTERVAL) :
		LsSolverI(),
		IterLsSolver(maxIterCount),
		absDxTol(absDxTol),
		relDxTol(relDxTol),
		errCheckInter(errCheckInter) {}

	std::vector<double> solve(const SparseMatrix& A, const std::vector<double>& b) override;

	bool solveLS(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x);

	std::vector<double> solve(const SparseMatrix& A,
						  const std::vector<double>& b,
						  std::vector<double> x0);
protected:
	const double absDxTol;
	const double relDxTol;

	double errTol = 0.0;
	const int errCheckInter;
	double relResTol = 1e-6;
};