#include "parallel_GS_solver.hpp"
#include <stdexcept>
#include "util/parallel_for.hpp"
#include <functional>

matrix ParallelGSSolver::solve(const std::vector<matrix>& As, const matrix& bs) {
   isInputCorrect(As, bs);

    matrix Xs = matrix(As.size());
    notСonvergedMatrixNumbers = std::vector<int>(As.size(), 0);

    parallel::run(
        threadCount,
        &ParallelGSSolver::worker,
        this,
        std::cref(As),
        std::cref(bs),
        std::cref(Xs)
    );


    return Xs;
}

void ParallelGSSolver::worker(const std::vector<matrix>& As, const matrix& bs, matrix& Xs, std::atomic<int>& curInd, std::size_t threadId) {
    int last_ind = As.size();
    while(true) {
        int ind = curInd.fetch_add(1, std::memory_order_relaxed);
        if(ind >= last_ind)
            break;
        if(!solveLS(As[ind], bs[ind], Xs[ind]))
            notСonvergedMatrixNumbers[ind] = 1;
    }
}