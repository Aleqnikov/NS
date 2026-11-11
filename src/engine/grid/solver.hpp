#pragma once

#include <vector>


#include <array>   // <-- добавить

struct SparseRow
{
    static constexpr std::size_t MAX_NNZ = 5;
    std::array<std::pair<std::size_t, double>, MAX_NNZ> entries{};
    std::size_t nnz = 0;

    void add(std::size_t col, double val)
    {
        // при корректном построении — всегда ≤ 5
        entries[nnz++] = {col, val};
    }
};

using SparseMatrix = std::vector<SparseRow>;

inline SparseMatrix make_sparse_matrix(std::size_t n)
{
    return SparseMatrix(n);   // nnz = 0 по умолчанию
}

std::vector<double> GS(const SparseMatrix& A, const std::vector<double>& b, double tol);





