#pragma once
#include <vector>
#include <stdexcept>

struct SparseRow
{
    static const size_t max_nnz = 5;
    std::vector<double> values;
    std::vector<size_t> indexes;

    SparseRow() {
        values.reserve(max_nnz);
        indexes.reserve(max_nnz);
    }

    void add(std::size_t col, double val)
    {
        //закоментить после первых проверок
        if(values.size() >= max_nnz)
            throw std::logic_error("SparceRow can contain no more than five elements");
        values.push_back(val);
        indexes.push_back(col);
    }
};

using SparseMatrix = std::vector<SparseRow>;
