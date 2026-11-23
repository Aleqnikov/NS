#pragma once
#include <vector>
#include "../../sparse_matrix.hpp"

namespace util
{
    std::vector<double> multiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vec) {
        std::vector<double> res(matrix.size());
        res.reserve(matrix.size());

        long double sum;
        for(int i = 0; i < matrix.size(); i++) {
            sum = 0;
            for(int k = 0; k < matrix[0].size(); k++) {
                sum += matrix[i][k] * vec[k];
            }
            res[i] = sum;
        }
        return res;
    }

    long double multiply_vector(const std::vector<double>& vec1, const std::vector<double>& vec2) {
        long double sum = 0;
        for(int i = 0; i < vec1.size(); i++) 
            sum += vec1[i] * vec2[i]; 
        return sum;
    }

    std::vector<double> multiply_sparse(const SparseMatrix& A, const std::vector<double>& x) {
        std::vector<double> res;
        res.reserve(A.size());

        long double sum;
        for(int i = 0; i < A.size(); i++) {
            sum = 0;
            for(int k = 0; k < A[i].indexes.size(); k++) {
                sum += A[i].values[k] * x[A[i].indexes[k]];
            }
            res.push_back(sum);
        }
        return res;
    }
    
}