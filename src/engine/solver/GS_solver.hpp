#pragma once
#include "iter_ls_solver.hpp"

/**
 * Условия применения:
 *  все матрицы передаваемые в солвер являются SPD (ну или хз, какие они там у нас сейчас???)
 */

class GsSolver : public IterLsSolver {
public:
    GsSolver(double tolerance = 0.01, double maxDx = 0.01, int maxIterCount = 100) : 
        LsSolverI(), IterLsSolver(tolerance, maxIterCount), maxDx(maxDx) {}

    matrix solve(const std::vector<matrix>& As, const matrix& bs) override;

    /**
     * В случае если вышли за максимальное число иттераций false,
     * После выхода из иттеративного цикла, если число иттерацией меньше максимального,
     * то в конце выполняется проверка невязки, 
     * если невязка < tolerance - true, 
     * если больше - false
     *
     * Пока не знаю как это использовать, но мб те, в которых не сошлись решать Холецким?
     */
    bool solveLS(const matrix& A, const std::vector<double>& b, std::vector<double>& x);

    const std::vector<int>& getNotCoveredMatrixNumbers() const { return notСonvergedMatrixNumbers; };

protected:
    void isInputCorrect(const std::vector<matrix>& As, const matrix& bs);

    const double maxDx;

    // матрица размера m (кол-во СЛАУ на решение), если матрица не 
    //сошлась, то значение по её индексу = 1. В ином случае 0
    std::vector<int> notСonvergedMatrixNumbers;
};