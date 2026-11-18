#include "CS_solver.hpp"
#include <iostream>
#include "util/util.hpp"
#include <cmath>
#include <limits>
#include <algorithm>

std::vector<double> CGSolver::solve(const SparseMatrix& A, const std::vector<double>& b) {
    if (A.size() == 0)
        throw std::logic_error("A matrix can't empty");
    
    std::vector<double> x(A.size(), 0);
    resTol = constants::calcErrorTolerance(A.size());

    if(!solveLS(A, b, x)) 
        std::cout << "The matrix did not converge, the final error is: " << calcBackwardError(A, b, x) << std::endl;
    return x;
}


std::vector<double> CGSolver::solve(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) {
    if (A.size() == 0)
        throw std::logic_error("A matrix can't empty");
    
    resTol = constants::calcErrorTolerance(A.size());

    if(!solveLS(A, b, x)) 
        std::cout << "The matrix did not converge, the final error is: " << calcBackwardError(A, b, x) << std::endl;
    return x;
}


bool CGSolver::solveLS(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x) {
    if(A.size() != x.size())
        x.resize(A.size(), 0);
    std::vector<double> r(b.size());
    {
        std::vector<double> Ax = util::multiply_sparse(A, x);
        for(int i = 0; i < b.size(); i ++) 
            r[i] = b[i] - Ax[i];
    }
    
    PriCondLC0 M;
    M.build(A);
    std::vector<double> z;
    M.precondition(r, z);

    std::vector<double> p = z;
    double rz_old_norm = util::multiply_vector(r, z);
    long double alph;
    double beta, rz_new_norm, pAp;


    double r0_norm = util::multiply_vector(r, r);
    if(std::abs(r0_norm) < std::numeric_limits<double>::epsilon())
        return true;
    double r0_norm_sqrt = std::sqrt(r0_norm);
    
    for(int k = 0; k < maxIterCount; k++) {
        std::vector<double> Ap = util::multiply_sparse(A, p);
        pAp =  util::multiply_vector(p, Ap);
        if(pAp <= 0) //может возникнуть только когда матрица не SPD
            return false;
        alph = rz_old_norm / pAp;
        for(int i = 0; i < A.size(); i++) {
            x[i] = x[i] + alph * p[i];
            r[i] = r[i] - alph * Ap[i];
        }
    
        //проверка сходимости
        if(std::sqrt(util::multiply_vector(r, r)) < r0_norm_sqrt * resTol)
            return true;

        M.precondition(r, z);
        rz_new_norm = util::multiply_vector(r, z);
        beta = rz_new_norm / rz_old_norm;
        rz_old_norm = rz_new_norm;
        for(int i = 0; i < A.size(); i++)
            p[i] = z[i] + beta * p[i];
    }
    return false;
}

void CGSolver::PriCondLC0::build(const SparseMatrix& A) {
    cols.clear();
    vals.clear();
    /**
     * на 1 больше тк храним границы отрезок. 
     * [a0, a1), [a1, a2), ... 
     * n+1 число
    */
    row_bounds.assign(A.size() + 1, 0);

    /** 
     * создаём контейнер для элемента матрицы A (строки). Т.к. порядок элементов (индексов) внутри A[n] не детерменирован,
     * а нам для разложения Холецкого (в такой реализации) нужно знать порядок. Сначала складываем всё в контейнер, затем 
     * сортируем по индексам (индексы повторяться не могут) 
    */
    std::vector<std::pair<int, double>> row;
    for(int i = 0; i < A.size(); i++) {
        row.clear();
        for(int k = 0; k < A[i].indexes.size(); k++) {
            if(A[i].indexes[k] <= i) //т.к. A - SPD
                row.push_back({A[i].indexes[k], A[i].values[k]});
        }
        std::sort(
            row.begin(), row.end(),
            [](const std::pair<int, double>& p1, const std::pair<int, double>& p2) -> bool {
                return p1.first < p2.first; //по возрастанию
            }
        );
        for(const auto& p : row) {
            cols.push_back(p.first);
            vals.push_back(p.second);
        }
        row_bounds[i + 1] = cols.size();
    }

    //Разложение
    int row_beg, row_end;
    for(int i = 0; i < A.size(); i++) {
        row_beg = row_bounds[i];
        row_end = row_bounds[i + 1];
        for(int p = row_beg; p < row_end; p++) {
            int j = cols[p];
            long double sum = vals[p];
            int p_i = row_beg;
            int p_j = row_bounds[j];
            int end_j = row_bounds[j+1] - 1;
            while(p_i < p && p_j < end_j) { //k<j (см. формулу разложения)
                if(cols[p_i] == cols[p_j]) { //нашли общий столбец (k). l_ik*l_jk
                    sum -= vals[p_i] * vals[p_j];
                    p_i++;
                    p_j++;
                } else if(cols[p_i] < cols[p_j]) {
                    p_i++;
                }else { //cols[p_i] > cols[p_j]
                    p_j++;
                }
            }
            //вычислили a_ij-sum_(k<j)(l_ik*l_jk)
            //в зависимости от i и j является частью или для внедиагонального или для диагонального (i==j)

            if(j < i) { //внедиагональный
                if(!std::isfinite(sum))
                    throw std::runtime_error("error in decomp!");
                vals[p] = sum / vals[row_bounds[j + 1] - 1]; 
                // -1 т.к. полуинтервалы. Когда мы попадаем сюда, соотвествующий диагональный элемент уже вычислен, т.к. j < i
            }else { //i==j. i < j - не бывает
                if(sum <= std::numeric_limits<double>::epsilon() || !std::isfinite(sum)) //возможно стоит заменить epsilon на >= 0
                    throw std::runtime_error("error in decomp!");
                vals[p] = std::sqrt(sum);
            }
        }
    }
    //итого мы построили нижнетреугольную матрицу L в разряженном формате
}

/**
 * решаем LL^T*z=r. Т.е.
 * L*y=r
 * L^T*z=y
*/
void CGSolver::PriCondLC0::precondition(const std::vector<double>& r, std::vector<double>& z) const {
    if(r.size() != row_bounds.size() - 1)
        throw std::logic_error("r size must equals A size");
    z = r;
    //L*y=r
    for(int i = 0; i < r.size(); i++) {
        long double sum = z[i];
        int diag = row_bounds[i + 1] - 1;
        for(int p = row_bounds[i]; p < diag; p++) 
            sum -= vals[p] * z[cols[p]];
        z[i] = sum / vals[diag];
    }
    //L^T*z=y
    for(int i = r.size() - 1; i >= 0; i--) {
        int diag = row_bounds[i+1] - 1;
        z[i] = z[i] / vals[diag];
        // т.к. отдельного представления L^T у нас нет, используем строки L как столбцы L^T
        // не считаем сумму по уже найденным компонентам при вычислении каждого z[i], а
        // после нахождения z[i] сразу вычитаем его вклад из правых частей ещё не решённых уравнений
        for(int p = row_bounds[i]; p < diag; p++) {
            z[cols[p]] -= vals[p] * z[i];
        }
    }
}   