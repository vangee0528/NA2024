#ifndef LAPACK_H
#define LAPACK_H

#include <vector>
#include <stdexcept>

// LAPACK 的三对角矩阵求解函数 dgtsv
extern "C" {
    void dgtsv_(int* n, int* nrhs, double* dl, double* d, double* du, double* b, int* ldb, int* info);
}

// LAPACK 的普通矩阵求解函数 dgesv
extern "C" {
    void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
}

// 解三对角线性方程组 Ax = b
std::vector<double> lapack_solve_t(const std::vector<double>& diag,
                                 const std::vector<double>& lambda,
                                 const std::vector<double>& mu,
                                 const std::vector<double>& b);

// 解一般线性方程组 Ax = b
std::vector<double> lapack_solve(int n,
                                   const std::vector<std::vector<double>>& A,
                                   const std::vector<double>& b);

#endif // LAPACK_H
