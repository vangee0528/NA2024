#include "lapack.h"
#include <vector>
#include <stdexcept>
#include <cstring>

// LAPACK 解三对角矩阵
std::vector<double> lapack_solve_t(const std::vector<double>& diag,
                                 const std::vector<double>& lambda,
                                 const std::vector<double>& mu,
                                 const std::vector<double>& b) {
    int n = diag.size();

    if (lambda.size() != n - 1 || mu.size() != n - 1 || b.size() != n) {
        throw std::invalid_argument("Input dimensions do not match!");
    }

    std::vector<double> dl = lambda;
    std::vector<double> d = diag;
    std::vector<double> du = mu;
    std::vector<double> B = b;

    int nrhs = 1;
    int ldb = n;
    int info;

    dgtsv_(&n, &nrhs, dl.data(), d.data(), du.data(), B.data(), &ldb, &info);

    if (info != 0) {
        throw std::runtime_error("LAPACK dgtsv failed with error code: " + std::to_string(info));
    }

    return B;
}

// LAPACK 解一般线性方程组
std::vector<double> lapack_solve(int n,
                                   const std::vector<std::vector<double>>& A,
                                   const std::vector<double>& b) {
    if (A.size() != n || b.size() != n) {
        throw std::invalid_argument("Matrix dimensions do not match!");
    }

    std::vector<double> a_flat(n * n);
    // 转为列优先顺序
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            a_flat[j * n + i] = A[i][j];
        }
    }

    std::vector<double> b_copy = b;
    std::vector<int> ipiv(n);
    int lda = n;
    int ldb = n;
    int nrhs = 1;
    int info;

    // 调用 LAPACK 函数
    dgesv_(&n, &nrhs, a_flat.data(), &lda, ipiv.data(), b_copy.data(), &ldb, &info);

    if (info > 0) {
        throw std::runtime_error("Matrix is singular and cannot be solved");
    } else if (info < 0) {
        throw std::runtime_error("Invalid argument passed to LAPACK dgesv");
    }

    return b_copy;
}
