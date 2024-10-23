#include "interpolation.h"
/* ====函数实现==== */

// POLYNOMIAL CLASS
Polynomial::Polynomial() {}

Polynomial::Polynomial(int degree, const std::vector<double>& coeffs) : coefficients(coeffs) {
    coefficients.resize(degree + 1, 0);
}

double Polynomial::evaluate(double x) const {
    double result = 0.0;
    for (int i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * pow(x, i);
    }
    return result;
}

Polynomial Polynomial::derivative() const {
    std::vector<double> newCoeffs(coefficients.size() - 1, 0.0);
    for (int i = 1; i < coefficients.size(); ++i) {
        newCoeffs[i - 1] = coefficients[i] * i;
    }
    return Polynomial(newCoeffs.size() - 1, newCoeffs);
}

void Polynomial::print() const {
    for (int i = coefficients.size() - 1; i >= 0; --i) {
        if (i < coefficients.size() - 1 && coefficients[i] >= 0) {
            std::cout << "+";
        }
        std::cout << std::fixed << std::setprecision(2) << coefficients[i] << "x^" << i << " ";
    }
    std::cout << std::endl;
}

Polynomial operator*(const Polynomial& p1, const Polynomial& p2) {
    std::vector<double> newCoeffs(p1.coefficients.size() + p2.coefficients.size() - 1, 0.0);
    for (int i = 0; i < p1.coefficients.size(); ++i) {
        for (int j = 0; j < p2.coefficients.size(); ++j) {
            newCoeffs[i + j] += p1.coefficients[i] * p2.coefficients[j];
        }
    }
    return Polynomial(newCoeffs.size() - 1, newCoeffs);
}

Polynomial operator+(const Polynomial& p1, const Polynomial& p2) {
    int maxSize = std::max(p1.coefficients.size(), p2.coefficients.size());
    std::vector<double> newCoeffs(maxSize, 0.0);

    for (int i = 0; i < maxSize; ++i) {
        if (i < p1.coefficients.size()) {
            newCoeffs[i] += p1.coefficients[i];
        }
        if (i < p2.coefficients.size()) {
            newCoeffs[i] += p2.coefficients[i];
        }
    }
    return Polynomial(newCoeffs.size() - 1, newCoeffs);
}

// NEWTON INTERPOLATOR CLASS
double NewtonInterpolator::dividedDifference(int i, int j) const {
    if (i == j) {
        return f_[i];
    }
    return (dividedDifference(i + 1, j) - dividedDifference(i, j - 1)) / (x_[j] - x_[i]);
}

NewtonInterpolator::NewtonInterpolator(const std::vector<double>& x, const std::vector<double>& f)
    : x_(x), f_(f), coeffs_(x.size()) {
    for (size_t i = 0; i < x.size(); ++i) {
        coeffs_[i] = dividedDifference(0, i);
    }
}

Polynomial NewtonInterpolator::interpolate() const {
    int n = x_.size();
    Polynomial result(0, {coeffs_[0]});
    Polynomial term(0, {1});

    for (int i = 1; i < n; ++i) {
        term = term * Polynomial(1, {-x_[i - 1], 1});
        result = result + term * Polynomial(0, {coeffs_[i]});
    }

    return result;
}

// HERMITE INTERPOLATOR CLASS
HermiteInterpolator::HermiteInterpolator(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& df)
    : x_(x), f_(f), df_(df), z_(2 * x.size()), q_(2 * x.size(), std::vector<double>(2 * x.size(), 0.0)) {
    computeDividedDifferences();
}

void HermiteInterpolator::computeDividedDifferences() {
    int n = x_.size();

    for (int i = 0; i < n; ++i) {
        z_[2 * i] = x_[i];
        z_[2 * i + 1] = x_[i];
        q_[2 * i][0] = f_[i];
        q_[2 * i + 1][0] = f_[i];
        q_[2 * i + 1][1] = df_[i];
        if (i != 0) {
            q_[2 * i][1] = (q_[2 * i][0] - q_[2 * i - 1][0]) / (z_[2 * i] - z_[2 * i - 1]);
        }
    }

    for (int i = 2; i < 2 * n; ++i) {
        for (int j = 2; j <= i; ++j) {
            q_[i][j] = (q_[i][j - 1] - q_[i - 1][j - 1]) / (z_[i] - z_[i - j]);
        }
    }
}

Polynomial HermiteInterpolator::interpolate() const {
    int n = x_.size();
    Polynomial result(0, {q_[0][0]});
    Polynomial term(0, {1});

    for (int i = 1; i < 2 * n; ++i) {
        term = term * Polynomial(1, {-z_[i - 1], 1});
        result = result + term * Polynomial(0, {q_[i][i]});
    }

    return result;
}

