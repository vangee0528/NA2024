#include "spline.h"
#include "lapack.h"

MathFunction :: MathFunction (double (*func)(double x)) {
    this -> function_ptr = func;
}

double MathFunction :: evaluate (double x) const {
    return (*function_ptr)(x);
}

// 使用系数构造多项式
Polynomial :: Polynomial (const std :: vector <double> &coef) : coefficients (coef) {}

// 使用 Newton 插值法构造多项式
Polynomial :: Polynomial (const std :: vector <double> &x_values, const std :: vector <double> &y_values) {
    if (x_values.size () != y_values.size ())
        throw "Size mismatch";
    int n = x_values.size ();
    std :: vector <double> coefficients (n, 0.0);
    std :: vector <double> dividedDiff (n, 0.0);
    for (int i = 0; i < n; ++ i)
        dividedDiff[i] = y_values[i];
    coefficients[0] = dividedDiff[0];
    for (int i = 1; i < n; ++ i) {
        for (int j = 0; j < n - i; ++ j) 
            dividedDiff[j] = (dividedDiff[j + 1] - dividedDiff[j]) / (x_values[j + i] - x_values[j]);
        coefficients[i] = dividedDiff[0];
    }
    std :: vector <double> standard_c (n, 0.0); 
    for (int i = 0; i < n; ++ i) {
        std :: vector <double> temp (1, coefficients[i]);
        for (int j = 0; j < i; ++ j) {
            temp.push_back(0); 
            for (int k = temp.size () - 1; k > 0; -- k)
                temp[k] = temp[k - 1] - x_values[j] * temp[k];
            temp[0] *= -x_values[j];
        }
        for (int j = 0; j < temp.size (); ++ j)
            standard_c[j] += temp[j];
    }
    this -> coefficients = standard_c;
}

// 输出多项式公式
void Polynomial :: print () const {
    for (int i = 0; i < coefficients.size (); ++ i)
        std :: cout << coefficients[i] << " ";
    std :: cout << std :: endl;
}

// 重载多项式的加法、减法、乘法
Polynomial Polynomial :: operator + (const Polynomial &other) const {
    std :: vector <double> coefficients (std :: max (this -> coefficients.size (), other.coefficients.size ()));
    for (int i = 0; i < coefficients.size (); ++ i) {
        if (i < this -> coefficients.size ())
            coefficients[i] += this -> coefficients[i];
        if (i < other.coefficients.size ())
            coefficients[i] += other.coefficients[i];
    }
    return Polynomial (coefficients);
}

Polynomial Polynomial :: operator - (const Polynomial &other) const {
    std :: vector <double> coefficients (std :: max (this -> coefficients.size (), other.coefficients.size ()));
    for (int i = 0; i < coefficients.size (); ++ i) {
        if (i < this -> coefficients.size ())
            coefficients[i] += this -> coefficients[i];
        if (i < other.coefficients.size ())
            coefficients[i] -= other.coefficients[i];
    }
    return Polynomial (coefficients);
}

Polynomial Polynomial :: operator * (const Polynomial &other) const {
    std :: vector <double> coefficients (this -> coefficients.size () + other.coefficients.size () - 1, 0.0);
    for (int i = 0; i < this -> coefficients.size (); ++ i)
        for (int j = 0; j < other.coefficients.size (); ++ j)
            coefficients[i + j] += this -> coefficients[i] * other.coefficients[j];
    return Polynomial (coefficients);
}

// 计算多项式在 x 点的值
double Polynomial :: evaluate (double x) const {
    double result = 0;
    for (int i = coefficients.size() - 1; i >= 0; -- i) {
        result *= x;
        result += coefficients[i];
    }
    return result;
}

// 返回多项式的导数
Polynomial Polynomial :: derivative () {
    std :: vector <double> coefficients (this -> coefficients.size () - 1);
    for (int i = 1; i < this -> coefficients.size (); ++ i)
        coefficients[i - 1] = this -> coefficients[i] * i;
    return Polynomial (coefficients);
}

// 通过各段的表达式和分段点构造分段多项式
PiecewisePolynomial :: PiecewisePolynomial (const std :: vector <Polynomial> &p, const std :: vector <double> &x) : polynomials (p), points (x) {}

// 计算分段多项式在 x 点的值
double PiecewisePolynomial :: evaluate(double x) const {
    if (x < this -> points[0] || x > this -> points.back ())
        throw "Out of range";
    for (int i = 0; i < this -> points.size () - 1; ++ i)
        if (x >= this -> points[i] && x <= this -> points[i + 1])
            return this -> polynomials[i].evaluate (x);
    throw "Out of range";
}

// 输出分段多项式公式
void PiecewisePolynomial :: print () const {
    for (int i = 0; i < polynomials.size (); ++ i) {
        std :: cout << points[i] << " " << points[i + 1] << std :: endl;
        polynomials[i].print ();
    }
}

// 返回分段多项式的导数
PiecewisePolynomial PiecewisePolynomial :: derivative () {
    std :: vector <Polynomial> polynomials (this -> polynomials.size ());
    for (int i = 0; i < polynomials.size (); ++ i)
        polynomials[i] = this -> polynomials[i].derivative ();
    return PiecewisePolynomial (polynomials, this -> points);
}


Curve :: Curve (int dim, const std :: vector <MathFunction *> &functions) {
    this -> dimensions = dim;
    if (functions.size () != dimensions)
        throw "Dimension mismatch";

    for (int i = 0; i < functions.size (); ++ i)
        this -> parametric_functions.push_back (functions[i]);
}

// 计算曲线在 t 点的值
std :: vector <double> Curve :: operator ()(double t) const {
    std :: vector <double> result (dimensions);
    for (int i = 0; i < dimensions; ++ i)
        result[i] = (*parametric_functions[i]).evaluate (t);
    return result;
}

std::vector<double> Spline::operator()(double t) const{
    std::vector<double> result(dimensions);
    for (int i = 0; i < dimensions; ++i)
        result[i] = (this ->segments[i]).evaluate(t);
    return result;
}
void Spline :: print () const {
    std :: cout << dimensions << " " << segments[0].segment_count() << std :: endl;
    for (int i = 0; i < dimensions; ++ i) {
        std :: cout << "Dimension " << i + 1 << std :: endl;
        (this ->segments[i]).print ();
    }
}

PiecewisePolynomial PPSpline::compute_spline_segments(SplineBoundaryCondition bc, const std::vector<double>& f, const std::vector<double>& t, double da, double db) {
    if (spline_order != 1 && spline_order != 3)
        throw "Order Error";
    if (f.size() != t.size())
        throw "Size mismatch";
    
    int N = f.size();
    std::vector<Polynomial> _c(N - 1);  // 存储每段多项式的系数

    if (spline_order == 1) {
        // 一阶样条函数（线性插值）
        for (int j = 0; j < N - 1; ++j) {
            std::vector<double> x(2), y(2);
            x[0] = t[j]; x[1] = t[j + 1];
            y[0] = f[j]; y[1] = f[j + 1];
            _c[j] = Polynomial(x, y);
        }
        return PiecewisePolynomial(_c, t);
    } else {
        if (bc == NO_CONDITION)
            throw "Boundary condition must be given";

        // 三阶样条函数，构造二阶差商表
        std::vector<double> dividedDiff(N - 1);  // 计算一阶差商
        for (int i = 0; i < N - 1; ++i)
            dividedDiff[i] = (f[i] - f[i + 1]) / (t[i] - t[i + 1]);

        std::vector<double> dividedDiff2(N - 2); // 计算二阶差商
        for (int i = 0; i < N - 2; ++i)
            dividedDiff2[i] = (dividedDiff[i + 1] - dividedDiff[i]) / (t[i + 2] - t[i]);

        // 构建 lambda 和 mu 数组
        std::vector<double> lambda(N - 2);
        std::vector<double> mu(N - 2);
        for (int i = 0; i < N - 2; ++i) {
            mu[i] = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
            lambda[i] = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
        }

        // 构造方程组右侧向量 B
        std::vector<double> B(N - 2);
        if (bc == CLAMPED) {
            for (int i = 0; i < N - 2; ++i)
                B[i] = 3 * mu[i] * dividedDiff[i + 1] + 3 * lambda[i] * dividedDiff[i];
            B[0] -= da * lambda[0];
            B[N - 3] -= db * mu[N - 3];
        } else if (bc == NATURAL_SPLINE) {
            for (int i = 0; i < N - 2; ++i)
                B[i] = 6 * dividedDiff2[i];
        } else if (bc == SECOND_DERIVATIVE_FIXED) {
            for (int i = 0; i < N - 2; ++i)
                B[i] = 6 * dividedDiff2[i];
            B[0] -= da * lambda[0];
            B[N - 3] -= db * mu[N - 3];
        }


        std::vector<double> B_copy = B;
        std::vector<double> diagonal(N - 2, 2.0); 
        std::vector<double> M = lapack_solve(diagonal, lambda, mu, B_copy);

        if (bc == CLAMPED) {
            std::vector<double> _M(N);
            _M[0] = da;
            for (int i = 0; i < N - 2; ++i)
                _M[i + 1] = M[i];
            _M[N - 1] = db;

            for (int i = 0; i < N - 1; ++i) {
                std::vector<double> coe(4);
                coe[0] = f[i];
                coe[1] = _M[i];
                coe[2] = (3 * dividedDiff[i] - 2 * _M[i] - _M[i + 1]) / (t[i + 1] - t[i]);
                coe[3] = (_M[i] + _M[i + 1] - 2 * dividedDiff[i]) / (t[i + 1] - t[i]) / (t[i + 1] - t[i]);

                Polynomial base({0}), power({1});
                for (int j = 0; j < 4; ++j) {
                    Polynomial temp({coe[j]});
                    temp = temp * power;
                    base = base + temp;
                    power = power * Polynomial({-t[i], 1});
                }
                _c[i] = base;
            }
        } else if (bc == SECOND_DERIVATIVE_FIXED || bc == NATURAL_SPLINE) {
            std::vector<double> _M(N);
            _M[0] = bc == NATURAL_SPLINE ? 0 : da;
            for (int i = 0; i < N - 2; ++i)
                _M[i + 1] = M[i];
            _M[N - 1] = bc == NATURAL_SPLINE ? 0 : db;

            for (int i = 0; i < N - 1; ++i) {
                std::vector<double> coe(4);
                coe[0] = f[i];
                coe[1] = dividedDiff[i] - (_M[i + 1] + 2 * _M[i]) / 6 * (t[i + 1] - t[i]);
                coe[2] = _M[i] / 2;
                coe[3] = (_M[i + 1] - _M[i]) / 6 / (t[i + 1] - t[i]);

                Polynomial base({0}), power({1});
                for (int j = 0; j < 4; ++j) {
                    Polynomial temp({coe[j]});
                    temp = temp * power;
                    base = base + temp;
                    power = power * Polynomial({-t[i], 1});
                }
                _c[i] = base;
            }
        }

        return PiecewisePolynomial(_c, t);
    }
}


//
PPSpline :: PPSpline (int dim, int order, const MathFunction &f, double a, double b, SplineBoundaryCondition bc, int N, double da, double db) : Spline (dim, order) {
    if (a >= b)
        throw "Invalid interval";
    if (N <= 1)
        throw "Invalid number of intervals";
    
    std :: vector <double> t (N), _f (N);
    for (int j = 0; j < N; ++ j) {
        t[j] = a + (b - a) * j / (N - 1);
        _f[j] = f.evaluate(t[j]);
    }
    this -> segments = { compute_spline_segments(bc, _f, t, da, db) };
}

PPSpline :: PPSpline (int dim, int order, const MathFunction &f, const std :: vector <double> &t, SplineBoundaryCondition bc, double da, double db) : Spline (dim, order) {
    if (dim != 1)
        throw "Dimension must be 1";
    if (t.size () <= 1)
        throw "Invalid number of intervals";
    
    std :: vector <double> _f (t.size ());
    for (int j = 0; j < t.size (); ++ j)
        _f[j] = f.evaluate(t[j]);
    this -> segments = { compute_spline_segments (bc, _f, t, da, db) };
}

PPSpline :: PPSpline (int dim, int order, const std :: vector <std :: vector <double> > &points, SplineBoundaryCondition bc) : Spline (dim, order) {
    if (points.size () <= 1)
        throw "Invalid number of intervals";
    if (points[0].size () != dim)
        throw "Dimension mismatch";
    
    std :: vector <double> t;
    // 使用累积弦长法来产生行止点
    t.push_back (0);
    for (int i = 1; i < points.size (); ++ i) {
        double dist = 0;
        static double totalDist = 0;
        for (int j = 0; j < points[i].size (); ++ j)
            dist += (points[i][j] - points[i - 1][j]) * (points[i][j] - points[i - 1][j]);
        dist = sqrt (dist);
        totalDist += dist;
        t.push_back (totalDist);
    }

    // 对每个维度分别产生样条曲线
    for (int i = 0; i < dim; ++ i) {
        std :: vector <double> _f;
        for (int j = 0; j < points.size (); ++ j)
            _f.push_back (points[j][i]);
        (this -> segments).push_back (compute_spline_segments (bc, _f, t, 0.0, 0.0));
    }
}

double BSpline :: evaluate_basis (int i, int k, double x) const {
    if (!k)
        return knot_vector[i - 1] < x && x <= knot_vector[i];
    else
        return (x - knot_vector[i - 1]) / (knot_vector[i + k - 1] - knot_vector[i - 1]) * evaluate_basis (i, k - 1, x) + (knot_vector[i + k] - x) / (knot_vector[i + k] - knot_vector[i]) * evaluate_basis (i + 1, k - 1, x);
}

double BSpline :: evaluate_basis_derivative (int i, int k, double x) const {
    return k * evaluate_basis (i, k - 1, x) / (knot_vector[i + k - 1] - knot_vector[i - 1]) - k * evaluate_basis (i + 1, k - 1, x) / (knot_vector[i + k] - knot_vector[i]);
}

double BSpline :: evaluate_basis_second_derivative (int i, int k, double x) const {
    return k * evaluate_basis_derivative (i, k - 1, x) / (knot_vector[i + k - 1] - knot_vector[i - 1]) - k * evaluate_basis_derivative (i + 1, k - 1, x) / (knot_vector[i + k] - knot_vector[i]);
}

PiecewisePolynomial BSpline :: compute_spline_segments (SplineBoundaryCondition bc, const std :: vector <double> &f, const std :: vector <double> &t, double da, double db) {
    if (dimensions > 3)
        throw "Dimension unsupported";
    (this -> knot_vector).clear ();
    for (int i = spline_order; i >= 1; -- i)
        (this -> knot_vector).push_back (t[0] - i);
    for (int i = 0; i < t.size (); ++ i)
        (this -> knot_vector).push_back (t[i]);
    for (int i = 1; i <= spline_order; ++ i)
        (this -> knot_vector).push_back (t[t.size () - 1] + i);
    int N = f.size ();
    std :: vector <std :: vector <double>> A (N + spline_order - 1, std :: vector <double> (N + spline_order - 1, 0.0));
    std :: vector <double> B (N + spline_order - 1, 0.0);
    if (spline_order == 3) {
        // 构造方程组
        for (int i = 0; i < N; ++ i) {
            for (int j = 0; j < spline_order; ++ j)
                A[i][i + j] = evaluate_basis (i + j + 1, spline_order, t[i]);
            B[i] = f[i];
        }
        // 默认使用自然边界条件
        for (int i = 0; i < spline_order; ++ i) {
            A[N][i] = evaluate_basis_second_derivative (i + 1, spline_order, t[0]);
            A[N + 1][N - 1 + i] = evaluate_basis_second_derivative (N + i, spline_order, t[N - 1]);
            B[N] = 0;
            B[N + 1] = 0;
        }
    } else if (spline_order == 2) {
        A[0][0] = evaluate_basis (1, spline_order, t[0]);
        A[0][1] = evaluate_basis (2, spline_order, t[0]);
        B[0] = f[0];
        for (int i = 0; i < N - 1; ++ i) {
            double m = (t[i] + t[i + 1]) / 2;
            A[i + 1][i] = evaluate_basis (i + 1, spline_order, m);
            A[i + 1][i + 1] = evaluate_basis (i + 2, spline_order, m);
            A[i + 1][i + 2] = evaluate_basis (i + 3, spline_order, m);
            B[i + 1] = f[i + 1];
        }
        A[N][N - 1] = evaluate_basis (N, spline_order, t[N - 1]);
        A[N][N] = evaluate_basis (N + 1, spline_order, t[N - 1]);
        B[N] = f[N - 1];
    }
    // 求解方程组
    std :: vector <double> coe;
    if (spline_order > 1)
        coe = GaussEliminate (N + spline_order - 1, A, B);
    else
        coe = f;
    // 构建返回结果的逐段多项式
    std :: vector <Polynomial> coefficients (N - 1);
    for (int i = 0; i < N - 1; ++ i) {
        // 算基函数太麻烦了，不如取几个点来插值
        std :: vector <double> x (spline_order + 1), y (spline_order + 1);
        for (int j = 0; j <= spline_order; ++ j) {
            x[j] = t[i] + (t[i + 1] - t[i]) * j / spline_order;
            y[j] = 0;
            for (int k = 0; k < N + spline_order - 1; ++ k)
                y[j] += coe[k] * evaluate_basis (k + 1, spline_order, x[j]);
        }
        coefficients[i] = Polynomial (x, y);
    }
    return PiecewisePolynomial(coefficients, t);
}

BSpline :: BSpline (int dim, int order, const MathFunction &f, double a, double b, int N) : Spline (dim, order) {
    if (a >= b)
        throw "Invalid interval";
    if (N <= 1)
        throw "Invalid number of intervals";
    
    std :: vector <double> t (N), _f (N);
    for (int j = 0; j < N; ++ j) {
        t[j] = a + (b - a) * j / (N - 1);
        _f[j] = f.evaluate(t[j]);
    }
    this -> segments = { compute_spline_segments (NO_CONDITION, _f, t, 0, 0) };
}

BSpline :: BSpline (int dim, int order, const MathFunction &f, const std :: vector <double> &t) : Spline (dim, order) {
    if (dim != 1)
        throw "Dimension must be 1";
    if (t.size () <= 1)
        throw "Invalid number of intervals";
    
    std :: vector <double> _f (t.size ());
    for (int j = 0; j < t.size (); ++ j)
        _f[j] = f.evaluate(t[j]);
    this -> segments = { compute_spline_segments (NO_CONDITION, _f, t, 0, 0) };
}

BSpline :: BSpline (int dim, int order, const std :: vector <std :: vector <double> > &points) : Spline (dim, order) {
    if (points.size () <= 1)
        throw "Invalid number of intervals";
    if (points[0].size () != dim)
        throw "Dimension mismatch";
    
    std :: vector <double> t;
    // 使用累积弦长法来产生行止点
    t.push_back (0);
    for (int i = 1; i < points.size (); ++ i) {
        double dist = 0;
        static double totalDist = 0;
        for (int j = 0; j < points[i].size (); ++ j)
            dist += (points[i][j] - points[i - 1][j]) * (points[i][j] - points[i - 1][j]);
        dist = sqrt (dist);
        totalDist += dist;
        t.push_back (totalDist);
    }

    // 对每个维度分别产生样条曲线
    for (int i = 0; i < dim; ++ i) {
        std :: vector <double> _f;
        for (int j = 0; j < points.size (); ++ j)
            _f.push_back (points[j][i]);
        (this -> segments).push_back (compute_spline_segments (NO_CONDITION, _f, t, 0, 0));
    }
}