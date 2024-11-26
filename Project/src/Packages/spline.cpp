#include "spline.h"
#include "lapack.h"

////////////////////////////////////////////////////////////////////

/* MathFunction 类 */

MathFunction :: MathFunction (double (*func)(double x)) {
    this -> function_ptr = func;
}

double MathFunction :: evaluate (double x) const {
    return (*function_ptr)(x);
}

////////////////////////////////////////////////////////////////////

/* Polynomial 类 */

// 使用系数构造多项式
Polynomial :: Polynomial (const std :: vector <double> &coef) : coefficients (coef) {}

// 使用 Newton 插值法构造多项式
Polynomial :: Polynomial (const std :: vector <double> &x_values, const std :: vector <double> &y_values) {
    if (x_values.size () != y_values.size ())
        throw "Newton interpolation : Size mismatch";
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

// 输出多项式公式
void Polynomial :: print () const {
    for (int i = 0; i < coefficients.size (); ++ i)
        std :: cout << coefficients[i] << " "; //从低到高输出系数
    std :: cout << std :: endl;
}


////////////////////////////////////////////////////////////////////

/* PiecewisePolynomial 类 */

// 计算分段多项式在 x 点的值
double PiecewisePolynomial :: evaluate(double x) const {
    if (x < this -> points[0] || x > this -> points.back ())
        throw "Evaluate PiecewisePolynomial : Out of range";

    for (int i = 0; i < this -> points.size () - 1; ++ i)
        if (x >= this -> points[i] && x <= this -> points[i + 1])
            return this -> polynomials[i].evaluate (x);

    throw "Evaluate PiecewisePolynomial : Out of range";
}

// 返回分段多项式的导数
PiecewisePolynomial PiecewisePolynomial :: derivative () {
    std :: vector <Polynomial> polynomials (this -> polynomials.size ());
    for (int i = 0; i < polynomials.size (); ++ i)
        polynomials[i] = this -> polynomials[i].derivative ();
    return PiecewisePolynomial (polynomials, this -> points);
}


// 打印分段多项式公式
void PiecewisePolynomial :: print () const {
    for (int i = 0; i < polynomials.size (); ++ i) {
        std :: cout << points[i] << "," << points[i + 1] << std :: endl;
        polynomials[i].print ();
    }
}


////////////////////////////////////////////////////////////////////

/* Curve 类 */

// // 使用参数方程构造曲线
// Curve :: Curve (int dim, const std :: vector <MathFunction *> &functions) {
//     this -> dimensions = dim;
//     if (functions.size () != dimensions)
//         throw "Curve : Dimension mismatch";

//     for (int i = 0; i < functions.size (); ++ i)
//         this -> parametric_functions.push_back (functions[i]);
// }

// // 计算曲线在 t 点的值
// std :: vector <double> Curve :: operator ()(double t) const {
//     std :: vector <double> result (dimensions);
//     for (int i = 0; i < dimensions; ++ i)
//         result[i] = (*parametric_functions[i]).evaluate (t);
//     return result;
// }

////////////////////////////////////////////////////////////////////

/* Spline 类 */

std::vector<double> Spline::operator()(double t) const{
    std::vector<double> result(dimensions);
    for (int i = 0; i < dimensions; ++i)
        result[i] = (this ->segments[i]).evaluate(t);
    return result;
}

// 输出曲线公式
void Spline :: print () const {
    for (int i = 0; i < dimensions; ++ i) {
        (this ->segments[i]).print ();
    }
}


////////////////////////////////////////////////////////////////////

/* PPSpline 类 */

// 计算分段样条函数
PiecewisePolynomial PPSpline::compute_spline_segments(SplineBoundaryCondition bc, const std::vector<double>& f, const std::vector<double>& t, double da, double db) {
    if (spline_order != 1 && spline_order != 3)
        throw "PPSpline : Order Error";
    if (f.size() != t.size())
        throw "PPSpline : Size mismatch";
    
    int num_points = f.size();
    std::vector<Polynomial> polynomials(num_points - 1);  

    if (spline_order == 1) {
        // 一阶样条函数
        for (int j = 0; j < num_points - 1; ++j) {
            std::vector<double> x(2), y(2);
            x[0] = t[j]; x[1] = t[j + 1];
            y[0] = f[j]; y[1] = f[j + 1];
            polynomials[j] = Polynomial(x, y);
        }
        return PiecewisePolynomial(polynomials, t);
    } else {
        if (bc == NO_CONDITION)
            throw "PPSpline : Boundary condition must be given";

        // 三阶样条函数
        std::vector<double> first_divided_diff(num_points - 1);  
        for (int i = 0; i < num_points - 1; ++i)
            first_divided_diff[i] = (f[i] - f[i + 1]) / (t[i] - t[i + 1]);

        std::vector<double> second_divided_diff(num_points - 2); 
        for (int i = 0; i < num_points - 2; ++i)
            second_divided_diff[i] = (first_divided_diff[i + 1] - first_divided_diff[i]) / (t[i + 2] - t[i]);

        // 构建 lambda 和 mu 数组
        std::vector<double> lambda(num_points - 2);
        std::vector<double> mu(num_points - 2);
        for (int i = 0; i < num_points - 2; ++i) {
            mu[i] = (t[i + 1] - t[i]) / (t[i + 2] - t[i]);
            lambda[i] = (t[i + 2] - t[i + 1]) / (t[i + 2] - t[i]);
        }

        // 构造方程组右侧向量 B
        std::vector<double> B(num_points - 2);
        if (bc == CLAMPED) {
            for (int i = 0; i < num_points - 2; ++i)
                B[i] = 3 * mu[i] * first_divided_diff[i + 1] + 3 * lambda[i] * first_divided_diff[i];
            B[0] -= da * lambda[0];
            B[num_points - 3] -= db * mu[num_points - 3];
        } else if (bc == NATURAL_SPLINE) {
            for (int i = 0; i < num_points - 2; ++i)
                B[i] = 6 * second_divided_diff[i];
        } else if (bc == SECOND_DERIVATIVE_FIXED) {
            for (int i = 0; i < num_points - 2; ++i)
                B[i] = 6 * second_divided_diff[i];
            B[0] -= da * lambda[0];
            B[num_points - 3] -= db * mu[num_points - 3];
        }

        std::vector<double> diagonal(num_points - 2, 2.0); 
        //对lambda和mu切片
        lambda.pop_back();
        mu.erase(mu.begin());

        std::vector<double> M = lapack_solve_t(diagonal, lambda, mu, B);

        if (bc == CLAMPED) {
            std::vector<double> M_extended(num_points);
            M_extended[0] = da;
            for (int i = 0; i < num_points - 2; ++i)
                M_extended[i + 1] = M[i];
            M_extended[num_points - 1] = db;

            for (int i = 0; i < num_points - 1; ++i) {
                std::vector<double> coefficients(4);
                coefficients[0] = f[i];
                coefficients[1] = M_extended[i];
                coefficients[2] = (3 * first_divided_diff[i] - 2 * M_extended[i] - M_extended[i + 1]) / (t[i + 1] - t[i]);
                coefficients[3] = (M_extended[i] + M_extended[i + 1] - 2 * first_divided_diff[i]) / (t[i + 1] - t[i]) / (t[i + 1] - t[i]);

                Polynomial base({0}), power({1});
                for (int j = 0; j < 4; ++j) {
                    Polynomial temp({coefficients[j]});
                    temp = temp * power;
                    base = base + temp;
                    power = power * Polynomial({-t[i], 1});
                }
                polynomials[i] = base;
            }
        } else if (bc == SECOND_DERIVATIVE_FIXED || bc == NATURAL_SPLINE) {
            std::vector<double> M_extended(num_points);
            M_extended[0] = bc == NATURAL_SPLINE ? 0 : da;
            for (int i = 0; i < num_points - 2; ++i)
                M_extended[i + 1] = M[i];
            M_extended[num_points - 1] = bc == NATURAL_SPLINE ? 0 : db;

            for (int i = 0; i < num_points - 1; ++i) {
                std::vector<double> coefficients(4);
                coefficients[0] = f[i];
                coefficients[1] = first_divided_diff[i] - (M_extended[i + 1] + 2 * M_extended[i]) / 6 * (t[i + 1] - t[i]);
                coefficients[2] = M_extended[i] / 2;
                coefficients[3] = (M_extended[i + 1] - M_extended[i]) / 6 / (t[i + 1] - t[i]);

                Polynomial base({0}), power({1});
                for (int j = 0; j < 4; ++j) {
                    Polynomial temp({coefficients[j]});
                    temp = temp * power;
                    base = base + temp;
                    power = power * Polynomial({-t[i], 1});
                }
                polynomials[i] = base;
            }
        }

        return PiecewisePolynomial(polynomials, t);
    }
}


PPSpline::PPSpline(int dim, int order, const MathFunction &f, double a, double b, SplineBoundaryCondition bc, int num_intervals, double da, double db) : Spline(dim, order) {
    if (a >= b)
        throw "PPSpline: Invalid interval";
    if (num_intervals <= 1)
        throw "PPSpline: Invalid number of intervals";
    
    std::vector<double> time_points(num_intervals), function_values(num_intervals);
    for (int j = 0; j < num_intervals; ++j) {
        time_points[j] = a + (b - a) * j / (num_intervals - 1);
        function_values[j] = f.evaluate(time_points[j]);
    }
    this->segments = { compute_spline_segments(bc, function_values, time_points, da, db) };
}

PPSpline::PPSpline(int dim, int order, const MathFunction &f, const std::vector<double> &time_points, SplineBoundaryCondition bc, double da, double db) : Spline(dim, order) {
    if (dim != 1)
        throw "PPSpline: Dimension must be 1";
    if (time_points.size() <= 1)
        throw "PPSpline: Invalid number of intervals";
    
    std::vector<double> function_values(time_points.size());
    for (int j = 0; j < time_points.size(); ++j)
        function_values[j] = f.evaluate(time_points[j]);
    this->segments = { compute_spline_segments(bc, function_values, time_points, da, db) };
}

PPSpline::PPSpline(int dim, int order, const std::vector<std::vector<double>> &points, SplineBoundaryCondition bc) : Spline(dim, order) {
    if (points.size() <= 1)
        throw "PPSpline: Invalid number of intervals";
    if (points[0].size() != dim)
        throw "PPSpline: Dimension mismatch";
    
    std::vector<double> time_points;

    time_points.push_back(0);
    for (int i = 1; i < points.size(); ++i) {
        double distance = 0;
        static double total_distance = 0;
        for (int j = 0; j < points[i].size(); ++j)
            distance += (points[i][j] - points[i - 1][j]) * (points[i][j] - points[i - 1][j]);
        distance = sqrt(distance);
        total_distance += distance;
        time_points.push_back(total_distance);
    }

    for (int i = 0; i < dim; ++i) {
        std::vector<double> function_values;
        for (int j = 0; j < points.size(); ++j)
            function_values.push_back(points[j][i]);
        this->segments.push_back(compute_spline_segments(bc, function_values, time_points, 0.0, 0.0));
    }
}


////////////////////////////////////////////////////////////////////

/* BSpline 类 */

// 计算 B 样条基函数 B_i^k 的值
double BSpline :: evaluate_basis (int i, int k, double x) const {
    if (!k)
        return knot_vector[i - 1] < x && x <= knot_vector[i];
    else
        return (x - knot_vector[i - 1]) / (knot_vector[i + k - 1] - knot_vector[i - 1]) * evaluate_basis (i, k - 1, x) + (knot_vector[i + k] - x) / (knot_vector[i + k] - knot_vector[i]) * evaluate_basis (i + 1, k - 1, x);
}

// 计算 B 样条基函数 B_i^k 的导数
double BSpline :: evaluate_basis_derivative (int i, int k, double x) const {
    return k * evaluate_basis (i, k - 1, x) / (knot_vector[i + k - 1] - knot_vector[i - 1]) - k * evaluate_basis (i + 1, k - 1, x) / (knot_vector[i + k] - knot_vector[i]);
}

// 计算 B 样条基函数 B_i^k 的二阶导数
double BSpline :: evaluate_basis_second_derivative (int i, int k, double x) const {
    return k * evaluate_basis_derivative (i, k - 1, x) / (knot_vector[i + k - 1] - knot_vector[i - 1]) - k * evaluate_basis_derivative (i + 1, k - 1, x) / (knot_vector[i + k] - knot_vector[i]);
}
// 计算分段样条函数
PiecewisePolynomial BSpline::compute_spline_segments(SplineBoundaryCondition bc, const std::vector<double> &function_values, const std::vector<double> &time_points, double da, double db) {
    if (dimensions > 3)
        throw "BSpline: Dimension unsupported";
    this->knot_vector.clear();
    for (int i = spline_order; i >= 1; --i)
        this->knot_vector.push_back(time_points[0] - i);
    for (int i = 0; i < time_points.size(); ++i)
        this->knot_vector.push_back(time_points[i]);
    for (int i = 1; i <= spline_order; ++i)
        this->knot_vector.push_back(time_points[time_points.size() - 1] + i);
    
    int num_points = function_values.size();
    std::vector<std::vector<double>> matrix_A(num_points + spline_order - 1, std::vector<double>(num_points + spline_order - 1, 0.0));
    std::vector<double> vector_B(num_points + spline_order - 1, 0.0);
    
    if (spline_order == 3) {
        // 构造方程组
        for (int i = 0; i < num_points; ++i) {
            for (int j = 0; j < spline_order; ++j)
                matrix_A[i][i + j] = evaluate_basis(i + j + 1, spline_order, time_points[i]);
            vector_B[i] = function_values[i];
        }
        // 默认使用自然边界条件
        for (int i = 0; i < spline_order; ++i) {
            matrix_A[num_points][i] = evaluate_basis_second_derivative(i + 1, spline_order, time_points[0]);
            matrix_A[num_points + 1][num_points - 1 + i] = evaluate_basis_second_derivative(num_points + i, spline_order, time_points[num_points - 1]);
            vector_B[num_points] = 0;
            vector_B[num_points + 1] = 0;
        }
    } else if (spline_order == 2) {
        matrix_A[0][0] = evaluate_basis(1, spline_order, time_points[0]);
        matrix_A[0][1] = evaluate_basis(2, spline_order, time_points[0]);
        vector_B[0] = function_values[0];
        for (int i = 0; i < num_points - 1; ++i) {
            double mid_point = (time_points[i] + time_points[i + 1]) / 2;
            matrix_A[i + 1][i] = evaluate_basis(i + 1, spline_order, mid_point);
            matrix_A[i + 1][i + 1] = evaluate_basis(i + 2, spline_order, mid_point);
            matrix_A[i + 1][i + 2] = evaluate_basis(i + 3, spline_order, mid_point);
            vector_B[i + 1] = function_values[i + 1];
        }
        matrix_A[num_points][num_points - 1] = evaluate_basis(num_points, spline_order, time_points[num_points - 1]);
        matrix_A[num_points][num_points] = evaluate_basis(num_points + 1, spline_order, time_points[num_points - 1]);
        vector_B[num_points] = function_values[num_points - 1];
    }

    std::vector<double> coefficients;
    if (spline_order > 1)
        coefficients = lapack_solve(num_points + spline_order - 1, matrix_A, vector_B);
    else
        coefficients = function_values;

    std::vector<Polynomial> polynomials(num_points - 1);
    for (int i = 0; i < num_points - 1; ++i) {
        std::vector<double> x_values(spline_order + 1), y_values(spline_order + 1);
        for (int j = 0; j <= spline_order; ++j) {
            x_values[j] = time_points[i] + (time_points[i + 1] - time_points[i]) * j / spline_order;
            y_values[j] = 0;
            for (int k = 0; k < num_points + spline_order - 1; ++k)
                y_values[j] += coefficients[k] * evaluate_basis(k + 1, spline_order, x_values[j]);
        }
        polynomials[i] = Polynomial(x_values, y_values);
    }
    return PiecewisePolynomial(polynomials, time_points);
}

// 通过均匀节点构造 B 样条
BSpline::BSpline(int dim, int order, const MathFunction &f, double a, double b, int num_intervals) : Spline(dim, order) {
    if (a >= b)
        throw "BSpline: Invalid interval";
    if (num_intervals <= 1)
        throw "BSpline: Invalid number of intervals";
    
    std::vector<double> time_points(num_intervals), function_values(num_intervals);
    for (int j = 0; j < num_intervals; ++j) {
        time_points[j] = a + (b - a) * j / (num_intervals - 1);
        function_values[j] = f.evaluate(time_points[j]);
    }
    this->segments = { compute_spline_segments(NO_CONDITION, function_values, time_points, 0, 0) };
}

// 通过不均匀节点构造 B 样条
BSpline::BSpline(int dim, int order, const MathFunction &f, const std::vector<double> &time_points) : Spline(dim, order) {
    if (dim != 1)
        throw "BSpline: Dimension must be 1";
    if (time_points.size() <= 1)
        throw "BSpline: Invalid number of intervals";
    
    std::vector<double> function_values(time_points.size());
    for (int j = 0; j < time_points.size(); ++j)
        function_values[j] = f.evaluate(time_points[j]);
    this->segments = { compute_spline_segments(NO_CONDITION, function_values, time_points, 0, 0) };
}

// 通过散点拟合构造 B 样条
BSpline::BSpline(int dim, int order, const std::vector<std::vector<double>> &points) : Spline(dim, order) {
    if (points.size() <= 1)
        throw "BSpline: Invalid number of intervals";
    if (points[0].size() != dim)
        throw "BSpline: Dimension mismatch";
    
    std::vector<double> time_points;
    
    time_points.push_back(0);
    double total_distance = 0;
    for (int i = 1; i < points.size(); ++i) {
        double distance = 0;
        for (int j = 0; j < points[i].size(); ++j)
            distance += (points[i][j] - points[i - 1][j]) * (points[i][j] - points[i - 1][j]);
        distance = sqrt(distance);
        total_distance += distance;
        time_points.push_back(total_distance);
    }

    for (int i = 0; i < dim; ++i) {
        std::vector<double> function_values;
        for (int j = 0; j < points.size(); ++j)
            function_values.push_back(points[j][i]);
        this->segments.push_back(compute_spline_segments(NO_CONDITION, function_values, time_points, 0, 0));
    }
}