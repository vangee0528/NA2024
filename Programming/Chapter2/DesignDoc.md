## 概述
本设计文档介绍了实现插值方法和 Bézier 曲线的设计思路及实现方法。本文档主要包括以下内容：
1. 插值类 (`NewtonInterpolator` 和 `HermiteInterpolator`) 设计。
2. Bézier 曲线生成类设计。
3. 各类的成员变量、成员函数说明及使用流程。

## 1. 设计思路

本项目包含三类算法：Newton 插值法、Hermite 插值法以及 Bézier 曲线生成。每种算法通过不同的类实现如下功能：

- **Newton 插值法与 Hermite 插值法**：分别由 `NewtonInterpolator` 和 `HermiteInterpolator` 类实现，通过输入数据构建实例并计算插值多项式。
- **Bézier 曲线**：由 `Bezier` 类实现，通过输入控制点生成 Bézier 曲线并输出参数方程到指定文件。

为了得到 Bézier 曲线的控制点，本设计参考了讲义中的算法 2.74，从目标曲线中选取标记点生成控制点。标记点顺序选取均匀分布以保证拟合效果，并按顺时针或逆时针顺序保证 Bézier 曲线连续性。

## 2. 类与其成员函数介绍

### 2.1 `Polynomial` 类

**描述**：表示一个多项式对象，支持评估与求导操作。

- **成员变量**：
  - `int degree`：多项式的次数。
  - `std::vector<double> coeffs`：多项式的系数。

- **成员函数**：
  - `double evaluate(double x) const`  
    - **描述**：计算多项式在 \( x \) 处的值。
  - `Polynomial derivative() const`  
    - **描述**：返回一个新的多项式，表示当前多项式的导数。

### 2.2 `NewtonInterpolator` 类

**描述**：实现 Newton 插值法。

- **成员变量**：
  - `std::vector<double> x`：插值点的自变量值。
  - `std::vector<double> f`：插值点的因变量值。

- **成员函数**：
  - `double dividedDifference(int i, int j) const`  
    - **描述**：递归计算 i 和 j 点之间的分差。
  - `Polynomial interpolate() const`  
    - **描述**：基于 Newton 插值法返回插值多项式。

### 2.3 `HermiteInterpolator` 类

**描述**：实现 Hermite 插值法。

- **成员变量**：
  - `std::vector<double> x`：插值点的自变量值。
  - `std::vector<double> f`：插值点的因变量值。
  - `std::vector<double> df`：插值点的导数值。

- **成员函数**：
  - `void computeDividedDifferences()`  
    - **描述**：计算 Hermite 插值所需的分差。
  - `Polynomial interpolate() const`  
    - **描述**：基于 Hermite 插值法生成插值多项式。

### 2.4 `Point` 结构

**描述**：表示二维坐标点。

- **成员变量**：
  - `double x_val`：X 坐标。
  - `double y_val`：Y 坐标。

- **成员函数**：
  - `Point operator+(const Point& other) const` 等  
    - **描述**：重载运算符，用于点的加减乘除。

### 2.5 `Bezier` 类

**描述**：生成 Bézier 曲线，输出至文件。

- **成员变量**：
  - `std::vector<Point> control_points`：贝塞尔曲线的控制点。

- **成员函数**：
  - `void printOut() const`  
    - **描述**：打印 Bézier 曲线参数方程。
  - `void FileOut(std::ofstream& outfile, double x) const`  
    - **描述**：将 Bézier 曲线的参数方程输出到文件。

## 3. 算法流程

### 3.1 Newton 插值法

1. **输入数据**：接收自变量数组 `x[]` 和因变量数组 `f[]`。
2. **构建类**：使用输入数据构建 `NewtonInterpolator` 实例。
3. **计算差商表**：在 `NewtonInterpolator` 类的构造函数中，调用 `dividedDifferences` 方法计算系数。
4. **返回插值多项式**：调用 `interpolate` 方法返回插值多项式 `Polynomial` 对象。

### 3.2 Hermite 插值法

1. **输入数据**：接收自变量数组 `x[]`，因变量数组 `f[]` 和导数数组 `df[]`。
2. **构建类**：使用输入数据构建 `HermiteInterpolator` 实例。
3. **计算差商表**：在 `HermiteInterpolator` 类的构造函数中，调用 `computeDividedDifferences` 计算差商。
4. **返回插值多项式**：调用 `interpolate` 方法返回 Hermite 插值多项式。

### 3.3 Bézier 曲线生成

1. **输入数据**：接收控制点数组 `control_points[]`。
2. **构建类**：使用控制点构建 `Bezier` 类的实例。
3. **输出参数方程**：调用 `FileOut` 方法将 Bézier 曲线的参数方程输出至文件。