# Chapter1 Programming Assignment
##### ==Author: Wanqi Chen (3220102895)==
## Introduction
This project is the first programming assignment of the course "Numerical Analysis" in ZJU. The project is to solve the following problems using C++:

## Document Structure
```
Programming/
└── Chapter1/
    ├── make.bat
    ├── Makefile
    ├── README.md
    ├── report.tex
    ├── picture/
    │   ├── ProblemB.png
    │   ├── ProblemC.png
    │   ├── ProblemD.png
    │   ├── ProblemE.png
    │   └── ProblemF.png
    ├── release/
    └── src/
        ├── EquationSolver.hpp
        ├── Function.hpp
        ├── ProblemB.cpp
        ├── ProblemC.cpp
        ├── ProblemD.cpp
        ├── ProblemE.cpp
        └── ProblemF.cpp

```


## Compile and Run
### Linux/macOS
Use the following command to **compile and run** the project:
```bash
make run
```

Use the following command to compile the TeX file:
```bash
make report
```

Use the following command to **clean** the project:
```bash
make clean
```

### Windows
In PowerShell, use the following command to **compile and run** the project:
```batch
.\make run     
```
In cmd, use the following command to **compile and run** the project:
```batch
make run
```

Use the following batch file to compile the TeX file:
```batch
.\make report
```

Use the following batch file to **clean** the project:
```batch
.\make clean
```

## Outputs and Results 
```
-----------PROBLEM B---------------

Testing f(x) = 1/x - tan(x) in [0, pi/2]:
Converged after 24 iterations.
Root: 0.860334
f(root) = 9.16679e-009
Root is correct.

--------------------------------

Testing f(x) = 1/x - 2^x in [0, 1]:
Converged after 24 iterations.
Root: 0.641186
f(root) = 4.43482e-008
Root is correct.

--------------------------------

Testing f(x) = 2^{-x} + e^x + 2cos(x) - 6 in [1, 3]:
Converged after 25 iterations.
Root: 1.82938
f(root) = -8.22642e-008
Root is correct.

--------------------------------

Testing f(x) = (x^3 + 4x^2 + 3x + 5) / (2x^3 - 9x^2 + 18x - 2) in [0, 4]:
Error: Maximum iterations reached without convergence.
Root: nan
f(root) = nan
Root is incorrect.

--------------------------------

Running ProblemC.exe...
-----------PROBLEM C---------------

Testing f(x) = x - tan(x) near 4.5:
Converged after 3 iterations.
Root near 4.5: 4.49341
f(root) = -3.69482e-012
Root is correct.

--------------------------------

Testing f(x) = x - tan(x) near 7.7:
Converged after 4 iterations.
Root near 7.7: 7.72525
f(root) = -4.51381e-011
Root is correct.

--------------------------------

Running ProblemD.exe...
-----------PROBLEM D---------------

Testing f(x) = sin(x/2) - 1 with initial value 0 and pi/2:
Converged after 15 iterations.
Root: 3.14093
f(root) = -5.45197e-008
Root is correct.

--------------------------------

Testing f(x) = sin(x/2) - 1 with another initial value 2pi and 3pi:
Converged after 0 iterations.
Root: 3.14159
f(root) = 0
Root is correct.

--------------------------------

Testing f(x) = e^x - tan(x) with initial value 1 and 1.4 :
Converged after 8 iterations.
Root: 1.30633
f(root) = -6.52234e-012
Root is correct.

--------------------------------

Testing f(x) = e^x - tan(x) with another initial value 2 and 3 :
Converged after 11 iterations.
Root: -6.28131
f(root) = -1.16272e-010
Root is correct.

--------------------------------

Testing f(x) = x^3 - 12x^2 + 3x + 1 with initial value 0 and -0.5:
Converged after 5 iterations.
Root: -0.188685
f(root) = -6.31696e-009
Root is correct.

--------------------------------

Testing f(x) = x^3 - 12x^2 + 3x + 1 with another initial value 1 and 2:
Converged after 6 iterations.
Root: 0.451543
f(root) = -1.51224e-008
Root is correct.

--------------------------------

Running ProblemE.exe...
-----------PROBLEM E---------------

Solving the equation with (Bisection Method)
Converged after 23 iterations.
Root: 0.166166 ft

--------------------------------

Solving the equation with (Newton Method)
Converged after 6 iterations.
Root: 0.166166 ft

--------------------------------

Solving the equation with (Secant Method)
Converged after 3 iterations.
Root: 0.166166 ft

--------------------------------

Running ProblemF.exe...
-----------PROBLEM F---------------

(a) Use the bisection method to find approximate values of alpha when l = 89 in., h = 49 in., D = 55 in., and beta_1 = 11.5 degrees:
Converged after 10 iterations.
Approximate root:32.959 degrees

--------------------------------

(b) Use Newton's method to find alpha with the initial guess 33 degrees when l,h,D, and beta_1 are as given above:
Converged after 2 iterations.
Root near 33 degrees: 32.9722
f(root) = 2.46072e-014
Root is correct.

--------------------------------

Use secant method with initial guesses far away from 33 to find alpha:

Converged after 9 iterations.
Root: -371.5
f(root) = -7.40393e-008
Root is correct.

--------------------------------
```




