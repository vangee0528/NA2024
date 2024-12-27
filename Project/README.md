## 项目文件结构
```
Project/
├── Makefile
├── README.md
├── compare.py
├── docs                    // 设计文档及项目报告
│   ├── design.md
│   ├── design.pdf
│   ├── report.md
│   ├── report.tex
│   └── 程序设计.png
├── figure                  // 使用python绘制的图片
│   ├── Sphere
│   ├── check
│   ├── jsontest
│   ├── problemA
│   ├── problemB
│   ├── problemC
│   ├── problemD
│   ├── problemE
│   └── test
├── logs                    // 编译日志
│   └── log.txt
├── main.cpp                // 测试要求的功能
├── output                  // cpp程序生成的输出
│   ├── Sphere
│   ├── check
│   ├── jsontest
│   ├── problemA
│   ├── problemB
│   ├── problemC
│   ├── problemD
│   ├── problemE
│   └── test
├── plot.py                 // 绘图程序
├── release                 // 编译结果
└── src
    ├── Json                // json解析测试
    │   ├── input
    │   │   ├── Example.json.txt
    │   │   ├── input_bspline.json
    │   │   └── input_ppspline.json
    │   ├── plot.py
    │   └── test.cpp    
    ├── Packages            // 项目依赖的头文件和实现
    │   ├── json.cpp
    │   ├── json.h
    │   ├── lapack.cpp
    │   ├── lapack.h
    │   ├── spline.cpp
    │   └── spline.h
    ├── ProblemA            // 问题A - E的测试程序
    │   ├── ProblemA.cpp
    │   └── plotA.py
    ├── ProblemB
    │   ├── ProblemB.cpp
    │   └── plotB.py
    ├── ProblemC
    │   ├── ProblemC.cpp
    │   └── plotC.py
    ├── ProblemD
    │   ├── ProblemD.cpp
    │   └── plotD.py
    ├── ProblemE
    │   ├── ProblemE.cpp
    │   └── plotE.py
    ├── Sphere              // 球面样条拟合测试
    │   ├── Sphere.cpp
    │   └── plotSphere.py
    └── test                // 更多函数测试
        ├── main.cpp
        └── plot.py
```

## 编译运行
```bash
make clean
make run
```

## 清理编译文件
```bash
make clean
make cleanlog
```

