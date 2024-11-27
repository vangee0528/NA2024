# 关闭输出命令
.SILENT:
all:
	g++ src/ProblemA/ProblemA.cpp src/Packages/spline.cpp src/Packages/lapack.cpp -o release/ProblemA.o -llapack -lblas
	release/ProblemA.o
	python3 src/ProblemA/plotA.py
	g++ src/ProblemC/ProblemC.cpp src/Packages/spline.cpp src/Packages/lapack.cpp -o release/ProblemC.o -llapack -lblas
	release/ProblemC.o
	python3 src/ProblemC/plotC.py
clean:
	rm -f release/*.o
	rm -f figure/problemA/*.png
	rm -f output/problemA/*.txt