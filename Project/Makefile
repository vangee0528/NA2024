# 关闭输出命令
.SILENT:
all:
	if [ -f logs/log.txt ]; then mv logs/log.txt logs/log_`date +%m_%d_%H:%M`.txt; fi
	echo "Make Time: `date`" > logs/log.txt
	g++ src/ProblemA/ProblemA.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/ProblemA.o -llapack -lblas
	release/ProblemA.o
	python3 src/ProblemA/plotA.py
	g++ src/ProblemB/ProblemB.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/ProblemB.o -llapack -lblas
	release/ProblemB.o
	python3 src/ProblemB/plotB.py
	g++ src/ProblemC/ProblemC.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/ProblemC.o -llapack -lblas
	release/ProblemC.o
	python3 src/ProblemC/plotC.py
	g++ src/ProblemD/ProblemD.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/ProblemD.o -llapack -lblas
	release/ProblemD.o
	python3 src/ProblemD/plotD.py
	g++ src/ProblemE/ProblemE.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/ProblemE.o -llapack -lblas
	release/ProblemE.o
	python3 src/ProblemE/plotE.py
	g++ src/Sphere/Sphere.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/Sphere.o -llapack -lblas
	release/Sphere.o
	python3 src/Sphere/plotSphere.py
	g++ main.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/main.o -llapack -lblas
	release/main.o
	g++ src/Json/test.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/test.o -llapack -lblas
	release/test.o
	python3 src/Json/plot.py
	g++ src/test/main.cpp src/Packages/spline.cpp src/Packages/lapack.cpp src/Packages/json.cpp -o release/test1.o -llapack -lblas
	release/test1.o
	python3 src/test/plot.py
	pdflatex -output-directory=docs docs/report.tex
clean:
# 向logs/log.txt中输出当前时间
	echo "Clean Time: `date`" >> logs/log.txt
	rm -f release/*.o
	rm -f figure/problemA/*.png
	rm -f output/problemA/*.txt
	rm -f figure/problemB/*.png
	rm -f output/problemB/*.txt
	rm -f figure/problemC/*.png
	rm -f output/problemC/*.txt
	rm -f figure/problemD/*.png
	rm -f output/problemD/*.txt
	rm -f figure/problemE/*.png
	rm -f output/problemE/*.txt
	rm -f figure/Sphere/*.png
	rm -f output/Sphere/*.txt
	rm -f output/check/*.txt
	rm -f figure/check/*.png
	rm -f figure/jsontest/*.png
	rm -f output/jsontest/*.txt
	rm -f output/test/*.txt
	rm -f figure/test/*.png

cleanlog:
	rm -f logs/*.txt