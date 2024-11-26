# 定义变量
CXX = g++
LDFLAGS = -llapack -lblas
SRC_DIR = src
RELEASE_DIR = release
OUTPUT_DIR = output
PLOT_DIR = src

# 定义目标文件
TARGETS = $(RELEASE_DIR)/ProblemA $(RELEASE_DIR)/ProblemB $(RELEASE_DIR)/ProblemC

# 定义源文件
SOURCES_A = $(SRC_DIR)/ProblemA/ProblemA.cpp $(SRC_DIR)/Packages/spline.cpp $(SRC_DIR)/Packages/lapack.cpp
SOURCES_B = $(SRC_DIR)/ProblemB/ProblemB.cpp $(SRC_DIR)/Packages/spline.cpp $(SRC_DIR)/Packages/lapack.cpp
SOURCES_C = $(SRC_DIR)/ProblemC/ProblemC.cpp $(SRC_DIR)/Packages/spline.cpp $(SRC_DIR)/Packages/lapack.cpp

# 定义依赖文件
DEPS = $(wildcard $(SRC_DIR)/Packages/*.h)

# 默认目标
all: clean $(TARGETS) run plot

# 清理目标
clean:
	rm -rf $(RELEASE_DIR)/*
	rm -rf $(OUTPUT_DIR)/*

# 编译目标
$(RELEASE_DIR)/ProblemA: $(SOURCES_A) $(DEPS)
	$(CXX)  -o $@ $(SOURCES_A) $(LDFLAGS)

$(RELEASE_DIR)/ProblemB: $(SOURCES_B) $(DEPS)
	$(CXX)  -o $@ $(SOURCES_B) $(LDFLAGS)

# $(RELEASE_DIR)/ProblemC: $(SOURCES_C) $(DEPS)
# 	$(CXX)  -o $@ $(SOURCES_C) $(LDFLAGS)

# 运行目标
run: $(TARGETS)
	@echo "Running ProblemA..."
	./$(RELEASE_DIR)/ProblemA
	@echo "Running ProblemB..."	
	./$(RELEASE_DIR)/ProblemB
	@echo "Running ProblemC..."
	./$(RELEASE_DIR)/ProblemC

# 绘图目标
plot:
	@echo "Running plotA.py..."
	python3 $(PLOT_DIR)/ProblemA/plotA.py
	@echo "Running plotB.py..."
	python3 $(PLOT_DIR)/ProblemB/plotB.py

# 提供一个输入通道
.PHONY: all clean run plot