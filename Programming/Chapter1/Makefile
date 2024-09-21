# Define variables
CXX = g++
CXXFLAGS = -Wall -Werror -std=c++11
TEX = pdflatex
SRC_DIR = .
BUILD_DIR = release
TARGET = $(BUILD_DIR)/all_programs
TEX_SRC = report.tex
TEX_OUT = report.pdf
CPP_FILES = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES = $(CPP_FILES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Phony targets
.PHONY: all clean run report

# Default target to compile and link all C++ programs
all: $(TARGET)

# Compile object files and link them
$(TARGET): $(OBJ_FILES) | $(BUILD_DIR)
	$(CXX) $(OBJ_FILES) -o $(TARGET)

# Compile C++ source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(TEX_OUT)

# Run all compiled programs
run: all
	$(TARGET)

# Compile LaTeX report
report:
	$(TEX) $(TEX_SRC)
