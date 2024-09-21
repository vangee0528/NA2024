# Define variables
CXX = g++
TEX = pdflatex
SRC_DIR = ./src
BUILD_DIR = release
PROBLEMS = ProblemB ProblemC ProblemD ProblemE ProblemF
TARGETS = $(patsubst %, $(BUILD_DIR)/%, $(PROBLEMS))
TEX_SRC = report.tex
TEX_OUT = report.pdf

# Phony targets
.PHONY: all clean run report

# Default target to compile and link all C++ programs
all: $(TARGETS)

# Compile each problem into its own executable
$(BUILD_DIR)/%: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $< -o $@

# Create build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(TEX_OUT) *.log *.aux

# Run all compiled programs sequentially
run: all
	@for target in $(TARGETS); do \
		echo "Running $$target...\n"; \
		$$target; \
	done

# Compile LaTeX report
report:
	$(TEX) $(TEX_SRC)
	$(TEX) $(TEX_SRC) # Run twice to ensure cross-references are correct
