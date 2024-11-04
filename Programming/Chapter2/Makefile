# Define variables
CXX = g++
TEX = pdflatex
SRC_DIR = ./src
BUILD_DIR = release
DATA_DIR = data
FIGURE_DIR = figure
PLOT_DIR = plot
PROBLEMS = ProblemB ProblemC ProblemD ProblemE ProblemF
TARGETS = $(patsubst %, $(BUILD_DIR)/%, $(PROBLEMS))
TEX_SRC = report.tex
TEX_OUT = report.pdf

# Phony targets
.PHONY: all clean run report

# Default target to compile and link all C++ programs
all: $(TARGETS)

# Compile each problem into its own executable
$(BUILD_DIR)/%: $(SRC_DIR)/%.cpp $(SRC_DIR)/interpolation.cpp | $(BUILD_DIR) $(DATA_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Create build directory and data directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(FIGURE_DIR):
	mkdir -p $(FIGURE_DIR)

$(DATA_DIR):
	mkdir -p $(DATA_DIR)

# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(DATA_DIR)/*.txt $(FIGURE_DIR)/*.png $(TEX_OUT) *.log *.aux *.out

# Run all compiled programs sequentially
run: all $(DATA_DIR) $(FIGURE_DIR)
	@echo "Running compiled programs..."
	@for target in $(TARGETS); do \
		echo "Running $$target..."; \
		./$$target; \
	done
	@echo "Running Python plot scripts..."
	python3 $(PLOT_DIR)/plot_ProblemB.py
	python3 $(PLOT_DIR)/plot_ProblemC.py
	python3 $(PLOT_DIR)/plot_ProblemD.py
	python3 $(PLOT_DIR)/plot_ProblemE.py
	python3 $(PLOT_DIR)/plot_ProblemF.py
	@echo "Process completed. Results saved in $(FIGURE_DIR)."

# Compile LaTeX report
report:
	$(TEX) $(TEX_SRC)
	$(TEX) $(TEX_SRC) # Run twice to ensure cross-references are correct
