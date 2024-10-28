# Define directories and files
SRC_DIR := src
DATA_DIR := data
RELEASE_DIR := release
FIGURE_DIR := figure
PLOT_DIR := plot
REPORT_FILE := report.tex

# Default target
all: run

# Run the program
run:
	@echo Cleaning old build files...
	@rm -f $(RELEASE_DIR)/*.exe
	@rm -f $(DATA_DIR)/*.txt
	@rm -f $(FIGURE_DIR)/*.png
	@rm -f *.aux
	@rm -f *.synctex.gz
	@rm -f *.log
	@rm -f *.out
	@rm -f *.pdf
	@echo Compiling C++ files...
	@for f in $(SRC_DIR)/Problem*.cpp; do \
		echo Compiling $$f... ; \
		g++ $$f src/interpolation.cpp -o $(RELEASE_DIR)/$$(basename $$f .cpp).exe ; \
	done
	@echo Running compiled programs...
	@for f in $(RELEASE_DIR)/*.exe; do \
		./$$f ; \
	done
	@echo Running Python plot scripts...
	@python $(PLOT_DIR)/plot_problemB.py
	@python $(PLOT_DIR)/plot_problemC.py
	@python $(PLOT_DIR)/plot_problemD.py
	@python $(PLOT_DIR)/plot_problemE.py
	@python $(PLOT_DIR)/plot_problemF.py
	@echo Process completed. Results saved in $(FIGURE_DIR).

# Clean up intermediate files
clean:
	@echo Cleaning up...
	@rm -f $(RELEASE_DIR)/*.exe
	@rm -f $(DATA_DIR)/*.txt
	@rm -f $(FIGURE_DIR)/*.png
	@echo Clean completed.

# Compile LaTeX report
report:
	@echo Compiling LaTeX report...
	@pdflatex $(REPORT_FILE)
	@echo Report compilation completed.
	@rm -f *.aux
	@rm -f *.log
	@rm -f *.out

# Phony targets
.PHONY: all run clean report