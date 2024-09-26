# Define variables
TEX = pdflatex
TEX_SRC = answer.tex
TEX_OUT = answer.pdf

# Phony targets
.PHONY: all clean  report

# Default target to compile and link all C++ programs
all: $(TARGETS)


# Clean up build files
clean:
	rm -rf $(BUILD_DIR) $(TEX_OUT) *.log *.aux


# Compile LaTeX report
report:
	$(TEX) $(TEX_SRC)
	$(TEX) $(TEX_SRC) # Run twice to ensure cross-references are correct
