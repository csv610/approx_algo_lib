# Makefile for Approximation Algorithms with Python - LaTeX Book

BOOK = book
TEX = $(BOOK).tex
PDF = $(BOOK).pdf

# Main targets
all: $(PDF)

$(PDF): $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)
	bibtex $(BOOK) || true
	pdflatex -interaction=nonstopmode $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)

# Quick compile (no bibtex)
quick: $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)

# Clean auxiliary files
clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot *.idx *.ind *.ilg *.synctex.gz *.fdb_latexmk *.fls

# Clean everything including PDF
distclean: clean
	rm -f $(PDF)

# Run Python demos
demo:
	python3 main.py

# Run specific chapter demo
ch1:
	python3 -m ch01_intro

ch2:
	python3 -m ch02_set_cover

ch3:
	python3 -m ch03_steiner_tsp

ch4:
	python3 -m ch04_multiway_kcut

ch5:
	python3 -m ch05_kcenter

ch12:
	python3 -m ch12_14_lp_algorithms

ch13:
	python3 -m ch13_knapsack

ch24:
	python3 -m ch24_facility_location

.PHONY: all quick clean distclean demo ch1 ch2 ch3 ch4 ch5 ch12 ch13 ch24