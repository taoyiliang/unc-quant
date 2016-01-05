SRC = CV_2014

all: tidy

tidy: CV_2014.tex

	pdflatex $(SRC).tex
	bibtex $(SRC)
	pdflatex $(SRC).tex
	pdflatex $(SRC).tex

.PHONY: clean

clean:
	@rm -f *~ *.aux *.bbl *.blg *.log *.out *.toc *.lot *.lof $(SRC).pdf
