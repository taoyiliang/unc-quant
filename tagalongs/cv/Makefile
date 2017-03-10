SRC = cv_2016

all: tidy

tidy: cv_2016.tex

	pdflatex $(SRC).tex
	# bibtex $(SRC)
	pdflatex $(SRC).tex
	pdflatex $(SRC).tex

.PHONY: clean

clean:
	@rm -f *~ *.aux *.bbl *.blg *.log *.out *.toc *.lot *.lof $(SRC).pdf
