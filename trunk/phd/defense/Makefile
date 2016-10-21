SRC = defense_present

all: tidy

tidy: defense_present.tex

	pdflatex $(SRC).tex
	pdflatex $(SRC).tex
	# pdflatex $(SRC).tex

.PHONY: clean

clean:
	@rm -f *~ *.aux *.bbl *.blg *.log *.out *.toc *.lot *.lof $(SRC).pdf
