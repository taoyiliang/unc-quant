SRC = 201606_transaction

all: tidy

tidy: 201606_transaction.tex

	pdflatex $(SRC).tex
	bibtex $(SRC)
	pdflatex $(SRC).tex
	pdflatex $(SRC).tex

.PHONY: clean

clean:
	@rm -f *~ *.aux *.bbl *.blg *.log *.out *.toc *.lot *.lof $(SRC).pdf
