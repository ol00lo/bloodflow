MAINFILE=main

all:
	pdflatex $(MAINFILE)
bib:
	bibtex ${MAINFILE}
clean:
	rm -f *.aux *.backup *.toc *.bbl *.blg *.log *.out *.brf
cleanall: clean
	rm -f *.pdf
