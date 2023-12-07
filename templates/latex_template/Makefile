PAPERFILENAME:=main
BIBFILE:=references.bib
PDF:=$(PAPERFILENAME).pdf

all: pdf

clean:
	latexmk -c -f $(PDF)
	rm $(PDF) $(PAPERFILENAME).bbl

pdf:
	xelatex $(PAPERFILENAME)
	bibtex $(PAPERFILENAME)
	xelatex $(PAPERFILENAME)
	xelatex $(PAPERFILENAME)


