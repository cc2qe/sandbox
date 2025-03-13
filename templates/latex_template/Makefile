PAPERFILENAME:=main
BIBFILE:=references.bib
PDF:=$(PAPERFILENAME).pdf
FLAGS=''

all: pdf

clean:
	latexmk -c -f $(PDF)
	rm $(PDF) $(PAPERFILENAME).bbl

pdf:
	xelatex $(FLAGS) $(PAPERFILENAME)
	if grep -q '\\citation{' $(PAPERFILENAME).aux; then bibtex $(PAPERFILENAME); fi
	xelatex $(FLAGS) $(PAPERFILENAME)
	xelatex $(FLAGS) $(PAPERFILENAME)


