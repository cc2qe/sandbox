PAPERFILENAME:=main
BIBFILE:=references.bib
PDF:=$(PAPERFILENAME).pdf
FLAGS=''

all: pdf

clean:
	latexmk -c -f $(PDF)
	rm -f $(PDF) $(PAPERFILENAME).bbl $(PAPERFILENAME).bcf \
	      $(PAPERFILENAME).blg $(PAPERFILENAME).run.xml \
	      $(PAPERFILENAME).aux $(PAPERFILENAME).log

pdf:
	xelatex $(FLAGS) $(PAPERFILENAME)
	biber $(PAPERFILENAME)
	xelatex $(FLAGS) $(PAPERFILENAME)
	xelatex $(FLAGS) $(PAPERFILENAME)

