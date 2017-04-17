# Makefile for creating an ATLAS LaTeX document
#------------------------------------------------------------------------------
# By default makes mydocument.pdf using target run_pdflatex
# Replace mydocument with your main filename or add another target set
# Replace BIBTEX = bibtex with BIBTEX = biber if you use biblatex and biber instead of bibtex
# Use "make clean" to cleanup.
# Use "make cleanpdf to delete $(BASENAME).pdf.
# "make cleanall" also deletes the PDF file $(BASENAME).pdf.

TEXLIVE  = 2015
LATEX    = latex
PDFLATEX = pdflatex
# BIBTEX   = bibtex
BIBTEX = biber
DVIPS    = dvips
DVIPDF   = dvipdf

BASENAME = susy2L2015

# Default target - make mydocument.pdf with pdflatex
default: run_pdflatex

lt:lt.pdf

.PHONY: new newtexmf draftcover preprintcover auxmat clean cleanpdf help

new:
	if [ $(TEXLIVE) -ge 2007 -a $(TEXLIVE) -lt 2100 ]; then \
	  sed s/atlas-document/$(BASENAME)/ template/atlas-document.tex | \
	    sed 's/texlive=20[0-9][0-9]/texlive=$(TEXLIVE)/' >$(BASENAME).tex; \
	else \
	  echo "Invalid value for TEXLIVE: $(TEXLIVE)"; \
	  sed s/atlas-document/$(BASENAME)/ template/atlas-document.tex >$(BASENAME).tex; \
	fi
	cp template/atlas-document-metadata.tex $(BASENAME)-metadata.tex
	touch $(BASENAME).bib
	touch $(BASENAME)-defs.sty

newtexmf:
	if [ $(TEXLIVE) -ge 2007 -a $(TEXLIVE) -lt 2100 ]; then \
	  sed s/atlas-document/$(BASENAME)/ template/atlas-document-texmf.tex | \
	    sed 's/texlive=20[0-9][0-9]/texlive=$(TEXLIVE)/' >$(BASENAME).tex; \
	else \
	  echo "Invalid value for TEXLIVE: $(TEXLIVE)"; \
	  sed s/atlas-document/$(BASENAME)/ template/atlas-document.tex >$(BASENAME).tex; \
	fi
	cp template/atlas-document-metadata.tex $(BASENAME)-metadata.tex
	touch $(BASENAME).bib
	touch $(BASENAME)-defs.sty

draftcover:
	if [ $(TEXLIVE) -ge 2007 -a $(TEXLIVE) -lt 2100 ]; then \
	  sed 's/texlive=20[0-9][0-9]/texlive=$(TEXLIVE)/' template/atlas-draft-cover.tex \
	    >$(BASENAME)-draft-cover.tex; \
	else \
	  echo "Invalid value for TEXLIVE: $(TEXLIVE)"; \
	  cp  template/$(BASENAME)-draft-cover.tex $(BASENAME)-draft-cover.tex; \
	fi
	
preprintcover:
	sed 's/texlive=20[0-9][0-9]/texlive=$(TEXLIVE)/' template/atlas-preprint-cover.tex >$(BASENAME)-preprint-cover.tex
	#cp template/atlas-preprint-cover.tex $(BASENAME)-preprint-cover.tex
	
auxmat:
	sed s/atlas-document/$(BASENAME)/ template/atlas-auxmat.tex | \
	sed 's/texlive=20[0-9][0-9]/texlive=$(TEXLIVE)/' >$(BASENAME)-auxmat.tex

run_pdflatex: $(BASENAME).pdf
	@echo "Made $<"

run_latex: $(BASENAME).dvi
	$(DVIPDF) -sPAPERSIZE=a4 -dPDFSETTINGS=/prepress ${BASENAME}
	@echo "Made $<"

# Standard Latex targets
# %.pdf:	%.tex *.bib
%.pdf:	%.tex *.bib chapters/*.tex
	$(PDFLATEX) $<
	-$(BIBTEX)  $(basename $<)
	$(PDFLATEX) $<
	$(PDFLATEX) $<

%.dvi:	%.tex 
	$(LATEX)    $<
	-$(BIBTEX)  $(basename $<)
	$(LATEX)    $<
	$(LATEX)    $<

%.bbl:	%.tex *.bib
	$(LATEX) $*
	$(BIBTEX) $*
	
help:
	@echo "To create a new document give the commands:"
	@echo "make new [BASENAME=mydocument] [TEXLIVE=YYYY]"
	@echo "make"
	@echo ""
	@echo "If atlaslatex is installed centrally, e.g. in ~/texmf:"
	@echo "make newtexmf [BASENAME=mydocument] [TEXLIVE=YYYY]"
	@echo ""
	@echo "If you need a standalone draft cover give the commands:"
	@echo "make draftcover [BASENAME=mydocument] [TEXLIVE=YYYY]"
	@echo "pdflatex mydocument-draft-cover"
	@echo ""
	@echo "If you need a standalone preprint cover give the commands:"
	@echo "make preprintcover [BASENAME=mydocument] [TEXLIVE=YYYY]"
	@echo "pdflatex mydocument-preprint-cover"
	@echo ""
	@echo "If you need a document for auxiliary material give the commands:"
	@echo "make auxmat [BASENAME=mydocument] [TEXLIVE=YYYY]"
	@echo "pdflatex mydocument-auxmat"
	@echo ""
	@echo "make clean    to clean auxiliary files (not output PDF)"
	@echo "make cleanpdf to clean output PDF files"
	@echo "make cleanall to clean all files"

%.ps:	%.dvi
	$(DVIPS) $< -o $@

clean:
	-rm *.dvi *.toc *.aux *.log *.out \
		*.bbl *.blg *.brf *.bcf *-blx.bib *.run.xml \
		*.cb *.ind *.idx *.ilg *.inx \
		*.synctex.gz *~ ~* spellTmp 
	
cleanpdf:
	-rm $(BASENAME).pdf 
	-rm $(BASENAME)-draft-cover.pdf $(BASENAME)-preprint-cover.pdf
	-rm $(BASENAME)-auxmat.pdf

cleanall: clean cleanpdf
