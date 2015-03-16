##############################################################################
# MSc HPC
# Performance Programming
# Coursework
# Exam no. B063622
#
# Makefile for LaTex report.
#


##############################################################################
# Running collection of things to clean up.

CLEANFILES :=


##############################################################################
# Build products.

REPORT := B063622-PP.pdf

CLEANFILES += $(REPORT)


##############################################################################
# Build system.

MF := Makefile


##############################################################################
# Image files.


##############################################################################
# LaTeX Compiler.

LTX       := pdflatex
LTXFLAGS  := -halt-on-error -shell-escape


##############################################################################
# LaTeX source.

LTXSRC := $(REPORT:.pdf=.tex)


##############################################################################
# LaTeX-generated files.

TOC     := $(REPORT:.pdf=.toc)
AUX     := $(REPORT:.pdf=.aux)
LTXLOG  := $(REPORT:.pdf=.log)

CLEANFILES += $(TOC) $(AUX) $(LTXLOG)


##############################################################################
# Primary target rules.

default: report

.PHONY: report clean

report: $(REPORT)

clean:
	rm -f $(CLEANFILES)


##############################################################################
# Build report from LaTeX source.

# $(REPORT): $(MF) $(LTXSRC) $(PDFS)
$(REPORT): $(MF) $(LTXSRC)
	# run pdflatex twice to get references right
	$(LTX) $(LTXFLAGS) $(LTXSRC)
	$(LTX) $(LTXFLAGS) $(LTXSRC)
