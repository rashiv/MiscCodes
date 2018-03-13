#!/bin/bash 

latex main.tex
bibtex main.aux
latex main.tex
pdflatex main.tex
rm -f main.dvi
