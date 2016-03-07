%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

%.pdf: %.md
	Rscript -e "rmarkdown::render('$^')"

all: 00-intro.pdf 01-msprot-rnaseq.pdf
	pdftk A=00-intro.pdf B=01-msprot-rnaseq.pdf C=thetatut.pdf cat A1-5 B C A6 output slides.pdf

