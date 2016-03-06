%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

01-msprot-rnaseq.pdf: 01-msprot-rnaseq.Rmd
	Rscript -e "rmarkdown::render('01-msprot-rnaseq.Rmd', output_file = 'mapping-slides.pdf', output_format = 'beamer_presentation')"

