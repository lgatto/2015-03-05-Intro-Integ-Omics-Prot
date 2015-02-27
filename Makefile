%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

all: caveats.md mapping.md mapping2.md mapping-rnaseq.md transfer-learning.md
