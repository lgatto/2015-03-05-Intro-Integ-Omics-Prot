%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

all: caveats.md mapping.md transfer-learning.md
