%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

all: caveats.md mapping.md mapping2.Rmd transfer-learning.md
