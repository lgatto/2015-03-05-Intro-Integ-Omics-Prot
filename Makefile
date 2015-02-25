%.md: %.Rmd
	Rscript -e "knitr::knit('$^')"

all: 02-caveats.md 03-mapping.md 04-transfer-learning.md
