02-caveats.md: 02-caveats.Rmd
	Rscript -e "knitr::knit('02-caveats.Rmd')"

all: 02-caveats.md
