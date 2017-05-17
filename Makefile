

RMDS=index.Rmd \
     AdvGRanges_Rtracklayer_Rsamtools.Rmd \
     data_structures.Rmd \
     strings_and_ranges.Rmd

HTMLS=$(patsubst %.Rmd,%.html,$(RMDS))


all : $(HTMLS)

%.html : %.Rmd
	Rscript -e 'rmarkdown::render("$<", "all")'

