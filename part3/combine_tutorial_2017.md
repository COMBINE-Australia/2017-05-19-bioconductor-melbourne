Melbourne Bioconductor - Part 3
================
Anna Qualieri
15th May 2017

-   [Chunks options](#chunks-options)
-   [Install packages](#install-packages)
-   [Load packages](#load-packages)
-   [Advanced `GenomicRanges`](#advanced-genomicranges)
    -   [Overlaps between two `GRanges` objects](#overlaps-between-two-granges-objects)
        -   [`findOverlaps()`](#findoverlaps)
        -   [`countOverlaps()`](#countoverlaps)
    -   [Nearest-methods in `GenomicRanges`](#nearest-methods-in-genomicranges)
        -   [`nearest()`](#nearest)
        -   [`distance()`](#distance)
        -   [`distanceToNearest()`](#distancetonearest)
    -   [`GRangesList`](#grangeslist)
    -   [Subsetting and looping over *GRanges* list](#subsetting-and-looping-over-granges-list)
        -   [Challenge](#challenge)
        -   [Solutions](#solutions)
-   [`Rtracklayer`](#rtracklayer)
    -   [GTF and GFF](#gtf-and-gff)
        -   [Challenge](#challenge-1)
        -   [Solution](#solution)
    -   [Wiggle (WIG) and bigWIG file formats for graphing tracks](#wiggle-wig-and-bigwig-file-formats-for-graphing-tracks)
    -   [`import` only a region of the bigWig file](#import-only-a-region-of-the-bigwig-file)
-   [`liftOver` from different genome releases](#liftover-from-different-genome-releases)
-   [`Rsamtools`](#rsamtools)
    -   [Input BAM into R](#input-bam-into-r)
    -   [Compute number of reads that falls within 100bp bins](#compute-number-of-reads-that-falls-within-100bp-bins)
        -   [Challenge](#challenge-2)
        -   [Solution](#solution-1)

Chunks options
==============

``` r
library(knitr)
knitr::opts_chunk$set(fig.width=8, 
    fig.height=6, echo=T, 
    warning=FALSE, message=FALSE,
    prompt=T,tidy=T,tidy.opts=list(width.cutoff=50),
    include=TRUE,cache=FALSE)
# To output a pdf or html document simply substitute github_document with pdf_document or html_document in the header of the file

# Working directory
dir <- "~/Documents/varie/Combine_tutorial"

# Extract R code
# this_rmd <- file.path(dir,"combine_tutorial_2017.Rmd")
# purl(this_rmd)
# purl(this_rmd, output = file.path(dir,"combine_tutorial_2017.R"), documentation = 2)
```

Install packages
================

``` r
> # Bioconductor packages
> source("https://bioconductor.org/biocLite.R")
> biocLite(c("GenomicRanges", "rtracklayer", "Rsamtools", 
+     "GenomicAlignments", "TxDb.Hsapiens.UCSC.hg19.knownGene"))
> 
> # CRAN packages
> install.packages("ggplot2")
```

Load packages
=============

``` r
> library("ggplot2")
> library("GenomicRanges")
> library("Rsamtools")
> library("GenomicAlignments")
> library("rtracklayer")
> library("TxDb.Hsapiens.UCSC.hg19.knownGene")
```

Advanced `GenomicRanges`
========================

Imagine you have done a ChIP-Seq experiment on *Sample 1* and your output is a set of ranges. Let's start by manually creating this very simple `GRanges` objects with ranges on *chr1* and *chr2*.

``` r
> # Sample 1
> gr_S1 <- GenomicRanges::GRanges(seqnames = Rle(c("chr1", 
+     "chr2"), c(3, 2)), ranges = IRanges::IRanges(start = c(5, 
+     8, 20, 8, 18), end = c(11, 15, 26, 16, 21), names = c(paste("Peak_", 
+     1:5, sep = ""))), strand = S4Vectors::Rle(strand(c("*")), 
+     c(5)))
> 
> gr_S1
```

    ## GRanges object with 5 ranges and 0 metadata columns:
    ##          seqnames    ranges strand
    ##             <Rle> <IRanges>  <Rle>
    ##   Peak_1     chr1  [ 5, 11]      *
    ##   Peak_2     chr1  [ 8, 15]      *
    ##   Peak_3     chr1  [20, 26]      *
    ##   Peak_4     chr2  [ 8, 16]      *
    ##   Peak_5     chr2  [18, 21]      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

A very common analysis to perform is to evaluate to what extent and where your ChIP-Seq peaks overlap with some **features**, such as genes, exons, other ChIP-Seq peaks, etc... As an example we create a simple gene annotation to use with the ranges created above.

``` r
> # Gene annotation
> genes <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(c("chr1", 
+     "chr2"), c(2, 2)), ranges = IRanges::IRanges(start = c(7, 
+     17, 7, 23), end = c(15, 23, 14, 26), names = c(paste("Gene_", 
+     1:4, sep = ""))), strand = S4Vectors::Rle(strand(c("+", 
+     "+")), c(2, 2)))
> 
> genes
```

    ## GRanges object with 4 ranges and 0 metadata columns:
    ##          seqnames    ranges strand
    ##             <Rle> <IRanges>  <Rle>
    ##   Gene_1     chr1  [ 7, 15]      +
    ##   Gene_2     chr1  [17, 23]      +
    ##   Gene_3     chr2  [ 7, 14]      +
    ##   Gene_4     chr2  [23, 26]      +
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

Figure 1 is a simple way of plotting ranges stored into *Granges* object using *geom\_rect* from *ggplot2*. More advanced plotting options are available but this simple strategy is sufficient for now and also, I love ggplot!

``` r
> gr <- as.data.frame(rbind(as.data.frame(gr_S1), as.data.frame(genes)))
> gr$rangeID <- c(names(gr_S1), names(genes))
> gr$Sample <- c(rep("Peaks", length(gr_S1)), rep("Genes", 
+     length(genes)))
> ggplot(data = gr, aes(xmin = start, xmax = end, ymin = 0, 
+     ymax = 1)) + geom_rect(aes(fill = rangeID), alpha = 0.4) + 
+     facet_wrap(~Sample + seqnames) + theme_bw() + labs(x = "Genomic position", 
+     y = "Ranges") + theme(axis.ticks.y = element_blank(), 
+     axis.text.y = element_blank())
```

![Graphic representation of the GRanges objects created above.](combine_tutorial_2017_files/figure-markdown_github/unnamed-chunk-6-1.png)

Overlaps between two `GRanges` objects
--------------------------------------

``` r
> # Combine overlapping ranges in the peak object
> gr_S1_reduced <- GenomicRanges::reduce(gr_S1)
> gr_S1_reduced
```

    ## GRanges object with 4 ranges and 0 metadata columns:
    ##       seqnames    ranges strand
    ##          <Rle> <IRanges>  <Rle>
    ##   [1]     chr1  [ 5, 15]      *
    ##   [2]     chr1  [20, 26]      *
    ##   [3]     chr2  [ 8, 16]      *
    ##   [4]     chr2  [18, 21]      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

``` r
> length(gr_S1)
```

    ## [1] 5

``` r
> length(gr_S1_reduced)
```

    ## [1] 4

One can also decide to keep every range distinct and evaluate the overlap for each one of them. For this analysis we will keep the overlapping ranges as disticnt regions.

### `findOverlaps()`

By default looks for overlaps of a minimum of 1bp between a `query` and a `subject`.

``` r
> `?`(GenomicRanges::findOverlaps)
> overlaps <- findOverlaps(query = gr_S1, subject = genes)
> overlaps
```

    ## Hits object with 4 hits and 0 metadata columns:
    ##       queryHits subjectHits
    ##       <integer>   <integer>
    ##   [1]         1           1
    ##   [2]         2           1
    ##   [3]         3           2
    ##   [4]         4           3
    ##   -------
    ##   queryLength: 5 / subjectLength: 4

``` r
> # Query the output
> queryHits(overlaps)
```

    ## [1] 1 2 3 4

``` r
> subjectHits(overlaps)
```

    ## [1] 1 1 2 3

``` r
> subjectLength(overlaps)
```

    ## [1] 4

``` r
> # You can allow for a gap to be ignored
> overlaps <- findOverlaps(gr_S1, genes, maxgap = 5)
> overlaps
```

    ## Hits object with 8 hits and 0 metadata columns:
    ##       queryHits subjectHits
    ##       <integer>   <integer>
    ##   [1]         1           1
    ##   [2]         2           1
    ##   [3]         2           2
    ##   [4]         3           1
    ##   [5]         3           2
    ##   [6]         4           3
    ##   [7]         5           3
    ##   [8]         5           4
    ##   -------
    ##   queryLength: 5 / subjectLength: 4

``` r
> # You can specify a minimum overlap
> overlaps <- findOverlaps(gr_S1, genes, minoverlap = 5)
> overlaps
```

    ## Hits object with 3 hits and 0 metadata columns:
    ##       queryHits subjectHits
    ##       <integer>   <integer>
    ##   [1]         1           1
    ##   [2]         2           1
    ##   [3]         4           3
    ##   -------
    ##   queryLength: 5 / subjectLength: 4

``` r
> # You can specify a minimum overlap
> overlaps <- findOverlaps(gr_S1, genes, type = "start")
> overlaps
```

    ## Hits object with 0 hits and 0 metadata columns:
    ##    queryHits subjectHits
    ##    <integer>   <integer>
    ##   -------
    ##   queryLength: 5 / subjectLength: 4

``` r
> overlaps <- findOverlaps(gr_S1, genes, type = "end")
> overlaps
```

    ## Hits object with 1 hit and 0 metadata columns:
    ##       queryHits subjectHits
    ##       <integer>   <integer>
    ##   [1]         2           1
    ##   -------
    ##   queryLength: 5 / subjectLength: 4

``` r
> overlaps <- findOverlaps(gr_S1, genes, type = "within")
> overlaps
```

    ## Hits object with 1 hit and 0 metadata columns:
    ##       queryHits subjectHits
    ##       <integer>   <integer>
    ##   [1]         2           1
    ##   -------
    ##   queryLength: 5 / subjectLength: 4

### `countOverlaps()`

``` r
> `?`(GenomicRanges::countOverlaps)
```

``` r
> N_overlaps <- GenomicRanges::countOverlaps(gr_S1, genes)
> N_overlaps
```

    ## Peak_1 Peak_2 Peak_3 Peak_4 Peak_5 
    ##      1      1      1      1      0

``` r
> # You can play around with the same options as
> # findOverlaps()
```

Nearest-methods in `GenomicRanges`
----------------------------------

### `nearest()`

``` r
> `?`(GenomicRanges::nearest)
```

It returns a vector of indeces referring to the nearest neighbour range in *subject* for every range in *x*. By default if one range overlaps with multiple genes then one overlap will be chosen at random:

``` r
> GenomicRanges::nearest(x = gr_S1, subject = genes)
```

    ## [1] 1 1 2 3 4

``` r
> GenomicRanges::nearest(x = gr_S1, subject = genes, 
+     select = "all")
```

    ## Hits object with 5 hits and 0 metadata columns:
    ##       queryHits subjectHits
    ##       <integer>   <integer>
    ##   [1]         1           1
    ##   [2]         2           1
    ##   [3]         3           2
    ##   [4]         4           3
    ##   [5]         5           4
    ##   -------
    ##   queryLength: 5 / subjectLength: 4

``` r
> GenomicRanges::nearest(gr_S1)
```

    ## [1] 2 1 2 5 4

``` r
> GenomicRanges::nearest(gr_S1, gr_S1)
```

    ## [1] 1 1 3 4 5

### `distance()`

``` r
> `?`(GenomicRanges::distance)
```

``` r
> GenomicRanges::distance(x = gr_S1[1], y = genes)
```

    ## [1]  0  5 NA NA

``` r
> GenomicRanges::distance(x = gr_S1[1:4], y = genes)
```

    ## [1]  0  1 NA  6

``` r
> GenomicRanges::distance(x = gr_S1, y = genes)
```

    ## [1]  0  1 NA  6 NA

`distance()` is a symmetric function which means that it requires x and y to have to have the same lenght and if one is shorter than the other one it will be recycled to match the length of the longest. Also, the distance between two consecutive blocks is 0 not 1 which affects the notion of overlaps. If `distance(x, y) == 0` then x and y can be either adjacent or overlapping ranges.

### `distanceToNearest()`

``` r
> `?`(GenomicRanges::distanceToNearest)
```

For every range in *x* it will return the index and the distance to its nearest neighbour in *subject*.

``` r
> GenomicRanges::distanceToNearest(x = gr_S1, subject = genes)
```

    ## Hits object with 5 hits and 1 metadata column:
    ##       queryHits subjectHits |  distance
    ##       <integer>   <integer> | <integer>
    ##   [1]         1           1 |         0
    ##   [2]         2           1 |         0
    ##   [3]         3           2 |         0
    ##   [4]         4           3 |         0
    ##   [5]         5           4 |         1
    ##   -------
    ##   queryLength: 5 / subjectLength: 4

`GRangesList`
-------------

`GRangesList` are lists of `GRanges` objects.

``` r
> # Sample 1
> gr_S1 <- GRanges(seqnames = Rle(c("chr1", "chr2"), 
+     c(3, 2)), ranges = IRanges(start = c(5, 8, 20, 
+     8, 18), end = c(11, 15, 26, 16, 21), names = c(paste("Peak_", 
+     1:5, sep = ""))), strand = Rle(strand(c("*")), 
+     c(5)))
> 
> gr_S1
```

    ## GRanges object with 5 ranges and 0 metadata columns:
    ##          seqnames    ranges strand
    ##             <Rle> <IRanges>  <Rle>
    ##   Peak_1     chr1  [ 5, 11]      *
    ##   Peak_2     chr1  [ 8, 15]      *
    ##   Peak_3     chr1  [20, 26]      *
    ##   Peak_4     chr2  [ 8, 16]      *
    ##   Peak_5     chr2  [18, 21]      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

``` r
> # Sample 2
> gr_S2 <- GRanges(seqnames = Rle(c("chr2", "chr3"), 
+     c(3, 5)), ranges = IRanges(start = c(1:8), width = 10, 
+     names = c(paste("Peak_", 1:8, sep = ""))), strand = Rle(strand(c("*")), 
+     c(8)))
> 
> gr_S2
```

    ## GRanges object with 8 ranges and 0 metadata columns:
    ##          seqnames    ranges strand
    ##             <Rle> <IRanges>  <Rle>
    ##   Peak_1     chr2   [1, 10]      *
    ##   Peak_2     chr2   [2, 11]      *
    ##   Peak_3     chr2   [3, 12]      *
    ##   Peak_4     chr3   [4, 13]      *
    ##   Peak_5     chr3   [5, 14]      *
    ##   Peak_6     chr3   [6, 15]      *
    ##   Peak_7     chr3   [7, 16]      *
    ##   Peak_8     chr3   [8, 17]      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

``` r
> # GRanges List
> 
> list_ranges <- GRangesList(Sample1 = gr_S1, Sample2 = gr_S2)
```

Many of the functions learnt for *GRanges* can also be applied to *GRangesList* objects even though the output will have to be interepreted accordingly:

``` r
> names(list_ranges)
```

    ## [1] "Sample1" "Sample2"

``` r
> length(list_ranges)
```

    ## [1] 2

``` r
> seqnames(list_ranges)
```

    ## RleList of length 2
    ## $Sample1
    ## factor-Rle of length 5 with 2 runs
    ##   Lengths:    3    2
    ##   Values : chr1 chr2
    ## Levels(3): chr1 chr2 chr3
    ## 
    ## $Sample2
    ## factor-Rle of length 8 with 2 runs
    ##   Lengths:    3    5
    ##   Values : chr2 chr3
    ## Levels(3): chr1 chr2 chr3

``` r
> strand(list_ranges)
```

    ## RleList of length 2
    ## $Sample1
    ## factor-Rle of length 5 with 1 run
    ##   Lengths: 5
    ##   Values : *
    ## Levels(3): + - *
    ## 
    ## $Sample2
    ## factor-Rle of length 8 with 1 run
    ##   Lengths: 8
    ##   Values : *
    ## Levels(3): + - *

``` r
> ranges(list_ranges)
```

    ## IRangesList of length 2
    ## $Sample1
    ## IRanges object with 5 ranges and 0 metadata columns:
    ##              start       end     width
    ##          <integer> <integer> <integer>
    ##   Peak_1         5        11         7
    ##   Peak_2         8        15         8
    ##   Peak_3        20        26         7
    ##   Peak_4         8        16         9
    ##   Peak_5        18        21         4
    ## 
    ## $Sample2
    ## IRanges object with 8 ranges and 0 metadata columns:
    ##              start       end     width
    ##          <integer> <integer> <integer>
    ##   Peak_1         1        10        10
    ##   Peak_2         2        11        10
    ##   Peak_3         3        12        10
    ##   Peak_4         4        13        10
    ##   Peak_5         5        14        10
    ##   Peak_6         6        15        10
    ##   Peak_7         7        16        10
    ##   Peak_8         8        17        10

``` r
> start(list_ranges)
```

    ## IntegerList of length 2
    ## [["Sample1"]] 5 8 20 8 18
    ## [["Sample2"]] 1 2 3 4 5 6 7 8

``` r
> end(list_ranges)
```

    ## IntegerList of length 2
    ## [["Sample1"]] 11 15 26 16 21
    ## [["Sample2"]] 10 11 12 13 14 15 16 17

``` r
> width(list_ranges)
```

    ## IntegerList of length 2
    ## [["Sample1"]] 7 8 7 9 4
    ## [["Sample2"]] 10 10 10 10 10 10 10 10

``` r
> unlist(list_ranges)
```

    ## GRanges object with 13 ranges and 0 metadata columns:
    ##                  seqnames    ranges strand
    ##                     <Rle> <IRanges>  <Rle>
    ##   Sample1.Peak_1     chr1  [ 5, 11]      *
    ##   Sample1.Peak_2     chr1  [ 8, 15]      *
    ##   Sample1.Peak_3     chr1  [20, 26]      *
    ##   Sample1.Peak_4     chr2  [ 8, 16]      *
    ##   Sample1.Peak_5     chr2  [18, 21]      *
    ##              ...      ...       ...    ...
    ##   Sample2.Peak_4     chr3   [4, 13]      *
    ##   Sample2.Peak_5     chr3   [5, 14]      *
    ##   Sample2.Peak_6     chr3   [6, 15]      *
    ##   Sample2.Peak_7     chr3   [7, 16]      *
    ##   Sample2.Peak_8     chr3   [8, 17]      *
    ##   -------
    ##   seqinfo: 3 sequences from an unspecified genome; no seqlengths

To get the number of ranges in every object of the list use *elementNROWS*:

``` r
> elementNROWS(list_ranges)
```

    ## Sample1 Sample2 
    ##       5       8

To append two *GRanges* lists simply use the R concatenate command *c*:

``` r
> # Sample 3
> gr_S3 <- GRanges(seqnames = Rle(c("chr1", "chr2"), 
+     c(3, 2)), ranges = IRanges(start = 20:24, width = 8, 
+     names = c(paste("Peak_", 1:5, sep = ""))), strand = Rle(strand(c("*")), 
+     c(5)))
> 
> gr_S1
```

    ## GRanges object with 5 ranges and 0 metadata columns:
    ##          seqnames    ranges strand
    ##             <Rle> <IRanges>  <Rle>
    ##   Peak_1     chr1  [ 5, 11]      *
    ##   Peak_2     chr1  [ 8, 15]      *
    ##   Peak_3     chr1  [20, 26]      *
    ##   Peak_4     chr2  [ 8, 16]      *
    ##   Peak_5     chr2  [18, 21]      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

``` r
> # Sample 4
> gr_S4 <- GRanges(seqnames = Rle(c("chr2", "chr3"), 
+     c(3, 5)), ranges = IRanges(start = 20:27, width = 10, 
+     names = c(paste("Peak_", 1:8, sep = ""))), strand = Rle(strand(c("*")), 
+     c(8)))
> 
> gr_S4
```

    ## GRanges object with 8 ranges and 0 metadata columns:
    ##          seqnames    ranges strand
    ##             <Rle> <IRanges>  <Rle>
    ##   Peak_1     chr2  [20, 29]      *
    ##   Peak_2     chr2  [21, 30]      *
    ##   Peak_3     chr2  [22, 31]      *
    ##   Peak_4     chr3  [23, 32]      *
    ##   Peak_5     chr3  [24, 33]      *
    ##   Peak_6     chr3  [25, 34]      *
    ##   Peak_7     chr3  [26, 35]      *
    ##   Peak_8     chr3  [27, 36]      *
    ##   -------
    ##   seqinfo: 2 sequences from an unspecified genome; no seqlengths

``` r
> # Second GRanges List
> list_ranges2 <- GRangesList(Sample3 = gr_S3, Sample4 = gr_S4)
```

``` r
> append_lists <- c(list_ranges, list_ranges2)
```

Subsetting and looping over *GRanges* list
------------------------------------------

Add an extra column to every GRanges in the GRangesList.

``` r
> addCols <- lapply(append_lists, function(x) {
+     elementMetadata(x) <- data.frame(NumberReads = rbinom(length(x), 
+         size = 100, prob = 0.5))
+     return(x)
+ })
> class(addCols)
```

    ## [1] "list"

``` r
> addCols <- GRangesList(addCols)
```

In many cases *GRanges* objects can be subsetted using the same rules that apply to normal lists or data frames in R:

``` r
> # As a list
> addCols[[1]]
```

    ## GRanges object with 5 ranges and 1 metadata column:
    ##          seqnames    ranges strand | NumberReads
    ##             <Rle> <IRanges>  <Rle> |   <integer>
    ##   Peak_1     chr1  [ 5, 11]      * |          57
    ##   Peak_2     chr1  [ 8, 15]      * |          50
    ##   Peak_3     chr1  [20, 26]      * |          56
    ##   Peak_4     chr2  [ 8, 16]      * |          47
    ##   Peak_5     chr2  [18, 21]      * |          50
    ##   -------
    ##   seqinfo: 3 sequences from an unspecified genome; no seqlengths

``` r
> addCols["Sample1"]
```

    ## GRangesList object of length 1:
    ## $Sample1 
    ## GRanges object with 5 ranges and 1 metadata column:
    ##          seqnames    ranges strand | NumberReads
    ##             <Rle> <IRanges>  <Rle> |   <integer>
    ##   Peak_1     chr1  [ 5, 11]      * |          57
    ##   Peak_2     chr1  [ 8, 15]      * |          50
    ##   Peak_3     chr1  [20, 26]      * |          56
    ##   Peak_4     chr2  [ 8, 16]      * |          47
    ##   Peak_5     chr2  [18, 21]      * |          50
    ## 
    ## -------
    ## seqinfo: 3 sequences from an unspecified genome; no seqlengths

``` r
> # As a data.frame object
> addCols[1, "NumberReads"]
```

    ## GRangesList object of length 1:
    ## $Sample1 
    ## GRanges object with 5 ranges and 1 metadata column:
    ##          seqnames    ranges strand | NumberReads
    ##             <Rle> <IRanges>  <Rle> |   <integer>
    ##   Peak_1     chr1  [ 5, 11]      * |          57
    ##   Peak_2     chr1  [ 8, 15]      * |          50
    ##   Peak_3     chr1  [20, 26]      * |          56
    ##   Peak_4     chr2  [ 8, 16]      * |          47
    ##   Peak_5     chr2  [18, 21]      * |          50
    ## 
    ## -------
    ## seqinfo: 3 sequences from an unspecified genome; no seqlengths

``` r
> addCols["Sample1", "NumberReads"]
```

    ## GRangesList object of length 1:
    ## $Sample1 
    ## GRanges object with 5 ranges and 1 metadata column:
    ##          seqnames    ranges strand | NumberReads
    ##             <Rle> <IRanges>  <Rle> |   <integer>
    ##   Peak_1     chr1  [ 5, 11]      * |          57
    ##   Peak_2     chr1  [ 8, 15]      * |          50
    ##   Peak_3     chr1  [20, 26]      * |          56
    ##   Peak_4     chr2  [ 8, 16]      * |          47
    ##   Peak_5     chr2  [18, 21]      * |          50
    ## 
    ## -------
    ## seqinfo: 3 sequences from an unspecified genome; no seqlengths

``` r
> lapply(addCols, length)
```

    ## $Sample1
    ## [1] 5
    ## 
    ## $Sample2
    ## [1] 8
    ## 
    ## $Sample3
    ## [1] 5
    ## 
    ## $Sample4
    ## [1] 8

``` r
> sapply(addCols, length)
```

    ## Sample1 Sample2 Sample3 Sample4 
    ##       5       8       5       8

``` r
> lapply(addCols, start)
```

    ## $Sample1
    ## [1]  5  8 20  8 18
    ## 
    ## $Sample2
    ## [1] 1 2 3 4 5 6 7 8
    ## 
    ## $Sample3
    ## [1] 20 21 22 23 24
    ## 
    ## $Sample4
    ## [1] 20 21 22 23 24 25 26 27

``` r
> # Concatenate and reduce the GRangesList
> Reduce(c, addCols)
```

    ## GRanges object with 26 ranges and 1 metadata column:
    ##          seqnames    ranges strand | NumberReads
    ##             <Rle> <IRanges>  <Rle> |   <integer>
    ##   Peak_1     chr1  [ 5, 11]      * |          57
    ##   Peak_2     chr1  [ 8, 15]      * |          50
    ##   Peak_3     chr1  [20, 26]      * |          56
    ##   Peak_4     chr2  [ 8, 16]      * |          47
    ##   Peak_5     chr2  [18, 21]      * |          50
    ##      ...      ...       ...    ... .         ...
    ##   Peak_4     chr3  [23, 32]      * |          56
    ##   Peak_5     chr3  [24, 33]      * |          61
    ##   Peak_6     chr3  [25, 34]      * |          49
    ##   Peak_7     chr3  [26, 35]      * |          56
    ##   Peak_8     chr3  [27, 36]      * |          46
    ##   -------
    ##   seqinfo: 3 sequences from an unspecified genome; no seqlengths

Similar considerations apply for the other functions explored above like *findOverlaps()*, *countOverlaps()*, etc...

#### Challenge

> Using the code below to download the `cpg_islands` object containing the CpG islands for the Human chr21 and the `txdb` gene annotation from `TxDb.Hsapiens.UCSC.hg19.knownGene`.

``` r
> # Download CpG islands from UCSC
> session <- browserSession("UCSC")
> query <- ucscTableQuery(session, "CpG Islands", GRangesForUCSCGenome("hg19", 
+     "chr21"))
> cpg_islands <- getTable(query)
> # Download genes
> txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
```

> 1.  What is the class of `cpg_islandsGR`? Convert `cpg_islandsGR` to a `GRanges` object.
> 2.  Using what you have previosly learnt, get the transcript for every gene from the `txdb` gene annotation and store them into the `transcriptHg19` object (suggestion: see `transcriptsBy`).
> 3.  Subset `transcriptHg19` to extract only the transcript on chr21.
> 4.  Find the nearest CpG island in `cpg_islandsGR` to every transcript in `transcriptHg19` as well as their distance.
> 5.  Optional: On average, every CpG island is close to how many genes?

#### Solutions

``` r
> # 1. Convert to GRanges
> cpg_islandsGR <- as(cpg_islands, "GRanges")
> 
> # 2. Get the transcript for every gene from the
> # `txdb`
> transcriptHg19 <- transcriptsBy(txdb, by = "gene")
> 
> # 3. Only select chr21
> transcriptHg19 <- unlist(transcriptHg19)
> chr21 <- transcriptHg19[seqnames(transcriptHg19) == 
+     "chr21", ]
> 
> # 4. Find the nearest CpG island in `cpg_islandsGR`
> # to every transcript in `transcriptHg19` as well
> # as their distance.
> chr21_cpg <- distanceToNearest(chr21, cpg_islandsGR)
> 
> # 5. On average, every CpG island is close to how
> # many genes?
> mean(table(subjectHits(chr21_cpg)))
```

    ## [1] 4.520833

`Rtracklayer`
=============

Now that the structure and manipulation of `GRanges` objects is clear, the `rtracklayer` package comes as a useful tool to import/export `GRanges` objects from and to different data formats commonly used in genomic analyses. `Rtracklayer` stands as an interface to mediate the crosstalk between R and genome browsers (UCSC built-in). Example of data files supported by this package are BED, bigWIG, wig, gff and gtf formats.

Another very important function of this package is to allow the visualisations of genomic annotation *tracks*. However, this part will not covered in the present tutorial. For more references see the [Rtracklayer Bioconductor page](http://bioconductor.org/packages/release/bioc/html/rtracklayer.html).

There is a very large number of different data formats used to store different sort of genomic data. In this tutororial only few of them will be covered as example. However, the UCSC genome growser offers a useful [FAQ page](https://genome.ucsc.edu/FAQ/FAQformat.html) where one can find the specification of all mandatory and optional fields for every data format.

GTF and GFF
-----------

The **GFF** and **GTF** are the preferred format used to store annotations and their exact specification is summarised in [GFF/GTF File Format](http://www.ensembl.org/info/website/upload/gff.html).

The function `rtracklayer::import()` is used to import all the supported data types into a `GRanges` object in R. The function recognises the format from the extension of the file but the argument `format` can be used to expliciclty define it. The function `rtracklayer::export()` works in the same way and it is used to export `GRanges` objects to files.

``` r
> # Import/Export with Rtrcaklayer
> rtracklayer::export(chr21, file.path(dir, "transcriptHg19.gff"))
> rtracklayer::export(chr21, file.path(dir, "transcriptHg19.gtf"))
> 
> importHg19 <- rtracklayer::import(file.path("transcriptHg19.gff"))
> class(importHg19)
```

    ## [1] "GRanges"
    ## attr(,"package")
    ## [1] "GenomicRanges"

``` r
> importHg19
```

    ## GRanges object with 868 ranges and 5 metadata columns:
    ##         seqnames               ranges strand |      source
    ##            <Rle>            <IRanges>  <Rle> |    <factor>
    ##     [1]    chr21 [47247755, 47256333]      - | rtracklayer
    ##     [2]    chr21 [31661463, 31661832]      - | rtracklayer
    ##     [3]    chr21 [ 9907189,  9968593]      - | rtracklayer
    ##     [4]    chr21 [ 9907189,  9968593]      - | rtracklayer
    ##     [5]    chr21 [ 9915250,  9968593]      - | rtracklayer
    ##     ...      ...                  ...    ... .         ...
    ##   [864]    chr21 [35012051, 35014160]      - | rtracklayer
    ##   [865]    chr21 [37536839, 37666572]      + | rtracklayer
    ##   [866]    chr21 [37537006, 37666572]      + | rtracklayer
    ##   [867]    chr21 [37617346, 37634210]      + | rtracklayer
    ##   [868]    chr21 [35736323, 35743440]      + | rtracklayer
    ##                     type     score     phase    group
    ##                 <factor> <numeric> <integer> <factor>
    ##     [1] sequence_feature      <NA>      <NA>    chr21
    ##     [2] sequence_feature      <NA>      <NA>    chr21
    ##     [3] sequence_feature      <NA>      <NA>    chr21
    ##     [4] sequence_feature      <NA>      <NA>    chr21
    ##     [5] sequence_feature      <NA>      <NA>    chr21
    ##     ...              ...       ...       ...      ...
    ##   [864] sequence_feature      <NA>      <NA>    chr21
    ##   [865] sequence_feature      <NA>      <NA>    chr21
    ##   [866] sequence_feature      <NA>      <NA>    chr21
    ##   [867] sequence_feature      <NA>      <NA>    chr21
    ##   [868] sequence_feature      <NA>      <NA>    chr21
    ##   -------
    ##   seqinfo: 93 sequences from hg19 genome

#### Challenge

> 1.  Create a simple `GRanges` object made out of 20 ranges and export it as a `gff` file.
> 2.  Read the file back into R.
> 3.  Print the number of ranges of the `GRanges` object just imported and create a histogram of the widths of the ranges.

#### Solution

``` r
> # 1. Create GRanges object made of 20 ranges
> exampleGR <- GRanges(seqnames = Rle("chr1", 20), IRanges(start = 1:20, 
+     width = rbinom(20, size = 100, prob = 0.6)))
> export(exampleGR, file.path(dir, "exampleGR.gff"))
> # 2. Read the file back into R.
> exampleGR_import <- import(file.path(dir, "exampleGR.gff"))
> # 3. Print the number of rows of the `GRanges`
> # object just imported and create a histogram of
> # the widths of the ranges.
> length(exampleGR_import)
```

    ## [1] 20

``` r
> hist(width(exampleGR_import))
```

![](combine_tutorial_2017_files/figure-markdown_github/unnamed-chunk-30-1.png)

Wiggle (WIG) and bigWIG file formats for graphing tracks
--------------------------------------------------------

The [WIG](http://www.ensembl.org/info/website/upload/wig.html#tracklines) data format is used for the disply of dense continuous data, like scores. Every line in the file is an intervals and *WIG* files only allow equally spaced intervals. For intervals of variable width one should use the [bedGraph](http://www.ensembl.org/info/website/upload/bed.html#bedGraph) format. There are two main data formats: the *variableStep* and the *fixedStep*. The *WIG* file is constituted by one or more blocks separated by declaration lines. Let's go through them with two examples.

-   **variableStep**

The following example is a *WIG* file with only one block of coordinates. The field *chrom* is required.

``` bash
variableStep chrom=chr2
300701  12.5
300702  12.5
300703  12.5
300704  12.5
300705  12.5
```

The same file can be defined as follows using the *span* argument.

``` bash
variableStep chrom=chr2 span=5
300701  12.5
```

-   **fixedStep**

This format allows a more compact way of storing intervals. In this situation the *span* would not change the dimension of the file.

``` bash
fixedStep chrom=chr3 start=400601 step=100
11
22
33
```

The [bigWig](https://genome.ucsc.edu/goldenpath/help/wiggle.html) file is created from a *WIG* file and in `R` this can be achieved through `rtracklayer::wigToBigWig()`. *bigWig* is the recommended format for very large data tracks due to the way in which data are compressed and stored in an indexed binary format. Usually whole-genome coverage vectors are stored as *bigWig*. Howver, the loss is negligible when dealing with very large amount of data.

``` r
> wig_path <- file.path(dir, "exampleWIG.wig")
> import_wig <- rtracklayer::import(con = wig_path, seqinfo = Seqinfo(genome = "mm10"))
> import_wig
```

    ## GRanges object with 16 ranges and 1 metadata column:
    ##        seqnames           ranges strand |     score
    ##           <Rle>        <IRanges>  <Rle> | <numeric>
    ##    [1]     chr2 [300701, 300701]      * |      12.5
    ##    [2]     chr2 [300702, 300702]      * |      12.5
    ##    [3]     chr2 [300703, 300703]      * |      12.5
    ##    [4]     chr2 [300704, 300704]      * |      12.5
    ##    [5]     chr2 [300705, 300705]      * |      12.5
    ##    ...      ...              ...    ... .       ...
    ##   [12]     chr3 [401501, 401501]      * |      10.5
    ##   [13]     chr4 [500001, 500001]      * |        10
    ##   [14]     chr4 [500201, 500201]      * |        20
    ##   [15]     chr4 [500401, 500401]      * |        20
    ##   [16]     chr4 [500601, 500601]      * |        20
    ##   -------
    ##   seqinfo: 66 sequences (1 circular) from mm10 genome

``` r
> `?`(GenomeInfoDb::Seqinfo)
> seqinfo(import_wig)
```

    ## Seqinfo object with 66 sequences (1 circular) from mm10 genome:
    ##   seqnames       seqlengths isCircular genome
    ##   chr1            195471971      FALSE   mm10
    ##   chr2            182113224      FALSE   mm10
    ##   chr3            160039680      FALSE   mm10
    ##   chr4            156508116      FALSE   mm10
    ##   chr5            151834684      FALSE   mm10
    ##   ...                   ...        ...    ...
    ##   chrUn_GL456392      23629      FALSE   mm10
    ##   chrUn_GL456393      55711      FALSE   mm10
    ##   chrUn_GL456394      24323      FALSE   mm10
    ##   chrUn_GL456396      21240      FALSE   mm10
    ##   chrUn_JH584304     114452      FALSE   mm10

Convert the *WIG* file into a *bigWig* file and import the *bigWig* file into R.

``` r
> # Output file
> bigwig_path <- file.path(dir, "example_wigToBigWig.bw")
> rtracklayer::wigToBigWig(wig_path, dest = bigwig_path, 
+     seqinfo = seqinfo(import_wig))
> import_bigwig <- rtracklayer::import(con = bigwig_path)
```

`import` only a region of the bigWig file
-----------------------------------------

The `import` functions allows to load into a R only a specific regions defined using a `GRanges` object.

``` r
> which_region <- GRanges(seqnames = "chr3", IRanges(start = 1, 
+     end = 5e+05))
> import_bigwig_region <- rtracklayer::import(con = bigwig_path, 
+     which = which_region)
```

`liftOver` from different genome releases
=========================================

``` r
> library(AnnotationHub)
> `?`(AnnotationHub)
> ahub <- AnnotationHub()
> 
> ahub.chain <- subset(ahub, rdataclass == "ChainFile" & 
+     species == "Homo sapiens")
> query(ahub.chain, c("hg18", "hg19"))
```

    ## AnnotationHub with 2 records
    ## # snapshotDate(): 2016-10-11 
    ## # $dataprovider: UCSC
    ## # $species: Homo sapiens
    ## # $rdataclass: ChainFile
    ## # additional mcols(): taxonomyid, genome, description,
    ## #   coordinate_1_based, maintainer, rdatadateadded, preparerclass,
    ## #   tags, sourceurl, sourcetype 
    ## # retrieve records with, e.g., 'object[["AH14149"]]' 
    ## 
    ##             title                   
    ##   AH14149 | hg19ToHg18.over.chain.gz
    ##   AH14220 | hg18ToHg19.over.chain.gz

``` r
> chain <- ahub.chain[ahub.chain$title == "hg19ToHg18.over.chain.gz"]
> chain <- chain[[1]]
> gr.hg18 <- rtracklayer::liftOver(import_bigwig, chain)
> gr.hg18
```

    ## GRangesList object of length 16:
    ## [[1]] 
    ## GRanges object with 1 range and 1 metadata column:
    ##       seqnames           ranges strand |     score
    ##          <Rle>        <IRanges>  <Rle> | <numeric>
    ##   [1]     chr2 [290701, 290701]      * |      12.5
    ## 
    ## [[2]] 
    ## GRanges object with 1 range and 1 metadata column:
    ##       seqnames           ranges strand | score
    ##   [1]     chr2 [290702, 290702]      * |  12.5
    ## 
    ## [[3]] 
    ## GRanges object with 1 range and 1 metadata column:
    ##       seqnames           ranges strand | score
    ##   [1]     chr2 [290703, 290703]      * |  12.5
    ## 
    ## ...
    ## <13 more elements>
    ## -------
    ## seqinfo: 3 sequences from an unspecified genome; no seqlengths

`Rsamtools`
===========

[Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html) is an R interface for [Samtools](http://samtools.sourceforge.net/) which offers a large variety of utilities to manipulate *SAM* and *BAM* files. The main purose of [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html) is to import *BAM* files into R. This will allow the user to then create objects which can be used for several types of downstream analyses. *BAM* files are the standard way of storing 'short' aligned (and unaligned) reads produced by any type of experiment (RNA-Seq, ChIP-Seq, Methylation data, etc). *BAM* files are binary versions of *SAM* files and can be indexed allowing to access only localised regions of the genome, similar to what we have seen for the *bigWig* file. Each read is stored in a *BAM* file with several other measures like base quality, position of the 5' end of the read, read name, etc .... [Rsamtools](http://bioconductor.org/packages/release/bioc/html/Rsamtools.html) allows the user to specify which parameters of the reads to load into R. A detailed description of the fields in the *BAM* and their names is available at <http://samtools.github.io/hts-specs/SAMv1.pdf>. It has to be remembered that *BAM* files can be very large and there are limits in how much can be done into R and for more complex analyses is preferred to work with [Samtools](http://samtools.sourceforge.net/) is prefered.

``` bash
> 
+ D00626:239:CAFMCANXX:4:1201:16222:93271 403  chr1    17018   3  38M177N62M      =       16965   -330 TGGCCCAGGTCTGGCACATAGAAGTAGTTCTCTGGGACCTGCTGTTCCAGCTGCTCT  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
+ RG:Z:9_CAFMCANXX_L004   NH:i:2  HI:i:2  NM:i:0  nM:i:1  AS:i:200
```

A very important field which is worth understanding is the *FLAG* of the read. This is a binary number which defines the type of read. In the example above it is the number *403*. For example, the *FLAG* defines whether the read is unmapped or mapped or with paired end (PE) reads if the first mate is mapped and the second one is not, and so on. It is not always trivial to decode a *FLAG* number and the [Broad Institute](https://www.broadinstitute.org/) offers a useful online tool to do that for you <https://broadinstitute.github.io/picard/explain-flags.html>.

Input BAM into R
----------------

The function used to load BAM files into R is `scanBam()`. The function `ScanBamParam()` is used to specify the fields and the region/s of the *BAM* file to load into R. The genomic region to load is defined through `GRanges`. Below is an example on how to load a few regions of 1000 bp width from a *BAM* file containing reads on chromosome 21 from an RNA-Seq experiment. We also require only *uniquely mapped* and *properly paired* reads to be loaded as well as we do not load PCR duplicates. The help page of `?scanBam` lists the name of the fields that can be accessed through the `what` argument in the `ScanBamParam()` function.

``` r
> # Path to the bamfile. In the same folder there has
> # to be the chr1.bam.bai file
> bamfile <- file.path(dir, "chr21.bam")
> # Define parameters
> which <- GRanges(seqnames = rep("chr21", 3), IRanges(start = c(9827000, 
+     16267000, 15890000), width = 10000))
> what <- c("rname", "strand", "pos", "qwidth", "seq", 
+     "isize")
> # ?scanBamFlag
> flag = scanBamFlag(isDuplicate = FALSE, isProperPair = TRUE, 
+     isPaired = TRUE)
> 
> param <- ScanBamParam(which = which, what = what, flag = flag)
> # load reads into R
> reads = scanBam(bamfile, param = param)
> 
> # Explore the reads
> names(reads)
```

    ## [1] "chr21:9827000-9836999"   "chr21:16267000-16276999"
    ## [3] "chr21:15890000-15899999"

``` r
> region1 <- reads$`chr21:9827000-9836999`
> names(region1)
```

    ## [1] "rname"  "strand" "pos"    "qwidth" "isize"  "seq"

``` r
> head(region1$rname)
```

    ## [1] chr21 chr21 chr21 chr21 chr21 chr21
    ## 25 Levels: chrM chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 ... chrY

``` r
> head(region1$pos)
```

    ## [1] 9826933 9826933 9826934 9826934 9826934 9826940

``` r
> head(region1$qwidth)
```

    ## [1] 100 100 100 100 100 100

``` r
> head(region1$seq)
```

    ##   A DNAStringSet instance of length 6
    ##     width seq
    ## [1]   100 CCGGGCCCGTCCTCGCGAGGCCCCCCGGCCG...TACCTACCTGGTTGATCCTGCCACTAGCATA
    ## [2]   100 CCGGGCCCGTCCTCGCGAGGCCCCCCGGCCG...TAACTACCTGGTTGATTCTGCCAGTAGCATA
    ## [3]   100 CGGGCCCGTCCTCGCGAGGCCCCCCGGCCGG...ACCTACCTGGTTGATCCTGCCAGTAGCATAT
    ## [4]   100 CGGGCCCGTCCTCGCGAGGCCCCCCGGCCGG...ACCTACCTGGTTGATCCTGCCAGTAGCATAT
    ## [5]   100 CGGGCCCGTCCTCGCGAGGCCCCCCGGCCGG...ACCTACCTGGTTGATCCTGCCAGTAGCATAT
    ## [6]   100 CGTCCTCGCGAGGCCCCCCGGCCGGCCGTCC...CTGGTTGATCCTGCCAGTAGCATTTGCTTGT

``` r
> head(region1$isize)
```

    ## [1] 146 144 162 162 145 158

`GenomicAlignments::readGAlignments` also loads *BAM* files into R directly into a `GRanges` object. The parameters are specified as in `RSamtools` with the argument `param`.

``` r
> readsGA <- GenomicAlignments::readGAlignments(bamfile, 
+     param = param)
```

Compute number of reads that falls within 100bp bins
----------------------------------------------------

It is often useful to summarise the number of reads falling within a bin to have an idea about the distribution of the reads across a region. There are several ways in which this can be accomplished and here two will be shown. Often, a read is counted in a bin if its 5' end falls within bin. This is the approach that we will consider.

``` r
> bins <- GRanges(seqnames = rep("chr21", 11), IRanges(start = seq(9827000, 
+     9827000 + 10000, by = 1000), width = 1000))
> # Create vector of 0-1 positions
> positions <- region1$pos
> vector_positions <- rep(0, max(positions))
> vector_positions[positions] <- 1
> # Views Object
> rangeViews <- Views(vector_positions, start = seq(min(positions), 
+     max(positions), length.out = 10), width = 60)
> # viewSums
> bin_sum <- viewSums(rangeViews)
> # Add extra column to the ranges object
> binned_counts <- ranges(rangeViews)
> values(binned_counts) <- data.frame(ReadCounts = bin_sum)
```

#### Challenge

> 1.  Define a new `ScanBamParam()` object that satisfies the following options: use the same ranges as above, load both main alignments and PCR duplicates, properly paired and uniquely mapped reads.
> 2.  Plot an histogram of the fragment sizes of every read (Suggestion: look for the `isize` fields)

#### Solution

``` r
> # 1.  Path to the bamfile. In the same folder there
> # has to be the chr1.bam.bai file
> bamfile <- file.path(dir, "chr21.bam")
> # Define parameters
> which <- GRanges(seqnames = rep("chr21", 3), IRanges(start = c(9827000, 
+     16267000, 15890000), width = 10000))
> what <- c("rname", "strand", "pos", "qwidth", "seq", 
+     "isize")
> # ?scanBamFlag
> flag = scanBamFlag(isDuplicate = NA, isProperPair = TRUE, 
+     isPaired = TRUE)
> 
> param <- ScanBamParam(which = which, what = what, flag = flag)
> # load reads into R
> reads = scanBam(bamfile, param = param)
> 
> # 2.
> hist(reads$`chr21:16267000-16276999`$isize)
```

![](combine_tutorial_2017_files/figure-markdown_github/unnamed-chunk-42-1.png)
