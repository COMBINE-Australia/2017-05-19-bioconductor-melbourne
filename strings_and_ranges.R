#
# Contents
# ========
#
# Installing Bioconductor packages
#
# Finding documentation
#
# IRanges
#
# GRanges
#
# DNAString and DNAStringSet
#
# BSgenome and TxDb packages
#
# Reading and writing files with rtracklayer
#


# ==============================================================================
# Installing Bioconductor packages

# This will install the packages used in this file.
source("http://www.bioconductor.org/biocLite.R")
biocLite(c(
    "GenomicRanges",
    "Biostrings", 
    "BSgenome", 
    "BSgenome.Scerevisiae.UCSC.sacCer3",
    "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"))



# Load packages

library(GenomicRanges)
library(Biostrings)
library(BSgenome)



# ==============================================================================
# Finding documentation

# Bioconductor packages usually have documentation in the form of "vignettes".
# These are also available on the Bioconductor website.
# https://www.bioconductor.org/help/

vignette()
vignette(package="Biostrings")
vignette("BiostringsQuickOverview")

browseVignettes()



# ==============================================================================
# IRanges
#
# Integer ranges, 1-based, from start to end inclusive.
#
# This is R, where everything is a vector, so there is no singular IRange, 
# only plural IRanges.

myiranges <- IRanges(start=c(5,20,25), end=c(10,30,40))

# Like all objects from Bioconductor, myiranges is an "S4" object.
# It has a class, which can be obtained with class(). This determines 
# - The functions it can be used with (methods).
# - What "slots" it has for storing data.
# Data stored in slots is accessible with @ (like $ for lists),
# but code that uses accessor functions will be more generic
# and less liable to break in future.

myiranges
class(myiranges)
methods(class="IRanges")
?"IRanges-class"

# Accessor functions
start(myiranges)
#(or myiranges@start, but start(myiranges) is preferable)
end(myiranges)
width(myiranges)

# IRanges supports many operations
resize(myiranges, 3, fix="start")
resize(myiranges, 3, fix="end")

# ... ok, this gets confusing, let me illustrate...

show_iranges <- function(obj) {
    for(i in seq_along(obj))
        cat(rep(" ", start(obj)[i]),
            rep("=", width(obj)[i]),
            "\n", sep="")
}

show_iranges( myiranges                              )

show_iranges( resize(myiranges, 3, fix="start")      )
show_iranges( resize(myiranges, 3, fix="end")        )

show_iranges( flank(myiranges, 2, start=TRUE)        )
show_iranges( flank(myiranges, 2, start=FALSE)       )

show_iranges( range(myiranges)                       )
show_iranges( reduce(myiranges)                      )
show_iranges( disjoin(myiranges)                     )

show_iranges( setdiff(range(myiranges), myiranges)   )
show_iranges( intersect(myiranges[2], myiranges[3])  )
show_iranges( union(myiranges[2], myiranges[3])      )


?"intra-range-methods"
?"inter-range-methods"
?"setops-methods"
?"findOverlaps-methods"
?"nearest-methods"

# Note:
# Loading tidyverse packages may clobber some of these generic functions.
# If necessary use BiocGenerics::setdiff, etc.


# Challenge
# ---------
#
# Without using the promoters() function,
# how would you get IRanges from -5 to +2 of the starts of myiranges?
#


# ==============================================================================
# GRanges
#
# To refer to a location in a genome we also need
# - sequence name (chromosome)
# - strand, + or -
#

mygranges <- GRanges(
    seqnames = c("chrII", "chrI", "chrI"),
    ranges = IRanges(start=c(5,20,25), end=c(10,30,40)),
    strand = c("+", "-", "+"))

mygranges
class(mygranges)
methods(class="GRanges")
?"GRanges-class"

seqnames(mygranges)
strand(mygranges)
ranges(mygranges)
start(mygranges)
end(mygranges)
as.data.frame(mygranges)


# GRanges are like vectors:
mygranges[2]
c(mygranges, mygranges)

# Like other vectors, elements may be named
names(mygranges) <- c("foo", "bar", "baz")
mygranges

# GRanges can have metadata columns, so they are also a little like data frames:
mygranges$wobble <- c(10, 20, 30)
mygranges
mcols(mygranges)
mygranges$wobble

# A handy way to create a GRanges
as("chrI:3-5:+", "GRanges")


# All the operations we saw for IRanges are available for GRanges.
# However most GRanges operations will take strand into account.

resize(mygranges, 3, fix="start")
resize(mygranges, 3, fix="end")


# ==============================================================================
# DNAString and DNAStringSet

myseq <- DNAString("gcgctgctggatgcgaccgcgcatgcgagcgcgacctatccggaa")

myseq
class(myseq)
methods(class="DNAString")
?"DNAString-class"

reverseComplement(myseq)
translate(myseq)
subseq(myseq, 3,5)

# DNAString objects behave like vectors.
myseq[3:5]
c(myseq, myseq)

# as() converts between types.
as(myseq, "character")
as("ACGT", "DNAString")

# We often want to work with a collection of DNA sequences
myset <- DNAStringSet(list(chrI=myseq, chrII=DNAString("ACGTACGTAA")))

myset
class(myset)
?"DNAStringSet-class"

lengths(myset)
seqinfo(myset)

# Since a DNAString is like a vector, a DNAStringSet is has to be like a list.
myset$chrII
# or myset[["chrII"]]
# or myset[[2]]

# Getting sequences with GRanges
getSeq(myset, mygranges)

getSeq(myset, as("chrI:1-3:+", "GRanges"))
getSeq(myset, as("chrI:1-3:-", "GRanges"))
# Performs reverse complement if strand is "-".

# GRanges can also have seqinfo, allowing various forms of error checking.


# Challenge
# ---------
#
# Reverse complement the following DNA sequence and then translate to an 
# amino acid sequence:
# 
# TTCCATTTCCAT
#



# ==============================================================================
# BSgenome and TxDb packages
#
# Genomes and genes for many model organisms are available as 
# Bioconductor packages
#
# BSgenome. ...  - Genomes as DNAStringSet-like objects
# TxDB. ...      - Genes, transcripts, exons, CDS as GRanges
# org. ...       - Translation of gene names, assignment to GO categories, etc
#
# Further data may be available using AnnotationHub.
#
# We will be using packages for yeast in this part of the workshop.
#

library(BSgenome.Scerevisiae.UCSC.sacCer3)
# A yeast genome object has been loaded.
# Actual sequence data is only loaded from disk as needed.
genome <- Scerevisiae

class(genome)
?"BSgenome-class"

seqinfo(genome)
genome$chrM


library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
# An object referring to a yeast transcriptome database has been loaded.
# Actual data is only loaded from disk as needed.
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

class(txdb)
?"TxDb-class"

# Genes have transcripts. 
# Transcripts have exons, and CDSs if they are protein coding.
#
# The GRange for a transcript spans all of its exons.
# The GRange for a gene spans all the exons of all of its transcripts.
# For example:
# =============================================== gene
#
# ========================  transcript1 of gene
# =======   ====  ========  exons of transcript1
#    ====   ====  ===       CDSs of transcript1
#
#         ======================================= transcript2 of gene
#         ======              ======   ========== exons of transcript2
#           ====              ======   ===        CDSs of transcript2
#
genes(txdb)
transcriptsBy(txdb, "gene")
exonsBy(txdb, "tx", use.names=TRUE)
cdsBy(txdb, "tx", use.names=TRUE)



# To illustrate the use of a TxDb, let us try to extract the start codons of
# yeast. Our first step is to obtain the locations of the coding sequences.

cds_list <- cdsBy(txdb, "tx", use.names=TRUE)

cds_list
class(cds_list)
lengths(cds_list)
table( lengths(cds_list) )
cds_list[ lengths(cds_list) >= 2 ]

seqinfo(cds_list)
genome(cds_list)
# GRanges and GRangesLists extracted from the txdb have associated seqinfo.
# Bioconductor will give an error if you try to use them with the wrong
# genome.

# The ...By functions return GRangesList objects.
# To flatten down to a GRanges, use unlist.
unlist(cds_list)
# Note that names will not be unique unless each element has length 1.


cds_ranges <- unlist( range(cds_list) )
start_codons <- resize(cds_ranges, 3, fix="start")
start_seqs <- getSeq(genome, start_codons)
table(start_seqs)

# Many Bioconductor types have a List version.
# This allows use of the split-apply-combine pattern. 
# List classes support most of the same functions as their base type. 
# If a function isn't supported, use lapply.
#
# Here
# split:   cdsBy  has performed the splitting for us.
# apply:   range  has a GRangesList version that calculates the range for each
#                 element individually.
# combine: unlist goes from GRangesList to GRanges, 
#                 taking names from the list names.
#
# Use split() to perform the splitting step manually.



# Challenge
# ---------
#
# In DNA, a stop codon can be TAG, TAA or TGA.
#
# Are stop codons included in the CDS sequence, or just after it?
#
# Does yeast have a preferred stop codon?
#
# Hint:
# Use resize(..., fix="end") and flank(..., start=FALSE)
#



# The above approach to extracting start codons has a potential flaw:
# An intron may lie within the start codon.

# An alternative way to get start codons.
cds_seqs <- extractTranscriptSeqs(genome, cds_list)
table( subseq(cds_seqs, 1, 3) )

# It really happens!
cds_list[ start_seqs == "AGC" ]
cds_list[ start_seqs == "AGT" ]



# ==============================================================================
# Reading and writing files with rtracklayer

library(rtracklayer)

export(cds_ranges, "cds_ranges.gff3")
export(start_codons, "start_codons.gff3")

import("cds_ranges.gff3")

# Examine these files in IGV
# https://software.broadinstitute.org/software/igv/download
#
# Load "sacCer3" genome (dropdown box in top left of window, may need to choose "More...")
# File/Load from file... to load GFF files.
# YBR111W-A has introns, check we have handled this correctly.


# rtracklayer can also read and write sequences in FASTA format.
# (alternatively use Biostrings::writeXStringSet, Biostrings::readDNAStringSet)
export(start_seqs, "start_codons.fasta")

import("start_codons.fasta", type="DNA")


