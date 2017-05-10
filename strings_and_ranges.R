
source("http://www.bioconductor.org/biocLite.R")
biocLite(c(
    "GenomicRanges",
    "Biostrings", 
    "BSgenome", 
    "BSgenome.Scerevisiae.UCSC.sacCer3",
    "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"))

library(GenomicRanges)
library(Biostrings)
library(BSgenome)

# ==============================================================================
# DNA sequences


# Challenge

# ==============================================================================
# Genomic ranges



# Challenge

# Operations on GRanges

# Challenge

# ==============================================================================
# Rle weirdness
#
# Bioconductor sometimes likes to use Run-Length Encoded (RLE) vectors.
# This saves memory when there are long runs of the same value in a vector.
# These *mostly* behave like normal R vectors.
#
# If not, cast with as(), or with as.logical(), as.numeric(), etc.

myrange <- as("chrI:1-5:+","GRanges")

if (strand(myrange) == "+") {
    print("Forward strand")
}

# The problem 
strand(myrange) == "+"

# Solution
if (as(strand(myrange) == "+", "logical")) {
    print("Forward strand")
}


# ==============================================================================
# Genomes and genes for model organisms are available as Bioconductor packages
#
# BSgenome. ...  - Genomes as DNAStringSet-like objects
# TxDB. ...      - Genes, transcripts, exons, CDS as GRanges
# org. ...       - Translation of gene names, assignment to GO categories
#
# Further data may be available using AnnotationHub.
#
# We will be using packages for yeast in this part of the workshop.


library(BSgenome.Scerevisiae.UCSC.sacCer3)
# A yeast genome object has been loaded.
# Actual data is only loaded from disk as needed.
genome <- Scerevisiae
seqinfo(genome)
genome$chrM

library(TxDb.Scerevisiae.UCSC.sacCer3.sgdGene)
# An object referring to a yeast transcriptome database has been loaded.
# Actual data is only loaded from disk as needed.
txdb <- TxDb.Scerevisiae.UCSC.sacCer3.sgdGene

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



cds_list <- cdsBy(txdb, "tx", use.names=TRUE)

cds_list
class(cds_list)
table( lengths(cds_list) )
cds_list[ lengths(cds_list) >= 2 ]

seqinfo(cds_list)
genome(cds_list)

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


library(rtracklayer)
export(cds_ranges, "cds_ranges.gff3")
export(start_codons, "start_codons.gff3")

# Examine these files in IGV
# https://software.broadinstitute.org/software/igv/download
#
# Load "sacCer3" genome (dropdown box in top left of window, may need to choose "More...")
# File/Load from file... to load GFF files.
# YBR111W-A has introns, check we have handled this correctly.

# See also rtracklayer::import() to load many file formats.
import("cds_ranges.gff3")



