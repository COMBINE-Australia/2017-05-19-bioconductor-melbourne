## ------------------------------------------------------------------------
library(knitr)
# To output a pdf or html document simply substitute github_document with pdf_document or html_document in the header of the file

# Working directory
dir <- "~/Documents/varie/Combine_tutorial"

## ----message = FALSE, eval = F-------------------------------------------
## # Bioconductor packages
## source("https://bioconductor.org/biocLite.R")
## biocLite(c("GenomicRanges","rtracklayer","Rsamtools","GenomicAlignments","TxDb.Hsapiens.UCSC.hg19.knownGene"))
## 
## # CRAN packages
## install.packages("ggplot2")

## ------------------------------------------------------------------------
library("ggplot2")
library("GenomicRanges")
library("Rsamtools")
library("GenomicAlignments")
library("rtracklayer")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")

## ------------------------------------------------------------------------
# Sample 1
gr_S1 <- GenomicRanges::GRanges(
seqnames = Rle(c("chr1", "chr2"), c(3, 2)),   
ranges = IRanges::IRanges(start = c(5,8,20,8,18), 
                 end = c(11,15,26,16,21), 
                 names = c(paste("Peak_",1:5,sep = ""))), 
strand = S4Vectors::Rle(strand(c("*")), c(5)))

gr_S1


## ------------------------------------------------------------------------
# Gene annotation
genes <- GenomicRanges::GRanges(
seqnames = S4Vectors::Rle(c("chr1", "chr2"), c(2, 2)),  
ranges = IRanges::IRanges(start = c(7,17,7,23), 
                 end = c(15,23,14,26), 
                 names = c(paste("Gene_",1:4,sep = ""))), 
strand = S4Vectors::Rle(strand(c("+","+")), c(2,2)))

genes

## ----fig.cap = "Graphic representation of the GRanges objects created above."----
gr <- as.data.frame(rbind(as.data.frame(gr_S1), as.data.frame(genes)))
gr$rangeID <- c(names(gr_S1),names(genes))
gr$Sample <- c(rep("Peaks",length(gr_S1)), rep("Genes",length(genes)))
ggplot(data = gr, aes(xmin = start, xmax = end, ymin = 0, ymax = 1)) + 
  geom_rect( aes(fill = rangeID),alpha = 0.4) + facet_wrap(~ Sample + seqnames) + theme_bw()

## ------------------------------------------------------------------------
# Combine overlapping ranges in the peak object
gr_S1_reduced <- GenomicRanges::reduce(gr_S1)
gr_S1_reduced
length(gr_S1)
length(gr_S1_reduced)

## ------------------------------------------------------------------------
?GenomicRanges::findOverlaps
overlaps <- findOverlaps(query = gr_S1, subject = genes)
overlaps
# Query the output
queryHits(overlaps)
subjectHits(overlaps)
subjectLength(overlaps)

## ------------------------------------------------------------------------
# You can allow for a maxgap to be ignored 
overlaps <- findOverlaps(gr_S1,genes, maxgap = 5)
overlaps
# You can specify a minimum overlap
overlaps <- findOverlaps(gr_S1,genes, minoverlap = 5)
overlaps
# You can specify a minimum overlap
overlaps <- findOverlaps(gr_S1,genes, type = "start")
overlaps
overlaps <- findOverlaps(gr_S1,genes, type = "end")
overlaps
overlaps <- findOverlaps(gr_S1,genes, type = "within")
overlaps

## ------------------------------------------------------------------------
?GenomicRanges::countOverlaps

## ------------------------------------------------------------------------
N_overlaps <- GenomicRanges::countOverlaps(gr_S1,genes)
N_overlaps
# You can play around with the same options as findOverlaps()

## ------------------------------------------------------------------------
?GenomicRanges::nearest

## ------------------------------------------------------------------------
GenomicRanges::nearest(x = gr_S1, subject = genes)
GenomicRanges::nearest(x = gr_S1,subject = genes, select = "all")
GenomicRanges::nearest(gr_S1)
GenomicRanges::nearest(gr_S1,gr_S1)

## ------------------------------------------------------------------------
?GenomicRanges::distance

## ------------------------------------------------------------------------
GenomicRanges::distance(x = gr_S1[1],y = genes)
GenomicRanges::distance(x = gr_S1[1:4],y = genes)
GenomicRanges::distance(x = gr_S1, y = genes)

## ------------------------------------------------------------------------
?GenomicRanges::distanceToNearest

## ------------------------------------------------------------------------
GenomicRanges::distanceToNearest(x = gr_S1, subject = genes)

## ------------------------------------------------------------------------
# Sample 1
gr_S1 <- GRanges(
seqnames = Rle(c("chr1", "chr2"), c(3, 2)),   
ranges = IRanges(start = c(5,8,20,8,18), 
                 end = c(11,15,26,16,21), 
                 names = c(paste("Peak_",1:5,sep = ""))), 
strand = Rle(strand(c("*")), c(5))) 

gr_S1

# Sample 2
gr_S2 <- GRanges(
seqnames = Rle(c("chr2", "chr3"), c(3, 5)),   
ranges = IRanges(start = c(1:8), 
                 width = 10, 
                 names = c(paste("Peak_",1:8,sep = ""))), 
strand = Rle(strand(c("*")), c(8)))

gr_S2

# GRanges List
list_ranges <- GRangesList(Sample1 = gr_S1, Sample2 = gr_S2)


## ------------------------------------------------------------------------
names(list_ranges)
length(list_ranges)
seqnames(list_ranges)
strand(list_ranges)
ranges(list_ranges)
start(list_ranges)
end(list_ranges)
width(list_ranges)
unlist(list_ranges)

## ------------------------------------------------------------------------
elementNROWS(list_ranges) 

## ------------------------------------------------------------------------
# Sample 3
gr_S3 <- GRanges(
seqnames = Rle(c("chr1", "chr2"), c(3, 2)),   
ranges = IRanges(start = 20:24, 
                 width = 8, 
                 names = c(paste("Peak_",1:5,sep = ""))), 
strand = Rle(strand(c("*")), c(5)))

gr_S1

# Sample 4
gr_S4 <- GRanges(
seqnames = Rle(c("chr2", "chr3"), c(3, 5)),   
ranges = IRanges(start = 20:27, 
                 width = 10, 
                 names = c(paste("Peak_",1:8,sep = ""))), 
strand = Rle(strand(c("*")), c(8)))

gr_S4

# Second GRanges List
list_ranges2 <- GRangesList(Sample3 = gr_S3, Sample4 = gr_S4)

## ------------------------------------------------------------------------
append_lists <- c(list_ranges,list_ranges2)

## ------------------------------------------------------------------------
addCols <- lapply(append_lists, function(x){
elementMetadata(x) <- data.frame(NumberReads = rbinom(length(x),size = 100, prob = 0.5))
return(x)
})
class(addCols)
addCols <- GRangesList(addCols)

## ------------------------------------------------------------------------
# As a list
addCols[[1]]
addCols["Sample1"]
# As a data.frame object
addCols[1,"NumberReads"]
addCols["Sample1","NumberReads"]

## ------------------------------------------------------------------------
lapply(addCols, length)
sapply(addCols, length)
lapply(addCols, start)

## ------------------------------------------------------------------------
# Concatenate and reduce the GRangesList
Reduce(c,addCols)

## ------------------------------------------------------------------------
# Download CpG islands from UCSC
session <- browserSession("UCSC")
query <- ucscTableQuery(session, "CpG Islands",GRangesForUCSCGenome("hg19", "chr21"))
cpg_islands <- getTable(query)
# Download genes 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

## ------------------------------------------------------------------------
# 1. Convert to GRanges
cpg_islandsGR <- as(cpg_islands,"GRanges")

# 2. Get the transcript for every gene from the `txdb`
transcriptHg19 <- transcriptsBy(txdb,by = "gene")

# 3. Only select chr21
transcriptHg19 <- unlist(transcriptHg19)
chr21 <- transcriptHg19[seqnames(transcriptHg19) == "chr21",]

# 4. Find the nearest CpG island in `cpg_islandsGR` to every transcript in `transcriptHg19` as well as their distance.
chr21_cpg <- distanceToNearest(chr21, cpg_islandsGR)

# 5. On average, every CpG island is close to how many genes?
mean(table(subjectHits(chr21_cpg)))


## ------------------------------------------------------------------------
# Import
rtracklayer::export(chr21, file.path(dir,"transcriptHg19.gff"))
rtracklayer::export(chr21, file.path(dir,"transcriptHg19.gtf"))

importHg19 <- rtracklayer::import(file.path("transcriptHg19.gff"))
class(importHg19)
importHg19

## ------------------------------------------------------------------------
# 1. Create  GRanges object made of 20 ranges
exampleGR <- GRanges(seqnames = Rle("chr1",20), IRanges(start = 1:20, width = rbinom(20,size = 100, prob = 0.6)))
export(exampleGR, "./exampleGR.gff")
# 2. Read the file back into R.
exampleGR_import <- import("./exampleGR.gff")
# 3. Print the number of rows of the `GRanges` object just imported and create a histogram of the widths of the ranges.
length(exampleGR_import)
hist(width(exampleGR_import))

## ----eval = FALSE, prompt = FALSE, engine="bash"-------------------------
## variableStep chrom=chr2
## 300701  12.5
## 300702  12.5
## 300703  12.5
## 300704  12.5
## 300705  12.5

## ----eval = FALSE, prompt = FALSE, engine="bash"-------------------------
## variableStep chrom=chr2 span=5
## 300701  12.5

## ----eval = FALSE, prompt = FALSE, engine="bash"-------------------------
## fixedStep chrom=chr3 start=400601 step=100
## 11
## 22
## 33

## ------------------------------------------------------------------------
wig_path <- file.path(dir, "exampleWIG.wig")
import_wig <- rtracklayer::import(con = wig_path,seqinfo = Seqinfo(genome="mm10"))
import_wig
?GenomeInfoDb::Seqinfo
seqinfo(import_wig)

## ------------------------------------------------------------------------
# Output file
bigwig_path <- file.path(dir,"example_wigToBigWig.bw")
rtracklayer::wigToBigWig(wig_path, dest = bigwig_path, seqinfo = seqinfo(import_wig))
import_bigwig <- rtracklayer::import(con = bigwig_path)

## ------------------------------------------------------------------------
which_region <- GRanges(seqnames = "chr3", IRanges(start = 1, end = 500000))
import_bigwig_region <- rtracklayer::import(con = bigwig_path, which = which_region)

## ------------------------------------------------------------------------
library(AnnotationHub)
?AnnotationHub
ahub <- AnnotationHub()

ahub.chain <- subset(ahub, rdataclass == "ChainFile" & species == "Homo sapiens")
query(ahub.chain, c("hg18", "hg19"))

chain <- ahub.chain[ahub.chain$title == "hg19ToHg18.over.chain.gz"]
chain <- chain[[1]]
gr.hg18 <- rtracklayer::liftOver(import_bigwig, chain)
gr.hg18


## ----eval = FALSE, engine="bash"-----------------------------------------
## 
## D00626:239:CAFMCANXX:4:1201:16222:93271 403  chr1    17018   3  38M177N62M      =       16965   -330 TGGCCCAGGTCTGGCACATAGAAGTAGTTCTCTGGGACCTGCTGTTCCAGCTGCTCT  GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
## RG:Z:9_CAFMCANXX_L004   NH:i:2  HI:i:2  NM:i:0  nM:i:1  AS:i:200
## 

## ------------------------------------------------------------------------
# Path to the bamfile. In the same folder there has to be the chr1.bam.bai file
bamfile <- file.path(dir,"chr21.bam")
# Define parameters
which <- GRanges(seqnames = rep("chr21",3), 
                 IRanges(start = c(9827000, 16267000, 15890000), 
                         width = 10000)) 
what <- c("rname", "strand", "pos", "qwidth", "seq", "isize")
# ?scanBamFlag
flag = scanBamFlag(isDuplicate = FALSE, isProperPair = TRUE, isPaired = TRUE)

param <- ScanBamParam(which = which, what = what, flag = flag)
# load reads into R
reads = scanBam(bamfile, param = param)

# Explore the reads
names(reads)
region1 <- reads$`chr21:9827000-9836999`
names(region1)
head(region1$rname)
head(region1$pos)
head(region1$qwidth)
head(region1$seq)
head(region1$isize)

## ------------------------------------------------------------------------
readsGA <- GenomicAlignments::readGAlignments(bamfile, param = param)

## ------------------------------------------------------------------------
bins <- GRanges(seqnames = rep("chr21",11), 
                 IRanges(start = seq(9827000,9827000+10000, by = 1000), 
                         width = 1000)) 
# Create vector of 0-1 positions
positions <- region1$pos
vector_positions <- rep(0,max(positions))
vector_positions[positions] <- 1
# Views Object
rangeViews <- Views(vector_positions, start = seq(min(positions),max(positions), length.out = 10), width = 60)
# viewSums
bin_sum <- viewSums(rangeViews)
# Add extra column to the ranges object
binned_counts <- ranges(rangeViews)
values(binned_counts) <- data.frame(ReadCounts = bin_sum)

## ------------------------------------------------------------------------
# 1.
# Path to the bamfile. In the same folder there has to be the chr1.bam.bai file
bamfile <- file.path(dir,"chr21.bam")
# Define parameters
which <- GRanges(seqnames = rep("chr21",3), 
                 IRanges(start = c(9827000, 16267000, 15890000), 
                         width = 10000)) 
what <- c("rname", "strand", "pos", "qwidth", "seq", "isize")
# ?scanBamFlag
flag = scanBamFlag(isDuplicate = NA, isProperPair = TRUE, isPaired = TRUE)

param <- ScanBamParam(which = which, what = what, flag = flag)
# load reads into R
reads = scanBam(bamfile, param = param)

# 2. 
hist(reads$`chr21:16267000-16276999`$isize)

