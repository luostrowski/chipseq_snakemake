library(GenomicFeatures)
library(ChIPseeker)
library(writexl)

# Get command line arguments
args <- commandArgs(TRUE)
gtf <- args[1]
data_source <- args[2]
organism <- gsub("_", " ", args[3])
infile <- args[4]

# set outdir prefix for exported files
outdir <- "data/annot_peaks/"

# get sample name from infile
tmp <- gsub("data/macs2/", "", infile)
sample_name <- gsub("_peaks.narrowPeak", "", tmp)

#---------------#
# make a TxDB object from the annotation gtf file
txdb <- makeTxDbFromGFF(gtf,
    format = "gtf",
    data_source,
    organism)

# read in the region file to be annotated
#   here, SIG diff peaks txt file from diffbind
mydat <- makeGRangesFromDataFrame(read.delim(infile, skip = 1L, header = FALSE),
    keep.extra.columns = TRUE,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3",
    ignore.strand = TRUE)

# annotate via Chipseeker
res <- annotatePeak(mydat,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    level = "transcript",
    assignGenomicAnnotation = TRUE,
    overlap = "TSS")

# output annotations
write.table(as.data.frame(res),
    paste0(outdir, sample_name, ".annot.txt"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE)

# anno pie
pdf(paste0(outdir, sample_name, ".anno_pie_chart.pdf"))
plotAnnoPie(res)
dev.off()

# vennpie
pdf(paste0(outdir, sample_name, ".venn_pie_chart.pdf"))
#par(mar=c(0.5,2,1,3))
vennpie(res)
dev.off()
