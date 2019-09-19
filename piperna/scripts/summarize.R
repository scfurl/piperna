#!/usr/bin/env Rscript
require(docopt)
'Summarize.

Usage:
   summarize.R [options]

Options:
   -r --runsheet Runsheet location [default: ./runsheet.csv].
   -o --output Runsheet location [default: ./summarizedExperiment.RDS].
' -> doc

opts <- docopt(doc)
require("GenomicAlignments")
require("BiocParallel")
require("GenomicFeatures")
require("Rsamtools")

register(MulticoreParam(4))

df <- read.csv(opts$r, stringsAsFactors=F)
gtffile <- df$gtf[1]
if(is.null(gtffile)){stop("GTF file not found: NULL input")}
if(!file.exists(gtffile)){stop(paste0("GTF file not found: ", gtffile))}
message(paste0("Loading gtf file: ", gtffile))
txdb <- makeTxDbFromGFF(gtffile, format="gtf")
ebg <- exonsBy(txdb, by="gene")
if("fastq1" %in% colnames(df)){
    if("fastq2" %in% colnames(df)){
        singleend=FALSE
    }
}
if("fastqs" %in% colnames(df)){
	singleend=TRUE
}
if(df$software[1]=="STAR"){
	filenames <- file.path(paste0(df$output, "Aligned.sortedByCoord.out.bam"))
	if(!all(file.exists(filenames))){stop("All Bam files not found")}
	bamfiles <- BamFileList(filenames, yieldSize=20000000)

	se <- summarizeOverlaps(features=ebg,
	                        reads=bamfiles,
	                        mode="Union",
	                        singleEnd=singleend,
	                        ignore.strand=TRUE)
	saveRDS(se, opts$o)
}