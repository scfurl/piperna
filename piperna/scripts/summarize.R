#!/usr/bin/env Rscript
require(docopt)
'Usage:
   summarize.R [-r runsheet -o <output>]

Options:
   -r Runsheet [default: ./runsheet.csv]
   -o Output file [default: ./SummarizedExperiment.RDS]

 ]' -> doc

opts <- docopt(doc)

require("GenomicAlignments")
require("BiocParallel")
require("GenomicFeatures")
require("Rsamtools")

register(MulticoreParam(4))

df <- read.csv(opts$r)
gtffile <- df$gtf[1]
message(paste0("Loading gtf file: ", gtffile))
txdb <- makeTxDbFromGFF(gtffile, format="gtf")
ebg <- exonsBy(txdb, by="gene")
if("fastq1" %in% colnames(dtf)){
    if("fastq2" %in% colnames(dtf)){
        singleend=FALSE
    }
}
if("fastqs" %in% colnames(dtf)){
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