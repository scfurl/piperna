#!/usr/bin/env Rscript
require(docopt)
'Summarize.

Usage:
   summarize.R [options]

Options:
   -r --runsheet Runsheet location [default: ./runsheet.csv].
   -t --threads Number of threads [default: 8].
   -y --yield_size Number of records to yield each time the file is read in millions [default: 8].
   -m --mode Summarize overlap mode [default: Union].
   -o --output Runsheet location [default: ./summarizedExperiment.RDS].
   -b --by Summarize by gene_id, gene_name, exons, transcripts.  This string should be found in gtf file [default: gene_id].
   -s --geneshort An alternative name for the "by" argument.  This string should be found in gtf file [default: gene_name].
' -> doc

opts <- docopt(doc)
require("GenomicAlignments")
require("Rsamtools")
require("rtracklayer")
require("BiocParallel")


registered()

df <- read.csv(opts$r, stringsAsFactors=F)
threads <- as.numeric(opts$t)
register(MulticoreParam(workers=threads))
if("fastq1" %in% colnames(df)){
    if("fastq2" %in% colnames(df)){
        singleend=FALSE
    }
}
if("fastqs" %in% colnames(df)){
	singleend=TRUE
}
if(df$software[1]=="STAR"){
	mode <- as.character(opts$m)
	gtffile <- df$gtf[1]
	yield_size <- as.numeric(opts$y)*1e6
	if(is.null(gtffile)){stop("GTF file not found: NULL input")}
	if(!file.exists(gtffile)){stop(paste0("GTF file not found: ", gtffile))}
	message(paste0("Loading gtf file: ", gtffile))
	message(paste0("Summarizing data by: ", opts$b))
	message(paste0("Adding an additonal annotation column (gene_short_name): ", opts$geneshort))
	# txdb <- makeTxDbFromGFF(gtffile, format="gtf")
	# ebg <- exonsBy(txdb, by="gene")
	gff0 <- import(gtffile)
	idx <- mcols(gff0)$type == "exon"
	genes<- split(gff0[idx], mcols(gff0[idx])[[opts$by]])
	filenames <- file.path(paste0(df$output, "/Aligned.sortedByCoord.out.bam"))
	df$filenames <- filenames
	if(!all(file.exists(filenames))){stop("All Bam files not found")}
	message("All Bam files in runsheet found")
	bamfiles <- BamFileList(filenames, yieldSize=yield_size)
	names(bamfiles) <- filenames
	message(paste0("Summarizing Overlap: Using ", threads, " threads"))
	se <- summarizeOverlaps(features=genes,
	                        reads=bamfiles,
	                        mode=mode,
	                        singleEnd=singleend,
	                        ignore.strand=TRUE)
	message(paste0("Saving data as: ", opts$o))
	mcols(se)[[opts$by]]=mcols(gff0)[[opts$by]][match(rownames(se), mcols(gff0)[[opts$by]])]
	#all(mcols(se)[[opts$by]]==rownames(se))
	mcols(se)[['gene_short_name']]<-mcols(gff0)[[opts$geneshort]][match(rownames(se), mcols(gff0)[[opts$by]])]
	colnames(se)<-df$filenames
	colData(se)<-DataFrame(df)
	#saveRDS(se, "SE_32cores.RDS")
	saveRDS(se, opts$o)

}
if(df$software[1]=="kallisto"){
	require("tximport")
	require("rhdf5")
	df$filenames <- file.path(df$output, "abundance.h5")
	filenames <-df$filenames
	names(filenames)<- basename(df$output)
	if(!all(file.exists(filenames))){stop("All h5 files not found")}
	message("All h5 files in runsheet found")
	message(paste0("Running tximport using ", threads, " threads"))
	txi.kallisto <- tximport(filenames, type = "kallisto", txOut = TRUE)
	saveRDS(txi.kallisto, opts$o)

}