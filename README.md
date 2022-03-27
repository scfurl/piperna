[![PyPI](https://img.shields.io/pypi/v/simplesam.svg?)](https://pypi.org/project/piperna/)
<!-- [![Build Status](https://travis-ci.org/mdshw5/simplesam.svg?branch=master)](https://travis-ci.org/mdshw5/simplesam) -->
[![Documentation Status](https://readthedocs.org/projects/piperna/badge/?version=latest)](https://piperna.readthedocs.io/en/latest/?badge=latest)

# piperna
==========

A python wrapper for processing of bulk RNA seq data

## Requirements

1. Python > 3.5 (piperna uses the 'six' package but will attempt to install if not already installed)
2. Computing cluster with PBS or SLURM
3. Modules installed for python, STAR or kallisto, 
4. R Modules for GenomicAlignments, rtracklayer, and Rsamtools if running SUMMARIZE

## Installation

Installation can probably be done correctly many different ways.  Here are the methods that have worked for us.  We recommend that piperna be installed with pipx.

**At SCRI do the following**
```bash
module load python
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install --include-deps --pip-args '--trusted-host pypi.org --trusted-host files.pythonhosted.org' piperna
```


**At the FHCRC do the following...**
```bash
module load Python/3.6.7-foss-2016b-fh1
python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install --include-deps piperna
```

You should then be able to test installation by calling piperna.  After running the folllowing, you should see the help screen displayed.

```bash
piperna
```



## Usage

```bash
piperna usage: A wrapper for running RNASeq Alignment [-h] [--bam_folder BAM_FOLDER]
                                              [--fastq_folder FASTQ_FOLDER]
                                              [--organized_by {folder,file}]
                                              [--genome_key GENOME_KEY]
                                              [--split_char SPLIT_CHAR]
                                              [--R1_char R1_CHAR]
                                              [--R2_char R2_CHAR] [--ext EXT]
                                              [--select SELECT]
                                              [--sample_flag SAMPLE_FLAG]
                                              [--runsheet RUNSHEET]
                                              [--typeofseq {single,pe}]
                                              [--software {STAR,kallisto}]
                                              [--output OUTPUT] [--debug]
                                              [--cluster {PBS,SLURM}]
                                              [--user USER]
                                              [--threads THREADS]
                                              [--gb_ram GB_RAM]
                                              [--additional_header ADDITIONAL_HEADER]
                                              [--mfl MFL] [--sfl SFL]
                                              [--count] [--install INSTALL]
                                              [--outSAMtype OUTSAMTYPE]
                                              [--addSTARstring ADDSTARSTRING]
                                              [--log_prefix LOG_PREFIX]
                                              [--flow_cell_folders FLOW_CELL_FOLDERS]
                                              [--verbose]
                                              {MAKERUNSHEET,ALIGN,SUMMARIZE,CONCATFASTQ,GENOMESFILE,ENVIRONSFILE,UNBAM}

positional arguments:
  {MAKERUNSHEET,ALIGN,SUMMARIZE,CONCATFASTQ,GENOMESFILE,ENVIRONSFILE,UNBAM}
                        a required string denoting segment of pipeline to run.
                        1) "MAKERUNSHEET" - to parse a folder of fastqs; 2)
                        "ALIGN" - to perform alignment using STAR or KALLISTO;
                        3) "SUMMARIZE" - to summarize and count reads, 4)
                        "CONCATFASTQ" - function to concatenate fastq files
                        -i.e. for SRA upload; 5) "GENOMESFILE" - print
                        location of and cat genomes.json file; 6)
                        "ENVIRONSFILE" - print location of and cat
                        environs.json file; 7) "UNBAM" - convert bam to fastq

optional arguments:
  -h, --help            show this help message and exit
  --bam_folder BAM_FOLDER, -bf BAM_FOLDER
                        For UNBAM only: Pathname of bam locations
  --fastq_folder FASTQ_FOLDER, -fq FASTQ_FOLDER
                        For MAKERUNSHEET only: Pathname of fastq folder
  --organized_by {folder,file}, -b {folder,file}
                        Option to specify how fastq or unbam folder is
                        organized
  --genome_key GENOME_KEY, -gk GENOME_KEY
                        For MAKERUNSHEET only: abbreviation to use "installed"
                        genomes in the runsheet (See README.md for more
                        details
  --split_char SPLIT_CHAR, -sc SPLIT_CHAR
                        Character by which to split the fastqfile name into
                        samples, OPTIONAL and for MAKERUNSHEET only
  --R1_char R1_CHAR, -r1c R1_CHAR
                        Character by which to split the fastqfile name into
                        read1, OPTIONAL and for MAKERUNSHEET only
  --R2_char R2_CHAR, -r2c R2_CHAR
                        Character by which to split the fastqfile name into
                        read2, OPTIONAL and for MAKERUNSHEET only
  --ext EXT, -e EXT     suffix of fastq files, OPTIONAL and for MAKERUNSHEET
                        only
  --select SELECT, -s SELECT
                        To only run the selected row in the runsheet, OPTIONAL
                        and for MAKERUNSHEET only
  --sample_flag SAMPLE_FLAG, -f SAMPLE_FLAG
                        FOR MAKERUNSHEET only string to identify samples of
                        interest in a fastq folder
  --runsheet RUNSHEET, -r RUNSHEET
                        tab-delim file with sample fields as defined in the
                        script. - or name of runsheet to save if using
                        MAKERUNSHEET
  --typeofseq {single,pe}, -t {single,pe}
                        Type of sequencing performed - REQUIRED for
                        MAKERUNSHEET, UNBAM and CONCATFASTQ
  --software {STAR,kallisto}, -so {STAR,kallisto}
                        To set desired software, required and used for
                        MAKERUNSHEET only
  --output OUTPUT, -o OUTPUT
                        To set output path, required for MAKERUNSHEET, UNBAM;
                        OPTIONAL for SUMMARIZE-- default for SUMMARIZE is
                        ./SummarizedExperiment.RDS
  --debug, -d           To print commands (For testing flow)
  --cluster {PBS,SLURM}, -c {PBS,SLURM}
                        Cluster software. OPTIONAL Currently supported: PBS
                        and SLURM
  --user USER, -u USER  user for submitting jobs - defaults to username.
                        OPTIONAL
  --threads THREADS, -th THREADS
                        To set number of cores
  --gb_ram GB_RAM, -gb GB_RAM
                        To set gb_ram
  --additional_header ADDITIONAL_HEADER, -ah ADDITIONAL_HEADER
                        Additional bash header lines
  --mfl MFL, -mf MFL    Mean fragment length (kallisto ONLY)
  --sfl SFL, -sf SFL    SD fragment length (kallisto ONLY)
  --count, -co          Run Count (STAR Only)
  --install INSTALL, -i INSTALL
                        FOR GENOMESFILE: location of file to install as a new
                        genomes.json file, existing genomes.json will be
                        erased
  --outSAMtype OUTSAMTYPE, -st OUTSAMTYPE
                        To define type of SAM/BAM output (STAR Only)
  --addSTARstring ADDSTARSTRING, -a ADDSTARSTRING
                        Additional STAR arguments to be run on all jobs in
                        runsheet (STAR Only)
  --log_prefix LOG_PREFIX, -l LOG_PREFIX
                        Prefix specifying log files for piperna output from
                        henipipe calls. OPTIONAL
  --flow_cell_folders FLOW_CELL_FOLDERS, -fc FLOW_CELL_FOLDERS
                        For CONCATFASTQ only: Comma-seprated location of
                        flowcell folders - i.e. as output from CellRanger
                        REQUIRED for CONCATFASTQ
  --verbose, -v         Run with some additional ouput - not much though...
                        OPTIONAL
```


## Runsheet

The runsheet is the brains of the piperna workflow.  You can make a runsheet using the MAKERUNSHEET command.  This command will parse a directory of fastq folder (specified using the -fq flag; fastq files should be organized in subfolders named by sample) and will find fastq mates (R1 and R2 - Currently only PE sequencing is supported).  Running piperna MAKERUNSHEET will find and pair these fastqs for you and populate the runsheet with genome index locations (see below) and output filenames with locations as specified using the -o flag.  Note that piperna output will default to the current working directory if no location is otherwise specified.  There is an option for selecting only folders that contain a specific string (using the -sf flag).  *After generation of a runsheet (csv file), you should take a look at it in Excel or Numbers to make sure things look okay...*  Here are the columns that you can include.  Order is irrelevant.  Column names (headers) exactly as written below are required.

### Example Runsheet 

**absolute pathnames are required for runsheets**

| sample | index | fastq1 | fastq2 | output |  software  |     gtf    |
|--------|-------|--------|--------|--------|------------|------------|
|  mys1  |  path |  path  |  path  |  path  |   STAR     |    path    |
|  mys2  |  path |  path  |  path  |  path  |   STAR     |    path    |


* 'sample' name of the sample REQUIRED.  
* 'index' location of the indexed fasta file REQUIRED.  
* 'fastq1' a tab seperated string of filenames denoting location of all R1 files for a sample REQUIRED if paired end.  
* 'fastq2' a tab seperated string of filenames denoting location of all R2 files for a sample REQUIRED if paired end.  
* 'fastqs' a tab seperated string of filenames can be used for single end reads REQUIRED if single end.  
* 'output' name of the location for the aligned and sorted bam file.  
* 'software' either 'STAR' or 'kallisto'.  REQUIRED
* 'gtf' a location for annotation file in gtf format.  REQUIRED for SUMMARIZE.  

## Genomes and adding genome locations

you should have a previously indexed (by the software package of your choosing) location of your genome accessible to piperna.  This location is referred to in piperna as the 'index'.

piperna provides an easy way to add these locations to your system for repeated use using the --genome_key (-gk) option during MAKERUNSHEET commands.  A file called genomes.json can be found in the 'data' directory of the piperna install folder.  This file can be edited to include those locations you want to regularly put in the runsheet.  The following shows an example of a genomes.json file.  The files "top level" is a name that can be used in the --genome_key field (-gk) during runsheet generation to populate the columns of the runsheet with locations associated with that genome_key.  The 'default' key will be used when no genome_key is specified.

```json
{
    "default": {
      "fasta": "/path/path/hg38/STAR_index",
      "gtf": "/path/path/hg38/hg38.gtf"
    },
    "default_kallisto": {
        "fasta": "/path/path/hg38/kallisto_index",
        "gtf": "/path/path/hg38/hg38.gtf"
    }
}
```

## Cluster specific parameters

Piperna could be modified to run on any cluster  Templates are provided for SLURM and PBS

Location of the environs file can be found by running:

```bash
piperna ENVIRONSFILE
```

As an example the SLURM script generation  environs is shown below.  Broadly it's parameters is divided into three groups: popen, script_lines, and resources.

popen - an entry that lists the command to open a submission job ('sbatch' in this case)
script_lines - provides the header of the bash script
resources - gives job-specific parameters

```json
{
"SLURM": {
    "popen" : ["sbatch"],
    "script_lines": {
        "1" : ["#!/bin/bash\n#SBATCH --output=outtmp\n#SBATCH --error=errtmp", ""],
        "2" : ["#SBATCH --job-name=<--0-->", "JOB_NAME"],
        "3" : ["#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=<--0-->","THREADS"],
        "4" : ["#SBATCH --mem-per-cpu=<--0-->000","RAM"],
        "5" : ["<--0-->","HEADER"],
        "6" : ["{\necho '<--1-->';<--0-->;<--1-->;} 2>&1 | tee <--2-->", "MODULES|COMMAND|TEMP_LOG_FILE"],
        "7" : ["sed -e 's@^@[PIPERNA-<--0-->] JOB: <--1-->:\t\t@' <--2--> >> <--3-->", "TIME|JOB_NAME|TEMP_LOG_FILE|LOG_FILE"],
        "8" : ["rm <--0-->\n", "TEMP_LOG_FILE"]
    },
        "resources" : {
            "PIPERNA_STAR": {
                "ram" : "8",
                "threads" : "4",
                "modules" : "module load STAR/2.7.6a-foss-2019b"
            },
            "PIPERNA_KALLISTO": {
                "ram" : "8",
                "threads" : "4",
                "modules" : "module load kallisto"
            },
            "PIPERNA_SUMMARIZE": {
                "ram" : "16",
                "threads" : "8",
                "modules" : "module load R"
            },
            "PIPERNA_UNBAM": {
                "ram" : "2",
                "threads" : "8",
                "modules" : "module load SAMtools\nmodule load BEDTools"
            }
        }
    }
}
```

## Doing a piperna run

Say your fastqs live within within subfolders of a folder 'fastq' in the folder 'data'.  So if you were to...
```bash
cd /data/fastq
ls
```
... you'd get a bunch of folders, each of which would be filled with fastqs.  Each folder name should correspond to a sample name.

**To run piperna, do the following...**
1. Make a new output directory 'piperna'.
2. Go into that directory and make a runsheet pointing to the fastq folder i.e. the folder level above.  (at the command line, piperna is cool with either relative or absolute pathnames; but as stated earlier, absolute pathnames are required for the runsheet.)
3.  Optionally you can only select directories of fastq files that contain in their name the string denoted using the -sf flag.
4. After inspecting and completing the runsheet, run ALIGN, NORM, and SEACR.  

```bash
cd ..
mkdir piperna
cd piperna
piperna MAKERUNSHEET -fq ../fastq -sf MySampleDirectoriesStartWithThisString -o .
piperna ALIGN -r runsheet.csv
piperna SUMMARIZE -r runsheet.csv -o SummarizedExperiment.RDS
```


## Interfacing with DESeq2

**After running piperna, the SummarizedExperiment.RDS file can be input directly into DESeq2 like this**

```R
se<-readRDS(file.path(DATA_DIR, "summarizedExperiment.RDS"))
colnames(se)<-colData(se)$sample
se$group<-c("GroupA", "GroupA", "GroupB", "GroupB", "GroupB", "GroupA" ) #as an example
dds<-DESeqDataSet(se, design=~group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds <- DESeq(dds, parallel = T)
```


## Acknowledgements

Written by Scott Furlan.