import sys
import argparse
import logging
import getpass
import os
from . import piperna
"""
sys.path.append('/home/sfurlan/.local/pipx/venvs/piperna/lib/python3.7/site-packages/piperna')
import piperna
"""

POLL_TIME = 5
LOG_PREFIX = '[piperna]: '

# Set up a basic logger
LOGGER = logging.getLogger('something')
myFormatter = logging.Formatter('%(asctime)s: %(message)s')
handler = logging.StreamHandler()
handler.setFormatter(myFormatter)
LOGGER.addHandler(handler)
LOGGER.setLevel(logging.DEBUG)
myFormatter._fmt = "[piperna]: " + myFormatter._fmt


def run_piperna(args=None):
    parser = argparse.ArgumentParser('A wrapper for running RNASeq Alignment')
    parser.add_argument('job', type=str, choices=['MAKERUNSHEET', 'ALIGN', 'SUMMARIZE', 'CONCATFASTQ', 'GENOMESFILE', 'ENVIRONSFILE','UNBAM'], help='a required string denoting segment of pipeline to run.  1) "MAKERUNSHEET" - to parse a folder of fastqs; 2) "ALIGN" - to perform alignment using STAR or KALLISTO; 3) "SUMMARIZE" - to summarize and count reads, 4) "CONCATFASTQ" - function to concatenate fastq files -i.e. for SRA upload; 5) "GENOMESFILE" - print location of and cat genomes.json file; 6) "ENVIRONSFILE" - print location of and cat environs.json file; 7) "UNBAM" - convert bam to fastq')
    parser.add_argument('--bam_folder', '-bf', type=str, help='For UNBAM only: Pathname of bam locations')
    parser.add_argument('--fastq_folder', '-fq', type=str, help='For MAKERUNSHEET only: Pathname of fastq folder')
    parser.add_argument('--organized_by', '-b', type=str, choices=['folder', 'file'], default='folder', help='Option to specify how fastq or unbam folder is organized')
    parser.add_argument('--genome_key', '-gk', default="default", type=str, help='For MAKERUNSHEET only: abbreviation to use "installed" genomes in the runsheet (See README.md for more details')
    parser.add_argument('--split_char', '-sc', type=str, default="_R1_", help='Character by which to split the fastqfile name into samples, OPTIONAL and for MAKERUNSHEET only')
    parser.add_argument('--R1_char', '-r1c', type=str, default="_R1_", help='Character by which to split the fastqfile name into read1, OPTIONAL and for MAKERUNSHEET only')
    parser.add_argument('--R2_char', '-r2c', type=str, default="_R2_", help='Character by which to split the fastqfile name into read2, OPTIONAL and for MAKERUNSHEET only')
    parser.add_argument('--ext', '-e', type=str, default=".fastq.gz", help='suffix of fastq files, OPTIONAL and for MAKERUNSHEET only')
    parser.add_argument('--select', '-s', type=str, default=None, help='To only run the selected row in the runsheet, OPTIONAL and for MAKERUNSHEET only')
    parser.add_argument('--sample_flag', '-f', type=str, default="", help='FOR MAKERUNSHEET only string to identify samples of interest in a fastq folder')
    parser.add_argument('--runsheet', '-r', type=str, default = "runsheet.csv", help='tab-delim file with sample fields as defined in the script. - or name of runsheet to save if using MAKERUNSHEET')
    parser.add_argument('--typeofseq', '-t', type=str, default = "pe", choices=['single', 'pe'], help= 'Type of sequencing performed - REQUIRED for MAKERUNSHEET, UNBAM and CONCATFASTQ')
    parser.add_argument('--software', '-so', type=str, choices=['STAR', 'kallisto'], default="STAR", help='To set desired software, required and used for MAKERUNSHEET only')
    parser.add_argument('--output', '-o', type=str, default=".", help='To set output path, required for MAKERUNSHEET, UNBAM; OPTIONAL for SUMMARIZE-- default for SUMMARIZE is ./SummarizedExperiment.RDS')
    parser.add_argument('--debug', '-d', action='store_true', help='To print commands (For testing flow)')
    parser.add_argument('--cluster', '-c', type=str, default='SLURM', choices=['PBS', 'SLURM'], help='Cluster software.  OPTIONAL Currently supported: PBS and SLURM')
    parser.add_argument('--user', '-u', type=str, default='sfurlan', help='user for submitting jobs - defaults to username.  OPTIONAL')
    parser.add_argument('--threads', '-th', type=int, default=4, help='To set number of cores')
    parser.add_argument('--gb_ram', '-gb', type=int, default=8, help='To set gb_ram')
    parser.add_argument('--additional_header', '-ah', type=str, default=None, help='Additional bash header lines')
    parser.add_argument('--mfl', '-mf', type=int, default=400, help='Mean fragment length (kallisto ONLY)')
    parser.add_argument('--sfl', '-sf', type=int, default=20, help='SD fragment length (kallisto ONLY)')
    parser.add_argument('--count', '-co', action='store_true', default=True, help='Run Count (STAR Only)')
    parser.add_argument('--install', '-i', type=str, default=None, help='FOR GENOMESFILE: location of file to install as a new genomes.json file, existing genomes.json will be erased')
    parser.add_argument('--outSAMtype', '-st', type=str, default='BAM SortedByCoordinate', help='To define type of SAM/BAM output (STAR Only)')
    parser.add_argument('--addSTARstring', '-a', type=str, default=None, help='Additional STAR arguments to be run on all jobs in runsheet (STAR Only)')
    parser.add_argument('--log_prefix', '-l', type=str, default='piperna.log', help='Prefix specifying log files for piperna output from henipipe calls. OPTIONAL')
    parser.add_argument('--flow_cell_folders', '-fc', type=str, default="", help='For CONCATFASTQ only: Comma-seprated location of flowcell folders - i.e. as output from CellRanger REQUIRED for CONCATFASTQ')
    parser.add_argument('--verbose', '-v', default=False, action='store_true', help='Run with some additional ouput - not much though... OPTIONAL')

    args = parser.parse_args()
    """
    call='piperna CONCATFASTQ -fc /archive/furlan_s/seq/cellranger/181015-NHPTreg/HHJJ7BGX5/outs/fastq_path,/archive/furlan_s/seq/cellranger/181015-NHPTreg/HWVFMBGX3/outs/fastq_path'
    call='piperna MAKERUNSHEET -o . -fq /fh/scratch/delete90/furlan_s/mullighan/fastqs'
    args = parser.parse_args(call.split(" ")[1:])

    #log
    """
    if args.job=="UNBAM":
        if os.path.isabs(args.bam_folder) is False:
            if args.bam_folder == ".":
                args.bam_folder = os.getcwd()
            else :
                args.bam_folder = os.path.abspath(args.bam_folder)
        if os.path.exists(args.bam_folder) is False:
            raise ValueError('Path: '+args.bam_folder+' not found')
        if os.path.isabs(args.output) is False:
            if args.output == ".":
                args.output = os.getcwd()
            else :
                args.output = os.path.abspath(args.output)
        if os.path.exists(args.output) is False:
            raise ValueError('Path: '+args.output+' not found')
        LOGGER.info("Parsing bam folder - "+args.bam_folder+" ...")
        unbam_runsheet = piperna.unbam_prep(folder=args.bam_folder, output=args.output, typeofseq=args.typeofseq, \
            organized_by=args.organized_by)
        if args.select is not None:
            unbam_runsheet = [unbam_runsheet[i-1] for i in list(piperna.parse_range_list(args.select))]
        unbam_jobs = piperna.unbam(runsheet_data = unbam_runsheet, user=args.user, \
                debug=args.debug, threads=args.threads, additional_header=args.additional_header, gb_ram=str(args.gb_ram), log=args.log_prefix, \
                cluster=args.cluster)
        unbam_jobs.run_job()
        exit()


    if args.job=="GENOMESFILE":
        _ROOT = os.path.abspath(os.path.dirname(__file__))
        if args.install is None:
            GENOMES_JSON = os.path.join(_ROOT, 'data', 'genomes.json')
            print("Showing contents of genomes.json file located at:\n"+GENOMES_JSON+"\n\n\n\n")
            f = open(GENOMES_JSON, 'r')
            file_contents = f.read()
            print (file_contents+"\n\n\n\n")
            f.close()
        if args.install is not None:
            from shutil import copyfile
            args.install = os.path.abspath(args.install)
            copyfile(args.install, os.path.join(_ROOT, 'data', 'genomes.json'))
        exit()

    if args.job=="ENVIRONSFILE":
        _ROOT = os.path.abspath(os.path.dirname(__file__))
        if args.install is None:
            ENVIRONS_JSON = os.path.join(_ROOT, 'data', 'environs.json')
            print("Showing contents of environs.json file located at:\n"+ENVIRONS_JSON+"\n\n\n\n")
            f = open(ENVIRONS_JSON, 'r')
            file_contents = f.read()
            print (file_contents+"\n\n\n\n")
            f.close()
        if args.install is not None:
            from shutil import copyfile
            args.install = os.path.abspath(args.install)
            copyfile(args.install, os.path.join(_ROOT, 'data', 'environs.json'))
        exit()

    #log
    if args.debug == False:
        LOGGER.info("Logging to %s... examine this file if samples fail." % args.log_prefix)

    #deal with user
    if args.user is None:
        args.user = getpass.getuser()

    #deal with paths
    if args.job=="MAKERUNSHEET":
        if os.path.isabs(args.fastq_folder) is False:
            if args.fastq_folder == ".":
                args.fastq_folder = os.getcwd()
            else :
                args.fastq_folder = os.path.abspath(args.fastq_folder)
        if os.path.exists(args.fastq_folder) is False:
            raise ValueError('Path: '+args.fastq_folder+' not found')
        if os.path.isabs(args.output) is False:
            if args.output == ".":
                args.output = os.getcwd()
            else :
                args.output = os.path.abspath(args.output)
        if os.path.exists(args.output) is False:
            raise ValueError('Path: '+args.output+' not found')

    if args.job != "MAKERUNSHEET" and args.job != "CONCATFASTQ":
        if os.path.exists(args.runsheet) is False:
            raise ValueError('Path: '+args.runsheet+' not found')

    if args.job=="MAKERUNSHEET":
        LOGGER.info("Parsing fastq folder - "+args.fastq_folder+" ...")
        piperna.make_runsheet(folder=args.fastq_folder, output=args.output, typeofseq=args.typeofseq, \
            genome_key=args.genome_key, sample_flag = args.sample_flag, strsplit= args.split_char, ext=args.ext,  r1_char=args.R1_char, \
            r2_char=args.R2_char, fname = args.runsheet, software=args.software, organized_by=args.organized_by)
        exit()

    if args.job=="CONCATFASTQ":
        args.flow_cell_folders = os.path.abspath(args.flow_cell_folders)
        args.output = os.path.abspath(args.output)
        LOGGER.info("Concatenating fastqs in the folder(s) - "+" and ".join(args.flow_cell_folders.split(","))+" ...")
        concatfastqjob = piperna.concatfastq(folders=args.flow_cell_folders, output=args.output, typeofseq=args.typeofseq, \
            cluster=args.cluster, threads = 1, log=args.log_prefix, user=args.user, debug=args.debug)
        concatfastqjob.run_job()
        exit()

    #parse and chech runsheet
    args.runsheet = os.path.abspath(args.runsheet)
    parsed_runsheet = list(piperna.parse_runsheet(args.runsheet))
    piperna.check_runsheet(args, parsed_runsheet, verbose=args.verbose)
    #parsed_runsheet = piperna.parse_runsheet(args.runsheet)
    if args.select is not None:
        parsed_runsheet = [parsed_runsheet[i-1] for i in list(piperna.parse_range_list(args.select))]


    if args.job == "ALIGN":
        if args.software == "STAR":
            Star = piperna.Star(runsheet_data = list(parsed_runsheet), user=args.user, \
                debug=args.debug, threads=args.threads, additional_header=args.additional_header, gb_ram=str(args.gb_ram), log=args.log_prefix, \
                count=args.count, out_sam_type=args.outSAMtype, \
                global_add_STAR_string=args.addSTARstring, cluster=args.cluster)
            Star.run_job()

        if args.software == "kallisto":
            kallisto = piperna.kallisto(runsheet_data = list(parsed_runsheet), user=args.user, \
                debug=args.debug, threads=args.threads, additional_header=args.additional_header, gb_ram=str(args.gb_ram), log=args.log_prefix, \
                mfl = args.mfl, sfl = args.sfl, cluster=args.cluster)
            kallisto.run_job()

    if args.job == "SUMMARIZE":
        if os.path.isabs(args.output) is False:
            if args.output == ".":
                args.output = os.getcwd()
            else :
                args.output = os.path.abspath(args.output)
        summarize = piperna.summarize(runsheet = args.runsheet, user=args.user, \
                debug=args.debug, threads=args.threads, additional_header=args.additional_header, log=args.log_prefix, gb_ram=str(args.gb_ram), \
                cluster=args.cluster, output = args.output)
        summarize.run_job()
