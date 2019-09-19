#!/usr/bin/python
# PBS cluster job submission in Python for RNASeq Alignment
#v 0.04

# version notes:
# 0.04 made the wrapper agnostic to alignment software; added make runsheet column
# 0.03 added options for STAR arguments to be passed through a runsheet column (localSTARstr) or from command line (globalSTARstr).

# RNASeq wrapper....

# runsheet columns
# 1=sample (str)
# 2=index, (path)
# 3=output, (path)
# 4=fastqs, (tab separated string) - single reads only
# 5=fastq1, read1 (tab separated string) - PE reads only
# 6=fastq2, read2 (tab separated string) - PE reads only
# 7=software
# 8=localSTARStr = additional args for STAR given as a string i.e. '--seedPerWindowNmax 30'


# programs
#STAR - "star"

import os
from subprocess import Popen, PIPE
import argparse
import time
import sys
import csv
import re
import string
import random
from itertools import chain, compress
import json

_ROOT = os.path.abspath(os.path.dirname(__file__))
#GENOMES_JSON = os.path.join('/Users/sfurla/Box Sync/PI_FurlanS/computation/develop/piperna/piperna/data/genomes.json')
GENOMES_JSON = os.path.join(_ROOT, 'data', 'genomes.json')
SUMMARIZE_SCRIPT = os.path.join(_ROOT, 'scripts', 'summarize.R')

class SampleFactory:
    def __init__(self, *args, **kwargs):
        self.user = kwargs.get('user')
        self.cluster = kwargs.get('cluster')
        self.runsheet_data = kwargs.get('runsheet_data')
        self.debug = kwargs.get('debug')
        self.log_name = kwargs.get('log')
    def __call__():
        pass
        # if self.debug == False:
        #   open(self.log, 'w')

    def id_generator(self, size=6, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    def get_runmode(self):
        mode = ""
        if "fastq1" in self.runsheet_data[0]:
            if "fastq2" in self.runsheet_data[0]:
                mode = "pe"
        if "fastqs" in self.runsheet_data[0]:
            mode = "single"
        return mode
    def generate_job(self):
        job_string=[]
        torun = len(self.runsheet_data)
        for i in range(torun):
            log_file = self.id_generator()
            job_name = self.job + "_" + self.runsheet_data[i]['sample']
            command = self.command[i]
            if self.cluster=="PBS":
                #to_append = '#!/bin/bash\n#PBS -N %s\n#PBS -l %s\n#PBS -j oe\n#PBS -o $PBS_O_WORKDIR/%s\n#PBS -A %s\ncd $PBS_O_WORKDIR\n%s' % (job_name, self.processor_line, self.log_name, self.user, command)
                #to_append = '#!/bin/bash\n#PBS -N %s\n#PBS -l %s\n#PBS -j oe\n#PBS -o $PBS_O_WORKDIR/%s\n#PBS -A %s\ncd $PBS_O_WORKDIR\n%s' % (job_name, self.processor_line, log_file, self.user, command)
                #to_append = "#!/bin/bash\n#PBS -N %s\n#PBS -l  %s\n#PBS -j oe\n#PBS -o $PBS_JOBDIR/%s\n#PBS -A %s\ncd $PBS_O_WORKDIR\n%s\nsed -e 's/^/[HENIPIPE] %s: /' $PBS_JOBDIR/%s >> %s\n" % (job_name, self.processor_line, log_file, self.user, command, job_name, log_file, self.log_name)
                to_append = "#!/bin/bash\n#PBS -N %s\n#PBS -l %s\n#PBS -j oe\n#PBS -o $PBS_O_WORKDIR/logtmp\n#PBS -A %s\ncd $PBS_O_WORKDIR\n{%s} 2>&1 | tee %s\nsed -e 's/^/[piperna] JOB: %s:\t\t/' %s >> %s\nrm %s\n" % (job_name, self.processor_line, self.user, command, log_file, job_name, log_file, self.log_name, log_file)
            if self.cluster=="SLURM":
                to_append = '#!/bin/bash\n#SBATCH --job-name=%s\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --mem-per-cpu=8000\n%s' % (job_name, command)
            job_string.append(to_append)
        return job_string

    def get_processor_line(self):
        gb = self.threads * 8
        return """select=1:mem=%sGB:ncpus=%s""" % (gb, self.threads)

    def run_job(self):
        for script in self.script:
            if self.cluster=="PBS":
                if self.debug==False:
                    # Open a pipe to the command.
                    proc = Popen('qsub', shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
                    if (sys.version_info > (3, 0)):
                        proc.stdin.write(script.encode('utf-8'))
                        out, err = proc.communicate()
                    else:
                        proc.stdin.write(script)
                        out, err = proc.communicate()
                # Print your job and the system response to the screen as it's submitted
                print(script)
                if self.debug==False:
                    print(out)
                    time.sleep(0.1)
            if self.cluster=="SLURM":
                if self.debug==False:
                    #open("temp.sh", 'w').write(script+'\n')
                    proc = Popen(['sbatch'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
                    out, err = proc.communicate(script)
                print(script)
                if self.debug==False:
                    print(out)
                    time.sleep(0.1)

class Star(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Star, self).__init__(*args, **kwargs)
        self.threads = kwargs.get('threads')
        self.job = "STAR_ALIGN"
        self.out_sam_type = kwargs.get('out_sam_type')
        self.count = kwargs.get('count')
        self.global_add_STAR_string = kwargs.get('global_add_STAR_string')
        self.runmode = self.get_runmode()
        self.processor_line = self.get_processor_line()
        self.command = self.Star_executable()
        self.script = self.generate_job()
    def __call__():
        pass

    def Star_executable(self):
        commandline=""
        command = []
        fastq_line = ""
        localSTARStr = ""
        #print("Runmode is " + self.runmode)
        for sample in self.runsheet_data:
            if self.runmode=="single":
                fastq_line = sample['fastqs'].replace('\t', ',')
            if self.runmode=="pe":
                fastq_line = sample['fastq1'].replace('\t', ',') + " " + sample['fastq2'].replace('\t', ',')
                #print ("Fastq line is '" + fastq_line + "'")
            modules = """\nmodule load STAR/2.5.3a\n"""
            if 'localSTARStr' in sample:
                localSTARStr = sample['localSTARStr']
            commandline = """\nSTAR --genomeDir %s --runThreadN %s --readFilesIn %s --outFileNamePrefix %s --readFilesCommand zcat --outSAMtype %s %s""" % (sample['index'], self.threads, fastq_line, sample['output'], self.out_sam_type, localSTARStr)
            if self.count:
                commandline = commandline + """ --quantMode GeneCounts"""
            commandline = modules + commandline + " " + self.global_add_STAR_string
            #print(commandline.__class__.__name__)
            command.append(commandline+"\n")
        return command

class kallisto(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(kallisto, self).__init__(*args, **kwargs)
        self.runmode = self.get_runmode()
        self.mfl = kwargs.get('mfl')
        self.sfl = kwargs.get('sfl')
        self.processor_line = self.get_processor_line()
        self.command = self.kallisto_executable()
        self.script = self.generate_job()
    def __call__():
        pass

    def kallisto_executable(self):
        commandline=""
        command = []
        for sample in self.runsheet_data:
            modules = """\nmodule load %s\n""" % sample['module']
            if self.runmode=="single":
                sample['single_info'] = """--single -l %s -s %s"""% (self.mfl, self.sfl)
            else:
                sample['single_info'] = ""
            commandline = """\nkallisto quant -i %s -o %s %s %s""" % (sample['index'], sample['output'], sample['single_info'], sample['fastqs'].replace("\t", " "))
            commandline = modules + commandline
            #print(commandline.__class__.__name__)
            command.append(commandline)
        return command

class summarize(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(summarize, self).__init__(*args, **kwargs)
        self.runsheet = kwargs.get('runsheet')
        self.output = kwargs.get('output')
        self.runsheet_data = [{"sample":"all_samples"}]
        self.job = "SUMMARIZE"
        self.processor_line = "select=1:mem=100gb:ncpus=1"
        self.command = self.summarize_executable()
        self.script = self.generate_job()


    def __call__():
        pass

    def summarize_executable(self):
        commandline=""
        command = []
        modules = """\nmodule load R\n"""
        commandline = """echo '\n[SUMMARIZE] Running SUMMARIZE... Output:\n'\nRscript %s -r %s -o %s\n""" % (SUMMARIZE_SCRIPT, self.runsheet, self.output)
        commandline = modules + commandline
        #print(commandline.__class__.__name__)
        command.append(commandline)
        return command

def convert_windows_newlines(file_name):
    """
    Helper function for converting windows newlines in a file as a preprocessing step in case users make samplesheet
    in excel.

    Args:
    file_name (str): name of file to convert to unix newlines.

    """
    # Read in the file in entirety first
    with open(file_name) as input:
        file_text = input.read()

    # Replace newlines produced by excel
    file_text = file_text.replace('\r', '\n')

    # Write the new output back to the file
    with open(file_name, 'w') as output:
        output.write(file_text)


def find_colnames(runsheet, header=True):
    """
    Helper function for getting headers for a runsheet.
    Args:
        runsheet (file): file object to runsheet.  CSV file with header is supported.
    Yields:
        Colnames of CSV
    """
    if header==True:
        with open(runsheet, 'r') as f:
            reader = csv.reader(f)
            return(next(reader))            # read header

def parse_runsheet(runsheet, header=True, colnames=None):
    """
    Helper function for parsing runsheets provided by user.
    Args:
        runsheet (file): file object to runsheet. 
    Yields:
        Each line of the run as a dict by column name.
        Columns: as listed in csv
    """
    if header:
        columns = find_colnames(runsheet, header=header)
    if header==False and colnames is not None:
        columns = colnames
    if header==False and colnames is None:
        raise ValueError('No header indicated in runsheet and no colnames provided')

    entries = []
    with open(runsheet, 'r') as f:
        reader = csv.reader(f)
        if header==True:
            next(reader)            # skip header
        for line in f:
                #print(line)
            entries = [entry.strip() for entry in line.strip().split(',')]
            if len(entries) != len(columns):
                raise ValueError('Sample sheet does not match expected columns. Expects: %s' % ','.join(columns))
            entry_dict = dict(zip(columns, entries))
            yield entry_dict

def find_fastq_mate(dir, sample_flag=None, full_name=True):
    fastqs=[]
    fastq1=[]
    fastq2=[]
    for file in os.listdir(dir):
        if file.endswith(".fastq.gz"):
            fastqs.extend([file])
            if "_R1_" in file:
                fastq1.extend([file])
            if "_R2_" in file:
                fastq2.extend([file])
    fastq1_mate=[]
    for fastq in fastq1:
        #check if present
        put_R2=re.sub('_R1_', '_R2_', fastq)
        try:
            fastq1_mate.extend([fastq2[fastq2.index(put_R2)]])
        except ValueError:
            raise ValueError("Could not find matching file: "+put_R2+" for: "+fastq+" in "+dir)
    if full_name:
            for i in range(len(fastq1)):
                fastq1[i] = os.path.join(dir, fastq1[i])
                fastq1_mate[i] = os.path.join(dir, fastq1_mate[i])
    keys=['directory_long', 'directory_short','fastq1', 'fastq2', 'has_fastq']
    #values=[dir, os.path.basename(dir), "\t".join(os.path.join(dir, fastq1)), "\t".join(os.path.join(dir, fastq1_mate)), len(fastq1)>0]
    values=[dir, os.path.basename(dir), "\t".join(fastq1), "\t".join(fastq1_mate), len(fastq1)>0]
    return(dict(zip(keys, values)))



def find_colnames(runsheet, header=True):
    """
    Helper function for getting headers for a runsheet.
    Args:
        runsheet (file): file object to runsheet.  CSV file with header is supported.
    Yields:
        Colnames of CSV
    """
    if header==True:
        with open(runsheet, 'r') as f:
            reader = csv.reader(f)
            return(next(reader))            # read header


def load_genomes(genomes_file):
    with open(genomes_file, "r") as read_file:
        genome_data = json.load(read_file)
    return genome_data

def make_runsheet(folder, sample_flag, genome_key, typeofseq, output=None, fasta=None, software="STAR"):
    #folder = '/active/furlan_s/Data/CNR/190801_CNRNotch/fastq/mini/fastq'
    #genome_key = "shivani_bulk"
    genome_data = load_genomes(GENOMES_JSON).get(genome_key)
    if output is None:
        output = os.path.join(os.getcwd())
    ddir=[x[0] for x in os.walk(folder)]
    dat=list(map(find_fastq_mate, ddir))
    good_dat = [i for i in dat if i.get('has_fastq') is True]
    good_dat = [i for i in good_dat if re.compile(r'.*'+sample_flag).search(i.get('directory_short'))]

    for i in good_dat:
        i.update({'sample': i.get('directory_short'), \
            'output': os.path.join(output, i.get('directory_short')), \
            'software': software,
            'index': genome_data.get('fasta'),
            'gtf': genome_data.get('gtf')})
    #print(good_dat)
    keys = good_dat[0].keys()
    with open(os.path.join(output, 'runsheet.csv'), 'w') as output_file:
        dict_writer = csv.DictWriter(output_file, fieldnames = keys, extrasaction='ignore')
        dict_writer.writeheader()
        dict_writer.writerows(good_dat)

def check_runsheet_parameter(runsheet, parameter, verbose=False):
    for i in runsheet:
        if i.get(parameter) is None:
            if verbose: print("no data for parameter "+parameter)
            return 0
        if i.get(parameter) is "":
            if vebose: print("header present, but no or incomplete data for "+parameter)
            return 0
    if verbose: print("runsheet_okay_for_parameter_"+parameter)
    return 1

def check_runsheet(args, runsheet, verbose=False):
    required_args = ["sample", "output", "software", "index"]
    has_data = []
    for i in required_args:
        has_data.append(check_runsheet_parameter(runsheet, i, verbose = verbose))
    has_data = dict(zip(required_args, has_data))
    listOfKeys = getKeysByValues(has_data, [0])
    error_string = "Data not found for the following required runsheet columns: "
    if len(listOfKeys) != 0:
        for key in listOfKeys:
            error_string = error_string + key + ', '
        raise ValueError(error_string)
    return

def getKeysByValues(dictOfElements, listOfValues):
    listOfKeys = list()
    listOfItems = dictOfElements.items()
    for item  in listOfItems:
        if item[1] in listOfValues:
            listOfKeys.append(item[0])
    return  listOfKeys 


