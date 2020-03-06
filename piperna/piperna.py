#!/usr/bin/python

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
import pdb
from datetime import datetime

_ROOT = os.path.abspath(os.path.dirname(__file__))
#GENOMES_JSON = os.path.join('/Users/sfurla/Box Sync/PI_FurlanS/computation/develop/piperna/piperna/data/genomes.json')
GENOMES_JSON = os.path.join(_ROOT, 'data', 'genomes.json')
SUMMARIZE_SCRIPT = os.path.join(_ROOT, 'scripts', 'summarize.R')
ENVIRONS_JSON = os.path.join(_ROOT, 'data', 'environs.json')

class SampleFactory:
    def __init__(self, *args, **kwargs):
        self.environs = environs(cluster = kwargs.get('cluster'), user = kwargs.get('user'), log = kwargs.get('log'), additional_header = kwargs.get('additional_header'),  threads = kwargs.get('threads'), gb_ram = kwargs.get('gb_ram'))
        #remove later
        self.user = kwargs.get('user')
        self.cluster = kwargs.get('cluster')
        self.threads = kwargs.get('threads')
        self.gb_ram = kwargs.get('gb_ram')
        self.cluster = kwargs.get('cluster')
        self.runsheet_data = kwargs.get('runsheet_data')
        self.debug = kwargs.get('debug')
        # self.log_name = kwargs.get('log')
    def __call__():
        pass

    def id_generator(self, size=10, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    def get_runmode(self):
        mode = ""
        if "fastq1" in self.runsheet_data[0]:
            if "fastq2" in self.runsheet_data[0]:
                mode = "pe"
        if "fastqs" in self.runsheet_data[0]:
            mode = "single"
        return mode

    def run_job(self):
        popen_command = self.environs.popen_command
        for script in self.bash_scripts:
            if self.debug==False:
                # Open a pipe to the command.
                proc = Popen(popen_command, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, close_fds=True)
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

class environs:
    def __init__(self, *args, **kwargs):
        self.cluster = kwargs.get('cluster')
        self.user = kwargs.get('user')
        self.log = kwargs.get('log')
        self.additional_header = kwargs.get('additional_header')
        if self.additional_header is None:
            self.additional_header = "#"
        self.environs_data = self.load_environs(ENVIRONS_JSON).get(self.cluster)
        self.popen_command = self.environs_data["popen"]
        self.threads = kwargs.get('threads')
        self.ram = kwargs.get('gb_ram')

    def id_generator(self, size=10, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    def load_environs(self, environs_file):
        with open(environs_file, "r") as read_file:
            data = json.load(read_file)
        return data

    def get_processor_line(self, *args, **kwargs):
        if kwargs.get('threads') is None:
            threads = self.environs_data["resources"][kwargs.get('job')]["threads"]
        else:
            threads = kwargs.get('threads')
        if kwargs.get('gb_ram') is None:
            ram = self.environs_data["resources"][kwargs.get('job')]["ram"]
        else:
            ram = kwargs.get('ram')
        return [threads, ram]

    def generate_job(self, commands, job):
        bash_scripts=[]
        torun = len(commands)
        threads,ram = self.get_processor_line(threads = self.threads, ram = self.ram, job = job)
        for i in range(torun):
            # bash_script = {  "PBS" : "#!/bin/bash\n#PBS -N %s\n#PBS -l %s\n#PBS -j oe\n#PBS -o $PBS_O_WORKDIR/logtmp\n#PBS -A %s\ncd $PBS_O_WORKDIR\n{%s} 2>&1 | tee %s\nsed -e 's/^/[HENIPIPE] JOB: %s:\t\t/' %s >> %s\nrm %s\n" % (job_name, self.processor_line, self.user, command, log_file, job_name, log_file, self.log_name, log_file),
                            # "SLURM" : "#!/bin/bash\n#SBATCH --job-name=%s\n#SBATCH --output=tmp\n#SBATCH --error=tmp\n#SBATCH --ntasks=1\n%s\n{%s} 2>&1 | tee %s\nsed -e 's/^/[HENIPIPE] JOB: %s:\t\t/' %s >> %s\nrm %s\n" % (job_name, self.processor_line, command, log_file, job_name, log_file, self.log_name, log_file)}
            bash_script = self.assemble_script(   LOG_FILE = self.log, \
                                    MODULES = self.environs_data["resources"][job]["modules"], \
                                    TEMP_LOG_FILE = self.id_generator(), \
                                    JOB_NAME = (job + "_" + commands[i][0]), \
                                    COMMAND = commands[i][1],
                                    RAM = ram,
                                    HEADER = self.additional_header,
                                    THREADS = threads,
                                    USER = self.user,
                                    TIME = str(datetime.now().strftime("%Y-%m-%d_%H:%M:%S")))
            bash_scripts.append(bash_script)
        return bash_scripts

    def assemble_script(self, *args, **kwargs):
        script_list=[]
        order=[]
        fn_args = kwargs
        for key, value in self.environs_data["script_lines"].items():
            script_list.append(value)
            order.append(key)
        script_list = [x for _,x in sorted(zip(order,script_list))]
        lines_unparsed = [x[0] for x in script_list]
        values_to_insert = [x[1].split("|") for x in script_list]
        lines_parsed = []
        global to_test
        to_test=[lines_unparsed, values_to_insert, fn_args]
        #pdb.set_trace()
        for i in range(len(lines_unparsed)):
            if values_to_insert[i][0] is "":
                lines_parsed.append(lines_unparsed[i])
            else:
                string=lines_unparsed[i]
                for j in range(len(values_to_insert[i])):
                    string = re.sub("<--{0}-->".format(j), str(fn_args.get(values_to_insert[i][j])), string)
                lines_parsed.append(string)
        return "\n".join(lines_parsed)


class Star(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(Star, self).__init__(*args, **kwargs)
        self.job = "PIPERNA_STAR"
        self.out_sam_type = kwargs.get('out_sam_type')
        self.count = kwargs.get('count')
        self.global_add_STAR_string = kwargs.get('global_add_STAR_string')
        if self.global_add_STAR_string is None:
            self.global_add_STAR_string = ""
        else:
            self.global_add_STAR_string = " "+self.global_add_STAR_string
        self.runmode = self.get_runmode()
        self.commands = self.Star_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
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
            if 'localSTARStr' in sample:
                localSTARStr = sample['localSTARStr']
            commandline = """\nSTAR --genomeDir %s --runThreadN %s --readFilesIn %s --outFileNamePrefix %s --readFilesCommand zcat --outSAMtype %s %s""" % (sample['index'], self.threads, fastq_line, sample['output'], self.out_sam_type, localSTARStr)
            if self.count:
                commandline = commandline + """ --quantMode GeneCounts"""
            commandline = commandline + " " + self.global_add_STAR_string
            #print(commandline.__class__.__name__)
            command.append([sample['sample'], commandline])
        return command

class kallisto(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(kallisto, self).__init__(*args, **kwargs)
        self.job = "PIPERNA_KALLISTO"
        self.runmode = self.get_runmode()
        self.mfl = kwargs.get('mfl')
        self.sfl = kwargs.get('sfl')
        self.processor_line = self.get_processor_line()
        self.commands = self.kallisto_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)
    def __call__():
        pass

    def kallisto_executable(self):
        commandline=""
        command = []
        for sample in self.runsheet_data:
            if self.runmode=="single":
                sample['single_info'] = """--single -l %s -s %s"""% (self.mfl, self.sfl)
            else:
                sample['single_info'] = ""
            commandline = """\nkallisto quant -i %s -o %s %s %s""" % (sample['index'], sample['output'], sample['single_info'], sample['fastqs'].replace("\t", " "))
            commandline = modules + commandline
            #print(commandline.__class__.__name__)
            command.append([sample['sample'], commandline])
        return command

class summarize(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(summarize, self).__init__(*args, **kwargs)
        self.runsheet = kwargs.get('runsheet')
        self.output = os.path.join(kwargs.get('output'), "summarizedExperiment.RDS")
        self.runsheet_data = [{"sample":"all_samples"}]
        self.job = "PIPERNA_SUMMARIZE"
        self.threads = kwargs.get('threads')
        self.commands = self.summarize_executable()
        self.bash_scripts = self.environs.generate_job(self.commands, self.job)


    def __call__():
        pass

    def summarize_executable(self):
        commandline=""
        command = []
        commandline = """\necho '\n[PIPERNA_SUMMARIZE] Running SUMMARIZE... Output:\n'\nRscript %s -r %s -o %s -t %s""" % (SUMMARIZE_SCRIPT, self.runsheet, self.output, self.threads)
        #print(commandline.__class__.__name__)
        command.append(["Samples", commandline])
        return command

class concatfastq(SampleFactory, object):
    def __init__(self, *args, **kwargs):
        super(concatfastq, self).__init__(*args, **kwargs)
        self.flowcell_folders = kwargs.get('folders')
        self.threads = 1
        self.job="CONCATFASTQ"
        self.output = kwargs.get('output')
        self.typeofseq = kwargs.get('typeofseq')
        self.processor_line = self.get_processor_line()
        self.runsheet_data = self.prep_concatfastq()
        self.commands = self.concatfastq_executable()
        self.script = self.environs.generate_job(self.commands, self.job)

    def __call__():
        pass

    def prep_concatfastq(self):
        flowcell_folders = self.flowcell_folders.split(",")
        #flowcell_folders = args.flow_cell_folders.split(",")
        if self.output is None:
            output = os.path.join(os.getcwd())
        files_by_handle={}
        typeofseq=self.typeofseq
        output = self.output
        #typeofseq="pe"
        for fc in flowcell_folders:
            #ddir=[x[0] for x in os.walk(fc)]
            fastq1s=[os.path.basename(i) for i in find_fastq_mate(fc).get("fastq1").split("\t")]
            handles = unique_handle(fastq1s)
            for i in handles:
                merged_fn = os.path.join(output, i+"_R1.fastq.gz")
                r = re.compile("^"+i+"*")
                filelist = [os.path.join(fc, j) for j in sorted(list(filter(r.match, fastq1s)))]
                if not all([os.path.exists(i) for i in filelist]):
                    raise ValueError("One of these files was not found: "+", ".join(filelist))
                if merged_fn in files_by_handle:
                    files_by_handle.get(merged_fn).extend(filelist)
                else:
                    files_by_handle.update({merged_fn: filelist})
            if typeofseq =="pe": 
                fastq2s=[os.path.basename(i) for i in find_fastq_mate(fc).get("fastq2").split("\t")]
                handles = unique_handle(fastq2s)
                for i in handles:
                    merged_fn = os.path.join(output, i+"_R2.fastq.gz")
                    r = re.compile("^"+i+"*")
                    filelist = [os.path.join(fc, j) for j in sorted(list(filter(r.match, fastq2s)))]
                    if not all([os.path.exists(i) for i in filelist]):
                        raise ValueError("One of these files was not found: "+", ".join(filelist))
                    if merged_fn in files_by_handle:
                        files_by_handle.get(merged_fn).extend(filelist)
                    else:
                        files_by_handle.update({merged_fn: filelist})
        runsheet_data=[]
        i=1
        for key, value in files_by_handle.items():
            runsheet_data.append({  "sample" : str(i),
                                    "key" : key,
                                    "file_list" : value})
            i = i+1
        return runsheet_data

    def concatfastq_executable(self):
        commandline=""
        command = []
        for item in self.runsheet_data:
            files_to_concat = " ".join(item.get("file_list"))
            commandline = """\ncat {0} > {1}\n""".format(files_to_concat, item.get("key"))
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

def unique_handle(fastqlist):
    #this function takes a list of fastq files and finds the "unique handle" ignoring lane and file number and groups and orders the fastqs based on this handle
    #the goal of this function is to find the set of files that ultimately can concatenated to create one sample (i.e. one unique handle)
    regpat="_L..._R._....fastq.gz"
    handles = sorted(list(set([os.path.basename(i) for i in [re.sub(regpat, "", i) for i in fastqlist]])))
    return handles

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

def find_fastqs_by_sample(dir, strsplit="_R1_", R1_flag="_R1_", R2_flag="_R2_", full_name=True):
    fastqs=[]
    fastq1=fq1=[]
    fastq2=fq2=[]
    sample=[]
    sample_dict=[]
    for file in os.listdir(dir):
        if file.endswith(".fastq.gz"):
            fastqs.extend([file])
            if "_R1_" in file:
                fastq1.extend([file])
            if "_R2_" in file:
                fastq2.extend([file])
    sample=[i.split(strsplit)[0] for i in fastq1]
    for samp in sample:
        R1matches=list(filter(lambda x: samp+R1_flag in x, fastq1))
        R1matches=[os.path.join(dir, i) for i in R1matches]
        fastq1s="\t".join(R1matches)
        R2matches=list(filter(lambda x: samp+R2_flag in x, fastq2))
        R2matches=[os.path.join(dir, i) for i in R2matches]
        fastq2s="\t".join(R2matches)
        sample_dict.extend([{"sample":samp,
            "fastq1":fastq1s,
            "fastq2":fastq2s}])
    return(sample_dict)

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

def make_runsheet(folder, sample_flag, genome_key, typeofseq, organized_by, output=None, fasta=None, software="STAR"):
    #folder = '/active/furlan_s/Data/CNR/190801_CNRNotch/fastq/mini/fastq'
    #genome_key = "shivani_bulk"
    genome_data = load_genomes(GENOMES_JSON).get(genome_key)
    if output is None:
        output = os.path.join(os.getcwd())
    if(organized_by=="folder" and typeofseq=="pe"):
        ddir=[x[0] for x in os.walk(folder)]
        dat=list(map(find_fastq_mate, ddir))
        good_dat = [i for i in dat if i.get('has_fastq') is True]
        good_dat = [i for i in good_dat if re.compile(r'.*'+sample_flag).search(i.get('directory_short'))]
        #good_dat = [i.update({"sample":i.get('directory_short')}) 
        for i in good_dat:
            i.update({"sample":i.get('directory_short')})
    if(organized_by=="file" and typeofseq=="pe"):
        good_dat=find_fastqs_by_sample(folder)
    for i in good_dat:
        i.update({'sample': i.get('sample'), \
            'output': os.path.join(output, i.get('sample')), \
            'software': software,
            'index': genome_data.get('fasta'),
            'gtf': genome_data.get('gtf')})
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


