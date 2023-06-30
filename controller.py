import os
import subprocess
import yaml
import time
from glob import glob
import pandas as pd
import numpy as np

'''
meant to be run in same directory as README.md
'''

# helper functions
def load_configs():
    '''
    reads in the configurations
    '''
    configs = yaml.safe_load(open('config.yaml', 'r'))
    return configs


def run_command(command, prelude):
    '''
    runs a command and outputs prelude message
    '''
    print(f'{prelude}...')
    proc = subprocess.Popen(command, shell=True)
    return proc


# creates directories if they do not already exist
def prep_folders():
    '''
    create folders and prepares them to be run
    '''
    # create folders
    create_folders = ['data/trimmed_data', 'data/reference', 'data/reference_STAR', 'data/mapped_data', 'qc_reports', 'qc_reports/fastqc_raw', 'qc_reports/summary_raw', 'qc_reports/adapter_detection/', 'qc_reports/cutadapt_trimmed', 'qc_reports/fastqc_trimmed', 'qc_reports/summary_trimmed', 'qc_reports/mapping_reads', 'qc_reports/summary_mapped', 'qc_reports/summary_counted', 'qc_reports/summary', 'outputs/']
    for folder in create_folders:
        if(not os.path.exists(folder)):
            os.mkdir(folder)


# qc raw data
def qc_raw_data(configs):
    '''
    perform fastqc then multiqc on raw fastq.gz files
    '''
    # run the qc co-routines in parallel
    procs = []
    for fn in glob('data/raw_data/*' + configs['fastq_suffix']):
        proc = run_command(f'fastqc {fn} -o qc_reports/fastqc_raw', 'FASTQC ON RAW DATA')
        procs.append(proc)
    # wait for them to all finish then run multiqc
    for proc in procs:
        proc.wait()
    run_command('multiqc qc_reports/fastqc_raw -o qc_reports/summary_raw', 'SUMMARIZING QC ON RAW DATA')

# detect adapters in raw data
def detect_adapters(configs):
    '''
    uses bbmerge.sh to detect adapters
    '''
    # define the read1 suffix
    read1_suffix = configs['read1_name'] + configs['fastq_suffix']
    
    # run the adapter finding
    procs = []
    for fastq in glob('data/raw_data/*' + read1_suffix):
        # define the fastq names
        raw_fastq1 = fastq
        raw_fastq2 = raw_fastq1.replace(configs['read1_name'], configs['read2_name'])
        # define the adapter file name
        adapters = 'qc_reports/adapter_detection/' + raw_fastq1.split('/')[-1].replace(read1_suffix, 'adapters.fa')
        stats = adapters.replace('adapters.fa', 'stats.txt')
        # run command
        proc = run_command(f'bbmerge.sh -Xmx4g in1={raw_fastq1} in2={raw_fastq2} outa={adapters}', 'FINDING ADAPTERS')
        procs.append(proc)
    for proc in procs:
        proc.wait()
    
    # run the adapter quantification
    procs = []
    for fastq in glob('data/raw_data/*' + read1_suffix):
        # define the fastq names
        raw_fastq1 = fastq
        raw_fastq2 = raw_fastq1.replace(configs['read1_name'], configs['read2_name'])
        # define the adapter file name
        adapters = 'qc_reports/adapter_detection/' + raw_fastq1.split('/')[-1].replace(read1_suffix, 'adapters.fa')
        stats = adapters.replace('adapters.fa', 'stats.txt')
        # run command
        proc = run_command(f'bbduk.sh -Xmx4g in1={raw_fastq1} in2={raw_fastq2} stats={stats} ref=data/raw_data/adapters.fa', 'QUANTIFYING KNOWN ADAPTERS')
        procs.append(proc)
    for proc in procs:
        proc.wait()


# trim data
def trim_raw_data(configs):
    '''
    filters min length >= 20 and quality >= 20
    '''
    # get adapter sequences
    forward = configs['read1']
    reverse = configs['read2']
    # get fastq sequences
    read1_suffix = configs['read1_name'] + configs['fastq_suffix']
    fastqs = glob('data/raw_data/*' + read1_suffix)
    # trim the fastqs
    procs = []
    for fastq in fastqs:
        # define the fastq file names
        raw_fastq1 = fastq
        trimmed_fastq1 = raw_fastq1.replace('raw_data/', 'trimmed_data/').replace(configs['read1_name'], configs['read1_name'] + 'TRIMMED')
        raw_fastq2 = raw_fastq1.replace(configs['read1_name'], configs['read2_name'])
        trimmed_fastq2 = raw_fastq2.replace('raw_data/', 'trimmed_data/').replace(configs['read2_name'], configs['read2_name'] + 'TRIMMED')
        log = raw_fastq1.replace('data/raw_data/', 'qc_reports/cutadapt_trimmed/').replace(configs['fastq_suffix'], '.cutadapt.out')
        # run commands, trim if needed
        if(configs['trim']=="true"):
            proc = run_command(f'cutadapt -a {forward} -A {reverse} -m 20 -q 20 -o {trimmed_fastq1} -p {trimmed_fastq2} {raw_fastq1} {raw_fastq2} > {log}', 'TRIMMING RAW DATA')
            procs.append(proc)
        else:
            proc = run_command(f'mv {raw_fastq1} {trimmed_fastq1}', 'MOVING READ1')
            procs.append(proc)
            proc = run_command(f'mv {raw_fastq2} {trimmed_fastq2}', 'MOVING READ2')
            procs.append(proc)
    for proc in procs:
        proc.wait()

# qc raw data
def qc_trimmed_data(configs):
    '''
    perform fastqc then multiqc on trimmed fastq.gz files
    '''
    # run the trimmed fastq qc co-routines in parallel
    procs = []
    for fn in glob('data/trimmed_data/*' + configs['fastq_suffix']):
        proc = run_command(f'fastqc {fn} -o qc_reports/fastqc_trimmed', 'FASTQC ON TRIMMED DATA')
        procs.append(proc)
    # wait for them to all finish then run multiqc
    for proc in procs:
        proc.wait()
    run_command('multiqc qc_reports/fastqc_trimmed -o qc_reports/summary_trimmed', 'SUMMARIZING QC ON TRIMMED DATA')

# prep for mapping
def prep_map_reference(configs):
    # get the link names
    gtf_link = configs['gtf']
    fasta_link = configs['fasta']
    genes_link = configs['genes']
    overhang = configs['read_length'] - 1
    # get the file names
    gtf_fn = gtf_link.split('/')[-1]
    fasta_fn = fasta_link.split('/')[-1]
    genes_fn = genes_link.split('/')[-2]
    # get the number of cores
    n_cores = configs['n_cores']
    # downloading the links and unpack them
    run_command(f'wget {gtf_link} -O data/reference/{gtf_fn}', 'DOWNLOADING REFERENCE GTF')
    run_command(f'wget {genes_link} -O data/reference/{genes_fn}', 'DOWNLOADING REFERENCE GENES')
    run_command(f'wget {fasta_link} -O data/reference/{fasta_fn}', 'DOWNLOADING REFERENCE FASTA')
    # wait for them to all finish downloading
    downloading = True
    while downloading:
        downloaded = True
        if not os.path.exists(f'data/reference/{fasta_fn}'):
            downloaded = False
        if not os.path.exists(f'data/reference/{gtf_fn}'):
            downloaded = False
        if not os.path.exists(f'data/reference/{genes_fn}'):
            downloaded = False
        if downloaded:
            downloading = False
        time.sleep(60)
    run_command(f'gzip -d data/reference/{gtf_fn}', 'UNPACKING REFERENCE GTF')
    run_command(f'tar -xzvf data/reference/{fasta_fn} -C data/reference/', 'UNPACKING REFERENCE FASTA')
    run_command(f'gzip -d data/reference/{genes_fn}', 'UNPACKING REFERENCE GENES')
    # index genome
    fasta_files = ' '.join(glob('data/reference/*.fa'))
    run_command(f'STAR --runThreadN {n_cores} --runMode genomeGenerate --genomeDir data/reference_STAR --genomeFastaFiles {fasta_files} --sjdbGTFfile data/reference/{gtf_fn[:-2]} --sjdbOverhang {overhang}', 'INDEXING GENOME')

# map trimmed data
def map_trimmed_data(configs_trimming, configs_mapping):
    '''
    maps data using star and refGene genome from ucsc hg38
    '''
    # get the number of cores
    n_cores = configs_mapping['n_cores']
    # retrieve the trimmed file names
    read1_suffix = configs_trimming['read1_name'] + 'TRIMMED' + configs_trimming['fastq_suffix']
    read1_files = sorted(glob('data/trimmed_data/*' + read1_suffix))
    # map each paired read
    procs = []
    for read1_file in read1_files:
        # get read2's file names and compute the mapping filename
        read2_file = read1_file.replace(configs_trimming['read1_name'], configs_trimming['read2_name'])
        prefix = read1_file.replace('trimmed_data', 'mapped_data').replace(read1_suffix,'')
        # run the alignment
        proc = run_command(f'STAR --runThreadN {n_cores} --genomeDir data/reference_STAR/ --readFilesIn {read1_file} {read2_file} --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN {n_cores} --outFileNamePrefix {prefix} --readFilesCommand gunzip -c --quantMode GeneCounts', 'MAPPING READS')
        procs.append(proc)
    for proc in procs:
        proc.wait()

# qc raw data
def qc_mapped_data(configs):
    '''
    perform multiqc on mapped outputs
    '''
    # retrieve filenames
    genes_fn = configs['genes'].split('/')[-2][:-3]
    # move the logs to the correct locations
    run_command('mv data/mapped_data/*.out qc_reports/mapping_reads/', 'MOVING STAR MAPPED LOGS')
    # loop through each bam file
    bam_files = glob('data/mapped_data/*_Aligned.sortedByCoord.out.bam')
    # mark the duplicates
    procs = []
    for in_bam in bam_files:
        out_bam = in_bam.replace('.out.bam', '.dupMarker.out.bam')
        metrics = in_bam.replace('data/mapped_data', 'qc_reports/mapping_reads').replace('Aligned.sortedByCoord.out.bam', 'marked_dup_metrics.txt')
        proc = run_command(f'java -Xmx20g -jar $PICARD MarkDuplicates -I {in_bam} -O {out_bam} -M {metrics}', 'MARKING DUPLICATES')
        procs.append(proc)
    for proc in procs:
        proc.wait()
    # this is why we use normal https://www.biostars.org/p/14283/
    # index the bam
    procs = []
    for in_bam in bam_files:
        proc = run_command(f'samtools index {in_bam}', 'INDEXING BAM')
        procs.append(proc)
    for proc in procs:
        proc.wait()
    # count the chromosomes
    procs = []
    for in_bam in bam_files:
        stats = in_bam.replace('data/mapped_data', 'qc_reports/mapping_reads').replace('Aligned.sortedByCoord.out.bam', 'chr_stats.txt')
        proc = run_command(f'samtools idxstats {in_bam} > {stats}', 'COUNTING CHROMOSOMES')
        procs.append(proc)
    for proc in procs:
        proc.wait()
    # infer the strand
    procs = []
    for in_bam in bam_files:
        strands = in_bam.replace('data/mapped_data', 'qc_reports/mapping_reads').replace('Aligned.sortedByCoord.out.bam', 'strand_inference.txt')
        proc = run_command(f'infer_experiment.py -r data/reference/{genes_fn} -i {in_bam} > {strands}', 'INFERRING STRAND')
        procs.append(proc)
    for proc in procs:
        proc.wait()
    # get the read distribution
    procs = []
    for in_bam in bam_files:
        distribution = in_bam.replace('data/mapped_data', 'qc_reports/mapping_reads').replace('Aligned.sortedByCoord.out.bam', 'read_distribution.txt')
        proc = run_command(f'read_distribution.py -i {in_bam} -r data/reference/{genes_fn} > {distribution}', 'GATHERING READ DISTRIBUTION')
        procs.append(proc)
    for proc in procs:
        proc.wait()
    # randomly gather genes
    lines = []
    with open(f'data/reference/{genes_fn}', 'rt') as f:
        for line in f:
            lines.append(line)
    np.random.seed(0)
    lines = np.random.choice(lines, size=2000, replace=False)
    # write those random genes down
    genes_fn = genes_fn.replace('.bed','_2k_genes.bed')
    with open('data/reference/hg38_RefSeq_2k_genes.bed', 'wt') as f:
        for line in lines:
            f.writelines(line)
    # get gene bed coverage
    bam_files = ','.join(bam_files)
    run_command(f'geneBody_coverage.py -i {bam_files} -r data/reference/{genes_fn} -o qc_reports/mapping_reads/rseqc', 'GATHERING GENE BODY COVERAGE')
    run_command('multiqc qc_reports/mapping_reads -o qc_reports/summary_mapped', 'SUMMARIZING QC ON MAPPED DATA')


def generate_count_matrix(configs):
    '''
    generates the count matrix per sample given the type
    '''
    count_files = glob('data/mapped_data/*_ReadsPerGene.out.tab')
    dfs = []
    for count_file in count_files:
        df = pd.read_table(count_file, header=None)
        # sense = read1-stranded bc read1+=transcript+=sense
        # and then read2+=transcript-=sense so all makes sense
        df.columns = ['GeneID', 'Unstranded', 'Sense-Stranded', 'Antisense-Stranded']
        # subset off the top four summary rows
        df = df.iloc[4:,][['GeneID', configs['type']]].set_index('GeneID')
        df.columns = [count_file.split('/')[-1].split('_ReadsPerGene.out.tab')[0]]
        dfs.append(df)
        del df
    df = pd.concat(dfs, axis=1).fillna(0)  # combine samples
    df.to_csv('outputs/raw_counts.csv')  # write it


def clean_up():
    '''
    remove irrelevant files and generates final qc summary
    '''
    dup_bams = ' '.join(glob('data/mapped_data/*_Aligned.sortedByCoord.dupMarker.out.bam'))
    splice_tabs = ' '.join(glob('data/mapped_data/*_SJ.out.tab'))
    run_command(f'rm {dup_bams} {splice_tabs}', 'REMOVING DUP AND SPLICE VALUES')
    run_command('multiqc qc_reports/fastqc_raw qc_reports/fastqc_trimmed qc_reports/cutadapt_trimmed/ qc_reports/mapping_reads/ -o qc_reports/summary', 'SUMMARIZING QC FOR ALL DATA')


def main():
    # prep
    configs = load_configs()
    # - screen --> just find adapters and qc the raw data
    # - full   --> run the entire thing end to end
    # - wrapup --> after defining the adapters
    run_type = configs['analysis']['type']
    # run
    # create folders
    prep_folders()
    if (run_type in ['full','screen']):
        # qc the raw sequences
        qc_raw_data(configs['trimming'])
        # retrieve the adapters
        detect_adapters(configs['trimming'])
    # if we just wanted to screen stop here
    if (run_type == 'screen'):
        return
    # otherwise keep going into trimming
    # trim the raw reads
    if(configs['trimming']['trim'] == 'true'):
        trim_raw_data(configs['trimming'])
        # qc the trimmed reads
        qc_trimmed_data(configs['trimming'])
    # if the genome is not indexed then indexe it
    if(configs['mapping']['index_genome'] == 'true'):
        prep_map_reference(configs['mapping'])
    # map the trimmed reads
    if(configs['mapping']['map'] == 'true'):
        map_trimmed_data(configs['trimming'], configs['mapping'])
    # qc the mapped data
    qc_mapped_data(configs['mapping'])
    # count up the reads
    generate_count_matrix(configs['counting'])
    # clean up the environment
    clean_up()

# run the program
main()
