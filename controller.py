import os, subprocess, yaml
from glob import glob

'''
meant to be run in same directory as readme.md
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
    try:
        print(f'{prelude}...')
        subprocess.check_output(command, shell=True)
        print('...FINISHED')
    except subprocess.CalledProcessError as e:
        print('ERROR --->>>')
        print(e.output.decode())
        print('<<<--- ERROR')


# cleans up all directories and creates them
def prep_folders():
    '''
    remove relevant directories and re-create them for clean environment
    '''
    # remove stuff
    dirty_folders = ['data/trimmed_data/*', 'qc_reports/*']  # TODO add on to this
    print('REMOVING DIRTY FILES...')
    for folder in dirty_folders:
        if(os.path.exists(folder[:-1])):
            run_command(f'rm -r {folder}', f'REMOVING DIRTY {folder}')
    print('...FINISHED')
    # create stuff
    create_folders = ['data/trimmed_data', 'data/reference', 'data/reference_STAR', 'qc_reports', 'qc_reports/fastqc_raw', 'qc_reports/summary_raw', 'qc_reports/adapter_detection/', 'qc_reports/trimming_reads', 'qc_reports/fastqc_trimmed', 'qc_reports/summary_trimmed', 'qc_reports/mapping_reads', 'qc_reports/summary_mapped']
    for folder in create_folders:
        if(not os.path.exists(folder)):
            os.mkdir(folder)


# qc raw data
def qc_raw_data():
    '''
    perform fastqc then multiqc on raw fastq.gz files
    '''
    run_command('fastqc data/raw_data/*.fastq.gz -o qc_reports/fastqc_raw', 'FASTQC ON RAW DATA')
    run_command('multiqc qc_reports/fastqc_raw -o qc_reports/summary_raw', 'SUMMARIZED QC ON RAW DATA')


# detect adapters in raw data
def detect_adapters():
    '''
    uses bbmerge.sh to detect adapters
    '''
    fastqs = glob('data/raw_data/*read1.fastq.gz')  # TODO change this to correct one
    for fastq in fastqs:
        # get data
        raw_fastq1 = fastq
        raw_fastq2 = raw_fastq1.replace('read1', 'read2')
        adapters = 'qc_reports/adapter_detection/' + raw_fastq1.split('/')[-1].replace('read1.fastq.gz', 'adapters.fa')
        # run command
        run_command(f'bbmerge.sh -Xmx4g in1={raw_fastq1} in2={raw_fastq2} outa={adapters}')


# trim data
def trim_raw_data(configs):
    '''
    filters min length >= 20 and quality >= 20
    '''
    # get adapters
    forward = configs['read1']
    reverse = configs['read2']
    fastqs = glob('data/raw_data/*read1.fastq.gz')  # TODO change this to correct one
    for fastq in fastqs:
        # get data
        raw_fastq1 = fastq
        trimmed_fastq1 = raw_fastq1.replace('raw_data/', 'trimmed_data/')
        raw_fastq2 = raw_fastq1.replace('read1', 'read2')  # TODO change this to correct one
        trimmed_fastq2 = raw_fastq2.replace('raw_data/', 'trimmed_data/')
        # run commands
        run_command(f'cutadapt -a {forward} -A {reverse} -o {trimmed_fastq1} -p {trimmed_fastq2} {raw_fastq1} {raw_fastq2} -m 20 -q 20 --info-file qc_reports/trimming_reads/cutadapt_outputs.txt', 'TRIMMING RAW DATA')

# qc raw data
def qc_trimmed_data():
    '''
    perform fastqc then multiqc on trimmed fastq.gz files
    '''
    run_command('fastqc data/trimmed_data/*.fastq.gz -o qc_reports/fastqc_trimmed', 'FASTQC ON TRIMMED DATA')
    run_command('multiqc qc_reports/fastqc_trimmed -o qc_reports/summary_trimmed', 'SUMMARIZED QC ON TRIMMED DATA')


# map trimmed data
def map_trimmed_data(configs):
    '''
    maps data using star and refGene genome from ucsc mm10
    '''
    # download reference
    gtf_link = configs['gtf']
    fasta_link = configs['fasta']
    run_command(f'wget {gtf_link} -O data/reference/mm10.refGene.gtf.gz', 'DOWNLOADING REFERENCE GTF')
    run_command(f'wget {fasta_link} -O data/reference/chromFa.tar.gz', 'DOWNLOADING REFERENCE FASTA')
    run_command(f'tar -xzvf data/reference/chromFa.tar.gz -C data/reference/', 'UNPACKING REFERENCE FASTA')
    run_command(f'rm data/reference/*random.fa data/reference/chrUn*.fa', 'CLEANING REFERENCE FASTA')
    # index genome
    fasta_files = ' '.join(glob('data/reference/*.fa'))
    run_command(f'STAR --runThreadN 10 --runMode genomeGenerate --genomeDir data/reference_STAR --genomeFastaFiles {fasta_files} -sjdbGTFfile data/reference/mm10.refGene.gtf.gz --sjdbOverhand 36', 'INDEXING GENOME')  # TODO CHANGE TO READ LENGTH - 1
    # map reads
    read1_files = ','.join(sorted(glob('data/trimmed_reads/*read1.fastq.gz')))
    read2_files = read1_files.replace('read', 'read2')
    run_command(f'STAR --runThreadN 10 --genomeDir data/reference_STAR/mapped_ --readFilesIn {read1_files} {read2_files} --outSAMtype BAM SortedByCoordinate --outFileNamePrefix data/mapped_data/ --readFilesCommand gunzip -c --outFilterIntronMotifs', 'MAPPING READS')


# qc raw data
def qc_mapped_data():
    '''
    perform multiqc on mapped outputs
    '''
    run_command('mv data/reference_STAR/*.out qc_reports/mapping_reads/', 'MOVE STAR MAPPED LOGS')
    run_command('multiqc qc_reports/mapping_reads -o qc_reports/summary_mapped', 'SUMMARIZED QC ON MAPPED DATA')


def main():
    configs = load_configs()
    prep_folders()  # prep environment
    qc_raw_data()  # qc on original data
    detect_adapters()  # try to find adapters
    if(configs['analysis']['full']=="true"):
        trim_raw_data(configs['trimming'])  # trim reads
        qc_trimmed_data() # qc on trimmed data
        map_trimmed_data(configs['mapping'])  # map reads
        qc_mapped_data()  # qc on mapped data

main()  # run program
























# comment

