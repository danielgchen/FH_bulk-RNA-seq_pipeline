import os, subprocess, yaml
from glob import glob
import pandas as pd
import numpy as np

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
        print(command)
        print(e.output.decode())
        print('<<<--- ERROR')


# cleans up all directories and creates them
def prep_folders(configs):
    '''
    remove relevant directories and re-create them for clean environment
    '''
    # remove stuff
    dirty_folders = ['data/trimmed_data', 'data/reference', 'data/reference_STAR', 'data/mapped_data', 'qc_reports/', 'outputs/']  # TODO add on to this
    if(configs['index_genome']!="true"):
        dirty_folders.remove('data/reference')
        dirty_folders.remove('data/reference_STAR')
    print('REMOVING DIRTY FILES...')
    for folder in dirty_folders:
        if(os.path.exists(folder)):
            run_command(f'rm -r {folder}', f'REMOVING DIRTY {folder}')
    print('...FINISHED')
    # create stuff
    create_folders = ['data/trimmed_data', 'data/reference', 'data/reference_STAR', 'data/mapped_data', 'qc_reports', 'qc_reports/fastqc_raw', 'qc_reports/summary_raw', 'qc_reports/adapter_detection/', 'qc_reports/cutadapt_trimmed', 'qc_reports/fastqc_trimmed', 'qc_reports/summary_trimmed', 'qc_reports/mapping_reads', 'qc_reports/summary_mapped', 'qc_reports/summary_counted', 'qc_reports/summary', 'outputs/']
    for folder in create_folders:
        if(not os.path.exists(folder)):
            os.mkdir(folder)


# qc raw data
def qc_raw_data():
    '''
    perform fastqc then multiqc on raw fastq.gz files
    '''
    run_command('fastqc data/raw_data/*.fastq.gz -o qc_reports/fastqc_raw', 'FASTQC ON RAW DATA')
    run_command('multiqc qc_reports/fastqc_raw -o qc_reports/summary_raw', 'SUMMARIZING QC ON RAW DATA')


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
        stats = adapters.replace('adapters.fa', 'stats.txt')
        # run command
        run_command(f'bbmerge.sh -Xmx4g in1={raw_fastq1} in2={raw_fastq2} outa={adapters}', 'FINDING ADAPTERS')
        # run command
        run_command(f'bbduk.sh -Xmx4g in1={raw_fastq1} in2={raw_fastq2} stats={stats} ref=data/raw_data/adapters.fa', 'QUANTIFYING KNOWN ADAPTERS')


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
        trimmed_fastq1 = raw_fastq1.replace('raw_data/', 'trimmed_data/').replace('read1', 'trimmed_read1')
        raw_fastq2 = raw_fastq1.replace('read1', 'read2')
        trimmed_fastq2 = raw_fastq2.replace('raw_data/', 'trimmed_data/').replace('read2', 'trimmed_read2')
        log = raw_fastq1.replace('data/raw_data/', 'qc_reports/cutadapt_trimmed').replace('.fastq.gz', '.cutadapt.out')
        # run commands
        run_command(f'cutadapt -a {forward} -A {reverse} -o {trimmed_fastq1} -p {trimmed_fastq2} {raw_fastq1} {raw_fastq2} -m 20 -q 20 > {log}', 'TRIMMING RAW DATA')


# qc raw data
def qc_trimmed_data():
    '''
    perform fastqc then multiqc on trimmed fastq.gz files
    '''
    run_command('fastqc data/trimmed_data/*.fastq.gz -o qc_reports/fastqc_trimmed', 'FASTQC ON TRIMMED DATA')
    run_command('multiqc qc_reports/fastqc_trimmed qc_reports/cutadapt_trimmed/ -o qc_reports/summary_trimmed', 'SUMMARIZING QC ON TRIMMED DATA')


# prep for mapping
def prep_map_reference(configs):
    # download reference
    gtf_link = configs['gtf']
    fasta_link = configs['fasta']
    genes_link = configs['genes']
    overhang = configs['read_length'] - 1
    run_command(f'wget {gtf_link} -O data/reference/mm10.refGene.gtf.gz', 'DOWNLOADING REFERENCE GTF')
    run_command(f'gzip -d data/reference/mm10.refGene.gtf.gz', 'UNPACKING REFERENCE GTF')
    run_command(f'wget {fasta_link} -O data/reference/chromFa.tar.gz', 'DOWNLOADING REFERENCE FASTA')
    run_command(f'tar -xzvf data/reference/chromFa.tar.gz -C data/reference/', 'UNPACKING REFERENCE FASTA')
    run_command(f'wget {genes_link} -O data/reference/mm10_RefSeq.bed.gz', 'DOWNLOADING REFERENCE GENES')
    run_command(f'gzip -d data/reference/mm10_RefSeq.bed.gz', 'UNPACKING REFERENCE GENES')
    # index genome
    fasta_files = ' '.join(glob('data/reference/*.fa'))
    run_command(f'STAR --runThreadN 10 --runMode genomeGenerate --genomeDir data/reference_STAR --genomeFastaFiles {fasta_files} --sjdbGTFfile data/reference/mm10.refGene.gtf --sjdbOverhang {overhang}', 'INDEXING GENOME')


# map trimmed data
def map_trimmed_data():
    '''
    maps data using star and refGene genome from ucsc mm10
    '''
    # map reads - assuming that each sample has one pair
    read1_files = sorted(glob('data/trimmed_data/*read1.fastq.gz'))
    for read1_file in read1_files:
        read2_file = read1_file.replace('read1', 'read2')
        prefix = read1_file.replace('trimmed_data', 'mapped_data').replace('read1.fastq.gz','')
        run_command(f'STAR --runThreadN 10 --genomeDir data/reference_STAR/ --readFilesIn {read1_file} {read2_file} --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 10 --outFileNamePrefix {prefix} --readFilesCommand gunzip -c --quantMode GeneCounts', 'MAPPING READS')


# qc raw data
def qc_mapped_data():
    '''
    perform multiqc on mapped outputs
    '''
    run_command('mv data/mapped_data/*.out qc_reports/mapping_reads/', 'MOVING STAR MAPPED LOGS')
    bam_files = glob('data/mapped_data/*_Aligned.sortedByCoord.out.bam')
    for in_bam in bam_files:
        out_bam = in_bam.replace('.out.bam', '.dupMarker.out.bam')
        metrics = in_bam.replace('data/mapped_data', 'qc_reports/mapping_reads').replace('Aligned.sortedByCoord.out.bam', 'marked_dup_metrics.txt')
        stats = metrics.replace('marked_dup_metrics.txt', 'samtools_idxstats.txt')
        strands = stats.replace('samtools_idxstats.txt', 'rseqc_strand_specificity.txt')
        distribution = strands.replace('strand_specificity.txt', 'distribution.txt')
        run_command(f'java -Xmx10g -jar $PICARD MarkDuplicates -I {in_bam} -O {out_bam} -M {metrics}', 'MARKING DUPLICATES')
        run_command(f'samtools index {in_bam}', 'INDEXING BAM')
        run_command(f'samtools idxstats {in_bam} > {stats}', 'COUNTING CHROMOSOMES')
        run_command(f'infer_experiment.py -r data/reference/mm10_RefSeq.bed -i {in_bam} > {strands}', 'INFERRING STRAND')
        run_command(f'read_distribution.py -i {in_bam} -r data/reference/mm10_RefSeq.bed > {distribution}', 'GATHERING READ DISTRIBUTION')
    # randomly gather genes
    lines = []
    with open('data/reference/mm10_RefSeq.bed', 'rt') as f:
        for line in f:
            lines.append(line)
    lines = np.random.choice(lines, size=2000, replace=False)
    with open('data/reference/mm10_RefSeq_2k_genes.bed', 'wt') as f:
        for line in lines:
            f.writelines(line)
    # get gene bed coverage
    bam_files = ','.join(bam_files)
    run_command(f'geneBody_coverage.py -i {bam_files} -r data/reference/mm10_RefSeq_2k_genes.bed -o qc_reports/mapping_reads/rseqc', 'GATHERING GENE BODY COVERAGE')
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
    # - screen  --> means just find adapters
    # - partial --> means map it but stop there
    # - full    --> means finish end-to-end
    # - wrapup    --> means continue from mapping
    run_type = configs['analysis']['type']
    # run
    if(run_type!='wrapup'):
        prep_folders(configs['mapping'])  # prep environment
        qc_raw_data()  # qc on original data
        detect_adapters()  # try to find adapters
        if(run_type!='screen'):
            trim_raw_data(configs['trimming'])  # trim reads
            qc_trimmed_data() # qc on trimmed data
            if(configs['mapping']['index_genome']=='true'):
                prep_map_reference(configs['mapping'])  # prep environment
            map_trimmed_data()  # map reads
            qc_mapped_data()  # qc on mapped data
    if(run_type in ['wrapup', 'full']):
        generate_count_matrix(configs['counting'])  # get raw counts
        clean_up()  # remove unnecessary bams


main()  # run program
























# comment
