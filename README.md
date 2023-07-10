# Install Packages
We run the following command in rhino (make sure that `Anaconda` is not installed).
1. **Install modules**
  - module load Python/3.7.4-foss-2019b-fh1
  - module load FastQC/0.11.9-Java-11
  - module load MultiQC/1.9-foss-2019b-Python-3.7.4
  - module load cutadapt/2.9-foss-2019b-Python-3.7.4
  - module load BBMap/38.91-GCC-10.2.0
  - module load STAR/2.7.7a-GCC-10.2.0
  - module load picard/2.25.0-Java-11
  - module load SAMtools/1.11-GCC-10.2.0
  - export PICARD=$EBROOTPICARD/picard.jar
2. Install pip packages
  - pip install RSeQC

    The output should look like: `Successfully installed RSeQC-4.0.0 bx-python-0.8.11 pyBigWig-0.3.18 pysam-0.16.0.1`

# Prepare Directory Environment
1. **Make Directories**
  - Create a directory for the project (`workdir`)
  - Inside `workdir` make a `data/raw_data` folder
  - `controller.py` and `config.yaml` will reside in `workdir`
2. **Download the raw data inside the data/raw_data folder**

# Run Main Pipeline
1. **Adjust `config.yaml`**
  - Usually this will need to be adjusted by species and first a `screen` run is run to discover adapters
  - Then a full or wrapup run is done with the known adapters and read length
2. **Run `controller.py`**
  - Ideally open a `tmux` window or use `nohup` because this command generally takes a while to run. Final QC reports can all be found in the `workdir/qc_reports` folder. The full summary is in `workdir/qc_reports/summary/multiqc.html`.
  - I prefer to run this with `nohup python controller.py > Log.log &` in order to save a log and run async
