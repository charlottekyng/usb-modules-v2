# Basic project set up

Assuming the new project is called PROJ.
```
PROJ_DIR=PROJ
mkdir $PROJ_DIR
cd $PROJ_DIR
```

From here on, we will assume that you are in this `$PROJ_DIR`.

Clone the code base
```
git clone https://github.com/charlottekyng/usb-modules-v2.git
```

### Setting up sample sheets

1. Sample sheet (required, default: `samples.txt`), consisting of a single column of sample names. 
    ```
    head samples.txt
    SAMPLE1T
    SAMPLE1N
    SAMPLE2T
    SAMPLE2N
    SAMPLE3T1
    SAMPLE3T2
    SAMPLE3N
    ```
1. Sample sets (required for `ANALYSIS_TYPE = SOMATIC`, default: `sample_sets.txt`, consisting of samples within sample sheet above, each row contains all samples from a given individual, with the germline sample last, tab or space delimited
    ```
    >head sample_sets.txt
    SAMPLE1T SAMPLE1N
    SAMPLE2T SAMPLE2N
    SAMPLE3T1 SAMPLE3T2 SAMPLE3N
    ```
1. Sample splits (required for one of these two scenarios, depending on your data, default: `samples.split.txt`)
    1. multiple sets of fastqs per sample, in which case the first column is the <SAMPLE_NAME>, the second column is <SAMPLE_NAME>\_<RUN_NAME>. See below on setting up the `unprocessed_fastq` directory; or 
        ```
        >head samples.split.txt
        SAMPLE1N SAMPLE1N_RUN1
        SAMPLE1N SAMPLE1N_RUN2
        SAMPLE1N SAMPLE1N_RUN3
        SAMPLE2N SAMPLE2N_RUN1
        SAMPLE3N SAMPLE3N_RUN1
        SAMPLE3N SAMPLE3N_RUN2
        ```	
    1. multiple bam files per samples that need to be merged. In which case, the first column is the <SAMPLE_NAME>, the second column is the <SAMPLE_NAME> of one of the bam files that need to be merged. See below on setting up the `unprocessed_bam` directory
        ```
        >head samples.split.txt
        SAMPLE1T SAMPLE1A
        SAMPLE1T SAMPLE1B
        SAMPLE2T SAMPLE2A
        SAMPLE2T SAMPLE2B
        ```

**NOTE**: sample names must be [A-Za-z0-9] (may include '-' but no other symbols allowed) that starts with an alphabet.

### Setting up data directories

There are several options in terms of data files:
1. If you start from FASTQs and you have a single or a single pair of fastqs per sample, then you put your files as `fastq/<SAMPLE_NAME>.1.fastq.gz` (and `fastq/<SAMPLE_NAME>.2.fastq.gz`). Then you are ready to run alignment.
    ```
    >ls fastq
    SAMPLE1T.1.fastq.gz SAMPLE1T.2.fastq.gz SAMPLE2T.1.fastq.gz SAMPLE2T.2.fastq.gz (...)
    ```
1. If you start from FASTQs and you have more than a single or more than a single pair of fastqs per sample, then you put your files as `unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.1.fastq.gz` (and `unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.2.fastq.gz`). Then you are ready to run alignment.
    ```
    >ls unprocessed_fasta
    SAMPLE1N_RUN1.1.fastq.gz SAMPLE1N_RUN1.2.fastq.gz SAMPLE1N_RUN2.1.fastq.gz SAMPLE1N_RUN2.2.fastq.gz (...)
    ```
1. If you start from BAMs (one bam per sample), you should put all your bams as `unprocessed_bam/<SAMPLE_NAME>.bam`. Then you do `make fix_rg` then you will have analysis-ready BAMs.
    ```
    >ls unprocessed_bam
    SAMPLE1A.bam SAMPLE1B.bam SAMPLE2A.bam SAMPLE2B.bam
    ```

### Setting up analysis parameters

You will need a project-level Makefile (`$PROJ_DIR/Makefile`). Note that this is different from the module-level Makefile (`${PROJ_DIR}/usb-modules-v2/Makefile`).
In its most basic form, it only needs one line
```
include usb-modules-v2/Makefile
```
Any additional parameters should go _before_ the line above.


This analysis pipeline is designed to be highly configurable. This also means that there are many possible combinations of parameters. 
The project-level Makefile is where user-configurable parameters are specified. 

You can specify as many parameters as required in your project-level `Makefile`, before the `include usb-modules-v2/Makefile` line.

Here are a few basic ones
```
# possible values: mm10, hg19, b37, hg19_ionref, b37_hbv_hcv, hg38
REF = b37

# possible values: ILLUMINA, IONTORRENT
SEQ_PLATFORM = ILLUMINA

# possible values: NONE (e.g. WGS), BAITS (bait-capture enrichment), PCR (amplicon-based enrichment), RNA (cDNA enrichment), CHIP
CAPTURE_METHOD = NONE

# e.g. HCC20160511, WXS
PANEL = NONE

# Single-end or paired-end, set to false if single-end
PAIRED_END = true

# possible values: SOMATIC,GERMLINE
ANALYSIS_TYPE = SOMATIC

include usb-modules-v2/Makefile
```
Most parameters are automatically set to the appropriate values if you set these above parameters correctly. 

Not all combinations are permissible. 'REF' and 'PANEL' combinations are valid if you can find a `genome_includes/<REF>.<PANEL>.inc` file.

Additional user-configurable parameters are defined (with default values) in the `usb-modules-v2/config.inc` file. 

---
# Executing the modules
This analysis pipeline is designed to be modular. The names of the modules are found in module-level Makefile (not project-level Makefile). To execute a nodule, you type
```
make <MODULE>
```
This will set the parameters you set up in the project-level Makefile, then it will go through the code to set the remaining parameters with default-values, then run your desired module. 
It is highly advisable to run this with either `nohup`, `screen` or `tmux`.

Here are some very common ones.

#### Alignment
For alignment, the following are implemented and tested.
```
make bwaaln
make bwamem
make hisat2
make star
```

#### QC
These will work for both Illumina and Ion Torrent sequencing, and will collect the appropriate metrics based on capture method and target panel.
```
make bam_metrics
make fastqc
make genotype
```

#### Germline variant calling
For Illumina, GATK v4, following the Best Practice guidelines is implemented and tested.
```
make gatk
```

For Ion Torrent, TVC is implemented but not well tested.
```
make tvc
```

#### Somatic variant calling
For Illumina, mutect (SNVs) and strelka (indels) are implemented and tested.
```
make mutect
make strelka
make mutation_summary
```

For Ion Torrent, TVC is implemented and tested.
```
make tvc_somatic
make mutation_summary
```

#### Somatic CNA detection
For Illumina, facets is implemented and tested.
```
make facets
```

For Ion Torrent, varscan is implemented and tested.
```
make varscan_cnv
```

#### ChIP-seq peak detection
MOSAICS is implemeted but not very well tested.
```
make mosaics
```

#### Others/ downstream tools
There are a lot more...

---

## Example use case 1: What to do when you get a set of Ion Torrent genomic data ##

```
PROJ_DIR=PROJ
mkdir $PROJ_DIR
cd $PROJ_DIR
```

From here on, we will assume that you are in this $PROJ_DIR.

Clone the code base
```
git clone https://github.com/charlottekyng/usb-modules-v2.git
```

Copy the Makefile_template to $PROJ (don't move or the file would disappear from the repo)
```
cp usb-modules-v2/Makefile_template Makefile
```

Rename the bam files to <sample_name>.bam and put them in $PROJ_DIR/unprocessed_bam
```
mkdir $PROJ_DIR/unprocessed_bam
cd $PROJ_DIR/unprocessed_bam
```

Make samples.txt
```
ls *bam | perl -p -e "s/\.bam//g;" > ../samples.txt
cd ..
```
Make sample_sets.txt. This file should be one patient per row. Each row should consist of the tumor samples, tab-delimited, followed by the matched normal sample as the last name on the row

Now fix read groups to ensure downstream processing do not fall over
```
make fix_rg
```

Generate some sequencing statistics
```
make bam_metrics
```

Genotype to make sure there are no mismatched samples
```
make genotype
```

Call somatic mutations
```
make tvc_somatic
```
Screen hotspots (for _TERT_ promoter) and for QC
```
make hotspot_screen
```

Make an Excel table of the mutations
```
make mutation_summary
```

Or do everything at once, if nothing falls over, this will do everything sequentially
```
make fix_rg bam_metrics genotype tvc_somatic hotspot_screen mutation_summary
```
