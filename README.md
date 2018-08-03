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
If you need to update the code based for PROJ
```
cd usb-modules-v2
git pull
```

### Setting up sample sheets

1. Sample sheet (required, default: `samples.txt`), consisting of a single column of sample names. 
    ```
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
    SAMPLE1T SAMPLE1N
    SAMPLE2T SAMPLE2N
    SAMPLE3T1 SAMPLE3T2 SAMPLE3N
    ```
1. Sample splits (required for one of these two scenarios, depending on your data, default: `samples.split.txt`)
    1. multiple sets of fastqs per sample, in which case the first column is the <SAMPLE_NAME>, the second column is <SAMPLE_NAME>\_<RUN_NAME>. This is usually required with Illumina sequencing, since the samples are usually multiplexed and sequenced over several lanes/runs. See below on setting up the `unprocessed_fastq` directory; or 
        ```
        SAMPLE1N SAMPLE1N_RUN1
        SAMPLE1N SAMPLE1N_RUN2
        SAMPLE1N SAMPLE1N_RUN3
        SAMPLE2N SAMPLE2N_RUN1
        SAMPLE3N SAMPLE3N_RUN1
        SAMPLE3N SAMPLE3N_RUN2
        ```	
    1. multiple bam files per samples that need to be merged. In which case, the first column is the <SAMPLE_NAME>, the second column is the <SAMPLE_NAME> of one of the bam files that need to be merged. This may be required when the a given sample was sequenced twice and you have 2 separately aligned BAMs. See below on setting up the `unprocessed_bam` directory
        ```
        SAMPLE1T SAMPLE1A
        SAMPLE1T SAMPLE1B
        SAMPLE2T SAMPLE2A
        SAMPLE2T SAMPLE2B
        ```
#### Important note on sample names

The preferred format is XXXnnn[TN]mm, where 
* XXX is a short project code (e.g. "HPU" for HCC plasma and urine.)
* nnn is the patient/individual identifier (e.g. 001)
* [TN] is tumor or normal
* mm is sample identifier (e.g. if there are more than one tumor sample).

Obviously, this format does not apply to all types of projects, and some variation is of course permissible.
Here are the rules
* all alphanumeric characters are allowed `[A-Za-z0-9]` 
* absolutely no space
* the only (tested) permissible symbol is `-`. Most other symbols are either known to break the pipeline or are untested.
* sample names should start with an alphabet [A-Za-z], although it is not know if the pipeline would actually fall over otherwise.

### Setting up data directories

There are several options in terms of data files:
1. If you start from FASTQs, you have a single or a single pair of fastqs per sample and you know your reads do not need trimming, then you put your files as `fastq/<SAMPLE_NAME>.1.fastq.gz` (and `fastq/<SAMPLE_NAME>.2.fastq.gz`). Then you are ready to run alignment.
    ```
    >ls fastq
    SAMPLE1T.1.fastq.gz SAMPLE1T.2.fastq.gz SAMPLE2T.1.fastq.gz SAMPLE2T.2.fastq.gz (...)
    ```
1. If you start from FASTQs, you have more than a single or more than a single pair of fastqs per sample, or your reads need trimming (e.g. adaptors) then you put your files as `unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.1.fastq.gz` or `unprocessed_fastq/<SAMPLE_NAME>.1.fastq.gz` (and `unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.2.fastq.gz` or `unprocessed_fastq/<SAMPLE_NAME>.2.fastq.gz`). Then you are ready to run alignment.
    ```
    >ls unprocessed_fastq
    SAMPLE1N_RUN1.1.fastq.gz SAMPLE1N_RUN1.2.fastq.gz SAMPLE1N_RUN2.1.fastq.gz SAMPLE1N_RUN2.2.fastq.gz (...)
    ```
    ```
    >ls unprocessed_fastq
    SAMPLE1T.1.fastq.gz SAMPLE1T.2.fastq.gz SAMPLE2T.1.fastq.gz SAMPLE2T.2.fastq.gz (...)
    ```
1. If you start from BAMs (one bam per sample), you should put all your bams as `unprocessed_bam/<SAMPLE_NAME>.bam`. Then you do `make fix_rg` then you will have analysis-ready BAMs.
    ```
    >ls unprocessed_bam
    SAMPLE1A.bam SAMPLE1B.bam SAMPLE2A.bam SAMPLE2B.bam
    ```

### Setting up analysis parameters

You will need a project-level Makefile (`$PROJ_DIR/Makefile`). 
Note that this is different from the module-level Makefile (`${PROJ_DIR}/usb-modules-v2/Makefile`).
In its most basic form, it only needs one line
```
include usb-modules-v2/Makefile
```
Any additional parameters should go _before_ the line above.

This analysis pipeline is designed to be highly configurable. 
This also means that there are many possible combinations of parameters. 
The project-level Makefile is where user-configurable parameters are specified. 

You can specify as many parameters as required in your project-level `Makefile`, 
before the `include usb-modules-v2/Makefile` line.

Here are the most basic ones and these should almost always be specified.
```
# possible values: mm10, b37, hg19_ionref, b37_hbv_hcv, hg38, b37_GRCm38 etc
REF = b37

# possible values: ILLUMINA, IONTORRENT
SEQ_PLATFORM = ILLUMINA

# possible values: NONE (e.g. WGS), BAITS (bait-capture enrichment), PCR (amplicon-based enrichment), RNA (cDNA enrichment), CHIP
CAPTURE_METHOD = NONE

# e.g. HCC20160511, WXS etc
PANEL = NONE

# Single-end or paired-end, set to false if single-end
PAIRED_END = true

# possible values: SOMATIC, GERMLINE
ANALYSIS_TYPE = SOMATIC

include usb-modules-v2/Makefile
```
Most parameters are automatically set to the basic appropriate values if you set these above parameters correctly. 

Not all combinations of REF and PANEL are permissible. 
(In the near future, valid combinations of 'REF' and 'PANEL' will be found as a `usb-modules-v2/genome_inc/<REF>/<PANEL>.inc` file.)

Additional user-configurable parameters are defined (with default values) in the `usb-modules-v2/config.inc` file. 
(In the near future, the parameters will be better documented in the config file.)

---
# Executing the modules
This analysis pipeline is designed to be modular. 
The names of the modules are found in module-level Makefile (not project-level Makefile). 
To execute a nodule, you type
```
make <MODULE>
```
This will set the parameters you set up in the project-level Makefile, 
then it will go through the code to set the remaining parameters with the appropriate default values, 
then run your desired module. 
It is highly advisable to run this with either `nohup`, or within `screen` or `tmux`.

Here are some very common modules. 
**Note:** Some of them have dependencies that are not well-documented 
(hopefully this will be improved in the future). 
The sequences in "Example recipes" section below are valid sequences.

#### Alignment
For genomic Illumina alignment, the following are implemented and tested.
```
make bwaaln
make bwamem
```
For transcriptomic Illumina sequencing, the following are implemented and tested.
```
make hisat2
make star
```

#### QC
The following will work for both Illumina and Ion Torrent sequencing, 
and will collect the appropriate metrics based on capture method and target panel.
```
make bam_metrics
make fastqc
make genotype
```

#### Germline variant calling
For Illumina, GATK v4 following the Best Practice guidelines is implemented and tested.
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
MOSAICS is implemented but not very well tested.
```
make mosaics
```

#### Others/ downstream tools
There are a lot more...

---

# If something falls over...

You should look in `$PROJ_DIR/log`. The log file will be named in the format `$PROJ_DIR/log/<module>.<date>.<attempt>.log`.
For example in `log/gatk.2018-08-03.2.log`, you should find lines like these.
```
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00245.variants.vcf.gz] Error 1
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00237.variants.vcf.gz] Error 1
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00303.variants.vcf.gz] Error 1
```

Each step in the pipeline spits out a log file. The individual log file for the 1st line in the example above is then
```
log/gatk.2018-08-03.2/gatk/intervals_gvcf/0030/ESBIPGRA00245.variants.vcf.gz.log
```

This is where you go looking for the error. 
Most errors are related to incorrect parameters, empty or invalid input files or running out of resources (time or memory).
If it is related to an invalid/empty input file, then go one step back in the pipeline to see 
if a previous step fell over without throwing and error (it happens).

If the reason is not obvious, try deleting any invalid/empty files, then re-run it. 
Sometimes there are transient system glitches and a simple re-run is enough to fix it.

---

# Example recipes
Assuming that you have set up the above correctly, 
here are some suggested recipes that are valid sequences.

#### Whole-exome sequencing on Illumina
```
make bwamem genotype fastqc bam_metrics facets mutect strelka mutation_summary
```
#### RNA-sequencing on Illumina
```
make star bam_metrics rsem
```
#### ChIP-seq on Illumina (75bp reads or longer)
```
make bwamem mosaics
```
#### ChIP-seq on Illumina (<75bp reads)
```
make bwaaln mosaics
```
#### Targeted panel sequencing on Ion Torrent (from bam files)
```
make fix_rg genotype bam_metrics tvc_somatic varscan_cnv hotspot_screen
```
#### Whole-genome sequencing on Illumina
```
make bwamem bam_metrics gatk
```

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
Make sample_sets.txt. This file should be one patient per row. 
Each row should consist of the tumor samples, tab-delimited, 
followed by the matched normal sample as the last name on the row

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
