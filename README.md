# Basic project set up

Assuming the new project is called PROJ
> PROJ_DIR=/scicore/home/terracci/GROUP/data/PROJ

Make a project directory under /scicore/home/terracci/GROUP/data
> mkdir $PROJ_DIR
> cd $PROJ_DIR

From here on, we will assume that you are in this $PROJ_DIR.

Clone the code base
> git clone https://github.com/charlottekyng/usb-modules-v2.git

You will need to Makefiles. The first is the one usb-modules-v2/Makefile (hereafter known as module-level Makefile), and a project-level Makefile (described below in **set-up** files).

There are the mandatory/optional **set-up** files:
1. Project-level Makefile (required, has to called Makefile), see **usb-modules-v2/Makefile_template** for an example and **usb-modules-v2/config** for options
1. sample sheet (required, default: samples.txt), consisting of a single column of sample names. 
1. sample sets (required for ANALYSIS_TYPE = SOMATIC, default: sample_sets.txt, consisting of samples within sample sheet above, each row contains all samples from a given individual, with the germline sample last, tab or space delimited
1. sample splits (required for 
  1. multiple sets of fastqs per sample, in which case the first column is the <SAMPLE_NAME>, the second column is <SAMPLE_NAME>_<RUN_NAME>. See below on setting up the unprocessed_fastq directory; or 
  1. multiple bam files per samples that need to be merged. In which case, the first column is the <SAMPLE_NAME>, the second column is the <SAMPLE_NAME> of one of the bam files that need to be merged

**NOTE**: sample names must be [A-Za-z0-9\-] that starts with an alphabet.

There are several options in terms of **data** files:
1. If you start from FASTQs and you have a single or a single pair of fastqs per sample, then you put your files as **fastq/<SAMPLE_NAME>.1.fastq.gz (and fastq/<SAMPLE_NAME>.2.fastq.gz**). Then you are ready to run alignment.
1. If you start from FASTQs and you have more than a single or more than a single pair of fastqs per sample, then you put your files as **unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.1.fastq.gz (and **unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.2.fastq.gz**). Then you are ready to run alignment.
1. If you start from BAMs (one bam per sample), you should put all your bams as **unprocessed_bam/<SAMPLE_NAME>.bam**. Then you do **make fix_rg** then you will have analysis-ready BAMs.

This analysis pipeline is designed to be modular and highly configurable. The names of the modules are found in module-level Makefile (not project-level Makefile). To execute a nodule, you type
> make <MODULE>

Most user-configurable parameters are in the config.inc file. You can specify as many parameters as required in your project-level Makefile.



## Example use case 1: What to do when you get a set of Ion Torrent genomic data ##

Assuming the new project is called PROJ
> PROJ_DIR=/scicore/home/terracci/GROUP/data/PROJ

Make a project directory under /scicore/home/terracci/GROUP/data
> mkdir $PROJ_DIR
> cd $PROJ_DIR

Get the code base
> git clone https://github.com/charlottekyng/usb-modules-v2.git

Copy the Makefile_template to $PROJ (don't move or the file would disappear from the repo)
> cp usb-modules-v2/Makefile_template Makefile

Rename the bam files to <sample_name>.bam and put them in $PROJ_DIR/unprocessed_bam
> mkdir $PROJ_DIR/unprocessed_bam
> cd $PROJ_DIR/unprocessed_bam

Make samples.txt
> ls *bam | perl -p -e "s/\.bam//g;" > ../samples.txt
> cd ..

Make sample_sets.txt. This file should be one patient per row. Each row should consist of the tumor samples, tab-delimited, followed by the matched normal sample as the last name on the row

Now fix read groups to ensure downstream processing do not fall over
> make fix_rg

Generate some sequencing statistics
> make bam_metrics

Genotype to make sure there are no mismatched samples
> make genotype

Call somatic mutations
> make tvc_somatic

Screen hotspots (for TERT) and for QC
> make hotspot_screen

Make an Excel table of the mutations
> make mutation_summary

Or do everything at once, if nothing falls over, this will do everything sequentially
> make fix_rg bam_metrics genotype tvc_somatic hotspot_screen mutation_summary
