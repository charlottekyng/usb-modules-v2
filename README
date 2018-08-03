###############################################################
## Basic project set up

## git clone https://github.com/charlottekyng/usb-modules-v2.git

## There are the mandatory/optional files:
## Makefile (required, has to called Makefile), see 
##      **usb-modules-v2/Makefile_template** for an example and
##      **usb-modules-v2/config** for options
## sample sheet (required, default: samples.txt), consisting of a single 
##      column of sample names. 
## NOTE: sample names must be [A-Za-z0-9\-] that starts with an alphabet. 
##     
## sample sets (required for ANALYSIS_TYPE = SOMATIC,
##      default: sample_sets.txt, consisting of samples within
##      sample sheet above, each row contains all samples from
##      a given individual, with the germline sample last,
##      tab or space delimited
## sample splits (required for 
##      1, multiple sets of fastqs per sample, in which case 
##           the first column is the <SAMPLE_NAME>, 
##           the second column is <SAMPLE_NAME>_<RUN_NAME>
##           See below on setting up the unprocessed_fastq directory
##      or 
##      2, multiple bam files per samples that need to be merged.
##           in which case,
##           the first column is the <SAMPLE_NAME>
##           the second column is the <SAMPLE_NAME> of one of the
##                bam files that need to be merged

## If starting from FASTQs and you have a single or a single pair
##      of fastqs, then you put your files as
##      **fastq/<SAMPLE_NAME>.1.fastq.gz fastq/<SAMPLE_NAME>.2.fastq.gz**
## If starting from FASTQs and you have more than a single or more
##      than a single pair of fastqs, then you put your files as
##      **unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.1.fastq.gz**
##      **unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.2.fastq.gz**
## Then you just run alignment (see below)

## If you start from BAMs, you should put all your bams as
##      **unprocessed_bam/<SAMPLE_NAME>.bam**
## Then you do **make fix_rg** then you will have
##     analysis-ready BAMs







###############################################################
## What to do when you get a set of Ion Torrent genomic data ##
###############################################################

# assuming the new project is called PROJ
PROJ_DIR=/scicore/home/terracci/GROUP/data/PROJ

# make a project directory under /scicore/home/terracci/GROUP/data
mkdir $PROJ_DIR
cd $PROJ_DIR

# get usb-modules-v2 code from github
git clone git@github.com:charlottekyng/usb-modules-v2.git; 

# copy the Makefile_template to $PROJ (don't move or the file would disappear from the repo)
cp usb-modules-v2/Makefile_template Makefile

# rename the bam files to <sample_name>.bam and put them in $PROJ_DIR/unprocessed_bam
mkdir $PROJ_DIR/unprocessed_bam
cd $PROJ_DIR/unprocessed_bam
# dump the files in the directory

# make samples.txt
ls *bam | perl -p -e "s/\.bam//g;" > ../samples.txt
cd ..

# make sample_sets.txt
# this file should be one patient per row
# each row should consist of the tumor samples, tab-delimited, 
# followed by the matched normal sample as the last name on the row

# now fix read groups to ensure downstream processing do not fall over
make fix_rg

# generate some sequencing statistics
make bam_metrics

# genotype to make sure there are no mismatched samples
make genotype

# call somatic mutations
make tvc_somatic

# screen hotspots (for TERT) and for QC
make hotspot_screen

# make an Excel table of the mutations
make mutation_summary

### or do everything at once, if nothing falls over, this will do everything sequentially
make fix_rg bam_metrics genotype tvc_somatic hotspot_screen mutation_summary
