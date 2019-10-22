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
If you need to update the code base for PROJ
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
* _ABSOLUTE NO_ white space or underscore (`_`)
* the only (tested) permissible symbol is `-`. Most other symbols are either known to break the pipeline or are untested.
* sample names should start with an alphabet [A-Za-z], although it is not know if the pipeline would actually fall over otherwise.

Adhere to these guidelines to avoid unnecessary troubleshooting.

### Setting up data directories

There are several options in terms of data files:
1. If you start from FASTQs, you have a single fastq or a single pair of fastqs per sample and you know your reads do not need trimming, then you put your files as `fastq/<SAMPLE_NAME>.1.fastq.gz` (and `fastq/<SAMPLE_NAME>.2.fastq.gz`). Then you are ready to run alignment. 
    ```
    >ls fastq/
    SAMPLE1T.1.fastq.gz SAMPLE1T.2.fastq.gz SAMPLE2T.1.fastq.gz SAMPLE2T.2.fastq.gz (...)
    ```
1. If you start from FASTQs, you have more than a single fastq or more than a single pair of fastqs per sample, or your reads need trimming (e.g. adaptors) then you put your files as `unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.1.fastq.gz` or `unprocessed_fastq/<SAMPLE_NAME>.1.fastq.gz` (and `unprocessed_fastq/<SAMPLE_NAME>_<RUN_NAME>.2.fastq.gz` or `unprocessed_fastq/<SAMPLE_NAME>.2.fastq.gz`). With this option, you will need the `samples.split.txt` file (see above). Then you are ready to run alignment.
    ```
    >ls unprocessed_fastq/
    SAMPLE1N_RUN1.1.fastq.gz SAMPLE1N_RUN1.2.fastq.gz SAMPLE1N_RUN2.1.fastq.gz SAMPLE1N_RUN2.2.fastq.gz (...)
    ```
    ```
    >ls unprocessed_fastq/
    SAMPLE1T.1.fastq.gz SAMPLE1T.2.fastq.gz SAMPLE2T.1.fastq.gz SAMPLE2T.2.fastq.gz (...)
    ```
1. If you start from BAMs (one bam per sample), you should put all your bams as `unprocessed_bam/<SAMPLE_NAME>.bam`. Then you do `make fix_rg` then you will have analysis-ready BAMs.
    ```
    >ls unprocessed_bam/
    SAMPLE1A.bam SAMPLE1B.bam SAMPLE2A.bam SAMPLE2B.bam
    ```
**Note**: for single-end FASTQs, use `.1.fastq.gz` (i.e. include the `.1`).

### Setting up analysis parameters

You will need a project-level Makefile (`${PROJ_DIR}/Makefile`). 
Note that this is different from the module-level Makefile (`${PROJ_DIR}/usb-modules-v2/Makefile`).
In its most basic form, it only needs one line
```
include usb-modules-v2/Makefile
```
This analysis pipeline is designed to be highly configurable. 
This also means that there are many possible combinations of parameters. 
The project-level Makefile is where user-configurable parameters are specified. 

You can specify as many parameters as required in your project-level `Makefile`, 
_before_ the `include usb-modules-v2/Makefile` line.

Here are the most basic ones and these should almost always be specified.
```
# example values:  b37, hg19_ionref, hg38 etc. Values permitted will have a `usb-modules-v2/genome_inc/$<REF>` directory
REF = b37

# possible values: [ILLUMINA|IONTORRENT]
SEQ_PLATFORM = ILLUMINA

# possible values: NONE (e.g. WGS), BAITS (bait-capture enrichment), PCR (amplicon-based enrichment), RNA (cDNA enrichment), CHIP
CAPTURE_METHOD = NONE

# example values: HCC20160511, WXS etc. Values permitted will have a `usb-modules-v2/genome_inc/$<REF>/$<PANEL>.inc` file
PANEL = NONE

# Single-end or paired-end, set to false if single-end [true|false]
PAIRED_END = true

# possible values: [SOMATIC|GERMLINE]
ANALYSIS_TYPE = SOMATIC

include usb-modules-v2/Makefile
```
Most parameters are automatically set to the basic appropriate values if you set these above parameters correctly, but there is a lot of room for customization.

Not all combinations of REF and PANEL are permissible. With the exception of `PANEL=NONE`, make sure your combination exists as a `usb-modules-v2/genome_inc/$<REF>/$<PANEL>.inc` file.

Additional user-configurable parameters are defined (with default values) in the `usb-modules-v2/config.inc` file. 

Some `Makefile` templates are provided in `usb-modules-v2/Makefile_templates/`. Please *copy* them to your project directory and do not remove them from `usb-modules-v2/`.
```
>ls usb-modules-v2/Makefile_templates/
Makefile_template_all_basic_options  Makefile_template_agilentallexonv6  Makefile_template_iontorrent_comprehensive_panel  Makefile_template_rnaseq_xenografts
>cp usb-modules-v2/Makefile_templates/Makefile_template_agilentallexonv6 Makefile
```

There are many, many possibilities to customize the analysis. The ones listed above are merely the most basic ones.

---
# Executing the modules
This analysis pipeline is designed to be modular. 
The names of the modules are found in module-level Makefile (i.e. `usb-modules-v2/Makefile`, not project-level Makefile). 
The ones listed below are only the most basic modules. There are many more advanced and/or less commonly used modules implemeted, particularly with re-processing bam/VCF files and downstream analyses.
To execute a nodule, you type
```
make <MODULE>
```
This will set the parameters you set up in the project-level `Makefile`, 
then it will go through the code to set the remaining parameters with the appropriate default values, 
then run your desired module. 
It is highly advisable to run this with either `nohup`, or within `screen` or `tmux`.

Here are some very common modules. 
**Note:** Some of them have dependencies that are not well-documented 
(hopefully this will be improved in the future). 
The sequences in "Example recipes" section below are valid sequences.

### Alignment

This runs the chosen aligner on FASTQ files, including preprocessing (e.g. adaptor trimming) and postprocessing (e.g. sorting, deduplication).

*Pre-requisites:* FASTQs in `fastq/` or `unprocessed_fastq/` (see 'Setting up data directories' above).

For genomic Illumina alignment, the following are implemented and tested.
```
make bwaaln     ### for reads < 75bp
make bwamem     ### for reads >= 75bp
```
For transcriptomic Illumina sequencing, the following are implemented and tested.
```
make hisat2
make star
```

### QC
Most of the following will work for both Illumina and Ion Torrent sequencing, unless otherwise specified.

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner.

```
make bam_metrics       # This should be done for every dataset
make fastqc            # This is occasionally useful for checking the quality of the sequencing (but isn't too useful most of the time)
make genotype          # This is useful for confirming that sample pairs/sets came from the correct patient
make facets_poolednorm # (Illumina only) This is useful for confirming that the germline samples are not contaminated with tumor cells
make facets            # (Illumina only) This is useful for checking tumor content
```

### Germline variant calling

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner.

For Illumina, GATK v4 following the Best Practice guidelines is implemented and tested. 
*IMPORTANT*: If you are using targeted sequencing (e.g. WES or any other kind of targeted panel, set `NUM_JOBS=3` or a similarly small number as one of the steps in GATK writes a huge number of files. This does not apply to WGS.
```
make gatk 
```

For Ion Torrent, TVC is implemented.
```
make tvc
```

### Somatic variant calling

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner.

For Illumina, mutect/2 (SNVs) and strelka/2 (indels) are implemented and tested. 

*Note*: It is generally advisable to run facets for CNAs before these. If you do not have facets results, you have to set `ANN_FACETS=false` in your `Makefile`.
```
make mutect              ## for b37 and variants
make strelka             ## for b37 and variants
make mutect2             ## for hg38 and variants
make strelka2            ## for hg38 and variants
make mutation_summary
```

For Ion Torrent, TVC (both SNVs and indels) is implemented and tested.
```
make tvc_somatic         ## old implementarion of PipeIT
make pipeit              ## use this
make mutation_summary
```

`mutation_summary` returns an Excel file for all mutations with simplified annotation. Here is a description of the fields included in the mutation summary:
* IMPACT: Impact of effect (HIGH=frameshift indel, truncating, splice site and similar; MODERATE=missense, in-frame indel and similar; LOW=synonymous and similar) Details: snpeff.sourceforge.net/VCFannotationformat_v1.0.pdf
* FA: variant allele fraction, meaning number of variant reads divided by total number of reads at the locus
* DP: depth, meaning total number of reads at the locus
* AD: ‘allelic depth’, meaning number of reference allele reads and number of variant allele reads, separated by a comma
* Hotspot_*, Hotspot3D_*: mutation hotspots according to Chang et al (Cancer Discov 2017, cancerhotspots.org) and Gao et al (Genome Med 2017, 3dhotspots.org)
* CancerGeneSets: list of cancer gene lists to which the gene belongs (see list below)
* facetsCF: clonal fraction of the copy number segment at the locus (experimental, to be deprecated)
* facetsTCN_EM: Total copy number at the locus
* facetsLCN_EM: Lesser (minor) copy number at the locus
* facetsLOHCall: loss of heterozygosity at the locus (i.e. LCN_EM=0)
* CHROM, POS, REF, ALT: chromosome position, reference and alternate (variant) alleles
* ID: dbSNP ID, COSMIC ID
* ExACnontcga_AC/AF: Population level allele count (AC) and allele frequency (AF).

*Important note re: multi-tumor cases and the FILTER field*: 
For a given patient, somatic mutations found in one (or more) sample by the mutation calling pipeline were additionally 
interrogated in all remaining tumor samples of this patient. (This "interrogation" is known to the pipeline as "sufam", FYI).
“interrogation” and “interrogation_absent” in FILTER indicate mutations that were not called in the given sample but were 
found to be supported by some sequencing reads (i.e. “interrogation”), or found not to be supported by any read 
(i.e. “interrogation_absent”), respectively. These additional ‘interrogated’ mutations have been included 1) to ensure we don't
falsely claim a mutation is not present in a particular tumor sample because of low sequencing depth and/or low but non-zero VAF (false negative) 
and 2) to aid clonal evolution analyses (e.g. PyClone wants counts for all mutations in all samples regardless whether the mutations
are in a sample or now).  This "interrogation" step is performed per-row in the sample_sets.txt, regardless whether the tumors 
are clonally related or not, and should be considered meaningless for tumors that were not clonally related. 
Keep in mind that ‘interrogation’ mutations are usually at low variant allele frequencies, or were initially filtered out 
for other quality/depth reasons, and should therefore be considered lower confidence. In general, mutations at low 
variant allele fraction are also generally of lower confidence (e.g. sequencing error). Thus, two tumors that did not 
share a single high confidence mutation, or a very small number of mutations but all of which were of low confidence or 
low variant allele fractions, should be carefully studied to determine clonal relatedness.
If you want to turn of this "interrogation" step, set `USE_SUFAM = false`.

Cancer gene sets:
* KANDOTH_127: 127 SMGs from Kandoth et al
* LAWRENCE_CANCER5000S: Cancer5000-S set from Lawrence et al
* TCGA_LIHC: SMGs from TCGA-LIHC
* TCGA_LIHC_extended: Extended SMGs from TCGA-LIHC
* SCHULZE_HCC: SMGs from Schulze et al
* FUJIMOTO_HCC: SMGs from Fujimoto et al
* CANCER_GENE_CENSUS_TIER1_V88: Cancer Gene Census Tier 1 v88
* CANCER_GENE_CENSUS_TIERS1_AND_2_V88: Cancer Gene Census Tiers 1 and 2 v88
* MARTINCORENA_ACC: ACC SMGs from Martincorena et al
* MARTINCORENA_BLCA: BLCA SMGs from Martincorena et al
* MARTINCORENA_BRCA: BRCA SMGs from Martincorena et al
* MARTINCORENA_CESC: CESC SMGs from Martincorena et al
* MARTINCORENA_COREAD: COREAD SMGs from Martincorena et al
* MARTINCORENA_ESCA: ESCA SMGs from Martincorena et al
* MARTINCORENA_GBM: GBM SMGs from Martincorena et al
* MARTINCORENA_HNSC: HNSC SMGs from Martincorena et al
* MARTINCORENA_KICH: KICH SMGs from Martincorena et al
* MARTINCORENA_KIRC: KIRC SMGs from Martincorena et al
* MARTINCORENA_KIRP: KIRP SMGs from Martincorena et al
* MARTINCORENA_LAML: LAML SMGs from Martincorena et al
* MARTINCORENA_LGG: LGG SMGs from Martincorena et al
* MARTINCORENA_LIHC: LIHC SMGs from Martincorena et al
* MARTINCORENA_LUAD: LUAD SMGs from Martincorena et al
* MARTINCORENA_LUSC: LUSC SMGs from Martincorena et al
* MARTINCORENA_MESO: MESO SMGs from Martincorena et al
* MARTINCORENA_OV: OV SMGs from Martincorena et al
* MARTINCORENA_PAAD: PAAD SMGs from Martincorena et al
* MARTINCORENA_PCPG: PCPG SMGs from Martincorena et al
* MARTINCORENA_PRAD: PRAD SMGs from Martincorena et al
* MARTINCORENA_SARC: SARC SMGs from Martincorena et al
* MARTINCORENA_SKCM: SKCM SMGs from Martincorena et al
* MARTINCORENA_STAD: STAD SMGs from Martincorena et al
* MARTINCORENA_TGCT: TGCT SMGs from Martincorena et al
* MARTINCORENA_THYM: THYM SMGs from Martincorena et al
* MARTINCORENA_UCEC: UCEC SMGs from Martincorena et al
* MARTINCORENA_UCS: UCS SMGs from Martincorena et al
* MARTINCORENA_PANCANCER: PANCAN SMGs from Martincorena et al
* BAILEY_ACC: ACC SMGs from Bailey et al
* BAILEY_BLCA: BLCA SMGs from Bailey et al
* BAILEY_BRCA: BRCA SMGs from Bailey et al
* BAILEY_CESC: CESC SMGs from Bailey et al
* BAILEY_CHOL: CHOL SMGs from Bailey et al
* BAILEY_COADREAD: COADREAD SMGs from Bailey et al
* BAILEY_DLBC: DLBL SMGs from Bailey et al
* BAILEY_ESCA: ESCA SMGs from Bailey et al
* BAILEY_GBM: GBM SMGs from Bailey et al
* BAILEY_HNSC: HNSC SMGs from Bailey et al
* BAILEY_KICH: KICH SMGs from Bailey et al
* BAILEY_KIRC: KIRC SMGs from Bailey et al
* BAILEY_KIRP: KIRP SMGs from Bailey et al
* BAILEY_LAML: LAML SMGs from Bailey et al
* BAILEY_LGG: LGG SMGs from Bailey et al
* BAILEY_LIHC: LIHC SMGs from Bailey et al
* BAILEY_LUAD: LUAD SMGs from Bailey et al
* BAILEY_LUSC: LUSC SMGs from Bailey et al
* BAILEY_MESO: MESO SMGs from Bailey et al
* BAILEY_OV: OV SMGs from Bailey et al
* BAILEY_PAAD: PAAD SMGs from Bailey et al
* BAILEY_PANCAN: PANCAN SMGs from Bailey et al
* BAILEY_PCPG: PCPG SMGs from Bailey et al
* BAILEY_PRAD: PRAD SMGs from Bailey et al
* BAILEY_SARC: SARC SMGs from Bailey et al
* BAILEY_SKCM: SKCM SMGs from Bailey et al
* BAILEY_STAD: STAD SMGs from Bailey et al
* BAILEY_THCA: THCA SMGs from Bailey et al
* BAILEY_THYM: THYM SMGs from Bailey et al
* BAILEY_UCEC: UCEC SMGs from Bailey et al
* BAILEY_UCS: UCS SMGs from Bailey et al
* BAILEY_UVM: UVM SMGs from Bailey et al

Sources: Kandoth et al (PMID 24132290), Lawrence et al (PMID 23770567), Schultze et al (PMID 25822088), Fujimoto et al (PMID 27064257), Martincorena et al (PMID 29056346), Bailey et al (PMID 30096302)
See https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations for TCGA code abbreviations

## If you have experiments other than matched tumor-normal pairs/sets from frozen samples...
The defaults of the pipeline are tuned towards tumor-normal pairs/sets from frozen samples. If you have other sample types, you might have to consider adding the following steps:
* for FFPE samples: formalin fixation artefacts are typically a strong enrichment of C>T/G>A variants at low VAF (<15%), 
frequently accounting for 90%+ of all somatic variants identified. You may have to consider performing an additional filter by removing C>T/G>A variants <10% or <15%, 
or supported by <5 reads. You may also consider remove all variants with <5 reads, in which case you can set `MIN_TUMOR_AD = 5` (sorry there is no parameter to do a blanket filtering for VAF for now).
* for tumors without a matched normal (e.g. cell lines, archival materials without matched normal): do as many of the following as you can
	* make a BAM file from a bunch of normals captured and sequenced the same way. To do this, you list these samples in `samples.poolednorm.txt`,
	make sure their BAM files are in the `bam/` directory, set `SAMTOOLS_DOWNSAMPLE_FACTOR` such that you are sampling roughly 1/n (where n is the number of normals,
	see `samtools view -s`), then run `make poolednorm_bam`, which generates a new sample called `poolednorm` - which you then use as your matched normal.
	Proceed to regular somatic mutations/CNA calling. If you don't have normals captured and sequenced in the same project, find the closest thing from our collection of data. Any normal is better than no normal.
	* make a panel of normal for filtering mutations. 


### Somatic CNA detection

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner.

For Illumina DNA sequencing, facets is implemented and tested. It is recommended that you run facets *before* 
somatic variant calling because the last step of variant annotation involves annotating copy number states to mutations. 
If you do not run facets before variant calling, you have to set `ANN_FACETS=false` in your `Makefile`.
```
make facets
```

`facets` returns figures and tables. The segment files (`*.cncf.txt`) have these columns:
* seg: segment number
* num.mark: number of SNPs in the segment
* nhet: number of heterozygous SNPs in the segment
* cnlr.median: log ratio
* mafR: log-odds-ratio summary of the segment
* segclust: segment cluster
* cnlr.median.clust: log ratio of segment cluster
* mafR.clust: mafR of segment cluser
* start/end: genomic position of segment start
* cf.em: cellular fraction of the segment
* tcn.em: total copy number
* lcn.em: minor (lesser) copy number
* clonal.cluster: clonal cluster

3 versions of gene-level copy number alterations have been generated
* GL_ASCNA: copy number status derived from total copy number. 
    * 2: ≥ploidy+4 (amplification)
    * 1: ≥ ploidy+1 (low-level copy number gain)
    * 0: ploidy (copy number neutral)
    * -1: < ploidy but not 0 (heterozygous loss)
    * -2: total copy number = 0 (homozygous deletion)

* GL_LRR: copy number status derived from log ratio
    * 2: amplification
    * 1: copy number gain
    * 0: copy number neutral
    * -1: copy number loss
    * -2: deep deletion
* tcn.em: total copy number (absolute copy number)

Additional files are 1) a list of amplifications and homozygous deletions derived from ASCNA (`*.ampdel.txt`) and 2) per-sample copy number profile figure (`*cncf.pdf`).

**Note**: FACETS profiling does not always work well. In particular samples in which >5% of the genome at copy number 0 (homozygous deletions) should be excluded from copy number analysis.

For Ion Torrent DNA sequencing, Varscan is implemented and tested.
```
make varscan_cnv
```

For Illumina RNA sequencing, cnvkit is implemeted (but currently not well tested or documented).
```
make cnvkit
```

### RNA-seq transcript quantification
RSEM is tested to be run after STAR alignment.
```
make rsem
```

`rsem` runs RSEM and performs some additional normalization for convenience.
* expected_count.results: expected counts from RSEM. 
* expected_count.results_coding: subset of protein coding genes from genes.expected_count.results_coding. 
* expected_count.results_coding_uq: expected_count.results_coding normalized by upper-quartile normalization according to the TCGA HCC paper (in this case to the fixed value of 1000). Values can be compared between samples, but NOT between genes. This is what TCGA used for clustering.
* RSEM-TPM/FPKM: transcripts-per-million and fragments-per-kilobase-million: see http://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/


### ChIP-seq peak detection
MOSAICS is implemented but not very well tested. In particular, it almost always falls over with paired-end data.

*Pre-requisites:* BAMs in `bam/` after alignment with an appropriate aligner (bwaaln or bwamem).
```
make mosaics
```

### Others/ downstream tools
There are a lot more... 

For exome analysis, there are a few things that are useful. These should work if you use them in the context of the suggested recipes below. Some of them may only work on the b37 genome.
You may run into errors if you run them outside of the context of in-house, standard data as they have complex (and cryptic) rules to obtain input files. 
```
make deconstruct_sigs # For mutational signatures, requires mutations
make lst              # For the detection of large-scale transitions, requires facets output
make msisensor        # For the detection of microsatellite instability, requires bam files only
make pyclone          # For clonality analysis, requires mutations and facets output
make absolute_seq     # For clonality analysis, requires mutations and facets output (see note below)
make pvacseq          # For the detection of neo-antigens, requires mutations (not well tested...)
```

_Note regarding absolute_seq_, it would attempt to run all 3 steps, choosing the default solutions. But, it does not always run to the end, as RunAbsolute error handing is not great.


### Note regarding sanity checks

At the moment, there are no checks in place to see if what you are attempting to run is a sensible thing to do given your parameters.

---

# Troubleshooting

You should look in `$PROJ_DIR/log/`. The log file will be named in the format `$PROJ_DIR/log/<module>.<date>.<attempt>.log`.
If you find the file `$PROJ_DIR/log/<module>.<date>.<attempt>.log`, but not the directory `$PROJ_DIR/log/<module>.<date>.<attempt>/`, then no job was submitted. 
If there was the log file and the log directory, then your jobs were submitted.

### If it falls over immediately... (jobs not submitted)

This usually happens because 1) files not found or pre-requisite not met, 2) there is a bug (or ten) in the code, or 3) there is a problem with your sample sheets. 
(1) and (2) will usually be met with a `'No rule found to make <file>'`, which means make could not locate the correct recipes or the required file/s.
(3) will usually be met with something like `'*** non-numeric second argument to 'wordlist' function: '-1''`.

1. Check that your file and directory names are correct, especially if you are attempting to run the pipeline on a fresh set of data.

1. Many modules have un/documented, obvious or not obvious, prerequisites, e.g. `mutect` requires aligned data, `mutation_summary` requires mutations. 
Check/compile the prerequiresites then try again.

1. If there are bugs in the code, you will want to see where it stops finding the correct recipes. 
You can try, e.g., `make --debug=i -nf usb-modules-v2/aligners/bwamemAligner.mk REF=b37 SEQ_PLATFORM=ILLUMINA (...) | less`
(the parameters in your project-level Makefile). This will produce a verbose dry-run of the files make is attempting to generate. 
Here you can see where it stops finding the recipes. This is very useful for debugging, and for educational purposes if you want to know what is being done.

1. Check your samples sheets. Make sure there are no stray spaces/tabs at the end of the lines. Make sure there are no blank lines (after the last samples).
    ````
    >cat -A sample_sets.txt
    SSA001T^ISSA001N$          # OK
    SSA002T^ISSA002N $         # stray space at the end
    SSA005T^ISSA005N^I$        # stray tab at the end
    SSA006T^I SSA006N$         # stray space in the middle
    $                          # remove blank lines
    ````

1. Check your Makefile. Check for stray symbols and white spaces. Check for correct panel/target bed files.
But errors from Makefile settings will more often lead to failed jobs (see below).

### If submitted jobs fail...

For example in `log/gatk.2018-08-03.2.log`, you should find lines that say "Error" like this.
```
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00245.variants.vcf.gz] Error 1
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00237.variants.vcf.gz] Error 1
make[1]: *** [gatk/intervals_gvcf/0030/ESBIPGRA00303.variants.vcf.gz] Error 1
```

Each step in the pipeline spits out a log file. The individual log file for the 1st line in the example above is then
```
log/gatk.2018-08-03.2/gatk/intervals_gvcf/0030/ESBIPGRA00245.variants.vcf.gz.log
```

Most errors are related to incorrect parameters, empty or invalid input files or running out of resources (time or memory).

1. Check your Makefile. Check the combination of genome and target panel/platform is valid (i.e. defined in `usb-modules-v2/genome_inc/`). 
Check for typos. Check the log to see if the parameters for individual steps are correct.

1. If it is related to an invalid/empty input file, then go one step back in the pipeline to see if a previous step fell over without throwing and error (it happens).

1. If the reason is not obvious, try deleting any invalid/empty files, then re-run it. Sometimes there are transient system glitches and a simple re-run is enough to fix it.

1. If it is related to resources, re-run it once or twice more. Some tools (e.g. GATK) occasionally get stuck for unknown reasons, 
or there were transient system glitches that cause something to be stuck. If the tool fails on the same sample several times, then tell Charlotte...

1. If your parameters seem to be correct but the commands in the log file are not correct, there could be a bug (or ten).
---

# Example recipes
Assuming that you have set up the above correctly, 
here are some suggested recipes that are valid sequences.

#### Whole-exome sequencing on Illumina
```
make bwamem genotype bam_metrics facets mutect(2) strelka(2) mutation_summary (deconstruct_sigs lst msisensor pyclone) [absolute_seq pvacseq]
```
Those in parentheses `()` will work, those in brackets `[]` may fall over.
#### RNA-sequencing on Illumina
```
make star bam_metrics rsem (cnvkit)
```
#### ChIP-seq on Illumina (75bp reads or longer)
```
make bwamem mosaics
```
#### ChIP-seq on Illumina (<75bp reads)
```
make bwaaln mosaics
```
#### Targeted panel sequencing on Ion Torrent (from bam files from the Torrent server)
```
make fix_rg genotype bam_metrics tvc_somatic varscan_cnv hotspot_screen mutation_summary
```
#### Whole-genome sequencing on Illumina for germline analysis
```
make bwamem bam_metrics gatk
```

## Example use case 1: What to do when you get a set of Ion Torrent genomic data ##

1. Go to the project directory
    ```
    PROJ_DIR=PROJ
    mkdir $PROJ_DIR
    cd $PROJ_DIR
    ```

1. Clone the code base
    ```
    git clone https://github.com/charlottekyng/usb-modules-v2.git
    ```

1. Copy the Makefile_template to $PROJ (don't move or the file would disappear from the repo)
    ```
    cp usb-modules-v2/Makefile_templates/Makefile_template_iontorrent_comprehensive_panel Makefile
    ```

1. Rename the bam files to <sample_name>.bam and put them in $PROJ_DIR/unprocessed_bam
    ```
    mkdir $PROJ_DIR/unprocessed_bam
    cd $PROJ_DIR/unprocessed_bam
    ```

1. Make samples.txt
    ```
    ls *bam | perl -p -e "s/\.bam//g;" > ../samples.txt
    cd ..
    ```

1. Make sample_sets.txt. This file should be one patient per row. 
Each row should consist of the tumor samples, tab-delimited,  followed by the matched normal sample as the last name on the row

1. Now fix read groups to ensure downstream processing do not fall over
    ```
    make fix_rg
    ```

1. Generate some sequencing statistics
    ```
    make bam_metrics
    ```

1. Genotype to make sure there are no mismatched samples
    ```
    make genotype
    ```

1. Call somatic mutations
    ```
    make tvc_somatic
    ```

1. Screen hotspots (for _TERT_ promoter) and for QC
    ```
    make hotspot_screen
    ```

1. Make an Excel table of the mutations
    ```
    make mutation_summary
    ```

Or do everything at once, if nothing falls over, this will do everything sequentially
```
make fix_rg bam_metrics genotype tvc_somatic hotspot_screen mutation_summary
```
