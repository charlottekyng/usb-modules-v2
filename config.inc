ifndef CONFIG_INC


############################################################
########## OPTIONS MOST LIKELY NEED TO CUSTOMIZED ##########
############################################################


# reference genome
# e.g.: hg19, b37, hg19_ionref, b37_hbv_hcv etc
REF ?= b37

# possible values: ILLUMINA, IONTORRENT
SEQ_PLATFORM ?= ILLUMINA

# possible values: NONE, BAITS, PCR, RNA, CHIP
CAPTURE_METHOD ?= NONE

# e.g: NONE, AGILENT_CLINICAL_EXOME, CCP, AGILENT_ALLEXON_MOUSE, HCC, WXS
# For WGS, RNA-seq and ChiP-seq, this should be NONE.
PANEL ?= NONE

# Single-end or paired-end, set to false if single-end
PAIRED_END ?= true

# for pair-end RNA-seq, [NONE|FIRST_READ_TRANSCRIPTION_STRAND|SECOND_READ_TRANSCRIPTION_STRAND]
STRAND_SPECIFICITY ?= NONE

# possible values: SOMATIC, GERMLINE
ANALYSIS_TYPE ?= SOMATIC

# if you anticipate a complex project requiring analyzing subsets of samples
# this is useful for distinguishing them
PROJECT_PREFIX ?=

INCLUDE_CHR_Y ?= true

# This parameter determines which aligned BAM files get put into the bam/ directory
# This is useful when multiple aligners are used and you want to specify the
# set that goes into bam/, since most downstream tools get the data from bam/
ifeq ($(strip $(CAPTURE_METHOD)),RNA)
PRIMARY_ALIGNER ?= star
else
PRIMARY_ALIGNER ?= bwamem
endif

# Specify the sample sheets
SAMPLE_SET_FILE ?= sample_sets.txt
SAMPLE_FILE ?= samples.txt
SAMPLE_SPLIT_FILE ?= samples.split.txt

# Normal samples to use for pooled_norm and panel_of_normals: 
# particularly useful for when samples are being added to specify a fixed set of normals to use.
# Otherwise, if one single new normal sample is added, everything re-runs.
SAMPLE_PON_FILE ?= samples.pon.txt
SAMPLE_POOLEDNORM_FILE ?= samples.poolednorm.txt


MERGE_SAMPLE_FILE ?= merge_samples.txt

############################################################
########### END MOST FREQUENTLY CUSTOMIZED PARAMS ##########
############################################################

# this grabs the parameters for the genome and target panels
include usb-modules-v2/genome_inc/$(REF).inc


#############################################
############### Pre-alignment ###############
#############################################

## These parameters are used when you align FASTQs to BAMs
## If MERGE_SPLIT_FASTQ=true, this happens before FASTQ_FILTER
## If set, <FASTQ_FILTER> is run before alignment

# If multiple sets of fastqs per sample, should they be merged before alignment
# For transcriptome, true. For genomic, false.
# There are performance advantages to setting this to false.
ifeq ($(CAPTURE_METHOD),RNA)
MERGE_SPLIT_FASTQ ?= true # merge multiple sets of fastqs from a given sample before alignment
endif
MERGE_SPLIT_FASTQ ?= false

# Pre-processing the FASTQs, trimming bad cycles or adaptors
# [<NULL>|trim|trimgalore|cutadapt]
# trim: uniforming trim reads down to TRIM_LENGTH (see below)
# trimgalore: trims adaptors, and optionally clip reads to CLIP_FASTQ_R1/2 base pairs
# cutadapt: trims adaptors, but more flexible options than trimgalore
ifeq ($(CAPTURE_METHOD),RNA)
FASTQ_FILTER ?= cutadapt # for fastq trimming. RNA-seq needs the adaptors trimmed off
endif
FASTQ_FILTER ?=

# Uniforming trim to this length if FASTQ_FILTER=trim
TRIM_LENGTH ?= 150

# setting the --clip_R1 and --clip_R2 parameters if FASTQ_FILTER=trimgalore
CLIP_FASTQ_R1 ?=
CLIP_FASTQ_R2 ?=

# options for FASTQ_FILTER=cutadapt
ADAPTOR_POS ?= 3prime # [3prime|5prime]
ADAPTOR_SEQ ?= AGATCGGAAGAGC
NEXTSEQ_TRIM ?= false # if sequencing was done with NextSeq, this may help.


#############################################
############### Post-alignment ###############
#############################################

# After alignment with the favourite aligner, 
# the BAM files are post-processed in the following order
# For RNA-seq, if mutation calls will eventually be needed, 
# then post-process, otherwise just dont, it takes a long time and creates hard-clipped
# reads that may be problematic. Use true with caution
POST_PROCESS_RNA_BAM ?= false [true|false]

# this runs Picard ReorderSam
BAM_REORDER ?= true 

# this runs Picard SortSam
BAM_SORT ?= true 

# this runs Picard AddOrReplaceReadGroups
BAM_FIX_RG ?= false 

# this filters out reads with BAM_FILTER_FLAGS
BAM_FILTER ?= true 
BAM_FILTER_FLAGS ?= 768

# this removes inter-chromosomal chimeric paired reads
BAM_INTRACHR ?= false 

# this runs GATK SplitNCigarReads, default to false except for RNA-seq
ifeq ($(CAPTURE_METHOD),RNA)
ifeq ($(POST_PROCESS_RNA_BAM),true)
BAM_SPLIT_N_TRIM ?= true 
endif
endif
BAM_SPLIT_N_TRIM ?= false 

# this splits bams into individual chroms betfore BAM_DUP_TYPE, 
# BAM_REALN and BAM_RECAL before remerging the BAMs.
# defaults to true except for RNA-seq, amplicon-based capture and CHIP
# setting this to false means BAM_DUP_TYPE, BAM_REALN and BAM_RECAL are done
# on a genome-wide level
ifeq ($(CAPTURE_METHOD),RNA)
SPLIT_CHR ?= false
endif
ifeq ($(CAPTURE_METHOD),PCR)
SPLIT_CHR ?= false
endif
ifeq ($(CAPTURE_METHOD),CHIP)
SPLIT_CHR ?= false
endif
SPLIT_CHR ?= true 


# This runs deduplication, GATK IndelRealigner and BaseRecalibrator
# deduplication: [none|markdup|rmdup]: markdup is Picard MarkDuplicates, rmdup is samtools rmdup
# defaults to markdup and true except for RNA-seq, amplicon-based capture and CHIP
ifeq ($(CAPTURE_METHOD),RNA)
BAM_DUP_TYPE ?= none
BAM_REALN ?= false
BAM_RECAL ?= false
endif
ifeq ($(CAPTURE_METHOD),PCR)
BAM_DUP_TYPE ?= none
BAM_REALN ?= false
BAM_RECAL ?= false
endif
ifeq ($(CAPTURE_METHOD),CHIP)
BAM_DUP_TYPE ?= none
BAM_REALN ?= false
BAM_RECAL ?= false
endif
BAM_DUP_TYPE ?= markdup 
BAM_REALN ?= true 
BAM_RECAL ?= true 
BAM_CHR2_BASE_RECAL ?= true

#########################################################
############### Additional BAM processing ###############
#########################################################

# if you have a set of aligned BAM files you want to put
# through (any step) of the post-alignment sequence above
# set this to true, and check and set the post-alignment
# parameters in the previous section accordingly, then
# put the bam files in the unprocessed_bam direcory
BAM_REPROCESS ?= false

BAM_PHRED64 ?= false # obsolete parameter
MERGE_SPLIT_BAMS ?= false  # internal parameter

# if you need to downsample BAM files, for example to
# make a pooled normal BAM file, then you set this.
# see the samtools view -s option
SAMTOOLS_DOWNSAMPLE_FACTOR ?= 1001.2

# For QC, genotyping is usually done. For performance reasons,
# we usually use chromosome 2 for genotyping (a large chromosome
# that is not usually gained/lost in cancer. For very small panel,
# this is not good. This is also specified in some genome_inc panel
# config files. Values in genome_inc configs overide the values here.
ifeq ($(CAPTURE_METHOD),PCR)
GENOTYPE_WITH_CHR2 = ?= false
endif
ifeq ($(CAPTURE_METHOD),CHIP)
GENOTYPE_WITH_CHR2 = ?= false
endif
GENOTYPE_WITH_CHR2 = ?= true



##################################
#### General MPS tools ###########
##################################

EXTRACT_TOOL ?= PICARD # obsolete

##################################
########## Aligners ##############
##################################

HISAT2_OPTS = --dta 

################################################
###### Somatic variant callers #################
################################################

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_GT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.5

################################################
###### Somatic variant filters #################
################################################

# after somatic variant callers, the calls go through
# a series of filters. If you don't want the filters
# set this to false
FILTER_VARIANTS ?= true

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
MIN_NORMAL_DEPTH ?= 5
MIN_TUMOR_DEPTH ?= 10
MIN_TUMOR_AD ?= 3
MIN_TN_AD_RATIO ?= 5.0
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
MIN_NORMAL_DEPTH ?= 10
MIN_TUMOR_DEPTH ?= 10
MIN_TUMOR_AD ?= 8
MIN_TN_AD_RATIO ?= 10.0
endif

MIN_MQ ?= 10
MIN_AF_SNP ?= 0.01
MIN_AF_INDEL ?= 0.02

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
USE_NFT ?= true
PON_VCF ?= mutect2/pon.mutect2.vcf
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
USE_NFT ?= false
PON_VCF ?= tvc/pon.tvc.vcf
endif


#################################################
###### Germline variant callers #################
#################################################

GATK_INTERVALS_COUNT = 200
GATK_HARD_FILTER_SNPS ?= false
GATK_HARD_FILTER_INDELS ?= false

VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL_SNPS = 99.5
VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL_INDELS = 99.0

## These are used if GATK_HARD_FILTER_SNPS ##
HAPCALL_SNP_MQ_THRESHOLD ?= 40.0
HAPCALL_SNP_QD_THRESHOLD ?= 2.0
HAPCALL_SNP_FS_THRESHOLD ?= 60.0
HAPCALL_SNP_HAP_SCORE_THRESHOLD ?= 13.0
HAPCALL_SNP_MQ_RANKSUM_THRESHOLD ?= -12.5
HAPCALL_SNP_READPOS_RANKSUM_THRESHOLD ?= -8.0

## These are used if GATK_HARD_FILTER_INDELS ##
HAPCALL_INDEL_INBREED_COEFF_THRESHOLD ?= -0.8
HAPCALL_INDEL_QD_THRESHOLD ?= 2.0
HAPCALL_INDEL_FS_THRESHOLD ?= 200.0
HAPCALL_INDEL_HAP_SCORE_THRESHOLD ?= 13.0
HAPCALL_INDEL_MQ_RANKSUM_THRESHOLD ?= -20.0
HAPCALL_INDEL_READPOS_RANKSUM_THRESHOLD ?= -8.0



##########################################
########## Variant Annotation ############
##########################################

ANNOTATE_VARIANTS ?= true
USE_SUFAM ?= true
ANN_FACETS ?= false

MUTATION_SUMMARY_FORMAT ?= EXCEL
INCLUDE_LNCRNA_IN_SUMMARY ?= false

SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-intergenic -noStats

EXAC_INFO_FIELDS ?= ExACnontcga_AC,ExACnontcga_AF
EXACNONPSYCH_INFO_FIELDS ?= ALL

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
CALLER_PREFIX ?= mutect strelka_indels
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
CALLER_PREFIX ?= tvc_snps tvc_indels
endif


##################################
########## Copy number ###########
##################################

FACETS_GATK_VARIANTS ?= false
FACETS_MINGC ?= 0
FACETS_MAXGC ?= 1

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
FACETS_SNP_PILEUP_MIN_DEPTH ?= 25
FACETS_SNP_PILEUP_MAX_DEPTH ?= 1000
FACETS_WINDOW_SIZE ?= 200
FACETS_CVAL1 ?= 200
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
FACETS_SNP_PILEUP_MIN_DEPTH ?= 50
FACETS_SNP_PILEUP_MAX_DEPTH ?= 5000
FACETS_WINDOW_SIZE ?= 10
FACETS_CVAL1 ?= 75
endif

FACETS_SNP_PILEUP_MINMAPQ ?= 1
FACETS_SNP_PILEUP_MINBASEQ ?= 13

FACETS_PRE_CVAL ?=
FACETS_CVAL2 ?=
FACETS_MIN_NHET ?= 15

FACETS_SPECIAL_CASES ?=
FACETS_PLOT_GENE_CN_OPTS = --sampleColumnPostFix '_EM'

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
CBS_SEG_SD ?= 2
CBS_SEG_SMOOTH ?= 10
CBS_SEG_ALPHA ?= 0.000001
CBS_TRIM ?= 0.025
CBS_CLEN ?= 10
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
CBS_SEG_SD ?= 1.5
CBS_SEG_SMOOTH ?= 15
CBS_SEG_ALPHA ?= 0.0000000001
CBS_TRIM ?= 0.1
CBS_CLEN ?= 5
endif

CBS_EXCL_N_OUTLIER_PC ?= 0.05
CBS_MIN_N_DEPTH ?= 0

CBS_MULTIPARAM_SEGMENT ?= false
CBS_SEG_SDS ?= 1.5 2 2.5 3
CBS_SEG_SMOOTHS ?= 5 10
CBS_SEG_ALPHAS ?= 0.01 0.0001 0.000001 0.0000000001

VARSCAN_GENE_CN_OPTS = $(FACETS_GENE_CN_OPTS)

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
VARSCAN_CNV_MAX_SEG_SIZE = 100
CBS_MPILEUP_OPTS = -q 1
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
VARSCAN_CNV_MAX_SEG_SIZE = 1000
CBS_MPILEUP_OPTS = -q 1 -d 20000
endif

##################################
############# GISTIC ############
##################################

GISTIC_THRESHOLD ?= 0.3
GISTIC_JS ?= 15

##################################
########## Gene expression #######
##################################

DESEQ_CONDITION ?= condition
DESEQ_REF_CONDITION ?= ref
DESEQ_ALT_CONDITION ?= alt
DESEQ_PHENO_FILE ?= pheno.txt


#################################
####### ChIP-seq ################
#################################

MOSAICS_FRAG_LEN ?= 200
MOSAICS_BIN_SIZE ?= 200
MOSAICS_NUM_CORES ?= 8
MOSAICS_PARALLEL ?= TRUE
MOSAICS_MAXGAP ?= 200
MOSAICS_MINSIZE ?= 50
MOSAICS_THRES ?= 10

##################################
########## Fusions ###############
##################################

STAR_FUSION_MAX_PROMISCUITY ?= 3
STAR_FUSION_MIN_JUNC_SUPP ?= 2
STAR_FUSION_MIN_ALT_PCT_JUNC ?= 5

##################################
####### CLONALITY ################
##################################

ABSOLUTE_SIGMA_P ?= 0
ABSOLUTE_MAX_SIGMA_H ?= 0.07
ABSOLUTE_MIN_PLOIDY ?= 0.9
ABSOLUTE_MAX_PLOIDY ?= 8
ABSOLUTE_DISEASE ?= BRCA
ABSOLUTE_PLATFORM ?= Illumina_WES
ABSOLUTE_MAX_SEG ?= 3500
ABSOLUTE_COPYNUMTYPE ?= total
ABSOLUTE_MAX_NEG_GENOME ?= 0.05
ABSOLUTE_MAX_NON_CLONAL ?= 0.2

PYCLONE_PRIOR ?= parental_copy_number

################################
###### Miscellaneous ###########
###############################

DECONSTRUCTSIGS_NUMITER ?= 100

endif
CONFIG_INC = true
