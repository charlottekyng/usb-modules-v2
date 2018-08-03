ifndef CONFIG_INC


############################################################
########## OPTIONS MOST LIKELY NEED TO CUSTOMIZED ##########
############################################################


# default reference:
# possible values: mm10, hg19, b37, hg19_ionref, b37_hbv_hcv
REF ?= b37

# possible values: ILLUMINA, IONTORRENT
SEQ_PLATFORM ?= ILLUMINA

# possible values: NONE, BAITS, PCR, RNA, CHIP
CAPTURE_METHOD ?= NONE

# possible values: NONE, AGILENT_CLINICAL_EXOME, CCP, AGILENT_ALLEXON_MOUSE, HCC, WXS
# For WGS, RNA-seq and ChiP-seq, this should be NONE.
PANEL ?= NONE

# Single-end or paired-end, set to false if single-end
PAIRED_END ?= true

# for RNA-seq, possible values NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND
STRAND_SPECIFICITY ?= NONE

# also for RNA-seq, if mutation calls will eventually be needed, 
# then post-process, otherwise just dont, it takes a long time and creates hard-clipped
# reads that may be problematic
POST_PROCESS_RNA_BAM ?= false

# possible values: SOMATIC, GERMLINE
ANALYSIS_TYPE ?= SOMATIC

# if you anticipate a complex project requiring analyzing subsets of samples
# this is useful for distinguishing them
PROJECT_PREFIX ?=

# Normal samples to use for pooled_norm and panel_of_normals: 
# particularly useful for when samples are being added to specify a fixed set of normals to use.
# Otherwise, if one single new normal sample is added, everything re-runs.
PANEL_OF_NORMAL_SAMPLES ?= $(NORMAL_SAMPLES)
POOLED_NORM_SAMPLES ?= $(NORMAL_SAMPLES)

INCLUDE_CHR_Y ?= true

# This parameter determines which aligned BAM files get put into the bam/ directory
# This is useful when multiple aligners are used and you want to specify the
# set that goes into bam/, since most downstream tools get the data from bam/
ifeq ($(strip $(CAPTURE_METHOD)),RNA)
PRIMARY_ALIGNER ?= star
else
PRIMARY_ALIGNER ?= bwamem
endif

############################################################
########### END MOST FREQUENTLY CUSTOMIZED PARAMS ##########
############################################################

include usb-modules-v2/genome_inc/$(REF).inc

#######################################
########### FASTQ/BAM PARAMS ##########
#######################################


## process FASTQs before alignment
TRIM_LENGTH ?= 150
CLIP_FASTQ_R1 ?=
CLIP_FASTQ_R2 ?=

ADAPTOR_POS ?= 3prime
ADAPTOR_SEQ ?= AGATCGGAAGAGC
NEXTSEQ_TRIM ?= false

## default FASTQ and BAM processing parameters based on application

ifeq ($(CAPTURE_METHOD),RNA)
FASTQ_FILTER ?= cutadapt
MERGE_SPLIT_FASTQ ?= true
BAM_DUP_TYPE ?= none
ifeq ($(POST_PROCESS_RNA_BAM),true)
BAM_SPLIT_N_TRIM ?= true
else
BAM_REALN ?= false
BAM_RECAL ?= false
SPLIT_CHR ?= false
endif
endif

ifeq ($(CAPTURE_METHOD),PCR)
BAM_DUP_TYPE ?= none
BAM_REALN ?= false
BAM_RECAL ?= false
BAM_CHR2_BASE_RECAL ?= false
GENOTYPE_WITH_CHR2 = ?= false
endif

ifeq ($(CAPTURE_METHOD),CHIP)
SPLIT_CHR ?= false
BAM_DUP_TYPE ?= none
BAM_REALN ?= none
BAM_RECAL ?= none
GENOTYPE_WITH_CHR2 = ?= false
endif

## default parameters if not set above
FASTQ_FILTER ?=
MERGE_SPLIT_FASTQ ?= false 
BAM_REPROCESS ?= false
BAM_PHRED64 ?= false
MERGE_SPLIT_BAMS ?= false  # merge processed split bams
BAM_REORDER ?= true
BAM_SORT ?= true
BAM_FIX_RG ?= false # fix the readgroups so that it corresponds with the bam filename
BAM_FILTER ?= true
BAM_INTRACHR ?= false
BAM_SPLIT_N_TRIM ?= false
SPLIT_CHR ?= true
BAM_DUP_TYPE ?= markdup
BAM_REALN ?= true
BAM_RECAL ?= true
BAM_CHR2_BASE_RECAL ?= true
GENOTYPE_WITH_CHR2 = ?= true
BAM_FILTER_FLAGS ?= 768
SAMTOOLS_DOWNSAMPLE_FACTOR ?= 1001.2


##################################
#### General MPS tools ###########
##################################

EXTRACT_TOOL ?= PICARD

##################################
########## Aligners ##############
##################################

HISAT2_OPTS = --dta 

################################################
########## Variant callers and filtering #######
################################################

FILTER_VARIANTS ?= true

### VARIANT FILTERS
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

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_GT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.5


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


ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
USE_NFT ?= true
PON_VCF ?= mutect2/pon.mutect2.vcf
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
USE_NFT ?= false
PON_VCF ?= tvc/pon.tvc.vcf
endif


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
