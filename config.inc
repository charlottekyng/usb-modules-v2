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
# possible values: NONE, AGILENT_CLINICAL_EXOME, IMPACTv5, CCP, CHPv2, IMPACT410, IMPACT341, AGILENT_ALLEXON_MOUSE, HCC, NOVARTIS_PC2, SEQCAP_EXOME_V3
PANEL ?= NONE
# Single-end or paired-end, set to false if single-end
PAIRED_END ?= true
# for RNA-seq, possible values NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND
STRAND_SPECIFICITY ?= NONE
# also for RNA-seq, if mutation calls will eventually be needed, 
# then post-process, otherwise just dont, it takes a long time and creates hard-clipped
# reads that may be problematic
POST_PROCESS_RNA_BAM ?= false


# possible values: CANCER, PSYCH
ANALYSIS_TYPE ?= CANCER 

PROJECT_PREFIX ?=

# Normal samples to use for pooled_norm and panel_of_normals: 
#particularly useful for when samples are being added...
PANEL_OF_NORMAL_SAMPLES ?= $(NORMAL_SAMPLES)
POOLED_NORM_SAMPLES ?= $(NORMAL_SAMPLES)

INCLUDE_CHR_Y ?= true

ifeq ($(strip $(CAPTURE_METHOD)),RNA)
PRIMARY_ALIGNER ?= star
else
PRIMARY_ALIGNER ?= bwamem
endif

############################################################
########### END MOST FREQUENTLY CUSTOMIZED PARAMS ##########
############################################################

# default to medium for exomes, other values are "low" for targeted panels and "high" for WGS
# possible values: low, medium, high, vhigh
include usb-modules-v2/genome_inc/$(REF).inc


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

FASTQ_FILTER ?=
MERGE_SPLIT_FASTQ ?= false 

BAM_REPROCESS ?= false
BAM_PHRED64 ?= false
MERGE_SPLIT_BAMS ?= false  # merge processed split bams
BAM_REORDER ?= true
BAM_SORT ?= true
BAM_FIX_RG ?= false # fix the readgroups so that it corresponds with the bam filename
BAM_FILTER ?= true
BAM_SPLIT_N_TRIM ?= false
SPLIT_CHR ?= true
BAM_DUP_TYPE ?= markdup
BAM_REALN ?= true
BAM_RECAL ?= true
BAM_CHR2_BASE_RECAL ?= true
GENOTYPE_WITH_CHR2 = ?= true
BAM_FILTER_FLAGS ?= 768
SAMTOOLS_DOWNSAMPLE_FACTOR ?= 1001.2




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

HRUN ?= false
HRUN_FILTER ?= 4

USE_SUFAM ?= true
ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
USE_NFT ?= true
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
USE_NFT ?= false
endif

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
PON_VCF ?= mutect2/pon.mutect2.vcf
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
PON_VCF ?= tvc/pon.tvc.vcf
endif

##################################
#### General MPS tools ###########
##################################

EXTRACT_TOOL ?= PICARD

##################################
########## Aligners ##############
##################################

BWA_NUM_CORES ?= 8
BWA_MEM_OPTS ?= -M
BWA_ALN_OPTS ?=

#TMAP_MODE ?= map3
#TMAP_OPTS =

TOPHAT_NUM_CORES ?= 4
NO_NOVEL_SPLICING ?= false
TOPHAT_OPTS = --keep-fasta-order --no-sort-bam -G $(GENES_GTF) -p ${TOPHAT_NUM_CORES} --tmp-dir $(TMPDIR)/$*
ifeq ($(BAM_PHRED64),true)
TOPHAT_OPTS += --solexa1.3-quals
endif
ifeq ($(NO_NOVEL_SPLICING),true)
TOPHAT_OPTS += --no-novel-juncs
endif

STAR_FUSION_MAX_PROMISCUITY ?= 3
STAR_FUSION_MIN_JUNC_SUPP ?= 2
STAR_FUSION_MIN_ALT_PCT_JUNC ?= 5

HISAT2_NUM_CORES ?= 8
HISAT2_OPTS = -x $(HISAT2_REF) --dta 
ifeq ($(STRAND_SPECIFICITY),FIRST_READ_TRANSCRIPTION_STRAND)
HISAT2_OPTS += --rna-strandness FR
endif
ifeq ($(STRAND_SPECIFICITY),SECOND_READ_TRANSCRIPTION_STRAND)
HISAT2_OPTS += --rna-strandness RF
endif
#ifeq ($(PHRED64),true)
#        HISAT_OPTS += --phred64
#endif

STRINGTIE_NUM_CORES ?= 4


##################################
########## Picard/GATK ###########
##################################


BAM_BASE_RECAL_OPTS = -knownSites $(DBSNP) \
	$(if $(findstring true,$(BAM_CHR2_BASE_RECAL)),-L $(word 2,$(CHROMOSOMES))) \
	$(if $(TARGETS_FILE_INTERVALS),-L $(TARGETS_FILE_INTERVALS))

ifneq ($(KNOWN_INDELS),)
BAM_REALN_OPTS = --knownAlleles $(KNOWN_INDELS)
BAM_REALN_TARGET_OPTS = --known $(KNOWN_INDELS)
endif

COVARIATES = -cov ReadGroupCovariate -cov QualityScoreCovariate -cov DinucCovariate -cov CycleCovariate

REPORT_STRATIFICATION := Filter

##################################
########## Variant callers #######
##################################

ANNOTATE_VARIANTS ?= true

MUTECT_MAX_ALT_IN_NORMAL ?= 500
MUTECT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.05
MUTECT_GT_MAX_ALT_IN_NORMAL_FRACTION ?= 0.5

STRELKA_NUM_CORES ?= 8

TVC_NUM_CORES ?= 4

VARSCAN_OPTS=--output-vcf 1 --mpileup 1 --min-var-freq $(MIN_AF_SNP)

GATK_HARD_FILTER_SNPS ?= true
GATK_POOL_SNP_RECAL ?= false

VARIANT_CALL_THRESHOLD = 30
HAPLOTYPE_CALLER_OPTS = --dbsnp $(DBSNP) -stand_call_conf $(VARIANT_CALL_THRESHOLD) -R $(REF_FASTA)

VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL = 99.0
VARIANT_RECAL_ANNOTATIONS = QD MQRankSum FS MQ ReadPosRankSum

## I don't think these are being used ##
HAPCALL_SNP_MQ_THRESHOLD ?= 40.0
HAPCALL_SNP_QD_THRESHOLD ?= 2.0
HAPCALL_SNP_FS_THRESHOLD ?= 60.0
HAPCALL_SNP_HAP_SCORE_THRESHOLD ?= 13.0
HAPCALL_SNP_MQ_RANKSUM_THRESHOLD ?= -12.5
HAPCALL_SNP_READPOS_RANKSUM_THRESHOLD ?= -8.0
HAPCALL_SNP_FILTERS := --filterName 'QD' --filterExpression 'QD < $(HAPCALL_SNP_QD_THRESHOLD)' \
--filterName 'MQ' --filterExpression 'MQ < $(HAPCALL_SNP_MQ_THRESHOLD)' \
--filterName 'FS' --filterExpression 'FS > $(HAPCALL_SNP_FS_THRESHOLD)' \
--filterName 'HapScore' --filterExpression 'HaplotypeScore > $(HAPCALL_SNP_HAP_SCORE_THRESHOLD)' \
--filterName 'MQRankSum' --filterExpression 'MQRankSum < $(HAPCALL_SNP_MQ_RANKSUM_THRESHOLD)' \
--filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < $(HAPCALL_SNP_READPOS_RANKSUM_THRESHOLD)' \
--filterName 'Depth' --filterExpression 'DP < 5'

HAPCALL_INDEL_INBREED_COEFF_THRESHOLD ?= -0.8
HAPCALL_INDEL_QD_THRESHOLD ?= 2.0
HAPCALL_INDEL_FS_THRESHOLD ?= 200.0
HAPCALL_INDEL_HAP_SCORE_THRESHOLD ?= 13.0
HAPCALL_INDEL_MQ_RANKSUM_THRESHOLD ?= -20.0
HAPCALL_INDEL_READPOS_RANKSUM_THRESHOLD ?= -8.0

INDEL_FILTERS = --filterName 'QD' --filterExpression 'QD < $(HAPCALL_INDEL_QD_THRESHOLD)' \
--filterName 'ReadPosRankSum' --filterExpression 'ReadPosRankSum < ($HAPCALL_INDEL_MQ_RANKSUM_THRESHOLD)' \
--filterName 'InbreedingCoeff' --filterExpression 'InbreedingCoeff < $(HAPCALL_INDEL_INBREED_COEFF_THRESHOLD)'  \
--filterName 'FS' --filterExpression 'FS > $(HAPCALL_INDEL_FS_THRESHOLD)' \
--filterName 'DP' --filterExpression 'DP < 5'


##################################
########## Annotation ############
##################################

MUTATION_SUMMARY_FORMAT ?= EXCEL
INCLUDE_LNCRNA_IN_SUMMARY ?= false


SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-intergenic

EXAC_INFO_FIELDS ?= AC_AFR,AC_AMR,AC_Adj,AC_EAS,AC_FIN,AC_Hemi,AC_Het,AC_Hom,AC_NFE,AC_OTH,AC_SAS,AN_AFR,AN_AMR,AN_Adj,AN_EAS,AN_FIN,AN_NFE,AN_OTH,AN_SAS,Hemi_AFR,Hemi_AMR,Hemi_EAS,Hemi_FIN,Hemi_NFE,Hemi_OTH,Hemi_SAS,Het_AFR,Het_AMR,Het_EAS,Het_FIN,Het_NFE,Het_OTH,Het_SAS,Hom_AFR,Hom_AMR,Hom_EAS,Hom_FIN,Hom_NFE,Hom_OTH,Hom_SAS,AC_MALE,AC_FEMALE,AN_MALE,AN_FEMALE
EXACNONPSYCH_INFO_FIELDS ?= ALL
#EXACNONPSYCH_INFO_FIELDS ?= ExACnonpsychAC_AFR,ExACnonpsychAC_AMR,ExACnonpsychAC_Adj,ExACnonpsychAC_EAS,ExACnonpsychAC_FIN,ExACnonpsychAC_Hemi,ExACnonpsychAC_Het,ExACnonpsychAC_Hom,ExACnonpsychAC_NFE,ExACnonpsychAC_OTH,ExACnonpsychAC_SAS,ExACnonpsychAN_AFR,ExACnonpsychAN_AMR,ExACnonpsychAN_Adj,ExACnonpsychAN_EAS,ExACnonpsychAN_FIN,ExACnonpsychAN_NFE,ExACnonpsychAN_OTH,ExACnonpsychAN_SAS,ExACnonpsychHemi_AFR,ExACnonpsychHemi_AMR,ExACnonpsychHemi_EAS,ExACnonpsychHemi_FIN,ExACnonpsychHemi_NFE,ExACnonpsychHemi_OTH,ExACnonpsychHemi_SAS,ExACnonpsychHet_AFR,ExACnonpsychHet_AMR,ExACnonpsychHet_EAS,ExACnonpsychHet_FIN,ExACnonpsychHet_NFE,ExACnonpsychHet_OTH,ExACnonpsychHet_SAS,ExACnonpsychHom_AFR,ExACnonpsychHom_AMR,ExACnonpsychHom_EAS,ExACnonpsychHom_FIN,ExACnonpsychHom_NFE,ExACnonpsychHom_OTH,ExACnonpsychHom_SAS,ExACnonpsychAC_MALE,ExACnonpsychAC_FEMALE,ExACnonpsychAN_MALE,ExACnonpsychAN_FEMALE
VCF_FIELDS = CHROM POS ID REF ALT FILTER QUAL
#NSFP_FIELDS ?= ExAC_Adj_AF 1000Gp3_AF ESP6500_AA_AF ESP6500_EA_AF \
#	clinvar_rs clinvar_clnsig \
#	Uniprot_acc_Polyphen2 Uniprot_id_Polyphen2 Polyphen2_HVAR_pred \
#	MutationAssessor_pred MutationTaster_pred \
#	PROVEAN_pred Interpro_domain
ANN_FIELDS = $(addprefix ANN[*].,ALLELE EFFECT IMPACT GENE GENEID FEATURE FEATUREID BIOTYPE RANK \
	HGVS_C HGVS_P CDNA_POS CDNA_LEN CDS_POS CDS_LEN AA_POS AA_LEN DISTANCE ERRORS)

#ANNOVAR = $(PERL) $(HOME)/share/usr/annovar/table_annovar.pl
#ANNOVAR_PROTOCOL ?= refGene$(,)cytoBand$(,)genomicSuperDups$(,)esp6500siv2_all$(,)1000g2014oct_all$(,)1000g2014oct_afr$(,)1000g2014oct_eas$(,)1000g2014oct_eur$(,)snp138$(,)ljb26_all
#ANNOVAR_OPERATION ?= g$(,)r$(,)r$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f
#ANNOVAR_OPTS = --dot2underline -remove -protocol $(ANNOVAR_PROTOCOL) -operation $(ANNOVAR_OPERATION) -nastring . -vcfinput -buildver $(ANNOVAR_REF)


#CHASM = $(RSCRIPT) modules/vcf_tools/chasmVcf.R
#CHASM_DIR = /ifs/opt/common/CHASM/CHASMDL.1.0.7
#CHASM_DIR = $(HOME)/share/usr/CHASM
#CHASM_PYTHON_ENV = $(HOME)/share/usr/anaconda-envs/pyenv27-chasm
#CHASM_CLASSIFIER ?= Breast

#FATHMM = $(MY_RSCRIPT) modules/vcf_tools/fathmmVcf.R
#FATHMM_DIR = $(HOME)/share/usr/fathmm
#FATHMM_PYTHON = $(HOME)/share/usr/bin/python
#FATHMM_PYTHONPATH = $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7
#FATHMM_OPTS = --genome $(REF) --ensemblTxdb $(ENSEMBL_TXDB) --ref $(REF_FASTA) --fathmmDir $(FATHMM_DIR) --python $(FATHMM_PYTHON) \
#	--mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) \
#	--mysqlPassword $(EMBL_MYSQLDB_PW) --mysqlDb (EMBL_MYSQLDB_DB)
#TRANSFIC = $(RSCRIPT) modules/vcf_tools/transficVcf.R
#TRANSFIC_PERL_SCRIPT = $(HOME)/share/usr/transfic/bin/transf_scores.pl

#PROVEAN = $(RSCRIPT) modules/vcf_tools/proveanVcf.R
#PROVEAN_OPTS = --genome $(REF) --aaTable $(AA_TABLE) --ensemblTxdb $(ENSEMBL_TXDB) --mysqlHost $(EMBL_MYSQLDB_HOST) \
#			   --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) \
#			   --mysqlDb $(EMBL_MYSQLDB_DB) --numThreads 8 --memPerThread 1G --queue $(QUEUE) --qsubPriority $(QSUB_PRIORITY)

#AA_TABLE = $(HOME)/share/reference/aa_table.tsv

#MUT_ASS = $(RSCRIPT) usb-modules-v2/vcf_tools/mutAssVcf.R
#MUTATION_TASTER = $(PYTHON) usb-modules-v2/vcf_tools/mutation_taster_vcf.py

#NON_SILENT_EFF = START_GAINED SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START
#NON_SILENT_CODING_EFF = START_GAINED START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START

#EMBL_MYSQLDB_HOST ?= 10.0.200.71
#EMBL_MYSQLDB_PORT ?= 38493
#EMBL_MYSQLDB_USER ?= embl
#EMBL_MYSQLDB_PW ?= embl
#EMBL_MYSQLDB_DB ?= homo_sapiens_core_75_37


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
FACETS_CVAL1 ?= 75 100 200
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
FACETS_SNP_PILEUP_MIN_DEPTH ?= 50
FACETS_SNP_PILEUP_MAX_DEPTH ?= 5000
FACETS_WINDOW_SIZE ?= 10
FACETS_CVAL1 ?= 75 100
endif

FACETS_SNP_PILEUP_MINMAPQ ?= 1
FACETS_SNP_PILEUP_MINBASEQ ?= 13

FACETS_PRE_CVAL ?=
FACETS_CVAL2 ?=
FACETS_MIN_NHET ?= 15

FACETS_SPECIAL_CASES ?=
FACETS_PLOT_GENE_CN_OPTS = --sampleColumnPostFix '_EM'

ANN_FACETS ?= false

CBS_SEG_SD ?= 2
CBS_SEG_SMOOTH ?= 10
CBS_SEG_ALPHA ?= 0.000001
CBS_TRIM ?= 0.025
CBS_CLEN ?= 10

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
# $(if $(TARGETS_FILE_INTERVALS),-l $(TARGETS_FILE_INTERVALS))

##################################
############# GISTIC ############

GISTIC_THRESHOLD ?= 0.3
GISTIC_JS ?= 15
GISTIC_OPTS = -genegistic 0 -smallmem 1 -maxseg 5000 -savegene 1 -saveseg 1 -savedata 0 -v 30
GISTIC_OPTS += -qvt 0.25 -conf 0.99 -broad 1 -brlen 0.5 -rx 0
GISTIC_OPTS += -ta $(GISTIC_THRESHOLD) -td $(GISTIC_THRESHOLD) -js $(GISTIC_JS)

##################################
########## Gene expression #######
##################################

DESEQ_CONDITION ?= condition
DESEQ_REF_CONDITION ?= ref
DESEQ_ALT_CONDITION ?= alt
DESEQ_PHENO_FILE ?= pheno.txt

RSEM_NUM_CORES ?= 4
RSEM_OPTIONS ?= --alignments --no-bam-output -p $(RSEM_NUM_CORES) --forward-prob 0

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

#DEFUSE ?= $(PERL) $(HOME)/share/usr/defuse-0.6.1/scripts/defuse.pl
#INTEGRATE ?= $(HOME)/share/usr/Integrate_0.2.0/bin/Integrate
#ONCOFUSE_JAR ?= $(JARDIR)/Oncofuse.jar
endif
CONFIG_INC = true
