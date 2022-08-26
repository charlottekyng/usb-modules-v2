include usb-modules-v2/config.inc

############################################
######## DEFINE SAMPLES/SETS ###############
############################################


ifndef MAKEFILE_INC

PROJECT_NAME = $(shell pwd | sed 's:.*/projects/::; s:.*/data/::; s:.*kinglab/::; s:/:_:g')

ifneq ($(wildcard $(SAMPLE_FILE)),)
  SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_FILE))
endif
ifneq ($(wildcard $(SAMPLE_PON_FILE)),)
  PANEL_OF_NORMAL_SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_PON_FILE))
endif
ifneq ($(wildcard $(SAMPLE_POOLEDNORM_FILE)),)
  POOLED_NORM_SAMPLES ?= $(shell sed '/^\#/d' $(SAMPLE_POOLEDNORM_FILE))
endif
ifneq ($(wildcard $(MERGE_SAMPLE_FILE)),)
  MERGE_SAMPLES ?= $(shell sed '/^\#/d' $(MERGE_SAMPLE_FILE))
endif


get_tumors = $(wordlist 1,$(shell expr $(words $1) - 1),$1)
get_normal = $(lastword $1)
get_space_sets = $(shell sed '/^\#/d; s/\s\+/ /g; s/\s\+$$//;' $(SAMPLE_SET_FILE) | sed -n '$(1)p')
get_underscore_sets = $(shell sed '/^\#/d; s/\s\+/_/g; s/\s\+$$//;' $(SAMPLE_SET_FILE) | sed -n '$(1)p')

ifneq ($(wildcard $(SAMPLE_SET_FILE)),)
  NUM_SETS := $(shell sed '/^\#/d' $(SAMPLE_SET_FILE) | wc -l | cut -f 1 -d' ')
  SETS_SEQ := $(shell seq 1 $(NUM_SETS))
  $(foreach i,$(SETS_SEQ), \
	  $(eval set.$i := $(call get_space_sets,$i)))
  # define set_lookup.$(sample)
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set_lookup.$(sample) := $i)))
  # define SAMPLE_SETS to contain sample_sample_normal
  $(foreach i,$(SETS_SEQ), \
	  $(eval SAMPLE_SETS += $(call get_underscore_sets,$i)))
  # tumor.SET => tumors
  $(foreach i,$(SETS_SEQ), \
	$(eval tumor.$(call get_underscore_sets,$i) := $(call get_tumors,$(set.$i))))
  # normal.SET => normal
  $(foreach i,$(SETS_SEQ), \
	$(eval normal.$(call get_underscore_sets,$i) := $(call get_normal,$(set.$i))))
  NORMAL_SAMPLES = $(filter-out poolednorm,$(foreach i,$(SETS_SEQ),$(call get_normal,$(set.$i))))
  TUMOR_SAMPLES = $(foreach i,$(SETS_SEQ),$(call get_tumors,$(set.$i)))
  SAMPLE_PAIRS = $(foreach set,$(SAMPLE_SETS),$(foreach tumor,$(tumor.$(set)),$(tumor)_$(normal.$(set))))
  # define tumor.pair and normal.pair to retrieve tumor/normal from pairs
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor.$(tumor)_$(normal.$(set)) := $(tumor))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval normal.$(tumor)_$(normal.$(set)) := $(normal.$(set)))))
  # tumor_normal.$(tumor) to retrieve tumor_normal
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(tumor) := $(tumor)_$(normal.$(set)))))
  $(foreach set,$(SAMPLE_SETS), \
	  $(foreach tumor,$(tumor.$(set)), \
		  $(eval tumor_normal.$(normal.$(set)) := $(tumor)_$(normal.$(set)))))
  # set.$(sample) to retrieve sample sets (underscores)
  $(foreach i,$(SETS_SEQ), \
	  $(foreach sample,$(set.$i), \
		  $(eval set.$(sample) := $(word $i,$(SAMPLE_SETS)))))
  UNMATCHED_SAMPLES = $(shell sed '/^\#/d' $(SAMPLE_FILE) $(SAMPLE_SET_FILE) | tr ' ' '\n' | sort | uniq -u)

  SAMPLE_SET_PAIRS = $(shell echo "$(SAMPLE_PAIRS) $(SAMPLE_SETS)" | tr ' ' '\n' | sort | uniq)
  $(foreach set,$(SAMPLE_SET_PAIRS), \
	$(eval samples.$(set) := $(shell echo "$(set)" | sed 's/_/ /g')))
endif

ifneq ($(wildcard $(SAMPLE_SPLIT_FILE)),)
  A = $(shell sed '/^\#/d' $(SAMPLE_SPLIT_FILE) | cut -f1)
  B = $(shell sed '/^\#/d' $(SAMPLE_SPLIT_FILE) | cut -f2)
  $(foreach i,$(shell seq 1 $(words $(A))),$(eval split.$(word $i,$(A)) += $(word $i,$(B))))
  UNSPLIT_SAMPLES = $(B)
  SPLIT_SAMPLES = $(shell sed '/^\#/d' $(SAMPLE_SPLIT_FILE) | cut -f1 | sort | uniq)
endif

ALL_SAMPLES = $(SAMPLE_PAIRS) $(SAMPLES)
COUNT_SAMPLES = $(shell expr `sed 's/\..*//; s:.*/::' <<< $1 | grep -o "_" | wc -l` + 1)
PANEL_OF_NORMAL_SAMPLES ?= $(NORMAL_SAMPLES)
POOLED_NORM_SAMPLES ?= $(NORMAL_SAMPLES)

#######################################################
########### DEFINE GENOME AND RESOURCES ###############
#######################################################

RESOURCE_REQ ?= medium
include usb-modules-v2/genome_inc/$(REF).inc
include usb-modules-v2/resources.inc

#####################################
### DEFINE BASIC FUNCTIONS ##########
#####################################

NO_RM ?= false
ifeq ($(NO_RM),true)
    RM := touch
    RMR = touch
else
    RM := rm -f
    RMR := rm -r
endif

# define $(,) and $( ) for escaping commas and spaces
, := ,
space :=
space +=
$(space) :=
$(space) +=

# $(call strip-suffix,filename)
strip-suffix = $(firstword $(subst ., ,$(1)))

LINK = ln -svf $(notdir $1) $(2) && touch $1
NOW := $(shell date +"%F")
MKDIR = mkdir -p -m 775
MKDIRS = $(MKDIR) $(LOGDIR)/$(@D) $(@D)
LOG = $(LOGDIR)/$(@).log
UMASK = 002

INIT = module purge; umask $(UMASK); $(MKDIRS); set -o pipefail;
PURGE_AND_LOAD = module purge && module load $1
#PURGE_AND_LOAD = export LMOD_DISABLE_SAME_NAME_AUTOSWAP='no' && module purge && module load $1

MD5 = md5sum $(@:.md5=) > $@
CHECK_MD5 = for i in {0..20}; do if md5sum -c $(filter %.md5,$^); then break; fi; sleep 5; done;
<M = $(<:.md5=)
^M = $(^:.md5=)
@M = $(@:.md5=)
<<M = $(word 2,$(^M))
<<<M = $(word 3,$(^M))

# $(call CHECK_VCF,vcf.file,target,command)
# check for valid header, ie non-empty vcf file
CHECK_VCF = if [ `grep -v '^\#' $1 | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi
CHECK_MAF = if [ `grep -v '^\#' $1 | wc -l` -eq 1 ] && [ `grep '^Hugo_Symbol' $1 | wc -l` -eq 1 ]; then cp $1 $2; else $3; fi
CHECK_VCF_CMD = if [ `grep -v '^\#' $1 | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then $2; else $3; fi
CHECK_VCF_HIGH_MODERATE_CMD = if [ `grep -v '^\#' $1 | egrep "HIGH|MODERATE" | wc -l` -eq 0 ] && [ `grep '^\#CHROM' $1 | wc -l` -eq 1 ]; then $2; else $3; fi

<< = $(word 2,$^)
<<< = $(word 3,$^)
<<<< = $(word 4,$^)
<<<<< = $(word 5,$^)

define R_INIT
dir.create('$(@D)', showWarnings = F, recursive = T)
qw <- function(x) unlist(strsplit(x, "[[:space:]]+"))
endef

##########################################
##### DEFINE THE MAIN RUN FUNCTION #######
##########################################

RUN_OPTS = $(if $(findstring false,$(USE_CLUSTER)),-l) -o $@ \
$(if $(CLUSTER_ENGINE),-g $(CLUSTER_ENGINE)) $(if $(EXCLUDE_NODE),-x $(EXCLUDE_NODE)) \
-S $(SHELL) --log_file $(LOG) $(if $(PROJECT_NAME),-p $(PROJECT_NAME)) -c 

## RUN is the main function!!
## arg1 = number of cores
## arg2 = memory per core
## arg3 = wall time

RUN = $(call PURGE_AND_LOAD,$(PYTHON_MODULE)); echo "\
$(call PURGE_AND_LOAD,$4); " $5 | python usb-modules-v2/scripts/run.py $(RUN_OPTS) -n $1 -m $2 -w $3

RUN += $(if $(findstring SLURM,$(CLUSTER_ENGINE)),\
$(if $(findstring $(QOS_2w),$3),--qos 2weeks,\
$(if $(findstring $(QOS_1w),$3),--qos 1week,\
$(if $(findstring $(QOS_1d),$3),--qos 1day,\
$(if $(findstring $(QOS_6h),$3),--qos 6hours,\
$(if $(findstring $(QOS_30min),$3),--qos 30min,infinite))))),)

###############################################################
####### DEFINE PATHS/NAMES FOR SCRIPTS, BINARIES, MODULES #####
###############################################################

USR_DIR ?= /scicore/home/pissal00/GROUP/usr_nobackup
JAR_DIR ?= $(USR_DIR)/local/lib/java
SCRIPTS_DIR ?= $(USR_DIR)/local/scripts
BIN_DIR ?= $(USR_DIR)/local/bin
MODULE_SCRIPTS_DIR ?= usb-modules-v2/scripts
LD_LIBRARY_PATH = $(USR_DIR)/local/MCR_R2014a/v83/runtime/glnxa64:$(USR_DIR)/local/MCR_R2014a/v83/bin/glnxa64:$(USR_DIR)/local/MCR_R2014a/v83/sys/os/glnxa64
IMG_DIR ?= /scicore/home/pissal00/GROUP/virtual_env/img
HOTSPOTS_DIR ?= $(USR_DIR)/local/cancer_hotspots

# use system $TEMPDIR on scicore (/scratch)
export TMPDIR := $(TEMPDIR)
export PATH := $(USR_DIR):$(PATH)
export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH)

SHELL ?= /bin/bash
UNZIP ?= /usr/bin/unzip
GZIP ?= gzip -f
GUNZIP ?= gunzip -c

PERL_MODULE ?= Perl/5.22.2-goolf-1.7.20
# R_MODULE ?= R/3.5.2-foss-2018b (breaks other modules that rely on local packages compiled for 3.5.0)
R_MODULE ?= R/3.5.0-goolf-1.7.20
R4_MODULE ?= R/4.0.3-foss-2018b
JAVA6_MODULE ?= Java/1.6.0_24
JAVA7_MODULE ?= Java/1.7.0_80
JAVA8_MODULE ?= Java/8u212b03
PYTHON_MODULE ?= Python/3.5.2-goolf-1.7.20
SAMTOOLS_MODULE ?= SAMtools/1.9-foss-2018b
BCFTOOLS_MODULE ?= BCFtools/1.5-goolf-1.7.20
BEDTOOLS_MODULE ?= BEDTools/2.26.0-goolf-1.7.20
VCFTOOLS_MODULE ?= VCFtools/0.1.14-goolf-1.7.20-Perl-5.22.2
BCL2FASTQ_MODULE ?= bcl2fastq/2.17.1.14-goolf-1.7.20
TRIM_GALORE_MODULE ?= Trim_Galore/0.4.4-goolf-1.7.20-Perl-5.22.2
CUTADAPT_MODULE ?= cutadapt/1.16-goolf-1.7.20-Python-2.7.11
HTSLIB_MODULE ?= HTSlib/1.3.1-goolf-1.7.20
TABIX_MODULE ?= $(HTSLIB_MODULE)
FASTQC_MODULE ?= FastQC/0.11.4-Java-1.7.0_21
BOWTIE2_MODULE ?= Bowtie2/2.3.2-goolf-1.7.20
BWA_MODULE ?= BWA/0.7.17-goolf-1.7.20
TOPHAT_MODULE ?= TopHat/2.0.6.Linux_x86_64
STAR_MODULE ?= STAR/2.7.3a-foss-2018b
STAR_FUSION_MODULE ?= STAR-Fusion/0.6.0-goolf-1.7.20-Perl-5.22.2
HISAT2_MODULE ?= HISAT2/2.0.4-goolf-1.7.20
STRINGTIE_MODULE ?= StringTie/1.3.3-goolf-1.7.20
SNP_EFF_MODULE ?= snpEff/4.3t-Java-1.8.0_92
CUFFLINKS_MODULE ?= Cufflinks/2.2.1-Linux_x86_64
STRINGTIE_MODULE ?= StringTie/1.3.3-goolf-1.7.20
RSEM_MODULE ?= RSEM/1.3.2-foss-2018b
OPENBLAS_MODULE ?= OpenBLAS/0.2.13-GCC-4.8.4-LAPACK-3.5.0
GATK36_MODULE ?= GATK/3.6-Java-1.8.0_92
GATK40_MODULE ?= GATK/4.0.7.0-Java-1.8.0_92
GATK42_MODULE ?= GATK/4.2.4.1-foss-2018b-Java-11
MUTECT_MODULE ?= MuTect/1.1.7-Java-1.7.0_80
PICARD_MODULE ?= picard/2.9.2
VARSCAN_MODULE ?= VarScan/2.4.1-Java-1.7.0_80
ANACONDA_MODULE ?= Anaconda3/5.0.1
#SINGULARITY_MODULE ?= Singularity/2.6.1 #deprecated on scicore. A newer singularity is available system-wide.
VEP_MODULE ?= ensembl-vep/92.0-goolf-1.7.20-Perl-5.22.2
MACS2_MODULE ?= MACS2/2.1.1.20160309-goolf-1.7.20-Python-2.7.11
MATPLOTLIB_MODULE ?= matplotlib/3.1.1-fosscuda-2019b-Python-3.7.4

PERL ?= perl
JAVA ?= java -Djava.io.tmpdir=$(TMPDIR)
PYTHON ?= python
PYTHON_ENV_ACTIVATE ?= ml load $(ANACONDA_MODULE) && source activate 
PYTHON_ENV_DEACTIVATE = source deactivate

SINGULARITY_EXEC ?= singularity exec --cleanenv --no-home -B /scicore/home/pissal00/GROUP
SINGULARITY_RUN ?= singularity run --cleanenv --no-home -B /scicore/home/pissal00/GROUP

R_LIBS ?= $(USR_DIR)/$(R_MODULE)/lib64/R/library
R4_LIBS ?= $(USR_DIR)/$(R4_MODULE)/lib64/R/library
R ?= export R_LIBS=$(R_LIBS); R
R4 ?= export R_LIBS=$(R4_LIBS); R
RSCRIPT ?= export R_LIBS=$(R_LIBS); Rscript
R4SCRIPT ?= export R_LIBS=$(R4_LIBS); Rscript
KNIT ?= $(RSCRIPT) $(MODULE_SCRIPTS_DIR)/knit.R
SWEAVE ?= $(RSCRIPT) $(MODULE_SCRIPTS_DIR)/Sweave.R
RBIND ?= $(RSCRIPT) $(MODULE_SCRIPTS_DIR)/rbind.R

SAMTOOLS ?= samtools
BCFTOOLS ?= bcftools
BEDTOOLS ?= bedtools
IGVTOOLS ?= $(USR_DIR)/local/IGV_2.12.0/igvtools
VCFTOOLS ?= vcftools
VCFTOOLS_MERGE ?= vcf-merge
BCL2FASTQ ?= bcl2fastq
TRIM_GALORE ?= trim_galore
CUTADAPT ?= cutadapt
BGZIP ?= bgzip
TABIX ?= tabix -p vcf
BAM_READCOUNT ?= $(BIN_DIR)/bam-readcount
BAMUTIL ?= $(USR_DIR)/local/bamUtil/bin/bam
VCF_SORT ?= $(PERL) $(SCRIPTS_DIR)/vcfsorter.pl
FASTQ_TRIMMER = $(PERL) usb-modules-v2/fastq_tools/trimFastq.pl
VCF_PASS ?= usb-modules-v2/vcf_tools/extract_pass_vcf.pl

ifeq ($(ADAPTOR_POS),3prime)
CUTADAPT_OPTS ?= -a $(ADAPTOR_SEQ) $(if $(findstring true,$(PAIRED_END)),-A $(ADAPTOR_SEQ))
endif
ifeq ($(ADAPTOR_POS),5prime)
CUTADAPT_OPTS ?= -g $(ADAPTOR_SEQ) $(if $(findstring true,$(PAIRED_END)),-G $(ADAPTOR_SEQ))
endif
ifeq ($(NEXTSEQ_TRIM),true)
CUTADAPT_OPTS += --nextseq-trim=20
endif
CUTADAPT_OPTS += -q 20,20 --trim-n -m 35 --max-n=3 

FASTQC ?= fastqc
FASTQC_SUMMARY_PLOT = $(RSCRIPT) usb-modules-v2/qc/fastqcSummaryPlot.R
PLOT_RNASEQ_METRICS = $(RSCRIPT) usb-modules-v2/qc/plotRnaseqMetrics.R

CLUSTER_VCF = $(R4SCRIPT) usb-modules-v2/qc/clusterVcf.R

BOWTIE2 ?= bowtie2
BWA_ALN ?= bwa aln
BWA_SAMSE ?= bwa samse
BWA_SAMPE ?= bwa sampe
BWA_MEM ?= bwa mem
TOPHAT ?= tophat
TOPHAT2 ?= tophat2
STAR ?= STAR
STAR_FUSION ?= STAR-Fusion
HISAT2 ?= hisat2
STRINGTIE ?= stringtie

#PICARD = java -Xmx$(RESOURCE_REQ_MEDIUM_MEM_JAVA) -jar \$$EBROOTPICARD/picard.jar
PICARD = $(JAVA) -Xmx$(RESOURCE_REQ_MEDIUM_MEM_JAVA) -jar $(USR_DIR)/local/picard_2.20.0/picard.jar $(1)
#PICARD += MAX_RECORDS_IN_RAM=4000000 TMP_DIR=$(TMPDIR)
#PICARD += REFERENCE_SEQUENCE=$(REF_FASTA) CREATE_INDEX=true
#PICARD += VALIDATION_STRINGENCY=LENIENT

GATK = $(JAVA) -Xmx$(2) -jar $(JAR_DIR)/GenomeAnalysisTK-3.6/GenomeAnalysisTK-3.6.jar 
GATK += -S LENIENT -T $(1) -filterNoBases
GATK37 = $(JAVA) -Xmx$(2) -jar $(JAR_DIR)/GenomeAnalysisTK-3.7/GenomeAnalysisTK-3.7.jar 
GATK37 += -S LENIENT -T $(1) -filterNoBases

GATK42 = gatk --java-options \"-Xmx$(2) -Xms$(2)\" $(1)
GATK4120 = $(USR_DIR)/local/gatk-4.1.2.0/gatk --java-options \"-Xmx$(2) -Xms$(2)\" $(1)
# gatk-4.1.3.0 is used in sufam (https://github.com/charlottekyng/usb-modules-v2/commit/f7f8ee6db4a33e48247812efc0430337af1d073c#commitcomment-36581870)
GATK4130 = $(USR_DIR)/local/gatk-4.1.3.0/gatk --java-options \"-Xmx$(2) -Xms$(2)\" $(1)
GATK4141 = $(USR_DIR)/local/gatk-4.1.4.1/gatk --java-options \"-Xmx$(2) -Xms$(2)\" $(1)
# gatk-4.1.5.0 is the latest that still has FilterByOrientationBias (see  https://gatk.broadinstitute.org/hc/en-us/community/posts/360060284132-FilterByOrientationBias-vanished-in-4-1-6-0)
GATK4150 = $(USR_DIR)/local/gatk-4.1.5.0/gatk --java-options \"-Xmx$(2) -Xms$(2)\" $(1)
GATK4241 = $(USR_DIR)/local/gatk-4.2.4.1/gatk --java-options \"-Xmx$(2) -Xms$(2)\" $(1)

FIX_GATK_VCF = $(PERL) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_gatk_vcf.pl
FIX_VARSCAN_VCF = $(PERL) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_varscan_vcf.pl
FIX_STRELKA_VCF = $(PERL) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_strelka_vcf.pl
FIX_TVC_VCF = $(PERL) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_tvc_vcf.pl
# These need pyhton3:
FIX_STRELKA2_VCF_INDELS = $(call PURGE_AND_LOAD,$(PYTHON_MODULE))&& $(PYTHON) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_strelka2_vcf_indels.py
FIX_STRELKA2_VCF_SNVS = $(call PURGE_AND_LOAD,$(PYTHON_MODULE))&& $(PYTHON) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_strelka2_vcf_snvs.py

#MUTECT ?= $(JAVA) -Xmx$(1) -jar $(JAR_DIR)/muTect-1.1.4.jar --analysis_type MuTect

## For the future
MUTECT ?= java -Xmx$(1) -jar \$$EBROOTMUTECT/mutect-1.1.7.jar --analysis_type MuTect
#MUTECT2 ?= java -Xmx$(2) -jar \$$EBROOTGATK/GenomeAnalysisTK.jar --analysis_type MuTect2

CONFIGURE_STRELKA = $(PERL) $(USR_DIR)/local/strelka_workflow-1.0.15/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = usb-modules-v2/variant_callers/somatic/strelka_config_bwa_default.ini

# Strelka2 runs with python2, which is the default on Scicore, and does not need to be loaded.
CONFIGURE_STRELKA2_SOMATIC = $(PYTHON) $(USR_DIR)/local/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py
STRELKA2_SOMATIC_CONFIG = usb-modules-v2/variant_callers/somatic/configureStrelkaSomaticWorkflow.py.ini
CONFIGURE_STRELKA2_GERMLINE = $(PYTHON) $(USR_DIR)/local/strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py
STRELKA2_GERMLINE_CONFIG = usb-modules-v2/variant_callers/configureStrelkaGermlineWorkflow.py.ini

PIPEIT_IMG ?= $(IMG_DIR)/PipeIT_2.0.1.img

TVC_ROOT_DIR ?= /PipeIT_resources/usr/local
TVC = $(SINGULARITY_EXEC) $(PIPEIT_IMG) $(TVC_ROOT_DIR)/bin/variant_caller_pipeline.py
TVC_SOMATIC_JSON = usb-modules-v2/variant_callers/somatic/TVC_somatic_parameters.json
TVC_SENSITIVE_JSON = usb-modules-v2/variant_callers/somatic/TVC_somatic_parameters_sensitive.json
TVC_MOTIF ?= $(TVC_ROOT_DIR)/share/TVC/sse/ampliseqexome_germline_p1_hiq_motifset.txt
TVC_NUM_CORES ?= 4

TS_IMG ?= $(IMG_DIR)/TS_v5.12.1.sp1_v1.0.1.img
FUSION_TS_FILTER = $(RSCRIPT) usb-modules-v2/sv_callers/fusionTSfilter.R
FUSION_TS_FILTER_MIN_TOTAL_READS ?= 500000
FUSION_TS_FILTER_MIN_OVERALL_E2E_PERCENT ?= 0.7
FUSION_TS_FILTER_MIN_BREAKPOINT_READS ?= 20
FUSION_TS_FILTER_MIN_BREAKPOINT_E2E_PERCENT ?= 0.6
FUSION_TS_FILTER_MAX_BREAKPOINT_E2E_STRANDBIAS ?= 0.2
FUSION_TS_FILTER_MIN_PARTNER_GENES_EXPRS ?= 20



VARSCAN_JAR ?= $(JAR_DIR)/varscan/VarScan.v2.4.3.jar
VARSCAN = $(JAVA) -Xmx$(RESOURCE_REQ_MEDIUM_MEM_JAVA) -jar $(VARSCAN_JAR)

### For the future
# VARSCAN ?= java -Xmx$(RESOURCE_REQ_LOW_MEM_JAVA) -jar \$$EBROOTVARSCAN/VarScan.v2.4.1.jar

#SNP_EFF_CONFIG ?= usb-modules-v2/vcf_tools/snpEff.config
SNP_EFF_DATA_DIR ?= /scicore/home/pissal00/GROUP/ref_nobackup/snpeff/4.3
SNP_EFF = $(JAVA) -Xmx$(1) -jar \$$EBROOTSNPEFF/snpEff.jar
SNP_EFF_OPTS = -i vcf -o vcf -configOption data.dir=$(SNP_EFF_DATA_DIR) -noStats $(SNP_EFF_FLAGS)
SNP_SIFT = $(JAVA) -Xmx$(1) -jar \$$EBROOTSNPEFF/SnpSift.jar
#SNP_SIFT_OPTS = -c $(SNP_EFF_CONFIG)
VCF_EFF_ONE_PER_LINE = $(PERL) \$$EBROOTSNPEFF/snpEff/scripts/vcfEffOnePerLine.pl
VCF_JOIN_EFF = $(PERL) $(MODULE_SCRIPTS_DIR)/joinEff.pl

VEP_CACHE_DIR ?= /scicore/home/pissal00/GROUP/ref_nobackup/vep/
VEP ?= vep --cache --offline --dir_cache $(VEP_CACHE_DIR) --no_stats --force_overwrite --pick --format vcf --vcf --symbol --terms SO
VEP_PLUGIN_DIR ?= $(USR_DIR)/downloads/VEP_plugins


ADD_GENE_LIST_ANNOTATION = $(RSCRIPT) usb-modules-v2/vcf_tools/addGeneListAnnotationToVcf.R
ANNOTATE_FACETS_VCF = $(RSCRIPT) usb-modules-v2/copy_number/annotateFacets2Vcf.R

FACETS_SNP_PILEUP = $(BIN_DIR)/snp-pileup
FACETS = $(R4SCRIPT) usb-modules-v2/copy_number/facetsRun.R
FACETS_GENE_CN = $(RSCRIPT) usb-modules-v2/copy_number/facetsGeneCN.R
FACETS_GENE_CN_PLOT = $(RSCRIPT) usb-modules-v2/copy_number/facetsGeneCNPlot.R
FACETS_CCF = $(RSCRIPT) usb-modules-v2/copy_number/facetsCCF.R

MUTATION_SUMMARY_RSCRIPT = $(RSCRIPT) usb-modules-v2/summary/mutation_summary_excel.R

SIGNATURES = $(if $(findstring NONE,$(CAPTURE_METHOD)), signatures.genome.cosmic.v3.may2019, $(if $(findstring BAITS,$(CAPTURE_METHOD)),signatures.exome.cosmic.v3.may2019))
DECONSTRUCTSIGS_TRI_COUNT_METHOD = default
DECONSTRUCTSIGS = $(RSCRIPT) usb-modules-v2/mut_sigs/deconstructSigs.R $(if $(findstring hg38,$(REF)), --hg38) $(if $(findstring NONE,$(CAPTURE_METHOD)), --signatures.ref $(SIGNATURES), $(if $(findstring BAITS,$(CAPTURE_METHOD)), --signatures.ref $(SIGNATURES))) --tri.counts.method $(DECONSTRUCTSIGS_TRI_COUNT_METHOD)
DECONSTRUCTSIGS_NUMCORES = 1

MUT_SIG_COSMIC_DIR = $(REF_DIR)/mut_sig_cosmic
MUT_SIG_COSMIC = signatures.exome.cosmic.v3.may2019
MUTATIONALPATTERNS = $(RSCRIPT) usb-modules-v2/mut_sigs/MutationalPatterns.R --signatures $(MUT_SIG_COSMIC_DIR)/$(MUT_SIG_COSMIC).rda $(if $(findstring hg38,$(REF)),--ref hg38, $(if $(findstring b37,$(REF)),--ref hg19))


PVACSEQ ?= pvacseq
PVACTOOLS_IMG ?= $(IMG_DIR)/pvactools.img
PVACSEQ_OPTS ?= -e 8,9,10,11 \
$(if $(findstring true,$(PVACSEQ_USE_LOCAL_IEDB)),--iedb-install-directory /MHC/,)
PVACSEQ_PREPARE_VCF ?= $(PERL) usb-modules-v2/neoepitopes/pvacseq_prepare_vcf.pl

LST = $(RSCRIPT) usb-modules-v2/mut_sigs/lst.R

PYCLONE ?= $(PYTHON_ENV_ACTIVATE) /scicore/home/pissal00/GROUP/usr_nobackup/local/conda/pyclonenv; PyClone
PYCLONE_ITER ?= 10000
PYCLONE_BURNIN ?= 1000
PYCLONE_SEED ?= 1234567
PYCLONE_MAKE_MUT_TXT ?= $(RSCRIPT) usb-modules-v2/clonality/pyclone_make_input.R
PYCLONE_MAKE_CLUSTER_TABLE ?= $(RSCRIPT) usb-modules-v2/clonality/pyclone_make_cluster_table.R
PYCLONE_DECONSTRUCTSIGS ?= $(RSCRIPT) usb-modules-v2/clonality/pyclone_deconstructSigs.R $(if $(findstring hg38,$(REF)), --hg38) $(if $(findstring NONE,$(CAPTURE_METHOD)), --signatures.ref $(SIGNATURES), $(if $(findstring BAITS,$(CAPTURE_METHOD)), --signatures.ref $(SIGNATURES)))
PYCLONE_POSTPY ?= $(RSCRIPT) usb-modules-v2/clonality/pyclone_postpy.R
PYCLONE_POSTPY_MAXSD ?= 0.3

EXPANDS_MAKE_INPUT ?= $(RSCRIPT) usb-modules-v2/clonality/expands_make_input.R
EXPANDS ?= $(RSCRIPT) usb-modules-v2/clonality/expandsRun.R
EXPANDS_NUM_CORES ?= 1

ABSOLUTE_MAKE_MUTS ?= $(RSCRIPT) usb-modules-v2/clonality/absolute_make_muts.R
ABSOLUTE_MAKE_SEGS ?= $(RSCRIPT) usb-modules-v2/clonality/absolute_make_segs.R
ABSOLUTE_STEP1 ?= $(RSCRIPT) usb-modules-v2/clonality/absolute_step1.R
ABSOLUTE_STEP2 ?= $(RSCRIPT) usb-modules-v2/clonality/absolute_step2.R
ABSOLUTE_STEP3 ?= $(RSCRIPT) usb-modules-v2/clonality/absolute_step3.R
ABSOLUTE_NUM_CORE ?= 1


CBS_SEGMENTCNV = $(RSCRIPT) usb-modules-v2/copy_number/varscanCNVSegment.R
VARSCAN_GENE_CN = $(RSCRIPT) usb-modules-v2/copy_number/varscanCNVGeneCN.R
VARSCAN_CONTROL_PROFILE_FROM_STAR ?= $(RSCRIPT) usb-modules-v2/copy_number/varscan_ekaryotyping_control_from_star.R
VARSCAN_RNA_LOGRATIO = $(RSCRIPT) usb-modules-v2/copy_number/varscan_input_from_RNAseq.R
VARSCAN_GENE_CN_PLOT ?= $(FACETS_GENE_CN_PLOT)
VARSCAN_GENE_CN_PLOT_OPTS ?= $(FACETS_GENE_CN_PLOT_OPTS)


CNVKIT_IMG ?= $(IMG_DIR)/cnvkit_0.9.5.img
CNVKIT ?= cnvkit.py
CNVKIT_GENE_CN ?= $(RSCRIPT) usb-modules-v2/copy_number/cnvkitGeneCN.R
CNVKIT_GENE_CN_PLOT ?= $(FACETS_GENE_CN_PLOT) 
CNVKIT_PROCESS_RSEM ?= $(RSCRIPT) usb-modules-v2/copy_number/cnvkit_process_rsem_input.R
CNVKIT_MEDIAN_CENTER ?= $(RSCRIPT) usb-modules-v2/copy_number/cnvkit_median_center.R
CNVKIT_SCALE_SEGMENTS ?= $(RSCRIPT) usb-modules-v2/copy_number/cnvkitScaleSegments.R

MUTSIGCV ?= $(USR_DIR)/local/MutSigCV_1.4/run_MutSigCV.sh
MCR ?= $(USR_DIR)/local/MCR_R2013a_glnxa64/v81/
MUTSIGCV_MAKE_INPUT ?= $(RSCRIPT) usb-modules-v2/siggenes/make_mutSigCV_input.R

YOUN_AND_SIMON ?= $(RSCRIPT) usb-modules-v2/siggenes/youn_and_simon.R
YOUN_AND_SIMON_MAKE_INPUT ?= $(RSCRIPT) usb-modules-v2/siggenes/youn_and_simon_make_input.R

GISTIC ?= $(USR_DIR)/local/gistic_2.0.23/gp_gistic2_from_seg
GISTIC_MAKE_MARKERS_FILE ?= $(RSCRIPT) usb-modules-v2/copy_number/gistic_make_markers.R
GISTIC_MAKE_SEG_FILE ?= $(RSCRIPT) usb-modules-v2/copy_number/gistic_make_seg_from_facets.R
GISTIC_MAKE_CNV_FILE ?= $(RSCRIPT) usb-modules-v2/copy_number/gistic_make_cnv.R

CUFFLINKS ?= cufflinks
CUFFCMP ?= currcompare

RSEM_CALC_EXPR ?= rsem-calculate-expression
RSEM_GEN_DATA_MATRIX ?= usb-modules-v2/rnaseq/rsem-generate-data-matrix.pl
RSEM_PROCCESS ?= usb-modules-v2/rnaseq/rsem_process.R

VIPER_RUN ?= $(RSCRIPT) usb-modules-v2/rnaseq/viper_run.R

STAR_PROCESS ?= usb-modules-v2/rnaseq/star_process.R

KALLISTO ?= $(USR_DIR)/local/kallisto_v0.46.1/kallisto
KALLISTO_NUM_CORES ?= 4
KALLISTO_OPTIONS ?= --fusion
# --bias: learns parameters for a model of sequences specific bias and corrects the abundances accordingly.
KALLISTO_BIAS ?= true
# --genomebam: Project pseudoalignments to genome sorted BAM file
KALLISTO_GENOMEBAM ?= true
# only used for single-end:
# kallisto manual: Typical Illumina libraries produce fragment lengths ranging from 180–200 bp but it’s best to determine this from a library quantification.
KALLISTO_FRAGMENT_LEN ?= 200

PIZZLY ?= $(USR_DIR)/local/pizzly_v0.37.3/pizzly
# k-mer size used in kallisto (kallisto default: 31)
PIZZLY_K ?= 31
PIZZLY_OPTIONS ?= 
PIZZLY_FLATTEN_JSON ?= $(USR_DIR)/local/pizzly_v0.37.3/scripts/flatten_json.py

MOSAICS_RUN ?= $(RSCRIPT) usb-modules-v2/chipseq/mosaics_run.R
MOSAICS_CONSTRUCTBINS ?= $(RSCRIPT) usb-modules-v2/chipseq/mosaics_constructBins.R

MACS2 ?= macs

MSISENSOR ?= $(BIN_DIR)/msisensor_Linux_x86_64

DETIN_DIR = /scicore/home/pissal00/GROUP/usr_nobackup/local/deTiN_2021-11-25/deTiN
DETIN ?= python $(DETIN_DIR)/deTiN.py

TCELLEXTRECT ?= $(RSCRIPT) usb-modules-v2/Tcell_quant/runTcellExTRECT.R

endif
MAKEFILE_INC = true
