#######################################################
############# DEFINE PATHS TO MAIN DIRS ###############
#######################################################

GROUP_DIR ?= /storage/research/dbmr_rubin_lab/pipeline
USR_DIR ?= $(GROUP_DIR)/usr
JAR_DIR ?=
SCRIPTS_DIR ?= $(USR_DIR)/local/scripts
BIN_DIR ?=
MODULE_SCRIPTS_DIR ?= usb-modules-v2/scripts
LD_LIBRARY_PATH = $(USR_DIR)/local/MCR_R2014a/v83/runtime/glnxa64:$(USR_DIR)/local/MCR_R2014a/v83/bin/glnxa64:$(USR_DIR)/local/MCR_R2014a/v83/sys/os/glnxa64
IMG_DIR ?= $(GROUP_DIR)/img
HOTSPOTS_DIR ?= $(USR_DIR)/local/cancer_hotspots

# export PATH := $(USR_DIR):$(PATH)
# export LD_LIBRARY_PATH := $(LD_LIBRARY_PATH)

include usb-modules-v2/config.inc

##########################################
##### DEFINE THE MAIN RUN FUNCTION #######
##########################################


RUN_OPTS = $(if $(findstring false,$(USE_CLUSTER)),-l) -o $@ \
$(if $(CLUSTER_ENGINE),-g $(CLUSTER_ENGINE)) $(if $(EXCLUDE_NODE),-x $(EXCLUDE_NODE)) \
$(if $(PARTITION),--partition $(PARTITION)) \
-S $(SHELL) --log_file $(LOG) $(if $(PROJECT_NAME),-p $(PROJECT_NAME)) -c 

## RUN is the main function!!
## arg1 = number of cores
## arg2 = memory per core
## arg3 = wall time

RUN = $(call PURGE_AND_LOAD,$(PYTHON_MODULE)); echo "\
$(call PURGE_AND_LOAD,$4); " $5 | python usb-modules-v2/scripts/run.py $(RUN_OPTS) -n $1 -m $2 -w $3

RUN += $(if $(findstring SLURM,$(CLUSTER_ENGINE)),\
$(if $(findstring $(QOS_2w),$3),--qos job_cpu_long,\
$(if $(findstring $(QOS_1w),$3),--qos job_cpu_long,\
$(if $(findstring $(QOS_1d),$3),--qos job_cpu,\
$(if $(findstring $(QOS_6h),$3),--qos job_cpu,\
$(if $(findstring $(QOS_30min),$3),--qos job_cpu,job_cpu_long))))),)


###############################################################
####### DEFINE PATHS/NAMES FOR SCRIPTS, BINARIES, MODULES #####
###############################################################

SHELL ?= /bin/bash
UNZIP ?= /usr/bin/unzip
GZIP ?= gzip -f
GUNZIP ?= gunzip -c

PERL_MODULE ?= Perl-bundle-CPAN/5.38.2-GCCcore-13.3.0
# will run R via singularity
# R_MODULE ?= R/3.5.0-goolf-1.7.20
# R4_MODULE ?= R/4.0.3-foss-2018b
# JAVA6_MODULE ?=
# JAVA7_MODULE ?=
# gatk-v3.7-0 won't work properly with java above 8. Current system default is 1.8.0, works fine:
JAVA8_MODULE ?=
PYTHON_MODULE ?=
SAMTOOLS_MODULE ?= SAMtools/1.13-GCC-10.3.0
SAMTOOLS_MODULE_GOOLF ?= SAMtools/1.13-GCC-10.3.0
BCFTOOLS_MODULE ?= BCFtools/1.12-GCC-10.3.0
BEDTOOLS_MODULE ?= BEDTools/2.30.0-GCC-10.3.0
VCFTOOLS_MODULE ?= VCFtools/0.1.16-GCC-10.3.0
BCL2FASTQ_MODULE ?=
TRIM_GALORE_MODULE ?=
CUTADAPT_MODULE ?= cutadapt/3.4-GCCcore-10.3.0
HTSLIB_MODULE ?= HTSlib/1.12-GCC-10.3.0
TABIX_MODULE ?= $(HTSLIB_MODULE)
FASTQC_MODULE ?= FastQC/0.11.9-Java-11
BOWTIE2_MODULE ?= Bowtie2/2.4.4-GCC-10.3.0
BWA_MODULE ?= BWA/0.7.17-GCC-10.3.0
TOPHAT_MODULE ?=
STAR_MODULE ?=
STAR_FUSION_MODULE ?=
HISAT2_MODULE ?= HISAT2/2.2.1-gompi-2021a
STRINGTIE_MODULE ?=
SNP_EFF_MODULE ?=
CUFFLINKS_MODULE ?=
RSEM_MODULE ?= RSEM/1.3.3-foss-2021a
OPENBLAS_MODULE ?= OpenBLAS/0.3.27-GCC-13.3.0
MUTECT_MODULE ?=
PICARD_MODULE ?=
VARSCAN_MODULE ?=
ANACONDA_MODULE ?=
VEP_MODULE ?=
MACS2_MODULE ?= MACS2/2.2.7.1-foss-2021a
MATPLOTLIB_MODULE ?= matplotlib/3.4.2-foss-2021a
PYCLONE_MODULE ?=
SIG_PROFILER_MODULE ?=

PERL ?= perl
PYTHON ?= python
PYTHON2 ?= python2

# shell $TMPDIR is set only on compute nodes (resolves to /scratch/local/<job id>), on login this is empty. Pass it properly
TMPDIR ?= \$$TMPDIR
JAVA ?= java -Djava.io.tmpdir=$(TMPDIR)

SINGULARITY_EXEC ?= singularity exec --cleanenv --no-home -B $(GROUP_DIR)/../
SINGULARITY_RUN ?= singularity run --cleanenv --no-home -B $(GROUP_DIR)/../

# on scicore there are some R modules that use v3, therefore the two R variables
# here everything runs in one version of R, but keeping the two variables for compatibility
R_LIBS ?= $(USR_DIR)/R_lib/Rocker_4.2.2
R4_LIBS ?= $(R_LIBS)
R ?= SINGULARITYENV_R_LIBS=$(R_LIBS) $(SINGULARITY_EXEC) $(IMG_DIR)/R_RStudio/Rocker_rstudio_4.2.2/Rocker_rstudio_4.2.2.sif R
R4 ?= $(R)
RSCRIPT ?= SINGULARITYENV_R_LIBS=$(R_LIBS) $(SINGULARITY_EXEC) $(IMG_DIR)/R_RStudio/Rocker_rstudio_4.2.2/Rocker_rstudio_4.2.2.sif Rscript
R4SCRIPT ?= $(RSCRIPT)

KNIT ?= $(RSCRIPT) $(MODULE_SCRIPTS_DIR)/knit.R
SWEAVE ?= $(RSCRIPT) $(MODULE_SCRIPTS_DIR)/Sweave.R
RBIND ?= $(RSCRIPT) $(MODULE_SCRIPTS_DIR)/rbind.R

SAMTOOLS ?= samtools
BCFTOOLS ?= bcftools
BEDTOOLS ?= bedtools
IGVTOOLS ?= $(USR_DIR)/local/IGV_2.16.1/igvtools
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
VCF_PASS ?= $(PERL) usb-modules-v2/vcf_tools/extract_pass_vcf.pl

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
STAR ?= $(USR_DIR)/local/STAR_2.7.11b/Linux_x86_64_static/STAR
STAR_CPU ?= 14
STAR_FUSION ?= STAR-Fusion
HISAT2 ?= hisat2
STRINGTIE ?= stringtie


PICARD = $(JAVA) -Xmx$(RESOURCE_REQ_MEDIUM_MEM_JAVA) -jar $(USR_DIR)/local/picard_2.25.2/picard.jar $(1) TMP_DIR=$(TMPDIR)

GATK = $(JAVA) -Xmx$(2) -jar $(USR_DIR)/local/gatk-3.7/GenomeAnalysisTK-3.7.jar 
GATK += -S LENIENT -T $(1) -filterNoBases

GATK42 = $(USR_DIR)/local/gatk-4.2.4.1/gatk --java-options \"-Xmx$(2) -Xms$(2)\" $(1)
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
FIX_STRELKA2_VCF_INDELS = $(PYTHON) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_strelka2_vcf_indels.py
FIX_STRELKA2_VCF_SNVS = $(PYTHON) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_strelka2_vcf_snvs.py

MUTECT ?=

CONFIGURE_STRELKA = $(PERL) $(USR_DIR)/local/strelka_workflow-1.0.15/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = usb-modules-v2/variant_callers/somatic/strelka_config_bwa_default.ini

STRELKA2_IMG ?= $(IMG_DIR)/strelka2/strelka2_v2.9.10.sif
CONFIGURE_STRELKA2_SOMATIC = $(SINGULARITY_EXEC) $(STRELKA2_IMG) configureStrelkaSomaticWorkflow.py
STRELKA2_SOMATIC_CONFIG = usb-modules-v2/variant_callers/somatic/configureStrelkaSomaticWorkflow.py.ini
CONFIGURE_STRELKA2_GERMLINE = $(SINGULARITY_EXEC) $(STRELKA2_IMG) configureStrelkaGermlineWorkflow.py
STRELKA2_GERMLINE_CONFIG = usb-modules-v2/variant_callers/configureStrelkaGermlineWorkflow.py.ini
STRELKA2 ?= $(SINGULARITY_EXEC) $(STRELKA2_IMG)

PIPEIT_IMG ?= $(IMG_DIR)/PipeIT_2.0.5.img

BAMIXCHECKER_IMG ?= $(IMG_DIR)/BAMixChecker/BAMixChecker_v1.0.1.sif

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

MUSE_IMG ?= $(IMG_DIR)/MuSE/MuSE_2.0.4.sif
MUSE = $(SINGULARITY_EXEC) $(MUSE_IMG) MuSE
FIX_MUSE_VCF = $(PERL) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_muse_vcf.pl

CAVEMAN_IMG ?= $(IMG_DIR)/CaVEMan/cgpCaVEManWrapper_1.18.3.sif
CAVEMAN ?= $(SINGULARITY_EXEC) $(CAVEMAN_IMG) caveman.pl
FIX_CAVEMAN_VCF = $(PYTHON) usb-modules-v2/variant_callers/fix_vcf_scripts/fix_caveman_vcf.py
MAKE_CAVEMAN_SEGS = $(RSCRIPT) usb-modules-v2/scripts/caveman_make_segs.R

MANTA_IMG ?= $(IMG_DIR)/Manta/manta_1.6.0.sif
MANTA ?= $(SINGULARITY_EXEC) $(MANTA_IMG)

DELLY_IMG ?= $(IMG_DIR)/Delly/delly_v1.2.6.sif
DELLY ?= $(SINGULARITY_EXEC) $(DELLY_IMG)

SVABA_IMG ?= $(IMG_DIR)/SvABA/svaba_1.2.0-ab39e3f.sif
SVABA ?= $(SINGULARITY_EXEC) $(SVABA_IMG)

### For the future
# VARSCAN ?= java -Xmx$(RESOURCE_REQ_LOW_MEM_JAVA) -jar \$$EBROOTVARSCAN/VarScan.v2.4.1.jar

#SNP_EFF_CONFIG ?= usb-modules-v2/vcf_tools/snpEff.config
SNP_EFF_DATA_DIR ?= $(REF_DIR)/snpeff/5.0

# needs java 21, available locally:
SNP_EFF = $(USR_DIR)/local/java/jdk-21.0.8/bin/java -Xmx$(1) -jar $(USR_DIR)/local/snpEff/v.5.2f/snpEff.jar
SNP_EFF_OPTS = -i vcf -o vcf -noStats $(SNP_EFF_FLAGS) -nodownload -dataDir $(SNP_EFF_DATA_DIR)
SNP_SIFT = $(USR_DIR)/local/java/jdk-21.0.8/bin/java -Xmx$(1) -jar $(USR_DIR)/local/snpEff/v.5.2f/SnpSift.jar
#SNP_SIFT_OPTS = -c $(SNP_EFF_CONFIG)
VCF_EFF_ONE_PER_LINE = $(PERL) $(USR_DIR)/local/snpEff/v.5.2f/scripts/_OLD/vcfEffOnePerLine.pl
VCF_JOIN_EFF = $(PERL) $(MODULE_SCRIPTS_DIR)/joinEff.pl

VEP_CACHE_DIR ?= $(REF_DIR)/vep/
VEP ?= vep --cache --offline --dir_cache $(VEP_CACHE_DIR) --no_stats --force_overwrite --pick --format vcf --vcf --symbol --terms SO
VEP_PLUGIN_DIR ?= $(USR_DIR)/downloads/VEP_plugins


ADD_GENE_LIST_ANNOTATION = $(RSCRIPT) usb-modules-v2/vcf_tools/addGeneListAnnotationToVcf.R
ANNOTATE_FACETS_VCF = $(RSCRIPT) usb-modules-v2/copy_number/annotateFacets2Vcf.R

FACETS_SNP_PILEUP_IMG ?= $(IMG_DIR)/snp-pileup/snp-pileup.v0.6.1.sif
FACETS_SNP_PILEUP = $(SINGULARITY_EXEC) $(FACETS_SNP_PILEUP_IMG) snp-pileup
FACETS = $(R4SCRIPT) usb-modules-v2/copy_number/facetsRun.R
FACETS_GENE_CN = $(RSCRIPT) usb-modules-v2/copy_number/facetsGeneCN.R
FACETS_GENE_CN_PLOT = $(RSCRIPT) usb-modules-v2/copy_number/facetsGeneCNPlot.R
FACETS_CCF = $(RSCRIPT) usb-modules-v2/copy_number/facetsCCF.R
FACETS_HETMARKFREQ = $(R4SCRIPT) usb-modules-v2/copy_number/facetsHetMarkFreq.R

MUTATION_SUMMARY_RSCRIPT = $(RSCRIPT) usb-modules-v2/summary/mutation_summary_excel.v2.R

SIGNATURES = $(if $(findstring NONE,$(CAPTURE_METHOD)), signatures.genome.cosmic.v3.may2019, $(if $(findstring BAITS,$(CAPTURE_METHOD)),signatures.exome.cosmic.v3.may2019))
DECONSTRUCTSIGS_TRI_COUNT_METHOD = default
DECONSTRUCTSIGS = $(RSCRIPT) usb-modules-v2/mut_sigs/deconstructSigs.R $(if $(findstring hg38,$(REF)), --hg38) $(if $(findstring NONE,$(CAPTURE_METHOD)), --signatures.ref $(SIGNATURES), $(if $(findstring BAITS,$(CAPTURE_METHOD)), --signatures.ref $(SIGNATURES))) --tri.counts.method $(DECONSTRUCTSIGS_TRI_COUNT_METHOD)
DECONSTRUCTSIGS_NUMCORES = 1

MUT_SIG_COSMIC_DIR = $(REF_DIR)/mut_sig_cosmic
MUT_SIG_COSMIC = signatures.exome.cosmic.v3.may2019
MUTATIONALPATTERNS = $(RSCRIPT) usb-modules-v2/mut_sigs/MutationalPatterns.R --signatures $(MUT_SIG_COSMIC_DIR)/$(MUT_SIG_COSMIC).rda $(if $(findstring hg38,$(REF)),--ref hg38, $(if $(findstring b37,$(REF)),--ref hg19))

SIG_PROFILER_ASSIGNMENT = $(PYTHON) usb-modules-v2/mut_sigs/run_SigProfilerAssignment.py

PVACSEQ ?= pvacseq
PVACTOOLS_IMG ?= $(IMG_DIR)/pvactools.img
PVACSEQ_OPTS ?= -e 8,9,10,11 \
$(if $(findstring true,$(PVACSEQ_USE_LOCAL_IEDB)),--iedb-install-directory /MHC/,)
PVACSEQ_PREPARE_VCF ?= $(PERL) usb-modules-v2/neoepitopes/pvacseq_prepare_vcf.pl

LST = $(RSCRIPT) usb-modules-v2/mut_sigs/lst.R

PYCLONE ?= PyClone
PYCLONE_ITER ?= 10000
PYCLONE_BURNIN ?= 1000
PYCLONE_SEED ?= 1234567
PYCLONE_MAKE_MUT_TXT ?= $(RSCRIPT) usb-modules-v2/clonality/pyclone_make_input.R
PYCLONE_MAKE_CLUSTER_TABLE ?= $(RSCRIPT) usb-modules-v2/clonality/pyclone_make_cluster_table.R
PYCLONE_DECONSTRUCTSIGS ?= $(RSCRIPT) usb-modules-v2/clonality/pyclone_deconstructSigs.R $(if $(findstring hg38,$(REF)), --hg38) $(if $(findstring NONE,$(CAPTURE_METHOD)), --signatures.ref $(SIGNATURES), $(if $(findstring BAITS,$(CAPTURE_METHOD)), --signatures.ref $(SIGNATURES))) --tri.counts.method $(DECONSTRUCTSIGS_TRI_COUNT_METHOD)
PYCLONE_POSTPY ?= $(RSCRIPT) usb-modules-v2/clonality/pyclone_postpy.R
PYCLONE_POSTPY_MAXSD ?= 0.3

EXPANDS_MAKE_INPUT ?= $(RSCRIPT) usb-modules-v2/clonality/expands_make_input.R
EXPANDS ?= $(RSCRIPT) usb-modules-v2/clonality/expandsRun.R
EXPANDS_NUM_CORES ?= 1

ABSOLUTE_IMG ?= $(IMG_DIR)/absolute/absolute_v1.0.6_fixed_ccf_posterior_grid.sif
ABSOLUTE_MAKE_MUTS ?= Rscript usb-modules-v2/clonality/absolute_make_muts.R
ABSOLUTE_MAKE_SEGS ?= Rscript usb-modules-v2/clonality/absolute_make_segs.R
ABSOLUTE_STEP1 ?= Rscript usb-modules-v2/clonality/absolute_step1.R
ABSOLUTE_STEP2 ?= Rscript usb-modules-v2/clonality/absolute_step2.R
ABSOLUTE_STEP3 ?= Rscript usb-modules-v2/clonality/absolute_step3.R
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

YOUN_AND_SIMON ?= $(RSCRIPT) usb-modules-v2/siggenes/youn_and_simon.R --resources_dir $(USR_DIR)/local/youn_and_simon
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

VIPER_RUN ?= $(R4SCRIPT) usb-modules-v2/rnaseq/viper_run.R

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

DETIN_DIR = $(GROUP_DIR)/usr_nobackup/local/deTiN_2021-11-25/deTiN
DETIN ?= python $(DETIN_DIR)/deTiN.py

TCELLEXTRECT_IMG ?= $(IMG_DIR)/TcellExTRECT/TcellExTRECT_v1.0.1.sif
TCELLEXTRECT ?= $(SINGULARITY_EXEC) $(TCELLEXTRECT_IMG) Rscript usb-modules-v2/Tcell_quant/runTcellExTRECT.R
