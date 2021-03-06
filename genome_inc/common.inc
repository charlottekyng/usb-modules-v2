REF_DIR = /scicore/home/pissal00/GROUP/ref_nobackup

GENOME_DIR = $(REF_DIR)/genomes
ANN_DIR = $(REF_DIR)/annotation_data/$(REF)
BED_DIR = $(REF_DIR)/target_beds/$(REF)


FACETS_TARGETS_INTERVALS ?= $(DBSNP_TARGETS_INTERVALS)

STAR_GENOME_DIR ?= $(REF_DIR)/genome_index/star_$(REF)
STAR_CTAT_DIR ?= $(ANN_DIR)/STAR-FusionFilter

RSEM_INDEX ?= $(REF_DIR)/genome_index/rsem_$(REF)/rsem_index

HISAT2_REF ?= $(REF_DIR)/genome_index/hisat2_$(REF)/genome_snp_tran

KALLISTO_INDEX ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).transcripts.shortname.k31.idx
KALLISTO_GTF ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.gtf
KALLISTO_CHR ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.chromosomes.len

PIZZLY_GTF_CACHE ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).primary_assembly.annotation.gtf.pizzly.cache
PIZZLY_FASTA ?= $(REF_DIR)/genome_index/kallisto/$(REF)/gencode.$(GENCODE_VER)/gencode.$(GENCODE_VER).transcripts.shortname.fa.gz
