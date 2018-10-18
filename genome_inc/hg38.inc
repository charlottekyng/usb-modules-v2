ifndef HG38_INC
REF = hg38
include usb-modules-v2/genome_inc/common.inc

REF_FASTA ?= $(GENOME_DIR)/Homo_sapiens_assembly38_gsa.fasta
REF_DICT = $(REF_FASTA:.fasta=.dict)

DBSNP ?= $(ANN_DIR)/dbsnp_138.b37.vcf.gz

ifeq ($(findstring NONE,$(PANEL)),NONE)
else
include usb-modules-v2/genome_inc/$(REF)/$(PANEL).inc
endif

RESOURCE_REQ ?= high
TARGETS_FILE_INTERVALS ?= $(REF_FASTA:.fasta=.bed)
DBSNP_TARGETS_INTERVALS ?= $(DBSNP)

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort -V | uniq)
else
CHROMOSOMES ?= $(shell seq 1 22) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif

### ANNOTATION SETS ####

CONTEST_POP_FILE ?= $(ANN_DIR)/hapmap_3.3_grch38_pop_stratified_af.vcf.gz

KNOWN_INDELS ?= $(ANN_DIR)/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
OMNI ?= $(ANN_DIR)/1000G_omni2.5.hg38.vcf.gz
HAPMAP ?= $(ANN_DIR)/hapmap_3.3.hg38.vcf.gzz
R1000G ?= $(ANN_DIR)/1000G_phase1.snps.high_confidence.hg38.vcf.gz

#CLINVAR = $(ANN_DIR)/clinvar.vcf
#EXAC_NONTCGA = $(ANN_DIR)/ExAC_nonTCGA.r0.3.1.sites.vep.edited.vcf
#EXAC_NONPSYCH ?= $(ANN_DIR)/ExAC.r0.3.nonpsych.sites.edited.vcf
#CADD_SNV ?= $(ANN_DIR)/CADD_whole_genome_SNVs_header.vcf.gz
#CADD_INDEL ?= $(ANN_DIR)/CADD_Indels_1000Gphase3_indels_header_sorted.vcf.gz
#COSMIC ?= $(ANN_DIR)/CosmicCodingMuts.sorted.vcf

#HAPLOTYPE_INSUF_BED = $(ANN_DIR)/haploinsufficient.genes.bed
#CANCER_GENE_CENSUS_BED = $(ANN_DIR)/cancer_gene_census.bed
#KANDOTH_BED = $(ANN_DIR)/kandoth127smg.bed
#LAWRENCE_BED = $(ANN_DIR)/lawrence_cancer5000s.bed
#DGD_BED = $(ANN_DIR)/dgd_Hsa_all_v71.bed
#CANCER_HOTSPOT_VCF = $(ANN_DIR)/publication_hotspots_inclnoncoding_incl3d.vcf

CENTROMERE_TABLE = $(ANN_DIR)/centromeres_hg38.txt

#DGV_FILE = $(ANN_DIR)/DGV_GRCh37_hg19_variants_2016-05-15.txt

GENCODE_GTF = $(ANN_DIR)/gencode.v27.primary_assembly.annotation.gtf
GENCODE_GENE_GTF = $(ANN_DIR)/gencode.v27.primary_assembly.annotation_gene.gtf
GENCODE_CODING_GENE_GTF = $(ANN_DIR)/gencode.v27.primary_assembly.annotation_protein_coding_gene.gtf
GENCODE_NONCODING_GENE_GTF = $(ANN_DIR)/gencode.v27.primary_assembly.annotation_gene.non_coding.gtf

RIBOSOMAL_INTERVALS ?= $(ANN_DIR)/ribosomal.interval_list
GENE_REF_FLAT = $(ANN_DIR)/refFlat.txt

### TOOL-SPECIFIC REF ###

#GISTIC_REF ?= /scicore/home/terracci/GROUP/usr_nobackup/local/gistic_2.0.23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat


SNP_EFF_GENOME = GRCh38.81
DB_NSFP ?= $(REF_DIR)/dbNSFP/GRCh38/dbNSFP3.2a.txt.gz

#MOSAICS_CHRFILE = $(GENOME_DIR)/human_g1k_v37_decoy.mosaics.chrfile

#MUTSIGCV_REF_DIR ?= $(REF_DIR)/MutSigCV
#MUTSIGCV_COVERAGE_REF ?= $(MUTSIGCV_REF_DIR)/coverage.exome_full192.txt
#MUTSIGCV_COV_REF ?= $(MUTSIGCV_REF_DIR)/gene.covariates.txt
#MUTSIGCV_DICT_REF ?= $(MUTSIGCV_REF_DIR)/mutation_type_dictionary_file.txt
#MUTSIGCV_SEQ_REF_DIR ?= $(MUTSIGCV_REF_DIR)/chr_files_hg19/

#YOUN_AND_SIMON_SEQ_DATA ?= $(ANN_DIR)/youn_and_simon_sequence_data.RData

MSISENSOR_REF ?= $(REF_DIR)/msi_reference/hg38.txt

endif
HG38_INC = true

