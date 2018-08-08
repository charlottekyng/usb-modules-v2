ifndef B37_INC
REF = b37
REF_DIR = /scicore/home/terracci/GROUP/ref_nobackup

GSA_REF_DIR ?= $(REF_DIR)/gsa_resource_bundle/2.8/b37
ANN_DIR ?= $(REF_DIR)/annotation_data/b37
BED_DIR ?= $(REF_DIR)/target_beds/b37

REF_FASTA ?= $(GSA_REF_DIR)/human_g1k_v37_decoy.fasta

DBSNP ?= $(GSA_REF_DIR)/dbsnp_138.b37.vcf

ifeq ($(findstring NONE,$(PANEL)),NONE)
else
include usb-modules-v2/genome_includes/$(REF).$(PANEL).inc
endif

RESOURCE_REQ ?= high
TARGETS_FILE_INTERVALS ?= $(GSA_REF_DIR)/human_g1k_v37_decoy.bed
DBSNP_TARGETS_INTERVALS ?= $(GSA_REF_DIR)/dbsnp_138.b37_snps_monoallelic.vcf
DBSNP_TARGETS_INTERVALS_COMMON ?= $(DBSNP_TARGETS_INTERVALS)

FACETS_TARGETS_INTERVALS ?= $(DBSNP_TARGETS_INTERVALS)


ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort -V | uniq)
else
CHROMOSOMES ?= $(shell seq 1 22) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif

### ANNOTATION SETS ####

CONTEST_POP_FILE ?= $(ANN_DIR)/hg19_population_stratified_af_hapmap_3.3.vcf

KNOWN_INDELS ?= $(GSA_REF_DIR)/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI ?= $(GSA_REF_DIR)/1000G_omni2.5.b37.vcf
HAPMAP ?= $(GSA_REF_DIR)/hapmap_3.3.b37.vcf
R1000G ?= $(GSA_REF_DIR)/1000G_phase1.snps.high_confidence.b37.vcf

CLINVAR = $(ANN_DIR)/clinvar.vcf
EXAC_NONTCGA = $(ANN_DIR)/ExAC_nonTCGA.r0.3.1.sites.vep.edited.vcf
COSMIC ?= $(ANN_DIR)/CosmicCodingMuts.sorted.vcf
EXAC_NONPSYCH = $(ANN_DIR)/ExAC.r0.3.nonpsych.sites.edited.vcf
CADD_SNV ?= $(ANN_DIR)/CADD_whole_genome_SNVs_header.vcf.gz
CADD_INDEL ?= $(ANN_DIR)/CADD_Indels_1000Gphase3_indels_header_sorted.vcf.gz

HAPLOTYPE_INSUF_BED = $(ANN_DIR)/haploinsufficient.genes.bed
CANCER_GENE_CENSUS_BED = $(ANN_DIR)/cancer_gene_census.bed
KANDOTH_BED = $(ANN_DIR)/kandoth127smg.bed
LAWRENCE_BED = $(ANN_DIR)/lawrence_cancer5000s.bed
DGD_BED = $(ANN_DIR)/dgd_Hsa_all_v71.bed
CANCER_HOTSPOT_VCF = $(ANN_DIR)/publication_hotspots_inclnoncoding_incl3d.vcf

CENTROMERE_TABLE = $(ANN_DIR)/centromeres_b37.txt

DGV_FILE = $(ANN_DIR)/DGV_GRCh37_hg19_variants_2016-05-15.txt
#KNOWN_CNVS = $(ANN_DIR)/GRCh37_hg19_variants_2013-07-23.txt
#REF_HOUSEKEEPING_GENE_BED = $(ANN_DIR)/share/reference/hg19.HouseKeepingGenes.bed
#ENSEMBL_TXDB = $(ANN_DIR)/hsapiens_ensembl_biomart.sqlite

GENCODE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_patch_hapl_scaff.gtf
GENCODE_GENE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_gsa_faidx.genes.gtf
GENCODE_CODING_GENE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_gsa_faidx.genes.protein_coding.gtf
GENCODE_NONCODING_GENE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_gsa_faidx.genes.non_coding.gtf

RIBOSOMAL_INTERVALS ?= $(ANN_DIR)/ribosomalRNA.interval_list
GENE_REF_FLAT = $(ANN_DIR)/refFlat.txt

#SNP_TXT = $(ANN_DIR)/b37_snp137.SingleDiNucl.1based.nochr.txt
MOSAICS_CHRFILE = $(GSA_REF_DIR)/human_g1k_v37_decoy.mosaics.chrfile

### TOOL-SPECIFIC REF ###

GISTIC_REF ?= /scicore/home/terracci/GROUP/usr_nobackup/local/gistic_2.0.23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat

STAR_GENOME_DIR = $(REF_DIR)/genome_index/star_b37
STAR_CTAT_DIR = $(ANN_DIR)/GRCh37_gencode_v19_CTAT_lib

RSEM_INDEX = $(REF_DIR)/genome_index/rsem_b37/rsem_index

HISAT2_REF = $(REF_DIR)/genome_index/hisat2_b37/genome_snp_tran

#CHIMERASCAN_REF ?= $(REF_DIR)/chimerascan_index
#BOWTIE_REF = $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37

#MAPSPLICE_REF_DIR = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/chr
#MAPSPLICE_REF_BASENAME = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/hg19

SNP_EFF_GENOME = GRCh37.75
DB_NSFP ?= $(REF_DIR)/dbNSFP/GRCh37/dbNSFP2.9.txt.gz

REF_DICT = $(REF_FASTA:.fasta=.dict)

MUTSIGCV_REF_DIR ?= $(REF_DIR)/MutSigCV
MUTSIGCV_COVERAGE_REF ?= $(MUTSIGCV_REF_DIR)/coverage.exome_full192.txt
MUTSIGCV_COV_REF ?= $(MUTSIGCV_REF_DIR)/gene.covariates.txt
MUTSIGCV_DICT_REF ?= $(MUTSIGCV_REF_DIR)/mutation_type_dictionary_file.txt
MUTSIGCV_SEQ_REF_DIR ?= $(MUTSIGCV_REF_DIR)/chr_files_hg19/

YOUN_AND_SIMON_SEQ_DATA ?= $(ANN_DIR)/youn_and_simon_sequence_data.RData

MSISENSOR_REF ?= $(REF_DIR)/msi_reference/b37.txt

endif
B37_INC = true

