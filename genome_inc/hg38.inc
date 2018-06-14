ifndef HG38_INC
REF = hg38
REF_DIR = /scicore/home/terracci/GROUP/ref_nobackup

GSA_REF_DIR ?= $(REF_DIR)/gsa_resource_bundle/hg38
ANN_DIR ?= $(REF_DIR)/annotation_data/hg38
BED_DIR ?= $(REF_DIR)/target_beds/hg38

REF_FASTA ?= $(GSA_REF_DIR)/Homo_sapiens_assembly38.fasta

DBSNP ?= $(GSA_REF_DIR)/dbsnp_146.hg38.vcf

RESOURCE_REQ ?= medium

$(info RESOURCE_REQ is $(RESOURCE_REQ))

ifeq ($(findstring bed,$(TARGETS_FILE_INTERVALS)),bed)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE_INTERVALS) | cut -f1 | sort -V | uniq)
else
CHROMOSOMES ?= chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX $(if $(findstring true,$(INCLUDE_CHR_Y)),chrY) chrM
endif

ifeq ($(PANEL),WXS)
TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/wxs.bed
TARGETS_FILE_INTERVALS ?= $(TARGETS_FILE_COVERED_INTERVALS)
TARGETS_FILE_MUTATIONS ?= $(TARGETS_FILE_COVERED_INTERVALS)
TARGETS_FILE_GENES ?= $(BED_DIR)/wxs.GeneSymbols.txt
DBSNP_TARGETS_INTERVALS_COMMON ?= $(BED_DIR)/dbsnp_146.hg38.common.wxs.vcf
RESOURCE_REQ ?= medium
SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -onlyProtein
endif



FACETS_TARGETS_INTERVALS ?= $(DBSNP_TARGETS_INTERVALS_COMMON)

### ANNOTATION SETS ####

CONTEST_REF_VCF ?= $(GSA_REF_DIR)/hapmap_3.3_grch38_pop_stratified_af.vcf

KNOWN_INDELS ?= $(GSA_REF_DIR)/Mills_and_1000G_gold_standard.indels.hg38.vcf
OMNI ?= $(GSA_REF_DIR)/1000G_omni2.5.hg38.vcf
HAPMAP ?= $(GSA_REF_DIR)/hapmap_3.3.hg38.vcf

#CLINVAR = $(ANN_DIR)/clinvar.vcf
#EXAC_NONTCGA = $(ANN_DIR)/ExAC_nonTCGA.r0.3.1.sites.vep.vcf
#COSMIC ?= $(ANN_DIR)/CosmicCodingMuts.sorted.vcf
#EXAC_NONPSYCH = $(ANN_DIR)/ExAC.r0.3.nonpsych.sites.edited.vcf
#CADD_SNV ?= $(ANN_DIR)/CADD_whole_genome_SNVs_header.vcf.gz
#CADD_INDEL ?= $(ANN_DIR)/CADD_Indels_1000Gphase3_indels_header_sorted.vcf.gz

#HAPLOTYPE_INSUF_BED = $(ANN_DIR)/haploinsufficient.genes.bed
#CANCER_GENE_CENSUS_BED = $(ANN_DIR)/cancer_gene_census.bed
#KANDOTH_BED = $(ANN_DIR)/kandoth127smg.bed
#LAWRENCE_BED = $(ANN_DIR)/lawrence_cancer5000s.bed
#DGD_BED = $(ANN_DIR)/dgd_Hsa_all_v71.bed
#CANCER_HOTSPOT_VCF = $(ANN_DIR)/publication_hotspots_inclnoncoding_incl3d.vcf

#CENTROMERE_TABLE = $(ANN_DIR)/centromeres_b37.txt

#DGV_FILE = $(ANN_DIR)/DGV_GRCh37_hg19_variants_2016-05-15.txt
#KNOWN_CNVS = $(ANN_DIR)/GRCh37_hg19_variants_2013-07-23.txt
#REF_HOUSEKEEPING_GENE_BED = $(ANN_DIR)/share/reference/hg19.HouseKeepingGenes.bed
#ENSEMBL_TXDB = $(ANN_DIR)/hsapiens_ensembl_biomart.sqlite

GENCODE_GTF = $(ANN_DIR)/gencode.v27.primary_assembly.annotation.gtf
GENCODE_GENE_GTF = $(ANN_DIR)/Homo_sapiens.GRCh37.82.chr_patch_hapl_scaff.genes.gtf

RIBOSOMAL_INTERVALS ?= $(ANN_DIR)/ribosomal.interval_list
GENE_REF_FLAT = $(ANN_DIR)/refFlat.txt

#SNP_TXT = $(ANN_DIR)/b37_snp137.SingleDiNucl.1based.nochr.txt

### TOOL-SPECIFIC REF ###

#GISTIC_REF ?= /scicore/home/terracci/GROUP/usr_nobackup/local/gistic_2.0.23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat

#STAR_GENOME_DIR = $(REF_DIR)/genome_index/star_b37
#STAR_CTAT_DIR = $(ANN_DIR)/GRCh37_gencode_v19_CTAT_lib

#RSEM_INDEX = $(REF_DIR)/genome_index/rsem_b37/rsem_index

#HISAT2_REF = $(REF_DIR)/genome_index/hisat2_b37/genome_snp_tran


SNP_EFF_GENOME = GRCh38.86
DB_NSFP ?= $(REF_DIR)/dbNSFP/GRCh38/dbNSFP3.2.txt.gz


REF_DICT = $(REF_FASTA:.fasta=.dict)

#MUTSIGCV_REF_DIR ?= $(REF_DIR)/MutSigCV
#MUTSIGCV_COVERAGE_REF ?= $(MUTSIGCV_REF_DIR)/coverage.exome_full192.txt
#MUTSIGCV_COV_REF ?= $(MUTSIGCV_REF_DIR)/gene.covariates.txt
#MUTSIGCV_DICT_REF ?= $(MUTSIGCV_REF_DIR)/mutation_type_dictionary_file.txt
#MUTSIGCV_SEQ_REF_DIR ?= $(MUTSIGCV_REF_DIR)/chr_files_hg19/

MSISENSOR_REF ?= $(REF_DIR)/msi_reference/hg38.txt

endif
B37_INC = true

