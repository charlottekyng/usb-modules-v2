TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed
TARGETS_FILE_INTERVALS ?= $(BED_DIR)/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed
TARGETS_FILE_MUTATIONS ?= $(BED_DIR)/SeqCap_EZ_Exome_v3_hg19_primary_targets.bed
DBSNP_TARGETS_INTERVALS ?= $(BED_DIR)/SeqCap_EZ_Exome_v3_hg19_primary_targets.dbsnp_138.vcf
RESOURCE_REQ ?= high
SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-downstream -no-intergenic -no-upstream -no-utr

