# Bed files and annotation files

TARGETS_FILE_COVERED_INTERVALS ?= $(BED_DIR)/S0276129_Covered_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_INTERVALS ?= $(BED_DIR)/S0276129_Padded_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_MUTATIONS ?= $(BED_DIR)/S0276129_Padded_noAnn_noHeader_GRCm37liftover.bed
TARGETS_FILE_GENES ?= $(BED_DIR)/S0276129_Targets.txt

DBSNP_TARGETS_INTERVALS ?= $(BED_DIR)/S0276129_Padded_noAnn_noHeader_GRCm37liftover.mgp.v5.merged.snps_all.dbSNP142_noGT.vcf
DBSNP_TARGETS_INTERVALS_COMMON ?= $(BED_DIR)/S0276129_Padded_noAnn_noHeader_GRCm37liftover.mgp.v5.merged.snps_all.dbSNP142_noGT.monoallelic.vcf

# Resources
RESOURCE_REQ ?= medium

# Annotation and reporting parameters
# SNP_EFF_FLAGS ?= -canon -no NEXT_PROT -no-downstream -no-intergenic -no-upstream -no-utr

