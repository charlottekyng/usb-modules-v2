# GATK variant caller module
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
# 

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/variantCaller.inc

LOGDIR ?= log/gatk.$(NOW)

VPATH ?= bam

VARIANT_TYPES ?= gatk_snps gatk_indels
PHONY += gatk

gatk : gatk_vcfs $(if $(findstring NONE,$(PANEL)),,gatk_tables)
gatk_vcfs : $(foreach type,$(VARIANT_TYPES),$(call VCFS,$(type)) $(addsuffix .idx,$(call VCFS,$(type))))
gatk_tables : $(foreach type,$(VARIANT_TYPES),$(call TABLES,$(type)))
gatk_reports : $(foreach type,$(VARIANT_TYPES),reports/$(type).dp_ft.grp)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

include usb-modules-v2/vcf_tools/vcftools.mk
include usb-modules-v2/variant_callers/gatk.mk

