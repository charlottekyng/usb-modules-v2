include usb-modules-v2/Makefile.inc

ifeq ($(findstring SOMATIC,$(ANALYSIS_TYPE)),SOMATIC)
USE_SUFAM = false
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc
else
include usb-modules-v2/variant_callers/variantCaller.inc
endif

EXT_NAMES ?=
FILTER_VARIANTS = false
ANNOTATE_VARIANTS ?= true



LOGDIR ?= log/ann_ext_vcf.$(NOW)

ann_ext_vcf : ext_vcfs ext_tables


ext_vcfs : $(foreach ext_name,$(EXT_NAMES),$(call MAKE_VCF_FILE_LIST,$(ext_name)))
ext_tables : $(foreach ext_name,$(EXT_NAMES),$(call MAKE_TABLE_FILE_LIST,$(ext_name)))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : ann_ext_vcf

include usb-modules-v2/vcf_tools/vcftools.mk
