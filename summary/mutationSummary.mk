include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/summary.$(NOW)

#.DELETE_ON_ERROR:
.SECONDARY: 
PHONY : mutation_summary

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
CALLER_PREFIX ?= mutect strelka_indels
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
CALLER_PREFIX ?= tvc_snps tvc_indels
endif


ifeq ($(findstring EXCEL,$(MUTATION_SUMMARY_FORMAT)),EXCEL)
mutation_summary : $(shell rm -f summary/mutation_summary$(PROJECT_PREFIX).xlsx) summary/mutation_summary$(PROJECT_PREFIX).xlsx
endif
ifeq ($(findstring TXT,$(MUTATION_SUMMARY_FORMAT)),TXT)
mutation_summary : $(shell rm -f summary/mutation_summary$(PROJECT_PREFIX).txt) summary/mutation_summary$(PROJECT_PREFIX).txt
endif

ALLTABLES_NS_SYNON_HS_LNC_PROM = $(foreach prefix,$(CALLER_PREFIX),alltables/allTN$(PROJECT_PREFIX).$(call SOMATIC_VCF_SUFFIXES,$(prefix)).tab.nonsynonymous_synonymous_hotspot_lincRNA_upstream.txt)
ALLTABLES_NS_SYNON_HS_LNC = $(foreach prefix,$(CALLER_PREFIX),alltables/allTN$(PROJECT_PREFIX).$(call SOMATIC_VCF_SUFFIXES,$(prefix)).tab.nonsynonymous_synonymous_hotspot_lincRNA.txt)
ALLTABLES_NS_SYNON_HS = $(foreach prefix,$(CALLER_PREFIX),alltables/allTN$(PROJECT_PREFIX).$(call SOMATIC_VCF_SUFFIXES,$(prefix)).tab.nonsynonymous_synonymous_hotspot.txt)

summary/mutation_summary$(PROJECT_PREFIX).xlsx : $(if $(findstring false,$(INCLUDE_LNCRNA_IN_SUMMARY)),$(ALLTABLES_NS_SYNON_HS),$(ALLTABLES_NS_SYNON_HS_LNC_PROM))
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(R_MODULE),"\
	$(MUTATION_SUMMARY_RSCRIPT) --outFile $@ $^")

summary/mutation_summary$(PROJECT_PREFIX).txt : $(if $(findstring false,$(INCLUDE_LNCRNA_IN_SUMMARY)),$(ALLTABLES_NS_SYNON_HS),$(ALLTABLES_NS_SYNON_HS_LNC_PROM))
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(R_MODULE),"\
	$(MUTATION_SUMMARY_RSCRIPT) --outputFormat TXT --outFile $@ $^")

#define mut_sum_txt
#summary/mutation_summary$(PROJECT_PREFIX)_$1.txt : $$(if $$(findstring false,$$(INCLUDE_LNCRNA_IN_SUMMARY)),alltables/allTN$(PROJECT_PREFIX).$$(call SOMATIC_VCF_SUFFIXES,$1).tab.nonsynonymous_synonymous_hotspot.txt,alltables/allTN.$$(call SOMATIC_VCF_SUFFIXES,$1).tab.nonsynonymous_synonymous_hotspot_lincRNA.txt)
#	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
#	$$(MUTATION_SUMMARY_RSCRIPT) --outFile $$@ --outputFormat TXT $$<")
#endef
#$(foreach prefix,$(CALLER_PREFIX),$(eval $(call mut_sum_txt,$(prefix))))
