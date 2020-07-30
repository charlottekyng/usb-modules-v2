include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/MutationalPatterns.$(NOW)

ifneq ($(words $(CALLER_PREFIX)),1)
  $(info CALLER_PREFIX contains more than one variant caller)
  $(info Choose only one by executing: make mutational_patterns CALLER_PREFIX=<variant caller>)
  $(info  )
  exit:
	val=1 && exit $${val}
endif

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: MutationalPatterns

MutationalPatterns : MutationalPatterns/all$(PROJECT_PREFIX).$(MSIGDB).RData

MutationalPatterns/all$(PROJECT_PREFIX).$(MSIGDB).RData : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(CALLER_PREFIX).*.hotspot.pass.vcf)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(MUTATIONALPATTERNS) --outPrefix all$(PROJECT_PREFIX) $^")

