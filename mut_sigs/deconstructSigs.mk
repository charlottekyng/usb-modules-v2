include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/deconstructSigs.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: deconstructSigs

deconstructSigs : deconstructSigs/all$(PROJECT_PREFIX).deconstructSigs.RData

deconstructSigs/all$(PROJECT_PREFIX).deconstructSigs.RData : summary/mutation_summary$(PROJECT_PREFIX).$(DOWNSTREAM_EFF_TYPES).txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
		$(DECONSTRUCTSIGS) --outPrefix=$(subst .RData,,$@) \
		--num_iter $(DECONSTRUCTSIGS_NUMITER) --num_cores $(DECONSTRUCTSIGS_NUMCORES) $<")

include usb-modules-v2/summary/mutationSummary.mk
