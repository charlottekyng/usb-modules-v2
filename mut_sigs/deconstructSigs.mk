include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/deconstructSigs.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: deconstructSigs

MUT_CALLER = mutect

ALLTABLE = $(lastword $(call SOMATIC_TABLES,$(MUT_CALLER)))
OUTPREFIX = $(notdir $(basename $(ALLTABLE)))

deconstructSigs : deconstructSigs/$(OUTPREFIX).deconstructSigs.RData

deconstructSigs/$(OUTPREFIX).deconstructSigs.RData : $(ALLTABLE)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
		$(DECONSTRUCTSIGS) --outPrefix=$(OUTPREFIX) $<")
