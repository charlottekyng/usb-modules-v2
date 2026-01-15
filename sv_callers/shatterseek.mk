include usb-modules-v2/Makefile.inc
include usb-modules-v2/sv_callers/conSV.mk


.DEFAULT_GOAL := shatterseek
LOGDIR ?= log/shatterseek.$(NOW)
PHONY += shatterseek

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

shatterseek : $(foreach pair,$(SAMPLE_PAIRS),shatterseek/$(tumor.$(pair))/$(tumor.$(pair)).shatterseekResults.csv)

define shatterseek
shatterseek/$1/$1.shatterseekResults.csv : facets/cncf/$1_$2.cncf.txt conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/$1.$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE).tsv
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(SHATTERSEEK) --sample $1 --CN $$(<) --SV $$(<<) --outFile $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call shatterseek,$(tumor.$(pair)),$(normal.$(pair)))) \
)