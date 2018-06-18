# run expands for determining tumor ploidy

include usb-modules-v2/Makefile.inc

EXPANDS_NUM_CORES ?= 4

LOGDIR = log/expands.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: expands

expands : $(foreach pair,$(SAMPLE_PAIRS),expands/input/$(pair).mutations.txt expands/input/$(pair).segments.txt)

expands/output/%.expands.RData : expands/input/%.mutations.txt expands/input/%.segments.txt
	$(MKDIR) expands/output; 
	$(call RUN,$(EXPANDS_NUM_CORES),$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(EXPANDS) --mutations $(<) --segs $$(<<) --nc $(EXPANDS_NUM_CORES) --outPrefix $(subst .RData,,$@)")

define expands_make_input
expands/input/$1_$2.mutations.txt : $$(foreach prefix,$$(CALLER_PREFIX),tables/$1_$2.$$(call SOMATIC_VCF_SUFFIXES,$$(prefix)).$$(TABLE_SUFFIX).txt)
        $$(MKDIR) expands/input; \
        $$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
        $$(RBIND) --tumorNormal $$^ > $$@.tmp1 && \
        $$(MUTATION_SUMMARY_RSCRIPT) --outFile $$@.tmp2 --outputFormat TXT $$@.tmp1 && \
        $$(EXPANDS_MAKE_INPUT) --outFile $$@ --type mutations $$@.tmp2 && \
        $$(RM) $$@.tmp1 $$@.tmp2")

expands/input/$1_$2.segments.txt : facets/cncf/$1_$2.cncf.txt
	$$(MKDIR) expands/input; \
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(EXPANDS_MAKE_INPUT) --outFile $$@ --type cna $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call expands_make_input,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/variants/somatic/
