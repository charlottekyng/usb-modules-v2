include usb-modules-v2/Makefile.inc

LOGDIR ?= log/rsem.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: rsem

rsem : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results) rsem/all$(PROJECT_PREFIX).genes.expected_count.results rsem/all$(PROJECT_PREFIX).genes.TPM.results rsem/all$(PROJECT_PREFIX).genes.FPKM.results rsem/all$(PROJECT_PREFIX).isoforms.expected_count.results rsem/all$(PROJECT_PREFIX).isoforms.TPM.results rsem/all$(PROJECT_PREFIX).isoforms.FPKM.results rsem/all$(PROJECT_PREFIX).genes.expected_count.results_coding_upper_quartiled

define rsem-calc-expression
rsem/$1.genes.results : star/$1.Aligned.toTranscriptome.out.bam 
	$$(call RUN,$$(RSEM_NUM_CORES),$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),$$(PERL_MODULE) $$(RSEM_MODULE),"\
	$$(RSEM_CALC_EXPR) $$(RSEM_OPTIONS) $$(if $$(findstring true,$$(PAIRED_END)),--paired-end) \
	$$< $$(RSEM_INDEX) $$(@D)/$1")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call rsem-calc-expression,$(sample))))

rsem/%.isoforms.results : rsem/%.genes.results
	

rsem/all$(PROJECT_PREFIX).genes.expected_count.results : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	$(RSEM_GEN_DATA_MATRIX) expected_count $^ > $@")

rsem/all$(PROJECT_PREFIX).genes.TPM.results : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	$(RSEM_GEN_DATA_MATRIX) TPM $^ > $@")

rsem/all$(PROJECT_PREFIX).genes.FPKM.results : $(foreach sample,$(SAMPLES),rsem/$(sample).genes.results)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	$(RSEM_GEN_DATA_MATRIX) FPKM $^ > $@")

rsem/all$(PROJECT_PREFIX).isoforms.expected_count.results : $(foreach sample,$(SAMPLES),rsem/$(sample).isoforms.results)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	$(RSEM_GEN_DATA_MATRIX) expected_count $^ > $@")

rsem/all$(PROJECT_PREFIX).isoforms.TPM.results : $(foreach sample,$(SAMPLES),rsem/$(sample).isoforms.results)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	$(RSEM_GEN_DATA_MATRIX) TPM $^ > $@")

rsem/all$(PROJECT_PREFIX).isoforms.FPKM.results : $(foreach sample,$(SAMPLES),rsem/$(sample).isoforms.results)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),$(PERL_MODULE),"\
	$(RSEM_GEN_DATA_MATRIX) FPKM $^ > $@")

rsem/all$(PROJECT_PREFIX).genes.expected_count.results_coding_upper_quartiled : rsem/all.genes.expected_count.results
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RSCRIPT) $(RSEM_NORM) --inputRSEMFile $< --gtf $(GENCODE_GTF) --outputFile $@ \
	--normalizationMethod \"uq\" --threshold_for_uq 1000")
	

include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/aligners/align.mk
include usb-modules-v2/aligners/starAligner.mk

