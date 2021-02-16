include usb-modules-v2/Makefile.inc

LOGDIR ?= log/kallisto.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: kallisto


kallisto : $(foreach sample,$(SAMPLES),kallisto/$(sample)/abundance.tsv) \
kallisto/all$(PROJECT_PREFIX).est_count.txt kallisto/all$(PROJECT_PREFIX).tpm.txt


define kallisto
kallisto/$1/abundance.tsv : fastq/$1.1.fastq.gz $(if $(findstring true,$(PAIRED_END)),fastq/$1.2.fastq.gz)
	$$(call RUN,4,3G,$$(RESOURCE_REQ_MEDIUM),,"\
	$$(KALLISTO) quant -i $$(KALLISTO_INDEX) $$(KALLISTO_OPTIONS) \
	$$(if $$(findstring FIRST_READ_TRANSCRIPTION_STRAND,$$(STRAND_SPECIFICITY)),--fr-stranded) \
	$$(if $$(findstring SECOND_READ_TRANSCRIPTION_STRAND,$$(STRAND_SPECIFICITY)),--rf-stranded) \
	--threads $$(KALLISTO_NUM_CORES) \
	$$(if $$(findstring true,$$(KALLISTO_BIAS)),--bias) \
	$$(if $$(findstring true,$$(KALLISTO_GENOMEBAM)),--genomebam --gtf $$(KALLISTO_GTF) --chromosomes $$(KALLISTO_CHR)) \
	-o kallisto/$1 \
	$$(if $$(findstring false,$$(PAIRED_END)),--single --fragment-length $$(KALLISTO_FRAGMENT_LEN)) \
	$$^ && \
	gzip kallisto/$1/fusion.txt")
endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call kallisto,$(sample))))


kallisto/all$(PROJECT_PREFIX).est_count.txt kallisto/all$(PROJECT_PREFIX).tpm.txt : $(foreach sample,$(SAMPLES),kallisto/$(sample)/abundance.tsv)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(RSCRIPT) usb-modules-v2/rnaseq/kallisto_merge.R --prefix kallisto/all$(PROJECT_PREFIX) $^")

include usb-modules-v2/fastq_tools/fastq.mk
