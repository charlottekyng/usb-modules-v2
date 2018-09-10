include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/pvacseq.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: pvacseq

pvacseq : $(foreach pair,$(SAMPLE_PAIRS),$(foreach prefix,$(CALLER_PREFIX),pvacseq/results/MHC_Class_I/$(pair).$(prefix).final.tsv))


pvacseq/input/%.vep.vcf : pvacseq/input/%.vcf
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(VEP_MODULE),"\
	$(VEP) --plugin Downstream --plugin Wildtype --dir_plugins $(VEP_PLUGIN_DIR) --assembly GRCh37 \
	--input_file $< --output_file $@")

define prepare-pvacseq
pvacseq/input/$1_$2.$3.vcf : vcf/$1_$2.$$(call VCF_SUFFIXES,$3).vcf
	$$(INIT) $$(PVACSEQ_PREPARE_VCF) < $$< > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach prefix,$(CALLER_PREFIX),\
		$(eval $(call prepare-pvacseq,$(tumor.$(pair)),$(normal.$(pair)),$(prefix)))))

define run-pvacseq
pvacseq/results/MHC_Class_I/$1_$2.$3.final.tsv : pvacseq/input/$1_$2.$3.vep.vcf pvacseq/optitype/$2.results.tsv
	$$(RM) pvacseq/tmp/$1_$2; \
	HLA_GENOTYPE=`tail -1 $$(<<) | cut -f2-7 | sed 's/^/HLA-/; s/\t/,HLA-/g'`; \
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
		$$(SINGULARITY_EXEC) $$(PVACTOOLS_IMG) $$(PVACSEQ) run $$< $1_$2 $$$$HLA_GENOTYPE $$(PVACSEQ_ALGO) \
		pvacseq/tmp/$1_$2 $$(PVACSEQ_OPTS) && \
		mv pvacseq/tmp/$1_$2/MHC_Class_I/$1_$2.final.tsv $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach prefix,$(CALLER_PREFIX),\
		$(eval $(call run-pvacseq,$(tumor.$(pair)),$(normal.$(pair)),$(prefix)))))

#define run-optitype
#pvacseq/optitype/$1.results.tsv : fastq/$1.1.fastq $$(if $$(PAIRED_END),fastq/$1.2.fastq,)
#	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
#		$$(SINGULARITY_RUN) $$(OPTITYPE_IMG) -i $^ $$(if NOT RNA,--dna,--rna) -o pvacseq/optitype/$1 && \
#		mv $(shell ls) && $(RMR) xx")
#endef
#$(foreach normal,$(NORMAL_SAMPLES),\
#	$(eval $(call run-optitype,$(normal.$(pair)))))



