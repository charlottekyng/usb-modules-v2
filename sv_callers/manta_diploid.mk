MUT_CALLER = manta_diploid

# Run manta in diploid mode (single normal sample, not single tumor)

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/manta_diploid.$(NOW)
PHONY += manta_diploid

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

manta_diploid : manta_diploid_vcfs
manta_diploid_vcfs : $(foreach sample,$(PANEL_OF_NORMAL_SAMPLES),manta/$(sample)/results/variants/diploidSV.vcf.gz)

define manta-normal
manta/$1/runWorkflow.py : bam/$1.bam
	$$(call RUN,1,1G,$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	$$(MANTA) configManta.py \
	--bam $$< \
	--referenceFasta $$(REF_FASTA) \
	$$(if $$(findstring BAITS,$$(CAPTURE_METHOD)),--exome,) \
	$$(if $$(findstring hg38,$$(REF)),--callRegions usb-modules-v2/resources/hg38_main_chr.bed.gz,) \
	--runDir $$(@D)")

# manta uses little RAM, 2G per cpu should be enough
manta/$1/results/variants/diploidSV.vcf.gz : manta/$1/runWorkflow.py
	$$(call RUN,$$(MANTA_NUM_CORES),$$(MANTA_MEM)G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(MANTA) $$< -j $$(MANTA_NUM_CORES) -g $$(MANTA_MEM)")

endef
$(foreach sample,$(PANEL_OF_NORMAL_SAMPLES),$(eval $(call manta-normal,$(sample))))

include usb-modules-v2/vcf_tools/vcftools.mk
