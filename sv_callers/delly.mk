MUT_CALLER = delly

# Run delly on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc

LOGDIR ?= log/delly.$(NOW)
PHONY += delly

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

delly : delly_vcfs

delly/sample_info.txt : $(SAMPLE_SET_FILE)
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),,"\
	grep -v '#' $(<) | sed -E -e 's/\t/\ttumor/g' -e 's/\$$/\tcontrol/' | sed 's/tumor/&\n/g' > $@")

delly_vcfs : delly/sample_info.txt $(foreach pair,$(SAMPLE_PAIRS),delly/$(tumor.$(pair))_$(normal.$(pair)).somatic.SVpass.vcf)

define delly-tumor-normal
# Main call
delly/$1_$2.bcf : bam/$1.bam bam/$2.bam
	$$(MKDIR) $$(@D)
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(DELLY) delly call \
	$$(if $$(findstring hg38,$$(REF)),-x $$(BED_DIR)/human.hg38.excl.bed) \
	$$(if $$(findstring b37,$$(REF)),-x $$(BED_DIR)/human.hg19.excl.bed) \
	$$(if $$(findstring GRCm38,$$(REF)),-x $$(BED_DIR)/mouse.mm10.excl.bed) \
	-g $$(REF_FASTA) \
	-o $$@ \
	$$^")

# Somatic pre-filtering
delly/$1_$2.pre.bcf : delly/$1_$2.bcf delly/sample_info.txt
	$$(call RUN,1,1G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(DELLY) delly filter \
	-f somatic \
	-o $$@ \
	-s $$(<<) \
	$$<")

# Genotype pre-filtered somatic sites across all normal samples.
# Do it with all normals in one go. If this becomes too slow with many normals modify this part to run per each normal and then merge.
# Eventually we might want to be able to pass external normals as well.
delly/$1_$2.geno.bcf : delly/$1_$2.pre.bcf bam/$1.bam bam/$2.bam $(foreach normal,$(PANEL_OF_NORMAL_SAMPLES),bam/$(normal).bam)
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(SINGULARITY_MODULE),"\
	$$(DELLY) delly call \
	$$(if $$(findstring hg38,$$(REF)),-x $$(BED_DIR)/human.hg38.excl.bed) \
	$$(if $$(findstring b37,$$(REF)),-x $$(BED_DIR)/human.hg19.excl.bed) \
	$$(if $$(findstring GRCm38,$$(REF)),-x $$(BED_DIR)/mouse.mm10.excl.bed) \
	-g $$(REF_FASTA) \
	-v $$< \
	-o $$@ \
	$$(filter-out $$<,$$^)")

# Post-filter
delly/$1_$2.somatic.bcf : delly/$1_$2.geno.bcf delly/sample_info.txt
	$$(call RUN,1,1G,$$(RESOURCE_REQ_VSHORT),$$(SINGULARITY_MODULE),"\
	$$(DELLY) delly filter \
	-f somatic \
	-o $$@ \
	-s $$(<<) \
	$$<")

delly/$1_$2.somatic.vcf : delly/$1_$2.somatic.bcf
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(BCFTOOLS_MODULE),"\
	$$(BCFTOOLS) view $$< > $$@")

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call delly-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk
