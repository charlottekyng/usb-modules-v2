MUT_CALLER = mutect2

#### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc
#include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/mutect2.$(NOW)

PHONY += mutect2 #mutect_vcfs mutect_tables ext_output #mut_report
..DUMMY := $(shell mkdir -p version; echo "$(MUTECT2) &> version/mutect2.txt")

mutect2 : mutect2_vcfs mutect2_tables #mutect2_vcfs_hotspotgt mutect2_tables_hotspotgt #ext_output
mutect2_vcfs : $(call MAKE_VCF_FILE_LIST,mutect2) $(addsuffix .idx,$(call MAKE_VCF_FILE_LIST,mutect2))
mutect2_tables : $(call MAKE_TABLE_FILE_LIST,mutect2)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

define mutect2-tumor-normal-chr
mutect2/chr_vcf/$1_$2.$3.mutect2%vcf : bam/$1%bam bam/$2%bam $(PON_VCF)
	$$(MKDIR) mutect2/chr_vcf; $$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_LONG),$$(JAVA8_MODULE),"\
	$$(call GATK,MuTect2,$$(RESOURCE_REQ_HIGH_MEM)) \
	--reference_sequence $$(REF_FASTA) --input_file:tumor $$< --input_file:normal $$(word 2,$$^) \
	--dbsnp $$(DBSNP) --cosmic $$(COSMIC) --intervals $3 \
	--max_alt_alleles_in_normal_count $(MUTECT_MAX_ALT_IN_NORMAL) \
	--max_alt_allele_in_normal_fraction $(MUTECT_MAX_ALT_IN_NORMAL_FRACTION) \
	--annotation TandemRepeatAnnotator --annotation OxoGReadCounts --normal_panel $$(word 3,$$^) \
	--out mutect2/chr_vcf/$1_$2.$3.mutect2.vcf")
endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
		$(eval $(call mutect2-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

# merge variants 
define mutect2-tumor-normal
vcf/$1_$2.mutect2.vcf : $$(foreach chr,$$(CHROMOSOMES),mutect2/chr_vcf/$1_$2.$$(chr).mutect2.vcf)
	$$(INIT) module load $$(PERL_MODULE); grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | \
	$$(VCF_SORT) $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call mutect2-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk
include usb-modules-v2/variant_callers/somatic/pon.mk
