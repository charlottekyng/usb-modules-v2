MUT_CALLER = strelka

# Run strelka on tumour-normal matched pairs

include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

##### DEFAULTS ######
STRELKA_NUM_CORES ?= 8


LOGDIR ?= log/strelka.$(NOW)
PHONY += strelka #strelka_vcfs strelka_tables

strelka : strelka_vcfs strelka_tables

VARIANT_TYPES := strelka_indels
#strelka_snps strelka_indels
strelka_vcfs : $(foreach type,$(VARIANT_TYPES),$(call MAKE_VCF_FILE_LIST,$(type)))
strelka_tables : $(foreach type,$(VARIANT_TYPES),$(call MAKE_TABLE_FILE_LIST,$(type)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY : $(PHONY)

define strelka-tumor-normal
# If hg38 *and* PAIRED_END then look for BAMs and REF in strelka_refs folder,
# elif hg38 then look for BAMs and REF also in strelka_refs folder,
# elif PAIRED_END then look for BAMs in bam_clipoverlap folder,
# elif any other case then look for BAMs in the bam folder.
strelka/$1_$2/Makefile : $(if $(and $(findstring hg38,$(REF)),$(findstring true,$(PAIRED_END))),strelka_refs,$(if $(findstring hg38,$(REF)),strelka_refs,$(if $(findstring true,$(PAIRED_END)),bam_clipoverlap,bam)))/$1.bam \
                         $(if $(and $(findstring hg38,$(REF)),$(findstring true,$(PAIRED_END))),strelka_refs,$(if $(findstring hg38,$(REF)),strelka_refs,$(if $(findstring true,$(PAIRED_END)),bam_clipoverlap,bam)))/$2.bam
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(PERL_MODULE),"\
	rm -rf $$(@D) &&\
	$$(CONFIGURE_STRELKA) --tumor=$$< --normal=$$(<<) \
	--ref=$$(if $$(findstring hg38,$(REF)),strelka_refs/$$(REF).fasta,$$(REF_FASTA)) --config=$$(STRELKA_CONFIG) --output-dir=$$(@D)")

# sed will revert any temporary contig name (strelka_refs module) back to the original in the final VCF results.
strelka/$1_$2/results/all.somatic.indels.vcf : strelka/$1_$2/Makefile
	$$(call RUN,$$(STRELKA_NUM_CORES),$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),,"\
	make -j $$(STRELKA_NUM_CORES) -C $$(<D) &&\
	sed -i -e 's/__asterisk__/\*/g' -e 's/__colon__/\:/g' strelka/$1_$2/results/all.somatic.indels.vcf ;\
	sed -i -e 's/__asterisk__/\*/g' -e 's/__colon__/\:/g' strelka/$1_$2/results/all.somatic.snvs.vcf ;\
	sed -i -e 's/__asterisk__/\*/g' -e 's/__colon__/\:/g' strelka/$1_$2/results/passed.somatic.indels.vcf ;\
	sed -i -e 's/__asterisk__/\*/g' -e 's/__colon__/\:/g' strelka/$1_$2/results/passed.somatic.snvs.vcf")

strelka/$1_$2/results/all.somatic.snvs.vcf : strelka/$1_$2/results/all.somatic.indels.vcf
	$$(INIT) rm -rf strelka/$1_$2/chromosomes

vcf/$1_$2.%.vcf : strelka/vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$<

strelka/vcf/$1_$2.strelka_snps.vcf : strelka/$1_$2/results/all.somatic.snvs.vcf
	$$(INIT) $$(FIX_STRELKA_VCF) $$< > $$@

strelka/vcf/$1_$2.strelka_indels.vcf : strelka/$1_$2/results/all.somatic.indels.vcf
	$$(INIT) $$(FIX_STRELKA_VCF) $$< > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call strelka-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include usb-modules-v2/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


