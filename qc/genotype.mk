
##### DEFAULTS ######
LOGDIR ?= log/genotype.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

ifeq ($(shell test $(words $(SAMPLES)) -gt 2; echo $$?),0)
all : genotype/all$(PROJECT_PREFIX).snps_filtered.sdp_ft.clust.png
endif


ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
define snps-chr
genotype/all$(PROJECT_PREFIX).snps.$1.vcf : $$(foreach sample,$$(SAMPLES),gatk/dbsnp/$$(sample).gatk_snps.vcf)
	$$(call RUN,1,$$(RESOURCE_REQ_VVHIGH_MEM),$$(RESOURCE_REQ_SHORT),$$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_VVHIGH_MEM_JAVA)) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) -o $$@ --genotypemergeoption UNSORTED \
	-L $1 -R $$(REF_FASTA)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call snps-chr,$(chr))))

ifeq ($(findstring true,$(GENOTYPE_WITH_CHR2)),true)
genotype/all$(PROJECT_PREFIX).snps.vcf : genotype/all$(PROJECT_PREFIX).snps.2.vcf
	$(INIT) ln -f $< $@
else
genotype/all$(PROJECT_PREFIX).snps.vcf : $(foreach chr,$(CHROMOSOMES),genotype/all$(PROJECT_PREFIX).snps.$(chr).vcf)
	$(call RUN,1,$(RESOURCE_REQ_VVHIGH_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
		$(call GATK,CombineVariants,$(RESOURCE_REQ_VVHIGH_MEM_JAVA)) \
		$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA) &&\
		$(RM) $^")
endif
endif #end illumina

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
genotype/all$(PROJECT_PREFIX).snps.vcf : $(foreach sample,$(SAMPLES),tvc/dbsnp/$(sample)/TSVC_variants.vcf)
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(JAVA8_MODULE),"\
	$(call GATK,CombineVariants,$(RESOURCE_REQ_HIGH_MEM_JAVA)) \
	$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")
endif

genotype/all$(PROJECT_PREFIX).snps_filtered.vcf : genotype/all$(PROJECT_PREFIX).snps.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@ && $(RM) $<

genotype/all$(PROJECT_PREFIX).%.clust.png : genotype/all$(PROJECT_PREFIX).%.vcf
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_MEDIUM),$(R_MODULE),"\
	$(CLUSTER_VCF) --outPrefix genotype/all$(PROJECT_PREFIX).$* \
	$(if $(findstring RNA,$(CAPTURE_METHOD)),--clusterType hetSameHom) $<")

include usb-modules-v2/vcf_tools/vcftools.mk
include usb-modules-v2/variant_callers/TVC.mk
include usb-modules-v2/variant_callers/gatk.mk
