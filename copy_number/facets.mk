include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/variantCaller.inc

LOGDIR ?= log/facets.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : facets

SNPPILEUP_SUFFIX = q$(FACETS_SNP_PILEUP_MINMAPQ)_Q$(FACETS_SNP_PILEUP_MINBASEQ)_d$(FACETS_SNP_PILEUP_MAX_DEPTH)_r$(FACETS_SNP_PILEUP_MIN_DEPTH)
FACETS_SUFFIX = $(SNPPILEUP_SUFFIX)_bin$(FACETS_WINDOW_SIZE)_mingc$(FACETS_MINGC)_maxgc$(FACETS_MAXGC)_nhet$(FACETS_MIN_NHET)_cval$(FACETS_CVAL1)

facets : facets/cncf/all$(PROJECT_PREFIX).summary.txt facets/cncf/all$(PROJECT_PREFIX).geneCN.GL_ASCNA.txt facets/cncf/all$(PROJECT_PREFIX).cncf.txt facets/cncf/all$(PROJECT_PREFIX).cncf.pdf.tar.gz

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
facets/base_pos/%.gatk.dbsnp.vcf : gatk/dbsnp/%.gatk_snps.vcf gatk/vcf/%.variants.vcf
	$(call RUN,1,$(RESOURCE_REQ_HIGH_MEM),$(RESOURCE_REQ_SHORT),$(JAVA8_MODULE),"\
	$(call GATK,CombineVariants,$(RESOURCE_REQ_HIGH_MEM)) \
	$(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")

define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2_$$(SNPPILEUP_SUFFIX).bc.gz : bam/$1.bam bam/$2.bam $$(if $$(findstring true,$$(FACETS_GATK_VARIANTS)),facets/base_pos/$1.gatk.dbsnp.vcf,$$(FACETS_TARGETS_INTERVALS))
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),,"\
	$$(FACETS_SNP_PILEUP) \
	-A -d $$(FACETS_SNP_PILEUP_MAX_DEPTH) -g -q $$(FACETS_SNP_PILEUP_MINMAPQ) \
	-Q $$(FACETS_SNP_PILEUP_MINBASEQ) -r $$(FACETS_SNP_PILEUP_MIN_DEPTH)$$(,)0 \
	$$(word 3,$$^) $$@ $$(word 2,$$^) $$(word 1,$$^)")
endef
$(info SAMPLE_PAIRS $(SAMPLE_PAIRS))
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call snp-pileup-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
define snp-pileup-tumor-normal
facets/snp_pileup/$1_$2.bc.gz : tvc/dbsnp/$2/TSVC_variants.vcf tvc/dbsnp/$1/TSVC_variants.vcf
	$$(call RUN,1,$$(RESOURCE_REQ_HIGH_MEM),$$(RESOURCE_REQ_SHORT),$$(SNP_EFF_MODULE) $$(JAVA8_MODULE),"\
	$$(call GATK,CombineVariants,$$(RESOURCE_REQ_HIGH_MEM)) \
	$$(foreach vcf,$$^,--variant $$(vcf) ) --genotypemergeoption UNSORTED -R $$(REF_FASTA) | \
	$$(SNP_SIFT) extractFields - CHROM POS REF ALT GEN[0].FRO GEN[0].FAO GEN[0].FXX GEN[0].FXX GEN[1].FRO GEN[1].FAO GEN[1].FXX GEN[1].FXX | \
	perl -p -e \"s/$$(,)[\w]+//g;\" | sed 's/^chr//g;' | sed 's/\t\t/\t0\t/g; s/\t$$$$/\t0/g; s/\t\t/\t0\t/g; s/\t\t/\t0\t/g;' | \
	perl -p -e \"s/\#CHROM.+$$$$/Chromosome$$(,)Position$$(,)Ref$$(,)Alt$$(,)File1R$$(,)File1A$$(,)File1E$$(,)File1D$$(,)File2R$$(,)File2A$$(,)File2E$$(,)File2D/g;\" | \
	grep -v '\.' | \
	sed 's/\t/$$(,)/g;' | gzip > $$@")
endef
endif

define facets-cval1-tumor-normal
facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).out : facets/snp_pileup/$1_$2_$$(SNPPILEUP_SUFFIX).bc.gz
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(FACETS) --minNDepth $$(FACETS_SNP_PILEUP_MIN_DEPTH) \
	--maxNDepth $$(FACETS_SNP_PILEUP_MAX_DEPTH) --snp_nbhd $$(FACETS_WINDOW_SIZE) \
	--minGC $$(FACETS_MINGC) --maxGC $$(FACETS_MAXGC) \
	--cval1 $$(FACETS_CVAL1) --genome $$(REF) --min_nhet $$(FACETS_MIN_NHET) \
	--outPrefix $$* $$< ")


facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).Rdata facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).cncf.txt facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).cncf.pdf facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).logR.pdf: facets/cncfTN/$1_$2_$$(FACETS_SUFFIX).out
	

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call facets-cval1-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


define facets-ln-files
ifdef FACETS_SPECIAL_CASES
FACETS_SUFFIX_FINAL = $$(shell grep $1_$2 $$(FACETS_SPECIAL_CASES) | cut -f2)
endif
FACETS_SUFFIX_FINAL ?= $$(FACETS_SUFFIX)


facets/cncf/$1_$2.out : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).out
	$$(INIT) ln -f $$(<) $$(@)

facets/cncf/$1_$2.Rdata : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).Rdata
	$$(INIT) ln -f $$(<) $$(@)

facets/cncf/$1_$2.cncf.txt : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).cncf.txt
	$$(INIT) ln -f $$(<) $$(@)

facets/cncf/$1_$2.cncf.pdf : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).cncf.pdf
	$$(INIT) ln -f $$(<) $$(@)

facets/cncf/$1_$2.logR.pdf : facets/cncfTN/$1_$2_$$(FACETS_SUFFIX_FINAL).logR.pdf
	$$(INIT) ln -f $$(<) $$(@)
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call facets-ln-files,$(tumor.$(pair)),$(normal.$(pair)))))
 
facets/cncf/all$(PROJECT_PREFIX).summary.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).out)
	$(INIT) \
	{ \
	cut -f2 -d' ' $< | tr '\n' '\t'; echo ""; \
	for metrics in $^; do \
		cut -f4 -d' ' $$metrics | tr '\n' '\t'; echo ""; \
	done;\
} >$@

facets/cncf/all$(PROJECT_PREFIX).geneCN.GL_ASCNA.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt) $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE),"\
	$(FACETS_GENE_CN) $(FACETS_GENE_CN_OPTS) --outFile $(@D)/all$(PROJECT_PREFIX).geneCN $(filter %.cncf.txt,$^)")

facets/cncf/all$(PROJECT_PREFIX).cncf.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.txt)
	$(INIT) head -1 $< | sed 's/^/TUMOR_NORMAL\t/' > $@; \
	for cncf in $^; do \
		samplename=`basename $$cncf | sed 's/\.cncf\.txt//;'`; \
		sed '/^chrom/d; s/^23/X/;' $$cncf | sed "s/^/$$samplename\t/" >> $@; \
	done

facets/cncf/all$(PROJECT_PREFIX).cncf.pdf.tar.gz : $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).cncf.pdf) $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).logR.pdf)
	$(INIT) tar -czf $@ $^

include usb-modules-v2/variant_callers/TVC.mk
include usb-modules-v2/bam_tools/processBam.mk
include usb-modules-v2/variant_callers/gatkVariantCaller.mk
