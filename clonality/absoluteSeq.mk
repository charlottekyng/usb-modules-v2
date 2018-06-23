include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/absoluteSeq.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: absolute
#absolute_rdata absolute_reviewed absolute_tables

whoami = $(shell whoami)

absolute : absolute/step3/reviewed/all$(PROJECT_PREFIX).$(whoami).ABSOLUTE.table.txt

define absolute_make_mutations
absolute/mutations/$1_$2.mutations.txt : $$(foreach prefix,$$(CALLER_PREFIX),tables/$1_$2.$$(call DOWMSTREAM_VCF_TABLE_SUFFIX,$$(prefix)).txt)
	$$(MKDIR) absolute/mutations; \
	sample=`basename $$@ .mutations.txt`; \
	$$(call RUN,1,$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(RBIND) --tumorNormal $$^ > $$@.tmp1 && \
	$$(MUTATION_SUMMARY_RSCRIPT) --outFile $$@.tmp2 --outputFormat TXT $$@.tmp1 && \
	$$(ABSOLUTE_MAKE_MUTS) --sample $$$$sample --outFile $$@ $$@.tmp2 && \
	$$(RM) $$@.tmp1 $$@.tmp2")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call absolute_make_mutations,$(tumor.$(pair)),$(normal.$(pair)))))

absolute/segments/%.seg.txt : facets/cncf/%.cncf.txt
	$(INIT) $(MKDIR) absolute/segments; \
	sample=`basename $$@ .seg.txt`; ml $(R_MODULE); \
	$(ABSOLUTE_MAKE_SEGS) --sample $$sample --outFile $@ $<

absolute/absolute_parameters.txt : $(foreach pair,$(SAMPLE_PAIRS),absolute/mutations/$(pair).mutations.txt absolute/segments/$(pair).seg.txt)
	$(INIT) echo "include patient sample seg.dat.fn sigma.p max.sigma.h min.ploidy max.ploidy \
	primary.disease platform sample.name results.dir max.as.seg.count copy_num_type \
	max.neg.genome max.non.clonal maf.fn min.mut.af output.fn.base" | perl -p -e "s/\s+/\t/g;" > $@ && \
	echo "" >> $@ && \
	for sample_pair in $(SAMPLE_PAIRS); do \
		echo "Y $$sample_pair $$sample_pair absolute/segments/$$sample_pair.seg.txt \
		$(ABSOLUTE_SIGMA_P) $(ABSOLUTE_MAX_SIGMA_H) $(ABSOLUTE_MIN_PLOIDY) $(ABSOLUTE_MAX_PLOIDY) \
		$(ABSOLUTE_DISEASE) $(ABSOLUTE_PLATFORM) $$sample_pair absolute/step1/$$sample_pair \
		$(ABSOLUTE_MAX_SEG) $(ABSOLUTE_COPYNUMTYPE) $(ABSOLUTE_MAX_NEG_GENOME) \
		$(ABSOLUTE_MAX_NON_CLONAL) absolute/mutations/$$sample_pair.mutations.txt 0 $$sample_pair" | \
		perl -p -e "s/\s+/\t/g;" >> $@ && echo "" >> $@; \
	done

define absolute_step1
absolute/step1/$1_$2/$1_$2.ABSOLUTE.RData : absolute/absolute_parameters.txt absolute/mutations/$1_$2.mutations.txt absolute/segments/$1_$2.seg.txt 
	$$(MKDIR) absolute/step1; $$(MKDIR) $$(@D); \
	sample=`basename $$@ .ABSOLUTE.RData`;\
	$$(call RUN,$$(ABSOLUTE_NUM_CORE),$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_SHORT),$$(R_MODULE),"\
	$$(ABSOLUTE_STEP1) --params $$< --numCores $$(ABSOLUTE_NUM_CORE) --sample $$$$sample")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call absolute_step1,$(tumor.$(pair)),$(normal.$(pair)))))

absolute/step2/all$(PROJECT_PREFIX).PP-modes.data.RData : $(foreach pair,$(SAMPLE_PAIRS),absolute/step1/$(pair)/$(pair).ABSOLUTE.RData)
	$(MKDIR) absolute/step2; \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(ABSOLUTE_STEP2) --outdir $(@D) --obj.name all$(PROJECT_PREFIX) $^")

absolute/step2/all$(PROJECT_PREFIX).PP-calls_tab.txt : absolute/step2/all$(PROJECT_PREFIX).PP-modes.data.RData
	

absolute/step3/reviewed/all$(PROJECT_PREFIX).$(whoami).ABSOLUTE.table.txt : absolute/step2/all$(PROJECT_PREFIX).PP-modes.data.RData absolute/step2/all$(PROJECT_PREFIX).PP-calls_tab.txt
	$(MKDIR) absolute/step3; \
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(R_MODULE),"\
	$(ABSOLUTE_STEP3) --obj.name all$(PROJECT_PREFIX) --analyst $(whoami) \
	--outdir absolute/step3 --modes.fn $(<) --pp.calls $(<<)")








