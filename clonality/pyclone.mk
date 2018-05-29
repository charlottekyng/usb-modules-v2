include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR ?= log/pyclone.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : pyclone

ifeq ($(findstring ILLUMINA,$(SEQ_PLATFORM)),ILLUMINA)
CALLER_PREFIX ?= mutect strelka_indels
endif
ifeq ($(findstring IONTORRENT,$(SEQ_PLATFORM)),IONTORRENT)
CALLER_PREFIX ?= tvc_snps tvc_indels
endif

ifeq ($(findstring false,$(INCLUDE_LNCRNA_IN_SUMMARY)),false)
TABLE_SUFFIX = tab.nonsynonymous_synonymous_hotspot
endif
ifeq ($(findstring true,$(INCLUDE_LNCRNA_IN_SUMMARY)),true)
TABLE_SUFFIX = tab.nonsynonymous_synonymous_hotspot_lincRNA_upstream
endif 

pyclone : $(foreach normal_sample,$(NORMAL_SAMPLES),pyclone/tables/$(normal_sample).cluster.txt)

pyclone/tables/%.cluster.txt : pyclone/configs/%.yaml pyclone/runs/%/alpha.tsv.bz2
	$(INIT) $(MKDIR) pyclone/tables; \
	$(PYCLONE) build_table --config_file $< --table_type cluster \
	--out_file $@ --burnin $(PYCLONE_BURNIN) && \
	$(PYCLONE) build_table --config_file $< --table_type loci \
	--out_file $(subst cluster,loci,$@) --burnin $(PYCLONE_BURNIN) && \
	$(PYTHON_ENV_DEACTIVATE)

pyclone/runs/%/alpha.tsv.bz2 : pyclone/configs/%.yaml 
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),,"\
	$(PYCLONE) run_analysis --config_file $< --seed $(PYCLONE_SEED) && $(PYTHON_ENV_DEACTIVATE)")

define pyclone_make_config
pyclone/configs/$1.yaml : $$(foreach tumor,$2,facets/cncf/$$(tumor)_$1.out pyclone/mutations/$$(tumor)_$1.mutations.yaml)
	$$(INIT) $$(MKDIR) pyclone/configs; \
	echo -n "num_iters: " > $$@; \
	echo $$(PYCLONE_ITER) >> $$@; \
	echo "base_measure_params:" >> $$@; \
	echo "  alpha: 1" >> $$@; \
	echo "  beta: 1" >> $$@; \
	echo "concentration:" >> $$@; \
	echo "  value: 1.0" >> $$@; \
	echo "  prior:" >> $$@; \
	echo "    shape: 1.0" >> $$@; \
	echo "    rate: 0.001" >> $$@; \
	echo "density: pyclone_beta_binomial" >> $$@; \
	echo "beta_binomial_precision_params:" >> $$@; \
	echo "  value: 1000" >> $$@; \
	echo "  prior:" >> $$@; \
	echo "    shape: 1.0" >> $$@; \
	echo "    rate: 0.0001" >> $$@; \
	echo "  proposal:" >> $$@; \
	echo "    precision: 0.01" >> $$@; \
	echo -n "working_dir: " >> $$@; echo `pwd` >> $$@; \
	echo -n "trace_dir: " >> $$@; echo "pyclone/runs/$1" >> $$@; \
	echo "samples:" >> $$@; \
	for cncf in $$(filter %.out,$$^); do \
		samplename=`basename $$$$cncf | cut -f1 -d'_'`; \
		tnname=`basename $$$$cncf | cut -f1 -d'.'`; \
		echo "  $$$$samplename:" >> $$@; \
		echo "    tumour_content: " >> $$@; \
		echo -n "      value: " >> $$@; \
		echo `grep Purity $$$$cncf | cut -f2 -d'=' | sed 's/NA/0.1/;'` >> $$@; \
		echo -n "    mutations_file: pyclone/mutations/" >> $$@; \
		echo "$$$${tnname}.mutations.yaml" >> $$@; \
		echo "    error_rate: 0.001"  >> $$@; \
	done
endef
$(foreach set,$(SAMPLE_SETS),\
	$(eval $(call pyclone_make_config,$(lastword $(subst _,$( ),$(set))),\
	$(wordlist 1,$(shell expr $(words $(subst _,$( ),$(set))) - 1),$(subst _,$( ),$(set))))))

pyclone/mutations/%.mutations.yaml : pyclone/mutations/%.mutations.txt
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_SHORT),,"\
	$(PYCLONE) build_mutations_file --prior $(PYCLONE_PRIOR) --in_file $< --out_file $@ \
	&& $(PYTHON_ENV_DEACTIVATE)")

define pyclone_make_mutations
pyclone/mutations/$1_$2.mutations.txt : $$(foreach prefix,$$(CALLER_PREFIX),tables/$1_$2.$$(call SOMATIC_VCF_SUFFIXES,$$(prefix)).$$(TABLE_SUFFIX).txt)
	$$(MKDIR) pyclone/mutations; \
	$$(call RUN,1,$$(RESOURCE_REQ_MEDIUM_MEM),$$(RESOURCE_REQ_VSHORT),$$(R_MODULE),"\
	$$(RBIND) --tumorNormal $$^ > $$@.tmp1 && \
	$$(MUTATION_SUMMARY_RSCRIPT) --outFile $$@.tmp2 --outputFormat TXT $$@.tmp1 && \
	$$(PYCLONE_MAKE_MUT_TXT) --outFile $$@ $$@.tmp2 && \
	$$(RM) $$@.tmp1 $$@.tmp2")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call pyclone_make_mutations,$(tumor.$(pair)),$(normal.$(pair)))))




