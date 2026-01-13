include usb-modules-v2/Makefile.inc
include usb-modules-v2/variant_callers/somatic/somaticVariantCaller.inc
include usb-modules-v2/sv_callers/conSV.mk

LOGDIR ?= log/sig_profiler_assignment.$(NOW)
CALLER_PREFIX ?= mutect2
ALLOWED_SIG_TYPES := SBS DBS ID CN SV
SIG_TYPE ?= SBS
.DEFAULT_GOAL := sig_profiler_assignment

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: sig_profiler_assignment

ifeq ($(SIG_TYPE), SBS)
	ifneq ($(words $(CALLER_PREFIX)),1)
		$(info CALLER_PREFIX contains more than one variant caller)
		$(info Choose only one by executing: make sig_profiler_assignment CALLER_PREFIX=<variant caller>)
		$(info  )
		$(error Stopping due to invalid CALLER_PREFIX)
	endif
endif

ifeq ($(filter $(SIG_TYPE),$(ALLOWED_SIG_TYPES)),)
	$(info SIG_TYPE must be one of: $(ALLOWED_SIG_TYPES).)
	$(info Choose the signature type by executing: make sig_profiler_assignment SIG_TYPE=DBS)
	$(info  )
	$(error Stopping due to invalid SIG_TYPE)
endif

sig_profiler_assignment : SigProfilerAssignment/$(SIG_TYPE)/output/JOB_METADATA_SPA.txt

ifeq ($(SIG_TYPE),SBS)
SigProfilerAssignment/SBS/output/JOB_METADATA_SPA.txt :$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(CALLER_PREFIX).*.hotspot.pass.vcf)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(SIG_PROFILER_MODULE),"\
	$(RM) -r SigProfilerAssignment/SBS/* && mkdir -p SigProfilerAssignment/SBS/input && mkdir -p SigProfilerAssignment/SBS/output && \
	ln $^ SigProfilerAssignment/SBS/input/ && \
	$(SIG_PROFILER_ASSIGNMENT) --samples SigProfilerAssignment/SBS/input --output SigProfilerAssignment/SBS/output \
	--cosmic_version $(SIG_PROFILER_COSMIC_VERSION) \
	$(if $(findstring True,$(SIG_PROFILER_COSMIC_EXOME)),--exome) \
	--genome_build $(SIG_PROFILER_COSMIC_GENOME) \
	$(if $(SIG_PROFILER_COSMIC_SIGNATURE_DB),--signature_database $(SIG_PROFILER_COSMIC_SIGNATURE_DB)) \
	$(if $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS),--signature_database $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS))")
endif

ifeq ($(SIG_TYPE),DBS)
SigProfilerAssignment/DBS/output/JOB_METADATA_SPA.txt :$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(CALLER_PREFIX).*.hotspot.pass.vcf)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(SIG_PROFILER_MODULE),"\
	$(RM) -r SigProfilerAssignment/DBS/* && mkdir -p SigProfilerAssignment/DBS/input && mkdir -p SigProfilerAssignment/DBS/output && \
	ln $^ SigProfilerAssignment/DBS/input/ && \
	$(SIG_PROFILER_ASSIGNMENT) --samples SigProfilerAssignment/DBS/input --output SigProfilerAssignment/DBS/output \
	--cosmic_version $(SIG_PROFILER_COSMIC_VERSION) \
	$(if $(findstring True,$(SIG_PROFILER_COSMIC_EXOME)),--exome) \
	--genome_build $(SIG_PROFILER_COSMIC_GENOME) \
	--context_type DINUC \
	$(if $(SIG_PROFILER_COSMIC_SIGNATURE_DB),--signature_database $(SIG_PROFILER_COSMIC_SIGNATURE_DB)) \
	$(if $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS),--signature_database $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS))")
endif

ifeq ($(SIG_TYPE),ID)
SigProfilerAssignment/ID/output/JOB_METADATA_SPA.txt :$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(CALLER_PREFIX).*.hotspot.pass.vcf)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(SIG_PROFILER_MODULE),"\
	$(RM) -r SigProfilerAssignment/ID/* && mkdir -p SigProfilerAssignment/ID/input && mkdir -p SigProfilerAssignment/ID/output && \
	ln $^ SigProfilerAssignment/ID/input/ && \
	$(SIG_PROFILER_ASSIGNMENT) --samples SigProfilerAssignment/ID/input --output SigProfilerAssignment/ID/output \
	--cosmic_version $(SIG_PROFILER_COSMIC_VERSION) \
	$(if $(findstring True,$(SIG_PROFILER_COSMIC_EXOME)),--exome) \
	--genome_build $(SIG_PROFILER_COSMIC_GENOME) \
	--context_type ID \
	$(if $(SIG_PROFILER_COSMIC_SIGNATURE_DB),--signature_database $(SIG_PROFILER_COSMIC_SIGNATURE_DB)) \
	$(if $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS),--signature_database $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS))")
endif

ifeq ($(SIG_TYPE),CN)

SigProfilerAssignment/CN/input/%.cncf.txt : facets/cncf/%.cncf.txt
	$(MKDIR) $(@D)
	$(INIT) name=$$(basename $< .cncf.txt | cut -d_ -f1); \
	awk -v name=$$name 'BEGIN{OFS="\t"} \
		NR==1{$$0="sample\t"$$0; print; next} \
		{print name$(,) $$0}' $< > $@

SigProfilerAssignment/CN/input/CNCF.txt : $(foreach pair,$(SAMPLE_PAIRS),SigProfilerAssignment/CN/input/$(pair).cncf.txt)
	$(INIT) awk 'FNR==1 && NR!=1 {next} {print}' $^ > $@

SigProfilerAssignment/CN/output/JOB_METADATA_SPA.txt : SigProfilerAssignment/CN/input/CNCF.txt
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(SIG_PROFILER_MODULE),"\
	mkdir -p SigProfilerAssignment/CN/output && \
	$(SIG_PROFILER_ASSIGNMENT) --samples $< --output SigProfilerAssignment/CN/output \
	--cosmic_version $(SIG_PROFILER_COSMIC_VERSION) \
	--genome_build $(SIG_PROFILER_COSMIC_GENOME) \
	--input_type "seg:FACETS" \
	$(if $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS),--signature_database $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS))")

endif

ifeq ($(SIG_TYPE),SV)
#could be implemented also for single callers (but practically I don't see the use case)
SigProfilerAssignment/SV/input/%.bedpe : conSV/conSV_$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))/%.$(MIN_CALLERS)_outof_$(words $(SV_CALLERS))_slope_$(SLOPE).tsv
	$(MKDIR) $(@D)
	$(INIT) \
	sample=$$(basename $< | cut -d. -f1); \
	awk -v sample=$$sample 'BEGIN{OFS="\t"} \
		NR==1 { \
			print "sample","chrom1","start1","end1","chrom2","start2","end2","svclass"; \
			for (i=1; i<=NF; i++) col[$$i]=i; \
			next \
		} \
		{ \
			print sample, \
			      $$col["CHROM"], \
			      $$col["START1"], \
			      $$col["END1"], \
			      $$col["CHROM2"], \
			      $$col["START2"], \
			      $$col["END2"], \
			      $$col["CONS_SVTYPE"] \
		}' $< > $@

SigProfilerAssignment/SV/output/SV.SV32.matrix.tsv : $(foreach pair,$(SAMPLE_PAIRS),SigProfilerAssignment/SV/input/$(tumor.$(pair)).bedpe)
	$(MKDIR) $(@D)
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(SIG_PROFILER_MODULE),"\
	$(SV_MATRIX_GENERATOR) \
	--input_dir SigProfilerAssignment/SV/input/ \
	--project SV \
	--output_dir SigProfilerAssignment/SV/output/")

SigProfilerAssignment/SV/output/JOB_METADATA_SPA.txt : SigProfilerAssignment/SV/output/SV.SV32.matrix.tsv
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_MEDIUM),$(SIG_PROFILER_MODULE),"\
	$(SIG_PROFILER_ASSIGNMENT) --samples $< --output SigProfilerAssignment/SV/output \
	--cosmic_version $(SIG_PROFILER_COSMIC_VERSION) \
	--context_type SV \
	--input_type matrix \
	--genome_build $(SIG_PROFILER_COSMIC_GENOME) \
	$(if $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS),--signature_database $(SIG_PROFILER_COSMIC_EXCLUDE_SIG_SUBGROUPS))")

endif