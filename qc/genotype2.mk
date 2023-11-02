##### DEFAULTS ######
LOGDIR ?= log/genotype2.$(NOW)

##### MAKE INCLUDES #####
include usb-modules-v2/Makefile.inc

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

ifeq ($(shell test $(words $(SAMPLES)) -gt 2; echo $$?),0)
all : genotype/BAMixChecker/Total_result.txt
endif

genotype/BAMixChecker/Total_result.txt : $(foreach sample,$(SAMPLES),bam/$(sample).bam)
	$(call RUN,6,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_VSHORT),,"\
	$(SINGULARITY_EXEC) $(BAMIXCHECKER_IMG) python /BAMixChecker-1.0.1/BAMixChecker.py -d bam \
	$(if $(findstring BAITS,$(CAPTURE_METHOD)),-b $(TARGETS_FILE_INTERVALS),if $(findstring PCR,$(CAPTURE_METHOD)),-b $(TARGETS_FILE_INTERVALS),) \
	-r $(REF_FASTA) \
	-o genotype --OFFFileNameMatching -p 6")
