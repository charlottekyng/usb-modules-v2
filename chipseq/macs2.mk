include usb-modules-v2/Makefile.inc

LOGDIR ?= log/macs2.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : macs2

macs2 : $(foreach pair,$(SAMPLE_PAIRS),macs2/$(pair)_summits.bed)

define macs2-callpeaks
macs2/$1_$2_summits.bed : bam/$1.bam bam/$2.bam
	$$(call RUN,$$(MACS2_NUM_CORES),$$(RESOURCE_REQ_LOW_MEM),$$(RESOURCE_REQ_SHORT),$$(MACS2_MODULE),"\
	$$(MACS2) callpeak -t $$(<) -c $$(<<) -f BAM -g hs -p $$(MACS2_PVALUE) -n $1_$2 \
	--tsize=$$(MACS2_READLENGTH) --outdir $$(@D)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call macs2-callpeaks,$(tumor.$(pair)),$(normal.$(pair)))))

macs2/%_peaks.narrowPeak : macs2/%_summits.bed
	

%.homer_ann.txt : %

	cut -f1-5 $< | sed 's/^/chr/' - > %.bed && |
