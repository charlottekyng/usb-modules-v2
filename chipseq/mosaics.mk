include usb-modules-v2/Makefile.inc

LOGDIR ?= log/mosaics.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : mosaics

mosaics : $(foreach pair,$(SAMPLE_PAIRS),mosaics/rdata_fragL$(MOSAICS_FRAG_LEN)_bin$(MOSAICS_BIN_SIZE)/$(pair).rdata mosaics/peaks_fragL$(MOSAICS_FRAG_LEN)_bin$(MOSAICS_BIN_SIZE)/$(pair).peakTFBS.annotated.bed)
#$(foreach sample,$(SAMPLES),mosaics/bin/$(sample).bam_fragL$(MOSAICS_FRAG_LEN)_bin$(MOSAICS_BIN_SIZE).txt)
#$(foreach pair,$(SAMPLE_PAIRS),mosaics/rdata/$(pair).rdata)

define mosaics-peaks
mosaics/rdata_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.rdata : mosaics/bin/$1.bam_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE).txt mosaics/bin/$2.bam_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE).txt
	$$(call RUN,$$(MOSAICS_NUM_CORES),$$(RESOURCE_REQ_VHIGH_MEM),$$(RESOURCE_REQ_SHORT),$$(R_MODULE) $$(BEDOOLS_MODULE),"\
	$$(MKDIR) mosaics mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE) \
	mosaics/rdata_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE) mosaics/plots_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE); \
	$$(MOSAICS_RUN) --parallel $$(MOSAICS_PARALLEL) --num_cores $$(MOSAICS_NUM_CORES) \
	--plotFileLoc mosaics/plots_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE) \
	--peakFileLoc mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE) \
	--rdataFileLoc mosaics/rdata_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE) \
	--maxgap $$(MOSAICS_MAXGAP) --minsize $$(MOSAICS_MINSIZE) --thres $$(MOSAICS_THRES) $$^ && \
	$$(BEDTOOLS) sort -i mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.peakTFBS.bed -faidx $$(REF_FASTA).fai > tmp && \
	mv tmp mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.peakTFBS.bed && \
	head -1 mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.peakTFBS.txt > tmp && \
	tail -n +2 mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.peakTFBS.txt | \
	$$(BEDTOOLS) sort -faidx $$(REF_FASTA).fai >> tmp && mv tmp mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.peakTFBS.txt")

mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.peakTFBS.bed : mosaics/rdata_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.rdata
	

mosaics/peaks_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.peakTFBS.txt : mosaics/rdata_fragL$$(MOSAICS_FRAG_LEN)_bin$$(MOSAICS_BIN_SIZE)/$1_$2.rdata
	
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call mosaics-peaks,$(tumor.$(pair)),$(normal.$(pair)))))

mosaics/%.annotated.bed : mosaics/%.bed
	$(call RUN,1,$(RESOURCE_REQ_LOW_MEM),$(RESOURCE_REQ_VSHORT),$(BEDTOOLS_MODULE),"\
	$(BEDTOOLS) closest -g chr_order -k 3 -t all -D b -a $< -b $(GENCODE_GENE_GTF) > $@")

mosaics/bin/%.bam_fragL$(MOSAICS_FRAG_LEN)_bin$(MOSAICS_BIN_SIZE).txt : bam/%.bam
	$(call RUN,1,$(RESOURCE_REQ_MEDIUM_MEM),$(RESOURCE_REQ_SHORT),$(R_MODULE) $(PERL_MODULE),"\
	$(MKDIR) mosaics mosaics/bin mosaics/wig; \
	$(MOSAICS_CONSTRUCTBINS) --fragLen $(MOSAICS_FRAG_LEN) --binSize $(MOSAICS_BIN_SIZE) --pet $(PAIRED_END) \
	--outfileLoc mosaics/bin $<")
