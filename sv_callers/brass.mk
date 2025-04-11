# run brass
MUT_CALLER = brass

include usb-modules-v2/Makefile.inc

LOGDIR = log/brass.$(NOW)
PHONY += brass

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

HIGH_DEPTH_BED = brass_files/merged/depth_mask.bed
VIRAL = $(IMG_DIR)/BRASS/ref/viral/2bit/viral.genomic.fa.2bit
MICROBIAL = $(IMG_DIR)/BRASS/ref/bacterial/2bit
CENTTEL = $(IMG_DIR)/BRASS/ref/centtel.bed
GCBINS = $(IMG_DIR)/BRASS/ref/gcBins.bed.gz.tbi
CYTOBAND = $(IMG_DIR)/BRASS/ref/cytoband.bed
GENOME_CACHE = $(GROUP_DIR)/ref/genomes/hg38/cache/vagrent.human.GRCh38.homo_sapiens_core_80_38.cache.gz
BRASS_PROTOCOL = WGS # WGS|WXS|RNA
SAMPSTAT = brass_files/sampstat

BRASS_OPTS = -g $(REF_FASTA) -s HUMAN -as 38 -pr $(BRASS_PROTOCOL) -gc $(GENOME_CACHE) -d $(HIGH_DEPTH_BED) -vi $(VIRAL) -mi $(MICROBIAL) -ct $(CENTTEL) -b $(GCBINS) -cb $(CYTOBAND)

brass : $(foreach pair,$(SAMPLE_PAIRS),brass/$(pair).brass_timestamp)

bam/%.bam.bas : bam/%.bam
	$(call RUN,1,6G,$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	$(BRASS) bam_stats \
	-i $< \
	-o $@ \
	-r $(REF_FASTA).fai")

brass/genome_window.bed :
	$(call RUN,1,2G,$(RESOURCE_REQ_MEDIUM),$(SINGULARITY_MODULE),"\
	$(BRASS) bedtools makewindows \
	-g $(REF_FASTA) \
	-w 10000 > $@")

define brass-tumor-normal
# Main call
brass/$1_$2.brass_input_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas
	$$(call RUN,2,6G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) brass.pl $$(BRASS_OPTS) \
	-ss $$(SAMPSTAT)/$2.txt \
	-t $$< \
	-n $$(<<) \
	-tn $1 \
	-nn $2 \
	-c 2 \
	-p input \
	-o brass/$1_$2 && touch $$@")

brass/$1_$2.brass_group_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_input_timestamp
	$$(call RUN,1,8G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) brass.pl $$(BRASS_OPTS) \
	-t $$< \
	-n $$(<<) \
	-tn $1 \
	-nn $2 \
	-ss $$(SAMPSTAT)/$2.txt \
	-p group \
	-o brass/$1_$2 && touch $$@")

brass/$1_$2.brass_filter_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_group_timestamp
	$$(call RUN,1,8G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) brass.pl $$(BRASS_OPTS) \
	-t $$< \
	-n $$(<<) \
	-tn $1 \
	-nn $2 \
	-ss $$(SAMPSTAT)/$2.txt \
	-p group \
	-o brass/$1_$2 && touch $$@")

brass/$1_$2.brass_assemble_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_filter_timestamp
	$$(call RUN,8,2G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) brass.pl $$(BRASS_OPTS) \
	-t $$< \
	-n $$(<<) \
	-tn $1 \
	-nn $2 \
	-ss $$(SAMPSTAT)/$2.txt \
	-c 8 \
	-p assemble \
	-o brass/$1_$2 && touch $$@")

brass/$1_$2.brass_grass_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_assemble_timestamp
	$$(call RUN,1,8G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) brass.pl $$(BRASS_OPTS) \
	-t $$< \
	-n $$(<<) \
	-tn $1 \
	-nn $2 \
	-ss $$(SAMPSTAT)/$2.txt \
	-p grass \
	-o brass/$1_$2 && touch $$@")

brass/$1_$2.brass_tabix_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_grass_timestamp
	$$(call RUN,1,2G,$$(RESOURCE_REQ_MEDIUM),$$(SINGULARITY_MODULE),"\
	$$(BRASS) brass.pl $$(BRASS_OPTS) \
	-t $$< \
	-n $$(<<) \
	-tn $1 \
	-nn $2 \
	-p tabix \
	-ss $$(SAMPSTAT)/$2.txt \
	-o brass/$1_$2 && touch $$@")

brass/$1_$2.brass_timestamp : brass/$1_$2.brass_tabix_timestamp
	$$(INIT) touch $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call brass-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

.PHONY: $(PHONY)